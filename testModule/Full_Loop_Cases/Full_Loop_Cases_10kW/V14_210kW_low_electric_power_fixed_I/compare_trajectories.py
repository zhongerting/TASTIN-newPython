"""Compare low-power trajectories from existing history_control.csv files."""

import argparse
import bisect
import csv
import json
import math
from pathlib import Path


TEMPERATURE_COLUMNS = {
    "core_outlet_T_K": "core_outlet_T",
    "hp_evaporator_temperature_mean_K": "hp_evaporator_temperature_mean",
    "hp_condenser_temperature_mean_K": "hp_condenser_temperature_mean",
    "collector_ring_wall_temperature_mean_K": "collector_ring_wall_temperature_mean",
    "radiator_fin_temperature_mean_K": "radiator_fin_temperature_mean",
}
REQUIRED_COLUMNS = {
    "elapsed_s",
    "electric_power_ratio",
    "thermal_power_setpoint_W",
    "flow_setpoint_kg_s",
    *TEMPERATURE_COLUMNS,
}


def _linear_slope(times, values):
    count = len(times)
    if count < 2:
        return None
    mean_t = sum(times) / count
    mean_y = sum(values) / count
    denominator = sum((time - mean_t) ** 2 for time in times)
    if denominator == 0.0:
        return None
    return sum(
        (time - mean_t) * (value - mean_y)
        for time, value in zip(times, values)
    ) / denominator


def _peak_window_slope(times, values, window_s):
    peak = None
    for start, start_time in enumerate(times):
        end = bisect.bisect_left(times, start_time + window_s, lo=start + 1)
        if end == len(times):
            break
        slope = _linear_slope(times[start : end + 1], values[start : end + 1])
        if slope is not None:
            peak = max(peak or 0.0, abs(slope))
    return peak


def _tv_over_net_change(values):
    net_change = abs(values[-1] - values[0])
    if net_change == 0.0:
        return None
    return sum(abs(right - left) for left, right in zip(values, values[1:])) / net_change


def _load_numeric_history(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")
        rows = []
        numeric_columns = REQUIRED_COLUMNS | {"tec_main_electric_power_W"}
        for line_number, raw in enumerate(reader, start=2):
            try:
                row = {
                    key: float(raw[key])
                    for key in numeric_columns
                    if key in raw and raw[key] != ""
                }
            except (TypeError, ValueError) as exc:
                raise ValueError(f"non-numeric history value at line {line_number}") from exc
            if not all(math.isfinite(value) for value in row.values()):
                raise ValueError(f"non-finite history value at line {line_number}")
            rows.append(row)
    if len(rows) < 2:
        raise ValueError("history must contain at least two records")
    times = [row["elapsed_s"] for row in rows]
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ValueError("elapsed_s must be strictly increasing")
    return rows


def analyze_history(history_path):
    history_path = Path(history_path)
    manifest_path = history_path.with_name("run_manifest.json")
    if not manifest_path.is_file():
        raise ValueError(f"missing run_manifest.json beside {history_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    trajectory = manifest["trajectory"]
    freeze_elapsed_s = float(trajectory["hold_before_ramp_s"]) + float(
        trajectory["ramp_duration_s"]
    )
    rows = _load_numeric_history(history_path)
    frozen = [row for row in rows if row["elapsed_s"] >= freeze_elapsed_s]
    if not frozen:
        raise ValueError("history ends before the configured freeze time")

    times = [row["elapsed_s"] for row in rows]
    ratios = [row["electric_power_ratio"] for row in rows]
    frozen_ratios = [row["electric_power_ratio"] for row in frozen]
    final_start_s = max(freeze_elapsed_s, times[-1] - 100.0)
    final_rows = [row for row in rows if row["elapsed_s"] >= final_start_s]
    final_times = [row["elapsed_s"] for row in final_rows]

    final_power_w = float(trajectory["final_power_w"])
    final_flow = float(trajectory["final_flow_kg_s"])
    setpoints_frozen = all(
        math.isclose(row["thermal_power_setpoint_W"], final_power_w, rel_tol=1e-10, abs_tol=1e-8)
        and math.isclose(row["flow_setpoint_kg_s"], final_flow, rel_tol=1e-10, abs_tol=1e-10)
        for row in frozen
    )

    result = {
        "history": str(history_path),
        "freeze_elapsed_s": freeze_elapsed_s,
        "frozen_duration_s": times[-1] - freeze_elapsed_s,
        "endpoint_power_compliant": 0.38 <= ratios[-1] <= 0.42,
        "frozen_power_band_compliant": all(0.35 <= ratio <= 0.45 for ratio in frozen_ratios),
        "setpoints_frozen": setpoints_frozen,
        "power_ratio_end": ratios[-1],
        "power_ratio_frozen_min": min(frozen_ratios),
        "power_ratio_frozen_max": max(frozen_ratios),
        "power_ratio_final_100s_slope_per_s": _linear_slope(
            final_times, [row["electric_power_ratio"] for row in final_rows]
        ),
        "power_ratio_tv_over_net_change": _tv_over_net_change(ratios),
        "power_ratio_peak_60s_slope_per_s": _peak_window_slope(times, ratios, 60.0),
        "core_outlet_T_max_K": max(row["core_outlet_T_K"] for row in rows),
        "core_outlet_T_final_100s_slope_K_per_s": _linear_slope(
            final_times, [row["core_outlet_T_K"] for row in final_rows]
        ),
    }

    if all("tec_main_electric_power_W" in row for row in rows):
        powers = [row["tec_main_electric_power_W"] for row in rows]
        result.update(
            electric_power_end_W=powers[-1],
            electric_power_final_100s_slope_W_per_s=_linear_slope(
                final_times,
                [row["tec_main_electric_power_W"] for row in final_rows],
            ),
            electric_power_tv_over_net_change=_tv_over_net_change(powers),
            electric_power_peak_60s_slope_W_per_s=_peak_window_slope(
                times, powers, 60.0
            ),
        )

    for column, label in TEMPERATURE_COLUMNS.items():
        values = [row[column] for row in rows]
        result[f"{label}_end_K"] = values[-1]
        result[f"{label}_delta_K"] = values[-1] - values[0]
        result[f"{label}_final_100s_slope_K_per_s"] = _linear_slope(
            final_times, [row[column] for row in final_rows]
        )
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histories", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, help="optional JSON output path")
    args = parser.parse_args(argv)
    results = [analyze_history(path) for path in args.histories]
    payload = json.dumps(results, indent=2, ensure_ascii=False, allow_nan=False)
    if args.output:
        args.output.write_text(payload + "\n", encoding="utf-8")
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
