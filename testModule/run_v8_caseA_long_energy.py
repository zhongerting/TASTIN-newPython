import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from audit_v7_caseA_global_storage_continue import (
    advance_interval,
    build_record,
    capture_storage_state,
    external_power_snapshot,
    finite_difference_storage_rates,
)
from run_v8_caseA_common import (
    DEFAULT_SOLID_ODE_METHOD,
    build_loaded_case,
    get_solid_ode_methods,
    json_default,
    parse_solid_ode_method,
    parse_v8_multipliers,
    passive_tec_source_totals,
)


DEFAULT_RESTART = "testModule/v8_caseA_migrated/v8_caseA_migrated_latest_restart.npz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V8 CaseA with append-only radial energy history.")
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v8_caseA_long_energy")
    parser.add_argument("--case-prefix", default="v8_caseA_long_energy")
    parser.add_argument("--duration", type=float, default=20000.0)
    parser.add_argument("--record-interval", type=float, default=200.0)
    parser.add_argument("--restart-interval", type=float, default=200.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument(
        "--solid-ode-method",
        type=parse_solid_ode_method,
        default=DEFAULT_SOLID_ODE_METHOD,
    )
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=False),
        default=parse_v8_multipliers("1,6,12,15,3"),
    )
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=True),
        default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True),
    )
    return parser.parse_args()


def read_existing_history_state(
    history_path: Path,
) -> Tuple[Optional[List[str]], Optional[float], Optional[float]]:
    if not history_path.exists():
        return None, None, None
    with history_path.open("r", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None, None, None
    first_absolute_time = float(rows[0]["absolute_time_s"])
    first_relative_time = float(rows[0]["relative_time_s"])
    return (
        list(rows[0].keys()),
        float(rows[-1]["absolute_time_s"]),
        first_absolute_time - first_relative_time,
    )


def append_row(history_path: Path, fieldnames: List[str], row: Dict[str, Any], *, write_header: bool) -> None:
    with history_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0:
        raise ValueError("duration must be positive.")
    if args.record_interval <= 0.0:
        raise ValueError("record-interval must be positive.")
    if abs(float(args.restart_interval) - float(args.record_interval)) > 1.0e-10:
        raise ValueError("V8 latest restart is written once per record; restart-interval must equal record-interval.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_energy_history.csv"
    summary_path = output_dir / f"{args.case_prefix}_latest_state.json"
    restart_path = output_dir / f"{args.case_prefix}_latest_restart.npz"

    build = build_loaded_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)
    existing_fields, last_recorded_time, recorded_origin_time = read_existing_history_state(history_path)
    if last_recorded_time is not None and last_recorded_time > start_time + 1.0e-8:
        raise ValueError(
            f"History already extends to {last_recorded_time}, beyond restart time {start_time}. "
            "Use the matching latest restart or a new output directory."
        )
    if last_recorded_time is not None and abs(last_recorded_time - start_time) > 1.0e-8:
        raise ValueError(
            f"History ends at {last_recorded_time}, but restart begins at {start_time}. "
            "Use matching history and restart files."
        )
    record_origin_time = start_time if recorded_origin_time is None else recorded_origin_time

    interval_start = start_time
    storage_before = capture_storage_state(build)
    power = external_power_snapshot(build)
    rows: List[Dict[str, Any]] = []
    print("=== V8 CaseA long radial energy audit ===", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)
    print(f"latest_restart={restart_path}", flush=True)
    print(f"coolant_material={build['coolant_material']}", flush=True)
    print(f"solid_ode_method={build['solid_ode_method']}", flush=True)
    print(f"wire_resistance_ohm={build['wire_resistance_ohm']}", flush=True)

    while interval_start < target_time - 1.0e-10:
        interval_end = min(interval_start + float(args.record_interval), target_time)
        advanced = advance_interval(build, interval_end, args, power)
        power = advanced["power"]
        storage_after = capture_storage_state(build)
        storage = finite_difference_storage_rates(
            build,
            storage_before,
            storage_after,
            interval_end - interval_start,
        )
        row = build_record(
            build,
            record_origin_time,
            interval_start,
            interval_end,
            advanced["integral"],
            storage,
            advanced["step_count"],
        )
        passive = passive_tec_source_totals(build)
        for key, value in passive.items():
            row[f"ring3_open_{key}"] = value
        if any(value != 0.0 for value in passive.values()):
            raise RuntimeError(f"Ring3_Open TEC sources are not zero: {passive}")

        fields = list(row.keys())
        if existing_fields is None:
            existing_fields = fields
        elif fields != existing_fields:
            raise ValueError("Existing V8 energy history columns do not match the current audit schema.")
        append_row(history_path, existing_fields, row, write_header=not history_path.exists())
        system.save_global_state(str(restart_path))
        rows.append(row)

        latest = {
            "restart_in": args.restart_in,
            "restart_out": str(restart_path),
            "history_csv": str(history_path),
            "start_time_s": start_time,
            "record_origin_time_s": record_origin_time,
            "end_time_s": float(system.global_time),
            "target_time_s": target_time,
            "record_interval_s": float(args.record_interval),
            "max_dt_s": float(args.max_dt),
            "inner_iter": int(args.inner_iter),
            "coolant_material": build["coolant_material"],
            "solid_ode_method": build["solid_ode_method"],
            "solid_ode_methods": get_solid_ode_methods(build),
            "wire_resistance_ohm": build["wire_resistance_ohm"],
            "ring_multipliers": build["ring_multipliers"],
            "tec_ring_multipliers": build["tec_ring_multipliers"],
            "latest": row,
        }
        with summary_path.open("w", encoding="utf-8") as f:
            json.dump(latest, f, indent=2, ensure_ascii=False, default=json_default)

        print(
            f"t_rel={row['relative_time_s']:.1f}s "
            f"Rmodel={row['global_residual_using_thermal_model_tec_heat_w']:.3f}W "
            f"Rterminal={row['global_residual_using_terminal_power_w']:.3f}W "
            f"Rfsc={row['fluid_solid_mapping_residual_w']:.3f}W "
            f"Qtec_delta={row['tec_thermal_model_minus_electric_count_w']:.6f}W",
            flush=True,
        )
        interval_start = interval_end
        storage_before = storage_after

    print("V8 long run completed.", flush=True)


if __name__ == "__main__":
    main()
