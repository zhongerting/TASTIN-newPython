from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network
from benke_validation import evaluate_benke_validation

CASE_DIR = Path(__file__).resolve().parent
DEFAULT_REGULATED_HE_K_VALUES = (0.073, 0.08, 0.087)
DEFAULT_WATER_H_VALUES = (528.0, 800.0, 1012.0)


def _as_float_list(values: Iterable[float]) -> list[float]:
    return [float(value) for value in values]


def _row_from_result(result, validation: dict) -> dict:
    config = result.config
    return {
        "regulated_he_effective_k_w_m_k": float(config.regulated_he_effective_k_w_m_k),
        "water_h_w_m2_k": float(config.water_h_w_m2_k),
        "active_zone_power_w": float(result.active_zone_power_w),
        "water_outlet_k": float(result.water_bulk_outlet_k),
        "water_delta_t_k": float(result.water_bulk_outlet_k - config.water_inlet_temperature_k),
        "energy_balance_error_w": float(result.energy_balance_error_w),
        "sleeve_outer_mean_k": float(result.sleeve_outer_temperature_k.mean()),
        "sleeve_outer_max_k": float(result.sleeve_outer_temperature_k.max()),
        "collector_inner_mean_k": float(result.collector_inner_temperature_k.mean()),
        "collector_inner_max_k": float(result.collector_inner_temperature_k.max()),
        "validation_status": validation["status"],
        "range_check_status": validation["range_check_status"],
    }


def _summary_from_rows(rows: list[dict]) -> dict:
    if not rows:
        raise ValueError("At least one scan row is required.")
    fields = (
        "water_outlet_k",
        "water_delta_t_k",
        "sleeve_outer_mean_k",
        "sleeve_outer_max_k",
        "collector_inner_mean_k",
        "collector_inner_max_k",
    )
    summary = {
        "grid_point_count": len(rows),
        "all_range_checks_passed": all(row["range_check_status"] == "passed" for row in rows),
        "all_energy_balance_abs_error_w_max": max(abs(row["energy_balance_error_w"]) for row in rows),
        "regulated_he_effective_k_values": sorted({row["regulated_he_effective_k_w_m_k"] for row in rows}),
        "water_h_values": sorted({row["water_h_w_m2_k"] for row in rows}),
    }
    for field in fields:
        values = [row[field] for row in rows]
        summary[f"{field}_min"] = min(values)
        summary[f"{field}_max"] = max(values)
    return summary


def scan_benke_parameter_envelope(
    regulated_he_k_values: Iterable[float] = DEFAULT_REGULATED_HE_K_VALUES,
    water_h_values: Iterable[float] = DEFAULT_WATER_H_VALUES,
    base_config: BenkeThermalNetworkConfig | None = None,
    experimental_data_dir: Path | None = None,
) -> tuple[list[dict], dict]:
    base_config = BenkeThermalNetworkConfig() if base_config is None else base_config
    experimental_data_dir = CASE_DIR / "experimental_data" if experimental_data_dir is None else Path(experimental_data_dir)
    rows: list[dict] = []
    for he_k in _as_float_list(regulated_he_k_values):
        for water_h in _as_float_list(water_h_values):
            config = BenkeThermalNetworkConfig(
                active_length_m=base_config.active_length_m,
                tisa_heated_length_m=base_config.tisa_heated_length_m,
                n_nodes=base_config.n_nodes,
                water_inlet_temperature_k=base_config.water_inlet_temperature_k,
                water_mass_flow_kg_s=base_config.water_mass_flow_kg_s,
                water_cp_j_kg_k=base_config.water_cp_j_kg_k,
                water_h_w_m2_k=water_h,
                regulated_he_effective_k_w_m_k=he_k,
                regulated_he_gap_m=base_config.regulated_he_gap_m,
                unregulated_he_effective_k_w_m_k=base_config.unregulated_he_effective_k_w_m_k,
                unregulated_he_gap_m=base_config.unregulated_he_gap_m,
                cs_gap_effective_k_w_m_k=base_config.cs_gap_effective_k_w_m_k,
                collector_k_w_m_k=base_config.collector_k_w_m_k,
                alumina_k_w_m_k=base_config.alumina_k_w_m_k,
                sleeve_k_w_m_k=base_config.sleeve_k_w_m_k,
                collector_inner_radius_m=base_config.collector_inner_radius_m,
                collector_wall_thickness_m=base_config.collector_wall_thickness_m,
                emitter_outer_radius_m=base_config.emitter_outer_radius_m,
                alumina_thickness_m=base_config.alumina_thickness_m,
                sleeve_outer_radius_m=base_config.sleeve_outer_radius_m,
                extra_resistance_k_per_w=base_config.extra_resistance_k_per_w,
            )
            result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, config)
            validation = evaluate_benke_validation(result, experimental_data_dir)
            rows.append(_row_from_result(result, validation))
    return rows, _summary_from_rows(rows)


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def create_run_dir(run_id: str | None = None, output_root: Path | None = None) -> Path:
    output_root = CASE_DIR / "runs" if output_root is None else Path(output_root)
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_benke_parameter_envelope"
    run_dir = output_root / run_id
    if run_dir.exists():
        raise FileExistsError(f"Run directory already exists: {run_dir}")
    (run_dir / "results").mkdir(parents=True)
    return run_dir


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Scan Benke thermal model over literature-bound parameter ranges.")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--regulated-he-k-values", nargs="*", type=float, default=list(DEFAULT_REGULATED_HE_K_VALUES))
    parser.add_argument("--water-h-values", nargs="*", type=float, default=list(DEFAULT_WATER_H_VALUES))
    args = parser.parse_args(argv)

    rows, summary = scan_benke_parameter_envelope(args.regulated_he_k_values, args.water_h_values)
    run_dir = create_run_dir(args.run_id, args.output_root)
    _write_csv(run_dir / "results" / "benke_parameter_envelope.csv", rows)
    with (run_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"run_dir={run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
