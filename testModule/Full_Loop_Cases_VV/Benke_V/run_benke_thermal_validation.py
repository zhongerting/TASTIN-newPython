from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Sequence

from benke_thermal_network import (
    BENKE_TYPICAL_CASE,
    BenkeThermalCase,
    BenkeThermalNetworkConfig,
    solve_benke_thermal_network,
)
from benke_report import build_markdown_report
from benke_validation import evaluate_benke_validation

CASE_DIR = Path(__file__).resolve().parent


def _csv_value(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.10g}"
    return value


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def create_run_dir(run_id: str | None = None, output_root: Path | None = None) -> Path:
    output_root = CASE_DIR / "runs" if output_root is None else Path(output_root)
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_benke_thermal_smoke"
    run_dir = output_root / run_id
    if run_dir.exists():
        raise FileExistsError(f"Run directory already exists: {run_dir}")
    for child in ("results", "logs", "input_snapshot"):
        (run_dir / child).mkdir(parents=True, exist_ok=False)
    for name in (
        "BENKE_THERMAL_CLOSURE_DATA.md",
        "BENKE_MODEL_DECISION.md",
        "BENKE_VALIDATION_STATUS.md",
        "BENKE_GOAL_COMPLETION_AUDIT.md",
        "benke_thermal_network.py",
        "run_benke_thermal_validation.py",
        "benke_validation.py",
        "benke_report.py",
    ):
        src = CASE_DIR / name
        if src.exists():
            (run_dir / "input_snapshot" / name).write_bytes(src.read_bytes())
    return run_dir


def run_case(config: BenkeThermalNetworkConfig, case: BenkeThermalCase = BENKE_TYPICAL_CASE) -> dict:
    result = solve_benke_thermal_network(case, config)
    axial_rows = []
    dz_m = config.active_length_m / config.n_nodes
    for idx in range(config.n_nodes):
        axial_rows.append(
            {
                "node": idx,
                "z_center_m": (idx + 0.5) * dz_m - 0.5 * config.active_length_m,
                "heat_source_w": float(result.heat_source_w[idx]),
                "water_bulk_k": float(result.water_bulk_temperature_k[idx]),
                "sleeve_outer_k": float(result.sleeve_outer_temperature_k[idx]),
                "collector_inner_k": float(result.collector_inner_temperature_k[idx]),
                "r_collector_to_water_k_per_w": float(result.r_collector_to_water_k_per_w[idx]),
                "r_sleeve_to_water_k_per_w": float(result.r_sleeve_to_water_k_per_w[idx]),
            }
        )
    tc_rows = []
    for idx, value in enumerate(result.sleeve_thermocouple_temperature_k):
        z_mm = result.sleeve_thermocouple_z_mm[idx]
        tc_rows.append(
            {
                "thermocouple_index": idx + 1,
                "z_mm": "" if z_mm is None else float(z_mm),
                "measurement_radius_m": result.sleeve_thermocouple_radius_m,
                "included_in_benke_average": bool(result.sleeve_thermocouple_included_in_benke_average[idx]),
                "sleeve_outer_k": float(value),
                "sleeve_thermocouple_k": float(value),
            }
        )
    return {
        "result": result,
        "axial_rows": axial_rows,
        "thermocouple_rows": tc_rows,
    }


def write_outputs(
    run_dir: Path,
    config: BenkeThermalNetworkConfig,
    case: BenkeThermalCase = BENKE_TYPICAL_CASE,
    experimental_data_dir: Path | None = None,
) -> dict:
    payload = run_case(config, case)
    result = payload["result"]
    _write_csv(run_dir / "results" / "axial_temperature_profile.csv", payload["axial_rows"])
    _write_csv(run_dir / "results" / "sleeve_thermocouple_12pt.csv", payload["thermocouple_rows"])
    validation = evaluate_benke_validation(result, CASE_DIR / "experimental_data" if experimental_data_dir is None else experimental_data_dir)
    summary = {
        "case": asdict(case),
        "config": asdict(config),
        "active_zone_power_w": result.active_zone_power_w,
        "water_bulk_outlet_k": result.water_bulk_outlet_k,
        "water_delta_t_k": result.water_bulk_outlet_k - config.water_inlet_temperature_k,
        "energy_balance_error_w": result.energy_balance_error_w,
        "collector_inner_mean_k": float(result.collector_inner_temperature_k.mean()),
        "collector_inner_max_k": float(result.collector_inner_temperature_k.max()),
        "sleeve_outer_mean_k": float(result.sleeve_outer_temperature_k.mean()),
        "sleeve_outer_max_k": float(result.sleeve_outer_temperature_k.max()),
        "regulated_he_effective_k_target_range_w_m_k": [0.073, 0.087],
        "water_h_reference_range_w_m2_k": [528.0, 1012.0],
        "validation": validation,
    }
    with (run_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    (run_dir / "validation_report.md").write_text(build_markdown_report(summary, run_dir=run_dir), encoding="utf-8")
    return summary


def append_process_log(run_dir: Path, summary: dict) -> None:
    log_path = CASE_DIR / "validation_process_log.md"
    now = datetime.now().isoformat(timespec="seconds")
    text = (
        f"\n## {now} - {run_dir.name}\n\n"
        f"Run directory: `{run_dir}`\n\n"
        f"Active-zone power: {summary['active_zone_power_w']:.6g} W\n\n"
        f"Water outlet temperature: {summary['water_bulk_outlet_k']:.6g} K\n\n"
        f"Energy balance error: {summary['energy_balance_error_w']:.6g} W\n\n"
        f"Sleeve outer mean/max: {summary['sleeve_outer_mean_k']:.6g} / {summary['sleeve_outer_max_k']:.6g} K\n\n"
        f"Validation status: {summary['validation']['status']} / range checks {summary['validation']['range_check_status']}\n\n"
        "Decision: Benke thermal-network validation baseline; compare with digitized Benke sleeve/water data when available.\n"
    )
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Benke single-cell TFE thermal-network smoke validation.")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--tisa-power-w", type=float, default=BENKE_TYPICAL_CASE.tisa_power_w)
    parser.add_argument("--regulated-he-pressure-torr", type=float, default=BENKE_TYPICAL_CASE.regulated_he_pressure_torr)
    parser.add_argument("--active-length-m", type=float, default=BenkeThermalNetworkConfig.active_length_m)
    parser.add_argument("--tisa-heated-length-m", type=float, default=BenkeThermalNetworkConfig.tisa_heated_length_m)
    parser.add_argument("--water-inlet-temperature-k", type=float, default=BenkeThermalNetworkConfig.water_inlet_temperature_k)
    parser.add_argument("--water-mass-flow-kg-s", type=float, default=BenkeThermalNetworkConfig.water_mass_flow_kg_s)
    parser.add_argument("--water-h-w-m2-k", type=float, default=BenkeThermalNetworkConfig.water_h_w_m2_k)
    parser.add_argument("--regulated-he-effective-k-w-m-k", type=float, default=BenkeThermalNetworkConfig.regulated_he_effective_k_w_m_k)
    parser.add_argument("--coolant-heat-fraction", type=float, default=BenkeThermalNetworkConfig.coolant_heat_fraction)
    parser.add_argument("--experimental-data-dir", type=Path, default=CASE_DIR / "experimental_data")
    args = parser.parse_args(argv)

    case = BenkeThermalCase(
        name="benke_cli_case",
        tisa_power_w=args.tisa_power_w,
        regulated_he_pressure_torr=args.regulated_he_pressure_torr,
    )
    config = BenkeThermalNetworkConfig(
        active_length_m=args.active_length_m,
        tisa_heated_length_m=args.tisa_heated_length_m,
        water_inlet_temperature_k=args.water_inlet_temperature_k,
        water_mass_flow_kg_s=args.water_mass_flow_kg_s,
        water_h_w_m2_k=args.water_h_w_m2_k,
        regulated_he_effective_k_w_m_k=args.regulated_he_effective_k_w_m_k,
        coolant_heat_fraction=args.coolant_heat_fraction,
    )
    run_dir = create_run_dir(args.run_id, args.output_root)
    summary = write_outputs(run_dir, config, case, experimental_data_dir=args.experimental_data_dir)
    append_process_log(run_dir, summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"run_dir={run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
