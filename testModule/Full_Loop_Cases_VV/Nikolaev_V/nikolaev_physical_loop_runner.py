from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nikolaev_source_data import SOURCE_CITATION, SOURCE_DOI, SOURCE_TITLE, TABLE2_OPERATING_POINTS
from nikolaev_physical_tfe_loop import PhysicalLoopConfig, result_to_dict, solve_physical_tfe_point, summarize_physical_results


def _csv_value(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.10g}"
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def _write_csv(path: Path, rows: Iterable[dict]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def _flat_result_row(result_dict: dict) -> dict:
    return {key: value for key, value in result_dict.items() if key != "iteration_history"}


def _history_rows(result_dicts: Sequence[dict]) -> list[dict]:
    rows = []
    for result in result_dicts:
        for item in result["iteration_history"]:
            rows.append({"voltage_v": result["voltage_v"], **item})
    return rows


def build_markdown_report(summary: dict, rows: list[dict], run_dir: Path) -> str:
    lines = [
        "# Nikolaev 1995 Physical Thermal-Hydraulic TFE Validation",
        "",
        f"Run directory: `{run_dir}`",
        "",
        "## Source",
        "",
        f"- Title: {SOURCE_TITLE}",
        f"- Citation: {SOURCE_CITATION}",
        f"- DOI: {SOURCE_DOI}",
        "",
        "## Model Boundary",
        "",
        "This path models an explicit heat-source to coolant flow process for one TFE. Heater power is distributed over the published heated length, ThermoCalc supplies electrical output plus electron cooling/heating and authoritative `joulePowerE/C`, and coolant outlet temperature is computed from `m_dot * cp * dT`. It does not prescribe collector temperature and does not set current from `I=P/V`.",
        "",
        "## Configuration",
        "",
        "| Parameter | Value |",
        "| --- | ---:|",
    ]
    for key, value in summary["config"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend([
        "",
        "## Table 2 Coupled Comparison",
        "",
        "| V | I exp A | I calc A | I err A | P calc W | eta calc % | Te exp K | Te mean K | Tc mean K | Tin K | Tout K | Qcool W | residual W | outer iter | closed? |",
        "| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| --- |",
    ])
    for row in rows:
        lines.append(
            "| {voltage_v:.3g} | {current_exp_a:.3f} | {current_calc_a:.3f} | {current_error_a:.3f} | {electric_power_calc_w:.3f} | {efficiency_calc_percent:.3f} | {emitter_temperature_exp_k:.1f} | {emitter_temperature_mean_k:.3f} | {collector_temperature_mean_k:.3f} | {coolant_inlet_temperature_k:.3f} | {coolant_outlet_temperature_k:.3f} | {coolant_heat_gain_w:.3f} | {thermal_balance_residual_w:.3e} | {outer_iterations} | {physical_loop_converged} |".format(**row)
        )
    lines.extend(["", "## Metrics", "", "| Metric | Value |", "| --- | ---:|"])
    for key, value in summary["metrics"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend([
        "",
        "## Interpretation",
        "",
        "This is a more physical reduced single-TFE model than the previous fixed-boundary thermal-resistance closure because it computes coolant heating explicitly. It is still not a full 2D TFEUnit/SystemManager solid-conduction run, so fitted coolant flow and heat-transfer coefficients remain validation parameters.",
    ])
    return "\n".join(lines) + "\n"


def run_table2(config: PhysicalLoopConfig):
    return [solve_physical_tfe_point(point, config) for point in TABLE2_OPERATING_POINTS]


def write_outputs(run_dir: Path, config: PhysicalLoopConfig) -> dict:
    results = run_table2(config)
    result_dicts = [result_to_dict(result) for result in results]
    rows = [_flat_result_row(item) for item in result_dicts]
    summary = {
        "source": {"title": SOURCE_TITLE, "citation": SOURCE_CITATION, "doi": SOURCE_DOI},
        "config": asdict(config),
        "metrics": summarize_physical_results(results),
        "status": "physical_thermal_hydraulic_tfe_validation",
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(run_dir / "results" / "table2_physical_tfe_comparison.csv", rows)
    history = _history_rows(result_dicts)
    if history:
        _write_csv(run_dir / "results" / "physical_loop_iteration_history.csv", history)
    with (run_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump({**summary, "table2": result_dicts}, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    (run_dir / "validation_report.md").write_text(build_markdown_report(summary, rows, run_dir), encoding="utf-8")
    snapshot_dir = run_dir / "input_snapshot"
    snapshot_dir.mkdir(exist_ok=True)
    for name in (
        "nikolaev_source_data.py",
        "nikolaev_thermocalc_model.py",
        "nikolaev_thermoelectric_closed_loop.py",
        "nikolaev_physical_tfe_loop.py",
        "nikolaev_physical_loop_runner.py",
        "PARAMETER_ADJUSTMENT_GUIDE.md",
    ):
        src = CASE_DIR / name
        if src.exists():
            shutil.copy2(src, snapshot_dir / name)
    return summary


def _default_run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S_nikolaev_physical_tfe")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Nikolaev 1995 physical thermal-hydraulic TFE validation.")
    parser.add_argument("--run-id", default=_default_run_id())
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--n-nodes", type=int, default=PhysicalLoopConfig.n_nodes)
    parser.add_argument("--coolant-inlet-temperature-k", type=float, default=PhysicalLoopConfig.coolant_inlet_temperature_k)
    parser.add_argument("--coolant-mass-flow-kg-s", type=float, default=PhysicalLoopConfig.coolant_mass_flow_kg_s)
    parser.add_argument("--coolant-heat-capacity-j-kg-k", type=float, default=PhysicalLoopConfig.coolant_heat_capacity_j_kg_k)
    parser.add_argument("--collector-convective-h-w-m2-k", type=float, default=PhysicalLoopConfig.collector_convective_h_w_m2_k)
    parser.add_argument("--emitter-to-collector-resistance-k-per-w", type=float, default=PhysicalLoopConfig.emitter_to_collector_resistance_k_per_w)
    parser.add_argument("--cesium-reservoir-temperature-k", type=float, default=PhysicalLoopConfig.cesium_reservoir_temperature_k)
    parser.add_argument("--wire-resistance-ohm", type=float, default=PhysicalLoopConfig.wire_resistance_ohm)
    parser.add_argument("--axial-shape-amplitude", type=float, default=PhysicalLoopConfig.axial_shape_amplitude)
    parser.add_argument("--enable-axial-conduction", dest="axial_conduction_enabled", action="store_true", default=PhysicalLoopConfig.axial_conduction_enabled)
    parser.add_argument("--disable-axial-conduction", dest="axial_conduction_enabled", action="store_false")
    parser.add_argument("--axial-conduction-smoothing", type=float, default=PhysicalLoopConfig.axial_conduction_smoothing)
    parser.add_argument("--axial-conduction-passes", type=int, default=PhysicalLoopConfig.axial_conduction_passes)
    parser.add_argument("--max-iterations", type=int, default=PhysicalLoopConfig.max_iterations)
    parser.add_argument("--relaxation", type=float, default=PhysicalLoopConfig.relaxation)
    parser.add_argument("--temperature-tolerance-k", type=float, default=PhysicalLoopConfig.temperature_tolerance_k)
    parser.add_argument("--i-guess-a", type=float, default=PhysicalLoopConfig.i_guess_a)
    args = parser.parse_args(argv)

    config = PhysicalLoopConfig(
        n_nodes=args.n_nodes,
        coolant_inlet_temperature_k=args.coolant_inlet_temperature_k,
        coolant_mass_flow_kg_s=args.coolant_mass_flow_kg_s,
        coolant_heat_capacity_j_kg_k=args.coolant_heat_capacity_j_kg_k,
        collector_convective_h_w_m2_k=args.collector_convective_h_w_m2_k,
        emitter_to_collector_resistance_k_per_w=args.emitter_to_collector_resistance_k_per_w,
        cesium_reservoir_temperature_k=args.cesium_reservoir_temperature_k,
        wire_resistance_ohm=args.wire_resistance_ohm,
        axial_shape_amplitude=args.axial_shape_amplitude,
        axial_conduction_enabled=args.axial_conduction_enabled,
        axial_conduction_smoothing=args.axial_conduction_smoothing,
        axial_conduction_passes=args.axial_conduction_passes,
        max_iterations=args.max_iterations,
        relaxation=args.relaxation,
        temperature_tolerance_k=args.temperature_tolerance_k,
        i_guess_a=args.i_guess_a,
    )
    run_dir = args.output_root / args.run_id
    summary = write_outputs(run_dir, config)
    print(json.dumps({"run_dir": str(run_dir), "metrics": summary["metrics"], "status": summary["status"]}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
