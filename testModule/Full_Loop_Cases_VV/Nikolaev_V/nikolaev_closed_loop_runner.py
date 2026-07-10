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

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nikolaev_source_data import SOURCE_CITATION, SOURCE_DOI, SOURCE_TITLE, TABLE2_OPERATING_POINTS
from nikolaev_thermocalc_model import NikolaevThermalNetworkConfig
from nikolaev_thermoelectric_closed_loop import (
    ClosedLoopConfig,
    result_to_dict,
    solve_closed_loop_point,
    summarize_closed_loop_results,
)


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
        "# Nikolaev 1995 Closed Thermoelectric Single-TFE Validation",
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
        "This run closes a local single-TFE thermal network around ThermoCalc. The literature geometry and heat input are fixed. Each outer iteration sends emitter/collector temperatures to ThermoCalc, reads electronic cooling/heating and authoritative `joulePowerE/C`, then relaxes the emitter-to-collector and collector-to-boundary thermal network. It does not set current from `I=P/V`.",
        "",
        "## Thermal Configuration",
        "",
        "| Parameter | Value |",
        "| --- | ---:|",
    ]
    for key, value in summary["thermal_config"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend(["", "## Closed-Loop Configuration", "", "| Parameter | Value |", "| --- | ---:|"])
    for key, value in summary["closed_loop_config"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend([
        "",
        "## Table 2 Comparison",
        "",
        "| V | I exp A | I calc A | I err A | P calc W | eta calc % | Te exp K | Te mean K | Tc mean K | e-cool W | c-heat W | Joule E/C W | outer iter | closed? |",
        "| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| --- |",
    ])
    for row in rows:
        lines.append(
            "| {voltage_v:.3g} | {current_exp_a:.3f} | {current_calc_a:.3f} | {current_error_a:.3f} | {electric_power_calc_w:.3f} | {efficiency_calc_percent:.3f} | {emitter_temperature_exp_k:.1f} | {emitter_temperature_mean_k:.3f} | {collector_temperature_mean_k:.3f} | {electron_cooling_power_w:.3f} | {collector_electron_heating_power_w:.3f} | {joule_power_emitter_w:.3f}/{joule_power_collector_w:.3f} | {outer_iterations} | {closed_loop_converged} |".format(**row)
        )
    lines.extend(["", "## Metrics", "", "| Metric | Value |", "| --- | ---:|"])
    for key, value in summary["metrics"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend([
        "",
        "## Interpretation",
        "",
        "This is a thermal-feedback validation step, not a final fitted V&V result. A physically acceptable parameter set must improve electrical agreement without pushing the emitter-temperature scale away from Table 2. Runs that fit current by overheating the emitter should be rejected.",
    ])
    return "\n".join(lines) + "\n"


def run_table2_closed_loop(thermal_config: NikolaevThermalNetworkConfig, closed_loop_config: ClosedLoopConfig):
    return [solve_closed_loop_point(point, thermal_config, closed_loop_config) for point in TABLE2_OPERATING_POINTS]


def write_outputs(run_dir: Path, thermal_config: NikolaevThermalNetworkConfig, closed_loop_config: ClosedLoopConfig) -> dict:
    results = run_table2_closed_loop(thermal_config, closed_loop_config)
    result_dicts = [result_to_dict(result) for result in results]
    rows = [_flat_result_row(item) for item in result_dicts]
    summary = {
        "source": {"title": SOURCE_TITLE, "citation": SOURCE_CITATION, "doi": SOURCE_DOI},
        "thermal_config": asdict(thermal_config),
        "closed_loop_config": asdict(closed_loop_config),
        "metrics": summarize_closed_loop_results(results),
        "status": "thermoelectric_closed_loop_table2_validation",
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(run_dir / "results" / "table2_closed_loop_comparison.csv", rows)
    history = _history_rows(result_dicts)
    if history:
        _write_csv(run_dir / "results" / "closed_loop_iteration_history.csv", history)
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
        "nikolaev_closed_loop_runner.py",
        "PARAMETER_ADJUSTMENT_GUIDE.md",
    ):
        src = CASE_DIR / name
        if src.exists():
            shutil.copy2(src, snapshot_dir / name)
    return summary


def _default_run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S_nikolaev_closed_loop")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Nikolaev 1995 thermoelectric closed-loop validation.")
    parser.add_argument("--run-id", default=_default_run_id())
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--n-nodes", type=int, default=NikolaevThermalNetworkConfig.n_nodes)
    parser.add_argument("--emitter-to-collector-resistance-k-per-w", type=float, default=NikolaevThermalNetworkConfig.emitter_to_collector_resistance_k_per_w)
    parser.add_argument("--axial-shape-amplitude", type=float, default=NikolaevThermalNetworkConfig.axial_shape_amplitude)
    parser.add_argument("--cesium-reservoir-temperature-k", type=float, default=NikolaevThermalNetworkConfig.cesium_reservoir_temperature_k)
    parser.add_argument("--wire-resistance-ohm", type=float, default=NikolaevThermalNetworkConfig.wire_resistance_ohm)
    parser.add_argument("--collector-to-boundary-resistance-k-per-w", type=float, default=ClosedLoopConfig.collector_to_boundary_resistance_k_per_w)
    parser.add_argument("--max-iterations", type=int, default=ClosedLoopConfig.max_iterations)
    parser.add_argument("--relaxation", type=float, default=ClosedLoopConfig.relaxation)
    parser.add_argument("--temperature-tolerance-k", type=float, default=ClosedLoopConfig.temperature_tolerance_k)
    parser.add_argument("--i-guess-a", type=float, default=ClosedLoopConfig.i_guess_a)
    args = parser.parse_args(argv)

    thermal_config = NikolaevThermalNetworkConfig(
        n_nodes=args.n_nodes,
        emitter_to_collector_resistance_k_per_w=args.emitter_to_collector_resistance_k_per_w,
        axial_shape_amplitude=args.axial_shape_amplitude,
        cesium_reservoir_temperature_k=args.cesium_reservoir_temperature_k,
        wire_resistance_ohm=args.wire_resistance_ohm,
    )
    closed_loop_config = ClosedLoopConfig(
        max_iterations=args.max_iterations,
        relaxation=args.relaxation,
        temperature_tolerance_k=args.temperature_tolerance_k,
        collector_to_boundary_resistance_k_per_w=args.collector_to_boundary_resistance_k_per_w,
        i_guess_a=args.i_guess_a,
    )
    run_dir = args.output_root / args.run_id
    summary = write_outputs(run_dir, thermal_config, closed_loop_config)
    print(json.dumps({"run_dir": str(run_dir), "metrics": summary["metrics"], "status": summary["status"]}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
