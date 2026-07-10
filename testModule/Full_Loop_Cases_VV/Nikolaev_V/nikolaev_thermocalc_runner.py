from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nikolaev_source_data import SOURCE_CITATION, SOURCE_DOI, SOURCE_TITLE, TABLE2_OPERATING_POINTS, Table2OperatingPoint
from nikolaev_thermocalc_model import NikolaevThermalNetworkConfig, build_thermocalc_case


@dataclass(frozen=True)
class ThermoCalcVoltageConfig:
    i_guess_a: float = 300.0


@dataclass(frozen=True)
class NikolaevThermoCalcResult:
    voltage_v: float
    thermal_power_kw: float
    emitter_temperature_exp_k: float
    emitter_temperature_mean_k: float
    emitter_temperature_min_k: float
    emitter_temperature_max_k: float
    collector_temperature_mean_k: float
    cesium_reservoir_temperature_k: float
    current_exp_a: float
    current_calc_a: float
    current_error_a: float
    electric_power_exp_w: float
    electric_power_calc_w: float
    electric_power_error_w: float
    efficiency_exp_percent: float
    efficiency_calc_percent: float
    efficiency_error_percent: float
    uout_v: float
    target_error_v: float
    converged: bool
    finite: bool
    iteration_count: int
    zero_emission_skipped: bool


def _default_model_factory(n_elements: int, n_nodes: int):
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    return ThermoCalcModel(n_elements=n_elements, n_nodes=n_nodes)


def apply_case_to_thermocalc(thermo_model, case_model, target_voltage_v: float, i_guess_a: float) -> None:
    arrays = case_model.arrays
    geometry = case_model.geometry
    thermo_model._input_data.dlE = np.asarray(arrays.dl_emitter_m, dtype=float)
    thermo_model._input_data.dlC = np.asarray(arrays.dl_collector_m, dtype=float)
    thermo_model._input_data.sideAreaE = np.asarray(arrays.side_area_emitter_m2, dtype=float)
    thermo_model._input_data.sideAreaC = np.asarray(arrays.side_area_collector_m2, dtype=float)
    thermo_model._input_data.crossAreaE = np.array([geometry.emitter_cross_area_m2], dtype=float)
    thermo_model._input_data.crossAreaC = np.array([geometry.collector_cross_area_m2], dtype=float)
    thermo_model._input_data.d_gap = np.array([geometry.interelectrode_gap_m * 1000.0], dtype=float)
    thermo_model._input_data.resistanceWire = np.full((1, 4), case_model.thermal_config.wire_resistance_ohm, dtype=float)
    thermo_model._input_data.wireU = np.array([[target_voltage_v, target_voltage_v, 0.0, 0.0]], dtype=float)
    thermo_model.setup_circuit_mode("fixed_u", float(target_voltage_v), I_guess=float(i_guess_a))
    thermo_model.set_temperatures(np.asarray(arrays.temitter_k, dtype=float), np.asarray(arrays.tcollector_k, dtype=float))
    thermo_model.set_tcs(np.asarray(arrays.tcs_k, dtype=float))


def scan_table2_point(
    point: Table2OperatingPoint,
    thermal_config: NikolaevThermalNetworkConfig = NikolaevThermalNetworkConfig(),
    voltage_config: ThermoCalcVoltageConfig = ThermoCalcVoltageConfig(),
    model_factory: Callable[[int, int], object] = _default_model_factory,
) -> NikolaevThermoCalcResult:
    case_model = build_thermocalc_case(point, thermal_config)
    thermo_model = model_factory(1, case_model.geometry.n_nodes)
    apply_case_to_thermocalc(thermo_model, case_model, point.voltage_v, voltage_config.i_guess_a)
    thermo_model.calculate(verbose=False)
    global_results = thermo_model.get_global_results() or {}
    current = float(global_results.get("Iout", math.nan))
    uout = float(global_results.get("Uout", point.voltage_v))
    power = current * uout
    finite = bool(np.isfinite([current, uout, power]).all())
    eta = 100.0 * power / (point.thermal_power_kw * 1000.0) if finite else math.nan
    return NikolaevThermoCalcResult(
        voltage_v=point.voltage_v,
        thermal_power_kw=point.thermal_power_kw,
        emitter_temperature_exp_k=point.emitter_temperature_k,
        emitter_temperature_mean_k=float(np.mean(case_model.arrays.temitter_k)),
        emitter_temperature_min_k=float(np.min(case_model.arrays.temitter_k)),
        emitter_temperature_max_k=float(np.max(case_model.arrays.temitter_k)),
        collector_temperature_mean_k=float(np.mean(case_model.arrays.tcollector_k)),
        cesium_reservoir_temperature_k=float(case_model.arrays.tcs_k[0, 0]),
        current_exp_a=point.current_a,
        current_calc_a=current,
        current_error_a=current - point.current_a if finite else math.nan,
        electric_power_exp_w=point.electric_power_w,
        electric_power_calc_w=power if finite else math.nan,
        electric_power_error_w=power - point.electric_power_w if finite else math.nan,
        efficiency_exp_percent=point.efficiency_percent,
        efficiency_calc_percent=eta,
        efficiency_error_percent=eta - point.efficiency_percent if finite else math.nan,
        uout_v=uout,
        target_error_v=abs(uout - point.voltage_v) if finite else math.nan,
        converged=bool(global_results.get("converged", False)),
        finite=finite,
        iteration_count=int(global_results.get("iteration_count", -1)),
        zero_emission_skipped=bool(global_results.get("zero_emission_skipped", False)),
    )


def run_table2_points(
    thermal_config: NikolaevThermalNetworkConfig = NikolaevThermalNetworkConfig(),
    voltage_config: ThermoCalcVoltageConfig = ThermoCalcVoltageConfig(),
    model_factory: Callable[[int, int], object] = _default_model_factory,
) -> list[NikolaevThermoCalcResult]:
    return [scan_table2_point(point, thermal_config, voltage_config, model_factory) for point in TABLE2_OPERATING_POINTS]


def summarize_results(results: Sequence[NikolaevThermoCalcResult]) -> dict:
    current_errors = np.asarray([r.current_error_a for r in results], dtype=float)
    power_errors = np.asarray([r.electric_power_error_w for r in results], dtype=float)
    eta_errors = np.asarray([r.efficiency_error_percent for r in results], dtype=float)
    return {
        "case_count": len(results),
        "finite_all": all(r.finite for r in results),
        "converged_all": all(r.converged for r in results),
        "current_mae_a": float(np.nanmean(np.abs(current_errors))),
        "current_max_abs_a": float(np.nanmax(np.abs(current_errors))),
        "electric_power_mae_w": float(np.nanmean(np.abs(power_errors))),
        "efficiency_mae_percent": float(np.nanmean(np.abs(eta_errors))),
    }


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


def build_markdown_report(summary: dict, rows: list[dict], run_dir: Path) -> str:
    lines = [
        "# Nikolaev 1995 ThermoCalc Single-TFE Electrical Validation",
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
        "This run fixes the Nikolaev geometry and heat input, builds emitter/collector temperature fields with a case-local thermal network, and obtains current from ThermoCalc fixed-voltage solves. It does not use `I=P/V` or force the Table 2 current.",
        "",
        "## Configuration",
        "",
        "| Parameter | Value |",
        "| --- | ---:|",
    ]
    for key, value in summary["thermal_config"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend([
        "",
        "## Table 2 Electrical Comparison",
        "",
        "| V | I exp A | I ThermoCalc A | I err A | P exp W | P calc W | eta exp % | eta calc % | Te mean K | converged |",
        "| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| --- |",
    ])
    for row in rows:
        lines.append(
            "| {voltage_v:.3g} | {current_exp_a:.3f} | {current_calc_a:.3f} | {current_error_a:.3f} | {electric_power_exp_w:.3f} | {electric_power_calc_w:.3f} | {efficiency_exp_percent:.3f} | {efficiency_calc_percent:.3f} | {emitter_temperature_mean_k:.3f} | {converged} |".format(**row)
        )
    lines.extend(["", "## Metrics", "", "| Metric | Value |", "| --- | ---:|"])
    for key, value in summary["metrics"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend([
        "",
        "## Interpretation",
        "",
        "This is the first non-table-reconstruction Nikolaev_V path. Current is produced by ThermoCalc from the fixed voltage, geometry, cesium temperature, and thermal-network temperature fields. Remaining discrepancies should be addressed by global adjustments to thermal resistance, cesium reservoir temperature, wire resistance, and then by digitizing Figure 4.",
    ])
    return "\n".join(lines) + "\n"


def write_outputs(run_dir: Path, thermal_config: NikolaevThermalNetworkConfig, voltage_config: ThermoCalcVoltageConfig) -> dict:
    results = run_table2_points(thermal_config, voltage_config)
    rows = [asdict(result) for result in results]
    summary = {
        "source": {"title": SOURCE_TITLE, "citation": SOURCE_CITATION, "doi": SOURCE_DOI},
        "thermal_config": asdict(thermal_config),
        "voltage_config": asdict(voltage_config),
        "metrics": summarize_results(results),
        "status": "thermocalc_table2_electrical_validation",
    }
    run_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(run_dir / "results" / "table2_thermocalc_electrical_comparison.csv", rows)
    with (run_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump({**summary, "table2": rows}, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    (run_dir / "validation_report.md").write_text(build_markdown_report(summary, rows, run_dir), encoding="utf-8")
    snapshot_dir = run_dir / "input_snapshot"
    snapshot_dir.mkdir(exist_ok=True)
    for name in ("nikolaev_source_data.py", "nikolaev_thermocalc_model.py", "nikolaev_thermocalc_runner.py", "PARAMETER_ADJUSTMENT_GUIDE.md"):
        src = CASE_DIR / name
        if src.exists():
            shutil.copy2(src, snapshot_dir / name)
    return summary


def _default_run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S_nikolaev_thermocalc")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Nikolaev 1995 ThermoCalc single-TFE electrical validation.")
    parser.add_argument("--run-id", default=_default_run_id())
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--n-nodes", type=int, default=NikolaevThermalNetworkConfig.n_nodes)
    parser.add_argument("--emitter-to-collector-resistance-k-per-w", type=float, default=NikolaevThermalNetworkConfig.emitter_to_collector_resistance_k_per_w)
    parser.add_argument("--axial-shape-amplitude", type=float, default=NikolaevThermalNetworkConfig.axial_shape_amplitude)
    parser.add_argument("--cesium-reservoir-temperature-k", type=float, default=NikolaevThermalNetworkConfig.cesium_reservoir_temperature_k)
    parser.add_argument("--wire-resistance-ohm", type=float, default=NikolaevThermalNetworkConfig.wire_resistance_ohm)
    parser.add_argument("--i-guess-a", type=float, default=ThermoCalcVoltageConfig.i_guess_a)
    args = parser.parse_args(argv)
    thermal_config = NikolaevThermalNetworkConfig(
        n_nodes=args.n_nodes,
        emitter_to_collector_resistance_k_per_w=args.emitter_to_collector_resistance_k_per_w,
        axial_shape_amplitude=args.axial_shape_amplitude,
        cesium_reservoir_temperature_k=args.cesium_reservoir_temperature_k,
        wire_resistance_ohm=args.wire_resistance_ohm,
    )
    voltage_config = ThermoCalcVoltageConfig(i_guess_a=args.i_guess_a)
    run_dir = args.output_root / args.run_id
    summary = write_outputs(run_dir, thermal_config, voltage_config)
    print(json.dumps({"run_dir": str(run_dir), "metrics": summary["metrics"], "status": summary["status"]}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
