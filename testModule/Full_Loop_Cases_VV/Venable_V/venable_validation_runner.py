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
from typing import Callable, Iterable, List, Sequence

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from venable_single_tfe_model import (
    DEFAULT_THERMAL_CLOSURE,
    AXIAL_PROFILE_COSINE,
    COLLECTOR_BOUNDARY_LINEAR,
    VALID_AXIAL_PROFILE_MODES,
    VALID_COLLECTOR_BOUNDARY_MODES,
    TCS_MODE_PLACEHOLDER,
    VALID_TCS_MODES,
    VALID_THERMAL_MODEL_MODES,
    VenableThermalClosure,
    build_case_model,
)
from venable_table71_data import TABLE71_CASES, VenableTable71Case


@dataclass(frozen=True)
class VoltageScanConfig:
    start_v: float = 0.05
    stop_v: float = 2.0
    points: int = 40
    i_guess_a: float = 150.0


@dataclass(frozen=True)
class CaseValidationResult:
    case_id: str
    q_az_w: float
    pcs_torr: float
    tcs_k: float
    p_out_exp_w: float
    eta_exp_percent: float
    u_max_calc_v: float
    i_at_max_calc_a: float
    p_out_calc_w: float
    eta_calc_percent: float
    p_out_abs_error_w: float
    p_out_rel_error: float
    eta_abs_error_percent: float
    converged_all: bool
    finite_all: bool
    scan_points: int
    best_iteration_count: int
    zero_emission_skipped_at_best: bool
    scan_records: List[dict]


def _voltage_values(config: VoltageScanConfig) -> np.ndarray:
    if config.points < 2:
        raise ValueError("Voltage scan requires at least two points.")
    if not math.isfinite(config.start_v) or not math.isfinite(config.stop_v):
        raise ValueError("Voltage scan bounds must be finite.")
    if config.stop_v <= config.start_v:
        raise ValueError("Voltage scan stop_v must be greater than start_v.")
    return np.linspace(float(config.start_v), float(config.stop_v), int(config.points))


def _default_model_factory(n_elements: int, n_nodes: int):
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    return ThermoCalcModel(n_elements=n_elements, n_nodes=n_nodes)


def _apply_case_model_to_thermocalc(thermo_model, case_model, target_voltage_v: float, i_guess_a: float) -> None:
    arrays = case_model.arrays
    geometry = case_model.geometry
    thermo_model._input_data.dlE = np.asarray(arrays.dl_emitter_m, dtype=float)
    thermo_model._input_data.dlC = np.asarray(arrays.dl_collector_m, dtype=float)
    thermo_model._input_data.sideAreaE = np.asarray(arrays.side_area_emitter_m2, dtype=float)
    thermo_model._input_data.sideAreaC = np.asarray(arrays.side_area_collector_m2, dtype=float)
    thermo_model._input_data.crossAreaE = np.array([geometry.emitter_cross_area_m2], dtype=float)
    thermo_model._input_data.crossAreaC = np.array([geometry.collector_cross_area_m2], dtype=float)
    thermo_model._input_data.d_gap = np.array([geometry.gap_m * 1000.0], dtype=float)
    thermo_model._input_data.wireU = np.array([[target_voltage_v, target_voltage_v, 0.0, 0.0]], dtype=float)
    thermo_model.setup_circuit_mode("fixed_u", target_voltage_v, I_guess=i_guess_a)
    thermo_model.set_temperatures(
        np.asarray(arrays.temitter_k, dtype=float),
        np.asarray(arrays.tcollector_k, dtype=float),
    )
    thermo_model.set_tcs(np.asarray(arrays.tcs_k, dtype=float))


def scan_case_fixed_voltage(
    case: VenableTable71Case,
    config: VoltageScanConfig,
    model_factory: Callable[[int, int], object] = _default_model_factory,
    thermal_closure: VenableThermalClosure = DEFAULT_THERMAL_CLOSURE,
    tcs_mode: str = TCS_MODE_PLACEHOLDER,
) -> CaseValidationResult:
    case_model = build_case_model(case, thermal_closure=thermal_closure, tcs_mode=tcs_mode)
    scan_records: List[dict] = []

    for voltage in _voltage_values(config):
        thermo_model = model_factory(1, case_model.geometry.n_nodes)
        _apply_case_model_to_thermocalc(
            thermo_model,
            case_model,
            target_voltage_v=float(voltage),
            i_guess_a=float(config.i_guess_a),
        )
        elapsed_ms = float(thermo_model.calculate(verbose=False))
        global_results = thermo_model.get_global_results() or {}
        current = float(global_results.get("Iout", math.nan))
        uout = float(global_results.get("Uout", voltage))
        power = current * uout
        finite = bool(np.isfinite([current, uout, power]).all())
        target_error_v = abs(uout - float(voltage)) if finite else math.nan
        target_tol_v = max(1.0e-6, 1.0e-4 * max(abs(float(voltage)), 1.0))
        target_matched = bool(finite and target_error_v <= target_tol_v)
        converged = bool(global_results.get("converged", False))
        scan_records.append(
            {
                "voltage_v": float(voltage),
                "uout_v": uout,
                "target_error_v": target_error_v,
                "target_matched": target_matched,
                "iout_a": current,
                "p_out_w": power if finite else math.nan,
                "elapsed_ms": elapsed_ms,
                "converged": converged,
                "iteration_count": int(global_results.get("iteration_count", -1)),
                "zero_emission_skipped": bool(global_results.get("zero_emission_skipped", False)),
                "finite": finite,
            }
        )

    finite_records = [row for row in scan_records if row["finite"] and row["target_matched"]]
    if not finite_records:
        best = scan_records[0]
    else:
        best = max(finite_records, key=lambda row: row["p_out_w"])

    p_calc = float(best["p_out_w"])
    eta_calc = 100.0 * p_calc / case.q_az_w if math.isfinite(p_calc) else math.nan
    p_abs = p_calc - case.p_out_exp_w if math.isfinite(p_calc) else math.nan
    p_rel = p_abs / case.p_out_exp_w if math.isfinite(p_abs) else math.nan
    eta_abs = eta_calc - case.eta_exp_percent if math.isfinite(eta_calc) else math.nan

    return CaseValidationResult(
        case_id=case.case_id,
        q_az_w=case.q_az_w,
        pcs_torr=case.pcs_torr,
        tcs_k=float(case_model.arrays.tcs_k[0, 0]),
        p_out_exp_w=case.p_out_exp_w,
        eta_exp_percent=case.eta_exp_percent,
        u_max_calc_v=float(best["uout_v"]),
        i_at_max_calc_a=float(best["iout_a"]),
        p_out_calc_w=p_calc,
        eta_calc_percent=eta_calc,
        p_out_abs_error_w=p_abs,
        p_out_rel_error=p_rel,
        eta_abs_error_percent=eta_abs,
        converged_all=all(bool(row["converged"]) and bool(row["target_matched"]) for row in scan_records),
        finite_all=all(bool(row["finite"]) for row in scan_records),
        scan_points=len(scan_records),
        best_iteration_count=int(best["iteration_count"]),
        zero_emission_skipped_at_best=bool(best["zero_emission_skipped"]),
        scan_records=scan_records,
    )


def run_validation_cases(
    cases: Sequence[VenableTable71Case],
    config: VoltageScanConfig,
    model_factory: Callable[[int, int], object] = _default_model_factory,
    thermal_closure: VenableThermalClosure = DEFAULT_THERMAL_CLOSURE,
    tcs_mode: str = TCS_MODE_PLACEHOLDER,
) -> List[CaseValidationResult]:
    return [
        scan_case_fixed_voltage(
            case,
            config,
            model_factory=model_factory,
            thermal_closure=thermal_closure,
            tcs_mode=tcs_mode,
        )
        for case in cases
    ]


def summarize_validation_results(results: Sequence[CaseValidationResult]) -> dict:
    rel_errors = np.asarray([result.p_out_rel_error for result in results], dtype=float)
    abs_rel = np.abs(rel_errors[np.isfinite(rel_errors)])
    p_abs = np.asarray([result.p_out_abs_error_w for result in results], dtype=float)
    return {
        "case_count": len(results),
        "finite_case_count": int(np.isfinite(rel_errors).sum()),
        "converged_all_cases": all(result.converged_all for result in results),
        "finite_all_scan_points": all(result.finite_all for result in results),
        "mean_abs_p_out_rel_error": float(np.mean(abs_rel)) if abs_rel.size else math.nan,
        "max_abs_p_out_rel_error": float(np.max(abs_rel)) if abs_rel.size else math.nan,
        "mean_p_out_abs_error_w": float(np.nanmean(p_abs)) if p_abs.size else math.nan,
        "max_abs_p_out_error_w": float(np.nanmax(np.abs(p_abs))) if p_abs.size else math.nan,
    }


def _csv_value(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        return f"{value:.10g}"
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def _result_row(result: CaseValidationResult) -> dict:
    row = asdict(result)
    row.pop("scan_records")
    return row


def _write_csv(path: Path, rows: Iterable[dict]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(value) for key, value in row.items()})


def create_run_dir(stage: str, run_id: str | None = None, output_root: Path | None = None) -> Path:
    output_root = CASE_DIR / "runs" if output_root is None else Path(output_root)
    if run_id is None:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + f"_{stage}"
    run_dir = output_root / run_id
    if run_dir.exists():
        raise FileExistsError(f"Run directory already exists: {run_dir}")
    for child in ("input_snapshot", "results", "logs", "plots"):
        (run_dir / child).mkdir(parents=True, exist_ok=False)
    return run_dir


def snapshot_inputs(run_dir: Path) -> None:
    names = [
        "goal.md",
        "PARAMETER_ADJUSTMENT_GUIDE.md",
        "venable_table71_data.py",
        "venable_single_tfe_model.py",
        "venable_thermal_network.py",
        "venable_validation_runner.py",
        "cases_table71.csv",
        "model_config_summary.json",
    ]
    target = run_dir / "input_snapshot"
    for name in names:
        src = CASE_DIR / name
        if src.exists():
            shutil.copy2(src, target / name)


def write_run_outputs(
    run_dir: Path,
    stage: str,
    config: VoltageScanConfig,
    thermal_closure: VenableThermalClosure,
    tcs_mode: str,
    results: Sequence[CaseValidationResult],
) -> dict:
    summary = summarize_validation_results(results)
    results_dir = run_dir / "results"
    _write_csv(results_dir / "validation_results.csv", [_result_row(result) for result in results])
    scan_rows = []
    for result in results:
        for row in result.scan_records:
            scan_rows.append({"case_id": result.case_id, **row})
    _write_csv(results_dir / "voltage_scan_records.csv", scan_rows)
    run_summary = {
        "stage": stage,
        "run_dir": str(run_dir),
        "scan_config": asdict(config),
        "thermal_closure": asdict(thermal_closure),
        "tcs_mode": tcs_mode,
        "summary": summary,
    }
    with (run_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(run_summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return run_summary

def append_process_log(run_id: str, stage: str, command: str, run_summary: dict) -> None:
    log_path = CASE_DIR / "validation_process_log.md"
    now = datetime.now().isoformat(timespec="seconds")
    summary = run_summary["summary"]
    text = (
        f"\n## {now} - {stage} - {run_id}\n\n"
        f"Command: `{command}`\n\n"
        f"Run directory: `{run_summary['run_dir']}`\n\n"
        f"Case count: {summary['case_count']}\n\n"
        f"Converged all cases: {summary['converged_all_cases']}\n\n"
        f"Mean abs relative power error: {summary['mean_abs_p_out_rel_error']}\n\n"
        f"Max abs relative power error: {summary['max_abs_p_out_rel_error']}\n\n"
        "Decision: baseline/diagnostic result recorded for follow-up analysis.\n"
    )
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text)


def _select_cases(stage: str, case_id: str | None) -> List[VenableTable71Case]:
    if case_id:
        matches = [case for case in TABLE71_CASES if case.case_id == case_id]
        if not matches:
            raise ValueError(f"Unknown case_id: {case_id}")
        return matches
    if stage == "single_case_smoke":
        return [TABLE71_CASES[-1]]
    return list(TABLE71_CASES)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Venable Table 7-1 single-TFE validation scans.")
    parser.add_argument("--stage", choices=("single_case_smoke", "baseline_14_cases"), default="single_case_smoke")
    parser.add_argument("--case-id", default=None)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-root", type=Path, default=CASE_DIR / "runs")
    parser.add_argument("--scan-start-v", type=float, default=0.05)
    parser.add_argument("--scan-stop-v", type=float, default=2.0)
    parser.add_argument("--scan-points", type=int, default=40)
    parser.add_argument("--i-guess-a", type=float, default=150.0)
    parser.add_argument("--thermal-model-mode", choices=VALID_THERMAL_MODEL_MODES, default=DEFAULT_THERMAL_CLOSURE.thermal_model_mode)
    parser.add_argument("--emitter-mean-min-k", type=float, default=DEFAULT_THERMAL_CLOSURE.emitter_mean_min_k)
    parser.add_argument("--emitter-mean-max-k", type=float, default=DEFAULT_THERMAL_CLOSURE.emitter_mean_max_k)
    parser.add_argument("--collector-mean-min-k", type=float, default=DEFAULT_THERMAL_CLOSURE.collector_mean_min_k)
    parser.add_argument("--collector-mean-max-k", type=float, default=DEFAULT_THERMAL_CLOSURE.collector_mean_max_k)
    parser.add_argument("--collector-boundary-mode", choices=VALID_COLLECTOR_BOUNDARY_MODES, default=COLLECTOR_BOUNDARY_LINEAR)
    parser.add_argument("--cooling-water-inlet-temperature-k", "--water-inlet-temperature-k", dest="cooling_water_inlet_temperature_k", type=float, default=DEFAULT_THERMAL_CLOSURE.cooling_water_inlet_temperature_k)
    parser.add_argument("--cooling-water-mass-flow-kg-s", "--water-mass-flow-kg-s", dest="cooling_water_mass_flow_kg_s", type=float, default=DEFAULT_THERMAL_CLOSURE.cooling_water_mass_flow_kg_s)
    parser.add_argument("--cooling-water-cp-j-kg-k", type=float, default=DEFAULT_THERMAL_CLOSURE.cooling_water_cp_j_kg_k)
    parser.add_argument("--coolant-heat-pickup-fraction", type=float, default=DEFAULT_THERMAL_CLOSURE.coolant_heat_pickup_fraction)
    parser.add_argument("--water-heat-transfer-coefficient-w-m2-k", "--water-h-w-m2-k", dest="water_heat_transfer_coefficient_w_m2_k", type=float, default=DEFAULT_THERMAL_CLOSURE.water_heat_transfer_coefficient_w_m2_k)
    parser.add_argument("--regulated-he-gap-effective-k-w-m-k", type=float, default=DEFAULT_THERMAL_CLOSURE.regulated_he_gap_effective_k_w_m_k)
    parser.add_argument("--regulated-he-gap-m", type=float, default=DEFAULT_THERMAL_CLOSURE.regulated_he_gap_m)
    parser.add_argument("--unregulated-he-gap-effective-k-w-m-k", type=float, default=DEFAULT_THERMAL_CLOSURE.unregulated_he_gap_effective_k_w_m_k)
    parser.add_argument("--unregulated-he-gap-m", type=float, default=DEFAULT_THERMAL_CLOSURE.unregulated_he_gap_m)
    parser.add_argument("--thermal-network-heat-pickup-fraction", type=float, default=DEFAULT_THERMAL_CLOSURE.thermal_network_heat_pickup_fraction)
    parser.add_argument("--cs-gap-effective-k-w-m-k", type=float, default=DEFAULT_THERMAL_CLOSURE.cs_gap_effective_k_w_m_k)
    parser.add_argument("--collector-extra-resistance-k-per-w", "--extra-thermal-resistance-k-w", dest="collector_extra_resistance_k_per_w", type=float, default=DEFAULT_THERMAL_CLOSURE.collector_extra_resistance_k_per_w)
    parser.add_argument("--axial-shape-amplitude", type=float, default=DEFAULT_THERMAL_CLOSURE.axial_shape_amplitude)
    parser.add_argument("--axial-profile-mode", choices=VALID_AXIAL_PROFILE_MODES, default=AXIAL_PROFILE_COSINE)
    parser.add_argument("--tisa-heated-length-m", type=float, default=DEFAULT_THERMAL_CLOSURE.tisa_heated_length_m)
    parser.add_argument("--emitter-quadratic-peak-k", type=float, default=DEFAULT_THERMAL_CLOSURE.emitter_quadratic_peak_k)
    parser.add_argument("--collector-quadratic-peak-k", type=float, default=DEFAULT_THERMAL_CLOSURE.collector_quadratic_peak_k)
    parser.add_argument("--tcs-mode", choices=VALID_TCS_MODES, default=TCS_MODE_PLACEHOLDER)
    args = parser.parse_args(argv)

    config = VoltageScanConfig(
        start_v=args.scan_start_v,
        stop_v=args.scan_stop_v,
        points=args.scan_points,
        i_guess_a=args.i_guess_a,
    )
    thermal_closure = VenableThermalClosure(
        thermal_model_mode=args.thermal_model_mode,
        emitter_mean_min_k=args.emitter_mean_min_k,
        emitter_mean_max_k=args.emitter_mean_max_k,
        collector_mean_min_k=args.collector_mean_min_k,
        collector_mean_max_k=args.collector_mean_max_k,
        collector_boundary_mode=args.collector_boundary_mode,
        cooling_water_inlet_temperature_k=args.cooling_water_inlet_temperature_k,
        cooling_water_mass_flow_kg_s=args.cooling_water_mass_flow_kg_s,
        cooling_water_cp_j_kg_k=args.cooling_water_cp_j_kg_k,
        coolant_heat_pickup_fraction=args.coolant_heat_pickup_fraction,
        water_heat_transfer_coefficient_w_m2_k=args.water_heat_transfer_coefficient_w_m2_k,
        regulated_he_gap_effective_k_w_m_k=args.regulated_he_gap_effective_k_w_m_k,
        regulated_he_gap_m=args.regulated_he_gap_m,
        unregulated_he_gap_effective_k_w_m_k=args.unregulated_he_gap_effective_k_w_m_k,
        unregulated_he_gap_m=args.unregulated_he_gap_m,
        thermal_network_heat_pickup_fraction=args.thermal_network_heat_pickup_fraction,
        cs_gap_effective_k_w_m_k=args.cs_gap_effective_k_w_m_k,
        collector_extra_resistance_k_per_w=args.collector_extra_resistance_k_per_w,
        axial_shape_amplitude=args.axial_shape_amplitude,
        axial_profile_mode=args.axial_profile_mode,
        tisa_heated_length_m=args.tisa_heated_length_m,
        emitter_quadratic_peak_k=args.emitter_quadratic_peak_k,
        collector_quadratic_peak_k=args.collector_quadratic_peak_k,
    )
    cases = _select_cases(args.stage, args.case_id)
    run_dir = create_run_dir(args.stage, run_id=args.run_id, output_root=args.output_root)
    snapshot_inputs(run_dir)
    results = run_validation_cases(cases, config, thermal_closure=thermal_closure, tcs_mode=args.tcs_mode)
    run_summary = write_run_outputs(run_dir, args.stage, config, thermal_closure, args.tcs_mode, results)
    append_process_log(
        run_id=run_dir.name,
        stage=args.stage,
        command=" ".join(sys.argv if argv is None else ["venable_validation_runner.py", *argv]),
        run_summary=run_summary,
    )
    print(json.dumps(run_summary["summary"], indent=2, ensure_ascii=False))
    print(f"run_dir={run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
