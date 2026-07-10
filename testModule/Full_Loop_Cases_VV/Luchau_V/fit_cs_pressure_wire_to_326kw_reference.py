from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Any
import sys

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases_VV.Luchau_V.luchau_single_tfe_model import (
    LuchauSingleTFEConfig,
    LuchauSingleTFEThermoCalcCoupler,
    build_luchau_single_tfe,
    configure_luchau_thermocalc,
    tcs_from_cesium_pressure,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_preheat_then_sweep import (
    DEFAULT_LOOKUP_DB,
    DEFAULT_LOOKUP_REGIONS,
    _ensure_lookup_loaded,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import json_default


TARGET_CSV = CASE_DIR / "reference_image_326kw_user_digitized_iv.csv"
THERMAL_RESTART = CASE_DIR / "runs" / "lookup_preheat_sweep_3260W_dt0p02_1000s" / "preheat_1000s_restart.npz"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Fit Cs pressure and scalar wire resistance against 3.26 kW I-V reference.")
    parser.add_argument("--pressure-list", type=str, required=True)
    parser.add_argument("--wire-list", type=str, required=True, help="Comma-separated scalar wire resistance values in ohm.")
    parser.add_argument("--duration-s", type=float, default=200.0)
    parser.add_argument("--dt-s", type=float, default=0.05)
    parser.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) // 2)))
    parser.add_argument("--output-label", type=str, default=None)
    parser.add_argument("--target-csv", type=Path, default=TARGET_CSV)
    parser.add_argument("--restart-in", type=Path, default=THERMAL_RESTART)
    parser.add_argument("--point-index-min", type=int, default=6)
    parser.add_argument("--point-index-max", type=int, default=20)
    return parser


def _parse_float_list(value: str, *, positive: bool = False, nonnegative: bool = False) -> list[float]:
    values = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not values:
        raise ValueError("list argument must contain at least one value.")
    for item in values:
        if positive and item <= 0.0:
            raise ValueError("values must be positive.")
        if nonnegative and item < 0.0:
            raise ValueError("values must be non-negative.")
    return values


def _load_targets(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "point_index": int(row["point_index"]),
                "voltage_v": float(row["voltage_v"]),
                "target_current_a": float(row["current_a"]),
            })
    return rows


def _temperature_vector(build: dict[str, Any]) -> np.ndarray:
    return np.concatenate(
        [np.asarray(solid.T, dtype=float).ravel() for solid in build["tfe"].solids.values() if hasattr(solid, "T")]
    )


def _temperature_stats(build: dict[str, Any]) -> dict[str, float]:
    tfe = build["tfe"]
    pellet = np.asarray(tfe.solids["pellet"].T, dtype=float)
    emitter = np.asarray(tfe.solids["emitter"].T, dtype=float)
    collector = np.asarray(tfe.solids["collector"].T, dtype=float)
    coolant = np.array([float(vol.T) for vol in build["channel"].volumes], dtype=float)
    return {
        "pellet_temperature_mean_k": float(np.mean(pellet)),
        "pellet_temperature_max_k": float(np.max(pellet)),
        "emitter_temperature_mean_k": float(np.mean(emitter)),
        "emitter_temperature_max_k": float(np.max(emitter)),
        "collector_temperature_mean_k": float(np.mean(collector)),
        "collector_temperature_max_k": float(np.max(collector)),
        "coolant_outlet_temperature_k": float(coolant[-1]),
    }


def _electric_summary(thermo_model: Any) -> dict[str, Any]:
    global_results = thermo_model.get_global_results() or {}
    res = thermo_model.get_tec_results(0) or {}
    j = np.asarray(res.get("J", []), dtype=float)
    ue = np.asarray(res.get("UE", []), dtype=float)
    uc = np.asarray(res.get("UC", []), dtype=float)
    ueuc = ue - uc if ue.size and uc.size else np.array([], dtype=float)
    current = _finite_float(global_results.get("Iout"))
    voltage = _finite_float(global_results.get("Uout"))
    return {
        "Iout_a": current,
        "Uout_v": voltage,
        "electric_power_w": None if current is None or voltage is None else float(current * voltage),
        "converged": bool(global_results.get("converged", False)),
        "iteration_count": global_results.get("iteration_count"),
        "zero_emission_skipped": bool(global_results.get("zero_emission_skipped", False)),
        "J_mean_A_cm2": None if j.size == 0 else float(np.mean(j)),
        "J_max_A_cm2": None if j.size == 0 else float(np.max(j)),
        "UE_minus_UC_mean_v": None if ueuc.size == 0 else float(np.mean(ueuc)),
    }


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _run_one_case(task: dict[str, Any]) -> dict[str, Any]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    pressure_torr = float(task["pressure_torr"])
    wire_resistance_ohm = float(task["wire_resistance_ohm"])
    voltage_v = float(task["voltage_v"])
    duration_s = float(task["duration_s"])
    dt_s = float(task["dt_s"])

    config = LuchauSingleTFEConfig(
        thermal_power_w=3260.0,
        target_voltage_v=voltage_v,
        heater_length_m=0.30,
        cesium_pressure_torr=pressure_torr,
        i_guess_a=max(50.0, float(task["target_current_a"])),
        wire_resistance_ohm=wire_resistance_ohm,
    )
    build = build_luchau_single_tfe(config)
    build["system"].load_global_state(str(task["restart_in"]))
    build["tfe"].clear_tec_sources()
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)

    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name=f"Luchau_fit_Pcs{pressure_torr:g}_Rw{wire_resistance_ohm:g}_U{voltage_v:g}",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    steps = int(round(duration_s / dt_s))
    previous = _temperature_vector(build)
    last_rate = math.nan
    wall_start = time.perf_counter()
    for _ in range(steps):
        build["system"].step(dt_s)
        current = _temperature_vector(build)
        last_rate = float(np.max(np.abs(current - previous)) / dt_s)
        previous = current
    wall_time_s = time.perf_counter() - wall_start

    electric = _electric_summary(thermo_model)
    iout = electric["Iout_a"]
    target = float(task["target_current_a"])
    error = None if iout is None else float(iout - target)
    return {
        "pressure_torr": pressure_torr,
        "wire_resistance_ohm": wire_resistance_ohm,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(pressure_torr),
        "point_index": int(task["point_index"]),
        "fit_selected": bool(task["fit_selected"]),
        "target_voltage_v": voltage_v,
        "target_current_a": target,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "wall_time_s": wall_time_s,
        "max_temperature_rate_k_s": last_rate,
        **_temperature_stats(build),
        **electric,
        "current_error_a": error,
        "current_abs_error_a": None if error is None else abs(error),
    }


def _metrics(rows: list[dict[str, Any]], *, selected_only: bool) -> dict[str, Any]:
    valid = [row for row in rows if row.get("Iout_a") is not None and bool(row.get("converged"))]
    if selected_only:
        valid = [row for row in valid if bool(row.get("fit_selected"))]
    total = [row for row in rows if (not selected_only or bool(row.get("fit_selected")))]
    errors = np.array([float(row["current_error_a"]) for row in valid], dtype=float)
    if errors.size == 0:
        return {"n_points": len(total), "n_converged": 0, "rmse_a": None, "mae_a": None, "max_abs_error_a": None, "mean_abs_rel_error": None, "signed_mean_error_a": None}
    abs_errors = np.abs(errors)
    targets = np.array([float(row["target_current_a"]) for row in valid], dtype=float)
    return {
        "n_points": len(total),
        "n_converged": len(valid),
        "rmse_a": float(np.sqrt(np.mean(errors**2))),
        "mae_a": float(np.mean(abs_errors)),
        "max_abs_error_a": float(np.max(abs_errors)),
        "mean_abs_rel_error": float(np.mean(abs_errors / np.maximum(np.abs(targets), 1.0))),
        "signed_mean_error_a": float(np.mean(errors)),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    pressures = _parse_float_list(args.pressure_list, positive=True)
    wires = _parse_float_list(args.wire_list, nonnegative=True)
    targets = _load_targets(Path(args.target_csv))
    label = args.output_label or f"fit_cs_wire_326kw_{int(time.time())}"
    output_dir = CASE_DIR / "runs" / label
    output_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for pressure in pressures:
        for wire in wires:
            for target in targets:
                point_index = int(target["point_index"])
                selected = int(args.point_index_min) <= point_index <= int(args.point_index_max)
                tasks.append({
                    **target,
                    "pressure_torr": float(pressure),
                    "wire_resistance_ohm": float(wire),
                    "fit_selected": selected,
                    "duration_s": float(args.duration_s),
                    "dt_s": float(args.dt_s),
                    "restart_in": str(args.restart_in),
                })

    wall_start = time.perf_counter()
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=int(args.workers)) as executor:
        futures = [executor.submit(_run_one_case, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(
                f"Pcs={row['pressure_torr']:.4g} torr Rw={row['wire_resistance_ohm']:.4g} ohm "
                f"U={row['target_voltage_v']:.5g} I={row['Iout_a']:.3f} "
                f"target={row['target_current_a']:.3f} err={row['current_error_a']:.3f}",
                flush=True,
            )
    rows.sort(key=lambda r: (float(r["pressure_torr"]), float(r["wire_resistance_ohm"]), int(r["point_index"])))

    metric_rows = []
    for pressure in pressures:
        for wire in wires:
            param_rows = [
                row for row in rows
                if abs(float(row["pressure_torr"]) - pressure) < 1.0e-12
                and abs(float(row["wire_resistance_ohm"]) - wire) < 1.0e-15
            ]
            all_m = _metrics(param_rows, selected_only=False)
            sel_m = _metrics(param_rows, selected_only=True)
            metric_rows.append({
                "pressure_torr": float(pressure),
                "wire_resistance_ohm": float(wire),
                **{f"all_{key}": value for key, value in all_m.items()},
                **{f"selected_{key}": value for key, value in sel_m.items()},
            })
    metric_rows.sort(key=lambda r: float("inf") if r["selected_rmse_a"] is None else float(r["selected_rmse_a"]))

    results_csv = output_dir / "fit_results_by_point.csv"
    metrics_csv = output_dir / "fit_metrics_by_params.csv"
    _write_csv(results_csv, rows)
    _write_csv(metrics_csv, metric_rows)
    summary = {
        "case": "center0p30_uniform_326kw_fit_cs_pressure_and_wire_to_user_digitized_iv",
        "target_csv": str(args.target_csv),
        "restart_in": str(args.restart_in),
        "output_dir": str(output_dir),
        "pressure_points_torr": pressures,
        "wire_resistance_points_ohm": wires,
        "selected_point_index_min": int(args.point_index_min),
        "selected_point_index_max": int(args.point_index_max),
        "duration_s": float(args.duration_s),
        "dt_s": float(args.dt_s),
        "workers": int(args.workers),
        "total_wall_time_s": time.perf_counter() - wall_start,
        "best_by_selected_rmse": metric_rows[0] if metric_rows else None,
        "metrics": metric_rows,
        "files": {"results_csv": str(results_csv), "metrics_csv": str(metrics_csv)},
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
