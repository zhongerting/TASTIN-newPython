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
    parser = argparse.ArgumentParser(description="Fit Cs pressure against user-digitized 3.26 kW I-V reference.")
    parser.add_argument("--pressure-list", type=str, required=True, help="Comma-separated Cs pressures in torr.")
    parser.add_argument("--duration-s", type=float, default=200.0)
    parser.add_argument("--dt-s", type=float, default=0.05)
    parser.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) // 2)))
    parser.add_argument("--output-label", type=str, default=None)
    parser.add_argument("--target-csv", type=Path, default=TARGET_CSV)
    parser.add_argument("--restart-in", type=Path, default=THERMAL_RESTART)
    return parser


def _parse_pressures(value: str) -> list[float]:
    pressures = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not pressures:
        raise ValueError("pressure-list must contain at least one pressure.")
    for pressure in pressures:
        if pressure <= 0.0:
            raise ValueError("all pressures must be positive.")
    return pressures


def _load_targets(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "point_index": int(row["point_index"]),
                    "voltage_v": float(row["voltage_v"]),
                    "target_current_a": float(row["current_a"]),
                }
            )
    if not rows:
        raise ValueError(f"no target rows found in {path}")
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
    emitter_surface = np.asarray(tfe.solids["emitter"].boundaries["right"].T_surface, dtype=float)
    collector_surface = np.asarray(tfe.solids["collector"].boundaries["left"].T_surface, dtype=float)
    return {
        "pellet_temperature_mean_k": float(np.mean(pellet)),
        "pellet_temperature_max_k": float(np.max(pellet)),
        "emitter_temperature_mean_k": float(np.mean(emitter)),
        "emitter_temperature_max_k": float(np.max(emitter)),
        "emitter_surface_temperature_mean_k": float(np.mean(emitter_surface)),
        "emitter_surface_temperature_max_k": float(np.max(emitter_surface)),
        "collector_temperature_mean_k": float(np.mean(collector)),
        "collector_temperature_max_k": float(np.max(collector)),
        "collector_surface_temperature_mean_k": float(np.mean(collector_surface)),
        "collector_surface_temperature_max_k": float(np.max(collector_surface)),
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
    voltage_v = float(task["voltage_v"])
    duration_s = float(task["duration_s"])
    dt_s = float(task["dt_s"])
    restart_in = Path(task["restart_in"])

    config = LuchauSingleTFEConfig(
        thermal_power_w=3260.0,
        target_voltage_v=voltage_v,
        heater_length_m=0.30,
        cesium_pressure_torr=pressure_torr,
        i_guess_a=max(50.0, float(task["target_current_a"])),
        wire_resistance_ohm=0.0,
    )
    build = build_luchau_single_tfe(config)
    build["system"].load_global_state(str(restart_in))
    build["tfe"].clear_tec_sources()
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)

    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name=f"Luchau_Center0p30_fit_Pcs{pressure_torr:g}_U{voltage_v:g}",
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
    stats = _temperature_stats(build)
    iout = electric["Iout_a"]
    target = float(task["target_current_a"])
    error = None if iout is None else float(iout - target)
    abs_error = None if error is None else abs(error)
    return {
        "pressure_torr": pressure_torr,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(pressure_torr),
        "point_index": int(task["point_index"]),
        "target_voltage_v": voltage_v,
        "target_current_a": target,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "wall_time_s": wall_time_s,
        "max_temperature_rate_k_s": last_rate,
        **stats,
        **electric,
        "current_error_a": error,
        "current_abs_error_a": abs_error,
    }


def _pressure_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    valid = [row for row in rows if row.get("Iout_a") is not None and bool(row.get("converged"))]
    errors = np.array([float(row["current_error_a"]) for row in valid], dtype=float)
    abs_errors = np.abs(errors)
    targets = np.array([float(row["target_current_a"]) for row in valid], dtype=float)
    rel = abs_errors / np.maximum(np.abs(targets), 1.0)
    return {
        "n_points": len(rows),
        "n_converged": len(valid),
        "rmse_a": None if len(errors) == 0 else float(np.sqrt(np.mean(errors**2))),
        "mae_a": None if len(errors) == 0 else float(np.mean(abs_errors)),
        "max_abs_error_a": None if len(errors) == 0 else float(np.max(abs_errors)),
        "mean_abs_rel_error": None if len(errors) == 0 else float(np.mean(rel)),
        "signed_mean_error_a": None if len(errors) == 0 else float(np.mean(errors)),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write empty CSV")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    pressures = _parse_pressures(args.pressure_list)
    targets = _load_targets(Path(args.target_csv))
    label = args.output_label or f"fit_cs_pressure_326kw_{int(time.time())}"
    output_dir = CASE_DIR / "runs" / label
    output_dir.mkdir(parents=True, exist_ok=True)

    tasks = []
    for pressure in pressures:
        for target in targets:
            tasks.append(
                {
                    **target,
                    "pressure_torr": float(pressure),
                    "duration_s": float(args.duration_s),
                    "dt_s": float(args.dt_s),
                    "restart_in": str(args.restart_in),
                }
            )

    wall_start = time.perf_counter()
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=int(args.workers)) as executor:
        futures = [executor.submit(_run_one_case, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(
                f"Pcs={row['pressure_torr']:.4g} torr U={row['target_voltage_v']:.5g} "
                f"I={row['Iout_a']:.3f} target={row['target_current_a']:.3f} "
                f"err={row['current_error_a']:.3f}",
                flush=True,
            )
    rows.sort(key=lambda r: (float(r["pressure_torr"]), int(r["point_index"])))

    metrics = []
    for pressure in pressures:
        pressure_rows = [row for row in rows if abs(float(row["pressure_torr"]) - pressure) < 1.0e-12]
        metrics.append({"pressure_torr": float(pressure), **_pressure_metrics(pressure_rows)})
    metrics.sort(key=lambda r: (float("inf") if r["rmse_a"] is None else float(r["rmse_a"])))

    results_csv = output_dir / "fit_results_by_point.csv"
    metrics_csv = output_dir / "fit_metrics_by_pressure.csv"
    _write_csv(results_csv, rows)
    _write_csv(metrics_csv, metrics)

    summary = {
        "case": "center0p30_uniform_326kw_fit_cesium_pressure_to_user_digitized_iv",
        "target_csv": str(args.target_csv),
        "restart_in": str(args.restart_in),
        "output_dir": str(output_dir),
        "pressure_points_torr": pressures,
        "duration_s": float(args.duration_s),
        "dt_s": float(args.dt_s),
        "workers": int(args.workers),
        "total_wall_time_s": time.perf_counter() - wall_start,
        "best": metrics[0] if metrics else None,
        "metrics": metrics,
        "files": {
            "results_csv": str(results_csv),
            "metrics_csv": str(metrics_csv),
        },
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
