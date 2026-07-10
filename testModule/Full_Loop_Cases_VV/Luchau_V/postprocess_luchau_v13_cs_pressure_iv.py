from __future__ import annotations

import csv
import json
from dataclasses import replace
from pathlib import Path
from typing import Any
import sys

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.common_core_builder import build_axial_power_profile
from testModule.Full_Loop_Cases_VV.Luchau_V.luchau_single_tfe_model import (
    LUCHAU_AXIAL_NODE_ALLOCATION,
    LuchauSingleTFEConfig,
    build_luchau_single_tfe,
    configure_luchau_thermocalc,
    tcs_from_cesium_pressure,
)
from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import json_default


RUN_DIR = CASE_DIR / "runs" / "lookup_v13_axial_power_3260W_cs_pressure_sweep_200s_dt0p05"
PRESSURE_POINTS_TORR = (0.4, 1.0, 2.0, 3.55, 5.0)
VOLTAGE_POINTS_V = (0.4, 0.6, 0.8, 1.0)
I_GUESSES_A = (50.0, 80.0, 120.0, 180.0, 250.0, 350.0, 550.0, 800.0)


def _build_loaded_case(config: LuchauSingleTFEConfig, restart_path: Path) -> dict[str, Any]:
    build = build_luchau_single_tfe(config)
    n_lower, n_active, n_upper = LUCHAU_AXIAL_NODE_ALLOCATION
    axial_profile = build_axial_power_profile(n_lower, n_active, n_upper)
    build["tfe"].axial_power_profile = axial_profile
    build["tfe"].solids["pellet"].set_axial_power_profile(axial_profile)
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    build["system"].load_global_state(str(restart_path))
    build["tfe"].axial_power_profile = axial_profile
    build["tfe"].solids["pellet"].set_axial_power_profile(axial_profile)
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    return build


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


def _surface_temperatures(build: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    tfe = build["tfe"]
    emitter_boundary = tfe.solids["emitter"].boundaries["right"]
    collector_boundary = tfe.solids["collector"].boundaries["left"]
    return (
        np.asarray(emitter_boundary.T_surface, dtype=float).reshape(1, tfe.mesh.n_axial),
        np.asarray(collector_boundary.T_surface, dtype=float).reshape(1, tfe.mesh.n_axial),
    )


def _electric_summary(model: Any) -> dict[str, Any]:
    global_results = model.get_global_results() or {}
    res = model.get_tec_results(0) or {}
    j = np.asarray(res.get("J", []), dtype=float)
    ue = np.asarray(res.get("UE", []), dtype=float)
    uc = np.asarray(res.get("UC", []), dtype=float)
    vd = np.asarray(res.get("Vd", []), dtype=float)
    current = _finite_float(global_results.get("Iout"))
    voltage = _finite_float(global_results.get("Uout"))
    ueuc = ue - uc if ue.size and uc.size else np.array([], dtype=float)
    return {
        "Iout_a": current,
        "Uout_v": voltage,
        "electric_power_w": None if current is None or voltage is None else float(current * voltage),
        "converged": bool(global_results.get("converged", False)),
        "iteration_count": global_results.get("iteration_count"),
        "zero_emission_skipped": bool(global_results.get("zero_emission_skipped", False)),
        "J_mean_A_cm2": None if j.size == 0 else float(np.mean(j)),
        "J_max_A_cm2": None if j.size == 0 else float(np.max(j)),
        "Vd_mean_v": None if vd.size == 0 else float(np.mean(vd)),
        "UE_minus_UC_mean_v": None if ueuc.size == 0 else float(np.mean(ueuc)),
        "UE_minus_UC_min_v": None if ueuc.size == 0 else float(np.min(ueuc)),
        "UE_minus_UC_max_v": None if ueuc.size == 0 else float(np.max(ueuc)),
    }


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    return number


def _calculate_with_retries(build: dict[str, Any], base_config: LuchauSingleTFEConfig, voltage: float) -> dict[str, Any]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    emitter_t, collector_t = _surface_temperatures(build)
    attempts: list[dict[str, Any]] = []
    for guess in I_GUESSES_A:
        config = replace(base_config, target_voltage_v=float(voltage), i_guess_a=float(guess))
        model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
        configure_luchau_thermocalc(model, build, config)
        model.set_temperatures(emitter_t, collector_t)
        model.calculate(verbose=False)
        result = _electric_summary(model)
        result["i_guess_a"] = float(guess)
        attempts.append(result)
        if result["converged"] and result["Iout_a"] is not None:
            return result
    best = attempts[-1]
    best["failed_attempts"] = attempts
    return best


def main() -> int:
    rows: list[dict[str, Any]] = []
    for pressure in PRESSURE_POINTS_TORR:
        config = LuchauSingleTFEConfig(
            thermal_power_w=3260.0,
            target_voltage_v=1.0,
            heater_length_m=0.30,
            cesium_pressure_torr=float(pressure),
            i_guess_a=185.0,
            wire_resistance_ohm=0.0,
        )
        restart_path = RUN_DIR / f"cs_{str(pressure).replace('.', 'p')}_final_restart.npz"
        build = _build_loaded_case(config, restart_path)
        stats = _temperature_stats(build)
        for voltage in VOLTAGE_POINTS_V:
            row = {
                "pressure_torr": float(pressure),
                "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(pressure)),
                "target_voltage_v": float(voltage),
                **stats,
                **_calculate_with_retries(build, config, float(voltage)),
            }
            rows.append(row)

    csv_path = RUN_DIR / "pressure_voltage_sweep_retry.csv"
    serializable_rows = [
        {key: value for key, value in row.items() if key != "failed_attempts"}
        for row in rows
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(serializable_rows[0].keys()))
        writer.writeheader()
        writer.writerows(serializable_rows)

    summary = {
        "case": "Luchau_V_single_TFE_v13_axial_power_cs_pressure_iv_retry",
        "source_run_dir": str(RUN_DIR),
        "pressure_points_torr": list(PRESSURE_POINTS_TORR),
        "voltage_points_v": list(VOLTAGE_POINTS_V),
        "i_guesses_a": list(I_GUESSES_A),
        "rows": rows,
        "csv": str(csv_path),
    }
    summary_path = RUN_DIR / "pressure_voltage_sweep_retry_summary.json"
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
