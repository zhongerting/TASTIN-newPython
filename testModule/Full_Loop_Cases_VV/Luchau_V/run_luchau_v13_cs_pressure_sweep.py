from __future__ import annotations

import csv
import json
import time
from dataclasses import asdict, replace
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


THERMAL_RESTART = (
    CASE_DIR
    / "runs"
    / "lookup_v13_axial_power_3260W_from0p30thermal_plus500s_sweep"
    / "v13_profile_thermal_500s_restart.npz"
)
OUTPUT_LABEL = "lookup_v13_axial_power_3260W_cs_pressure_sweep_200s_dt0p05"
PRESSURE_POINTS_TORR = (0.4, 1.0, 2.0, 3.55, 5.0)
VOLTAGE_POINTS_V = (0.4, 0.6, 0.8, 1.0)


def _build_case(config: LuchauSingleTFEConfig) -> dict[str, Any]:
    build = build_luchau_single_tfe(config)
    n_lower, n_active, n_upper = LUCHAU_AXIAL_NODE_ALLOCATION
    axial_profile = build_axial_power_profile(n_lower, n_active, n_upper)
    build["tfe"].axial_power_profile = axial_profile
    build["tfe"].solids["pellet"].set_axial_power_profile(axial_profile)
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    build["system"].load_global_state(str(THERMAL_RESTART))
    build["tfe"].clear_tec_sources()
    build["tfe"].axial_power_profile = axial_profile
    build["tfe"].solids["pellet"].set_axial_power_profile(axial_profile)
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)
    build["axial_power_profile"] = axial_profile
    return build


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


def _run_pressure_case(pressure_torr: float, output_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    config = LuchauSingleTFEConfig(
        thermal_power_w=3260.0,
        target_voltage_v=1.0,
        heater_length_m=0.30,
        cesium_pressure_torr=float(pressure_torr),
        i_guess_a=185.0,
        wire_resistance_ohm=0.0,
    )
    build = _build_case(config)
    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name=f"Luchau_SingleTFE_ThermoCalc_Cs{pressure_torr:g}",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    dt_s = 0.05
    duration_s = 200.0
    steps = int(round(duration_s / dt_s))
    history: list[dict[str, Any]] = []
    previous = _temperature_vector(build)
    wall_start = time.perf_counter()
    for step in range(steps):
        build["system"].step(dt_s)
        current = _temperature_vector(build)
        rate = float(np.max(np.abs(current - previous)) / dt_s)
        previous = current
        if step == 0 or (step + 1) % 100 == 0 or step + 1 == steps:
            row = {
                "pressure_torr": float(pressure_torr),
                "step": int(step + 1),
                "elapsed_s": float((step + 1) * dt_s),
                "absolute_time_s": float(build["system"].global_time),
                "max_temperature_rate_k_s": rate,
                **_temperature_stats(build),
                **_electric_summary(thermo_model),
            }
            history.append(row)
    wall_time_s = time.perf_counter() - wall_start

    iv_rows: list[dict[str, Any]] = []
    for voltage in VOLTAGE_POINTS_V:
        v_config = replace(config, target_voltage_v=float(voltage))
        v_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
        configure_luchau_thermocalc(v_model, build, v_config)
        v_model.calculate(verbose=False)
        row = {
            "pressure_torr": float(pressure_torr),
            "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(pressure_torr)),
            "target_voltage_v": float(voltage),
            **_temperature_stats(build),
            **_electric_summary(v_model),
        }
        iv_rows.append(row)

    restart = output_dir / f"cs_{str(pressure_torr).replace('.', 'p')}_final_restart.npz"
    build["system"].save_global_state(str(restart))
    pressure_summary = {
        "pressure_torr": float(pressure_torr),
        "config": asdict(config),
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "wall_time_s": wall_time_s,
        "restart": str(restart),
        "final_at_1p00v": history[-1],
        "iv_rows": iv_rows,
    }
    return history, iv_rows, pressure_summary


def main() -> int:
    output_dir = CASE_DIR / "runs" / OUTPUT_LABEL
    output_dir.mkdir(parents=True, exist_ok=True)
    lookup_info = _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)

    all_history: list[dict[str, Any]] = []
    all_iv: list[dict[str, Any]] = []
    pressure_summaries: list[dict[str, Any]] = []
    total_start = time.perf_counter()
    for pressure in PRESSURE_POINTS_TORR:
        history, iv_rows, pressure_summary = _run_pressure_case(float(pressure), output_dir)
        all_history.extend(history)
        all_iv.extend(iv_rows)
        pressure_summaries.append(pressure_summary)
    total_wall_s = time.perf_counter() - total_start

    history_csv = output_dir / "pressure_relaxation_history.csv"
    with history_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_history[0].keys()))
        writer.writeheader()
        writer.writerows(all_history)

    iv_csv = output_dir / "pressure_voltage_sweep.csv"
    with iv_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_iv[0].keys()))
        writer.writeheader()
        writer.writerows(all_iv)

    summary = {
        "case": "Luchau_V_single_TFE_v13_axial_power_cs_pressure_sweep",
        "thermal_restart": str(THERMAL_RESTART),
        "output_dir": str(output_dir),
        "thermal_power_w": 3260.0,
        "target_voltage_for_relaxation_v": 1.0,
        "pressure_points_torr": list(PRESSURE_POINTS_TORR),
        "voltage_points_v": list(VOLTAGE_POINTS_V),
        "lookup": lookup_info,
        "duration_per_pressure_s": 200.0,
        "dt_s": 0.05,
        "total_wall_time_s": total_wall_s,
        "pressure_summaries": pressure_summaries,
        "files": {
            "history_csv": str(history_csv),
            "iv_csv": str(iv_csv),
        },
    }
    summary_json = output_dir / "summary.json"
    with summary_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
