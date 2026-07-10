from __future__ import annotations

import csv
import json
import time
from dataclasses import asdict
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
OUTPUT_LABEL = "lookup_v13_axial_power_3260W_cs3p55_1p00V_plus200s_dt0p05"


def _surface_temperature(solid: Any, boundary_name: str, fallback: np.ndarray) -> np.ndarray:
    boundary = getattr(solid, "boundaries", {}).get(boundary_name)
    if boundary is not None and hasattr(boundary, "T_surface"):
        return np.asarray(boundary.T_surface, dtype=float).ravel()
    return np.asarray(fallback, dtype=float).ravel()


def _tec_profile_row(build: dict[str, Any], thermo_model: Any) -> list[dict[str, float]]:
    tfe = build["tfe"]
    node_lengths = np.asarray(build["node_lengths_m"], dtype=float)
    z = np.cumsum(node_lengths) - 0.5 * node_lengths
    emitter_surface_t = _surface_temperature(tfe.solids["emitter"], "right", tfe.solids["emitter"].T)
    collector_surface_t = _surface_temperature(tfe.solids["collector"], "left", tfe.solids["collector"].T)
    res = thermo_model.get_tec_results(0) or {}
    j_cm2 = np.asarray(res.get("J", np.zeros_like(z)), dtype=float)
    ue = np.asarray(res.get("UE", np.zeros_like(z)), dtype=float)
    uc = np.asarray(res.get("UC", np.zeros_like(z)), dtype=float)
    vd = np.asarray(res.get("Vd", np.zeros_like(z)), dtype=float)
    return [
        {
            "axial_index": int(i),
            "axial_center_m": float(z[i]),
            "current_density_a_cm2": float(j_cm2[i]),
            "current_density_a_m2": float(j_cm2[i] * 1.0e4),
            "emitter_potential_ue_v": float(ue[i]),
            "collector_potential_uc_v": float(uc[i]),
            "ue_minus_uc_v": float(ue[i] - uc[i]),
            "vd_v": float(vd[i]),
            "emitter_surface_temperature_k": float(emitter_surface_t[i]),
            "collector_surface_temperature_k": float(collector_surface_t[i]),
        }
        for i in range(len(z))
    ]


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


def main() -> int:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    output_dir = CASE_DIR / "runs" / OUTPUT_LABEL
    output_dir.mkdir(parents=True, exist_ok=True)

    config = LuchauSingleTFEConfig(
        thermal_power_w=3260.0,
        target_voltage_v=1.0,
        heater_length_m=0.30,
        cesium_pressure_torr=3.55,
        i_guess_a=185.0,
        wire_resistance_ohm=0.0,
    )
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

    lookup_info = _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)
    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name="Luchau_SingleTFE_ThermoCalc_Cs3p55",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    dt_s = 0.05
    duration_s = 200.0
    steps = int(round(duration_s / dt_s))
    history: list[dict[str, Any]] = []
    start_wall = time.perf_counter()

    previous_t = np.concatenate(
        [np.asarray(solid.T, dtype=float).ravel() for solid in build["tfe"].solids.values() if hasattr(solid, "T")]
    )
    for step in range(steps):
        build["system"].step(dt_s)
        current_t = np.concatenate(
            [np.asarray(solid.T, dtype=float).ravel() for solid in build["tfe"].solids.values() if hasattr(solid, "T")]
        )
        rate = float(np.max(np.abs(current_t - previous_t)) / dt_s)
        previous_t = current_t
        if step == 0 or (step + 1) % 100 == 0 or step + 1 == steps:
            global_results = thermo_coupler.last_global_results or {}
            stats = _temperature_stats(build)
            history.append(
                {
                    "step": int(step + 1),
                    "elapsed_s": float((step + 1) * dt_s),
                    "absolute_time_s": float(build["system"].global_time),
                    "Iout_a": float(global_results.get("Iout", np.nan)),
                    "Uout_v": float(global_results.get("Uout", np.nan)),
                    "electric_power_w": float(global_results.get("Iout", np.nan))
                    * float(global_results.get("Uout", np.nan)),
                    "max_temperature_rate_k_s": rate,
                    **stats,
                }
            )

    wall_time_s = time.perf_counter() - start_wall
    final_global = thermo_coupler.last_global_results or {}
    profile_rows = _tec_profile_row(build, thermo_model)
    final_restart = output_dir / "final_restart.npz"
    build["system"].save_global_state(str(final_restart))

    history_csv = output_dir / "history.csv"
    with history_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)

    profile_csv = output_dir / "tec_axial_profile_final.csv"
    with profile_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(profile_rows[0].keys()))
        writer.writeheader()
        writer.writerows(profile_rows)

    j = np.array([row["current_density_a_cm2"] for row in profile_rows], dtype=float)
    ueuc = np.array([row["ue_minus_uc_v"] for row in profile_rows], dtype=float)
    summary = {
        "case": "Luchau_V_single_TFE_v13_axial_power_cs3p55_fixed_u_1p00_continuation",
        "config": asdict(config),
        "thermal_restart": str(THERMAL_RESTART),
        "output_dir": str(output_dir),
        "lookup": lookup_info,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "wall_time_s": wall_time_s,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(config.cesium_pressure_torr)),
        "axial_power_profile_sum": float(np.sum(axial_profile)),
        "fission_gas_gap": {
            "h_eq_w_m2_k": 5678.0,
            "k_gas_w_m_k": 0.8517,
        },
        "final": {
            **_temperature_stats(build),
            "Iout_a": float(final_global.get("Iout", np.nan)),
            "Uout_v": float(final_global.get("Uout", np.nan)),
            "electric_power_w": float(final_global.get("Iout", np.nan)) * float(final_global.get("Uout", np.nan)),
            "converged": bool(final_global.get("converged", False)),
            "iteration_count": final_global.get("iteration_count"),
            "zero_emission_skipped": bool(final_global.get("zero_emission_skipped", False)),
            "J_mean_A_cm2": float(np.mean(j)),
            "J_max_A_cm2": float(np.max(j)),
            "UE_minus_UC_mean_v": float(np.mean(ueuc)),
            "UE_minus_UC_min_v": float(np.min(ueuc)),
            "UE_minus_UC_max_v": float(np.max(ueuc)),
        },
        "files": {
            "history_csv": str(history_csv),
            "profile_csv": str(profile_csv),
            "restart": str(final_restart),
        },
    }
    summary_json = output_dir / "summary.json"
    with summary_json.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
