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


THERMAL_RESTART = (
    CASE_DIR
    / "runs"
    / "lookup_preheat_sweep_3260W_dt0p02_1000s"
    / "preheat_1000s_restart.npz"
)
OUTPUT_LABEL = "lookup_center0p30_uniform_3260W_cs3p55_1p00V_plus200s_dt0p05"


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
    vd = np.asarray(res.get("Vd", []), dtype=float)
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


def _surface_temperature(solid: Any, boundary_name: str) -> np.ndarray:
    boundary = getattr(solid, "boundaries", {}).get(boundary_name)
    if boundary is not None and hasattr(boundary, "T_surface"):
        return np.asarray(boundary.T_surface, dtype=float).ravel()
    return np.asarray(solid.T, dtype=float).ravel()


def _write_axial_profile(build: dict[str, Any], thermo_model: Any, path: Path) -> None:
    tfe = build["tfe"]
    node_lengths = np.asarray(build["node_lengths_m"], dtype=float)
    z = np.cumsum(node_lengths) - 0.5 * node_lengths
    res = thermo_model.get_tec_results(0) or {}
    j = np.asarray(res.get("J", np.zeros_like(z)), dtype=float)
    ue = np.asarray(res.get("UE", np.zeros_like(z)), dtype=float)
    uc = np.asarray(res.get("UC", np.zeros_like(z)), dtype=float)
    vd = np.asarray(res.get("Vd", np.zeros_like(z)), dtype=float)
    emitter_surface = _surface_temperature(tfe.solids["emitter"], "right")
    collector_surface = _surface_temperature(tfe.solids["collector"], "left")
    rows = []
    for i in range(len(z)):
        rows.append(
            {
                "axial_index": i,
                "axial_center_m": float(z[i]),
                "current_density_a_cm2": float(j[i]),
                "current_density_a_m2": float(j[i] * 1.0e4),
                "emitter_potential_ue_v": float(ue[i]),
                "collector_potential_uc_v": float(uc[i]),
                "ue_minus_uc_v": float(ue[i] - uc[i]),
                "vd_v": float(vd[i]),
                "emitter_surface_temperature_k": float(emitter_surface[i]),
                "collector_surface_temperature_k": float(collector_surface[i]),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


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
    build["system"].load_global_state(str(THERMAL_RESTART))
    build["tfe"].clear_tec_sources()
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)

    lookup_info = _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)
    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name="Luchau_Center0p30_Cs3p55_1V",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    duration_s = 200.0
    dt_s = 0.05
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
            history.append(
                {
                    "step": int(step + 1),
                    "elapsed_s": float((step + 1) * dt_s),
                    "absolute_time_s": float(build["system"].global_time),
                    "max_temperature_rate_k_s": rate,
                    **_temperature_stats(build),
                    **_electric_summary(thermo_model),
                }
            )
    wall_time_s = time.perf_counter() - wall_start

    final_restart = output_dir / "final_restart.npz"
    build["system"].save_global_state(str(final_restart))
    history_csv = output_dir / "history.csv"
    with history_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    profile_csv = output_dir / "tec_axial_profile_final.csv"
    _write_axial_profile(build, thermo_model, profile_csv)

    summary = {
        "case": "Luchau_V_single_TFE_center0p30_uniform_cs3p55_fixed_u_1p00",
        "config": asdict(config),
        "thermal_restart": str(THERMAL_RESTART),
        "output_dir": str(output_dir),
        "lookup": lookup_info,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "wall_time_s": wall_time_s,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(config.cesium_pressure_torr)),
        "heater_profile": "center_0p30m_uniform",
        "fission_gas_gap": {
            "h_eq_w_m2_k": 5678.0,
            "k_gas_w_m_k": 0.8517,
        },
        "final": {
            **_temperature_stats(build),
            **_electric_summary(thermo_model),
        },
        "files": {
            "history_csv": str(history_csv),
            "profile_csv": str(profile_csv),
            "restart": str(final_restart),
        },
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
