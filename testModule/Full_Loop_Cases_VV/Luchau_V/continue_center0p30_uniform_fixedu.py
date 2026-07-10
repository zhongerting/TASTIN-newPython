from __future__ import annotations

import argparse
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Continue center-0.3m uniform Luchau single TFE with coupled ThermoCalc fixed_u.")
    parser.add_argument("--restart-in", type=Path, required=True)
    parser.add_argument("--output-label", type=str, required=True)
    parser.add_argument("--duration-s", type=float, required=True)
    parser.add_argument("--dt-s", type=float, default=0.05)
    parser.add_argument("--voltage-v", type=float, default=1.0)
    parser.add_argument("--cesium-pressure-torr", type=float, default=3.55)
    parser.add_argument("--thermal-power-w", type=float, default=3260.0)
    return parser


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


def run_case(args: argparse.Namespace) -> dict[str, Any]:
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    config = LuchauSingleTFEConfig(
        thermal_power_w=float(args.thermal_power_w),
        target_voltage_v=float(args.voltage_v),
        heater_length_m=0.30,
        cesium_pressure_torr=float(args.cesium_pressure_torr),
        i_guess_a=185.0,
        wire_resistance_ohm=0.0,
    )
    build = build_luchau_single_tfe(config)
    build["system"].load_global_state(str(args.restart_in))
    build["tfe"].clear_tec_sources()
    build["tfe"].update_neutronic_power(float(config.thermal_power_w), alpha=1.0)

    lookup = _ensure_lookup_loaded(DEFAULT_LOOKUP_DB, DEFAULT_LOOKUP_REGIONS)
    thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
    configure_luchau_thermocalc(thermo_model, build, config)
    thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
        name=f"Luchau_Center0p30_Cs{config.cesium_pressure_torr:g}_U{config.target_voltage_v:g}",
        tfe=build["tfe"],
        thermo_model=thermo_model,
        alpha_tec=1.0,
    )
    build["system"].add_component(thermo_coupler)

    output_dir = CASE_DIR / "runs" / str(args.output_label)
    output_dir.mkdir(parents=True, exist_ok=True)
    duration_s = float(args.duration_s)
    dt_s = float(args.dt_s)
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

    restart_out = output_dir / "final_restart.npz"
    build["system"].save_global_state(str(restart_out))
    history_csv = output_dir / "history.csv"
    with history_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    summary = {
        "case": "Luchau_V_single_TFE_center0p30_uniform_fixed_u_continuation",
        "config": asdict(config),
        "restart_in": str(args.restart_in),
        "output_dir": str(output_dir),
        "lookup": lookup,
        "duration_s": duration_s,
        "dt_s": dt_s,
        "steps": steps,
        "wall_time_s": wall_time_s,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(config.cesium_pressure_torr)),
        "final": history[-1],
        "files": {
            "history_csv": str(history_csv),
            "restart": str(restart_out),
        },
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return summary


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    print(json.dumps(run_case(args), indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
