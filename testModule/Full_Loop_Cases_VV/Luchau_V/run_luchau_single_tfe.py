from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict
from datetime import datetime
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build and optionally run the Luchau_V single-TFE fixed-voltage case.")
    parser.add_argument("--thermal-power-w", type=float, required=True)
    parser.add_argument("--target-voltage-v", type=float, required=True)
    parser.add_argument("--duration-s", type=float, default=0.0, help="Transient duration, or maximum duration when --steady is used.")
    parser.add_argument("--dt-s", type=float, default=0.1)
    parser.add_argument("--steady", action="store_true", help="Stop early when the temperature-rate criterion is continuously satisfied.")
    parser.add_argument("--steady-dtemp-k-s", type=float, default=1.0e-3)
    parser.add_argument("--steady-window-steps", type=int, default=20)
    parser.add_argument("--output-label", type=str, default=None)
    parser.add_argument("--skip-thermocalc-calc", action="store_true")
    parser.add_argument("--skip-transient", action="store_true")
    return parser


def json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def run_case(args: argparse.Namespace) -> dict[str, Any]:
    config = LuchauSingleTFEConfig(
        thermal_power_w=float(args.thermal_power_w),
        target_voltage_v=float(args.target_voltage_v),
    )
    build = build_luchau_single_tfe(config)

    thermo_coupler = None
    if not bool(args.skip_thermocalc_calc):
        from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

        thermo_model = ThermoCalcModel(n_elements=1, n_nodes=build["tfe"].mesh.n_axial)
        configure_luchau_thermocalc(thermo_model, build, config)
        thermo_coupler = LuchauSingleTFEThermoCalcCoupler(
            name="Luchau_SingleTFE_ThermoCalc",
            tfe=build["tfe"],
            thermo_model=thermo_model,
            alpha_tec=1.0,
        )
        build["thermo_coupler"] = thermo_coupler
        build["thermo_model"] = thermo_model
        build["system"].add_component(thermo_coupler)
        if bool(args.skip_transient) or float(args.duration_s) <= 0.0:
            thermo_coupler.sync_thermo_electric()

    transient_result: dict[str, Any] | None = None
    if not bool(args.skip_transient) and float(args.duration_s) > 0.0:
        if bool(args.steady):
            transient_result = _advance_system_until_steady(
                build["system"],
                build,
                float(args.duration_s),
                float(args.dt_s),
                float(args.steady_dtemp_k_s),
                int(args.steady_window_steps),
            )
        else:
            transient_result = _advance_system(build["system"], float(args.duration_s), float(args.dt_s))

    thermo_results = thermo_coupler.last_global_results if thermo_coupler is not None else None
    summary = _build_summary(config, build, thermo_results, transient_result)
    output_dir = _write_summary(summary, args.output_label)
    summary["output_dir"] = str(output_dir)
    return summary


def _advance_system(system, duration_s: float, dt_s: float) -> dict[str, Any]:
    if duration_s < 0.0:
        raise ValueError("duration_s must be non-negative.")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    elapsed = 0.0
    steps = 0
    while elapsed < duration_s - 1.0e-12:
        dt = min(dt_s, duration_s - elapsed)
        system.step(dt)
        elapsed += dt
        steps += 1
    return {"steady_reached": False, "steps": steps, "final_time_s": float(system.global_time), "max_temperature_rate_k_s": None}


def _advance_system_until_steady(
    system,
    build: dict[str, Any],
    duration_s: float,
    dt_s: float,
    steady_dtemp_k_s: float,
    steady_window_steps: int,
) -> dict[str, Any]:
    if duration_s <= 0.0:
        raise ValueError("duration_s must be positive for steady advancement.")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    if steady_dtemp_k_s <= 0.0:
        raise ValueError("steady_dtemp_k_s must be positive.")
    if int(steady_window_steps) < 1:
        raise ValueError("steady_window_steps must be positive.")

    elapsed = 0.0
    steps = 0
    satisfied = 0
    last_rate: float | None = None
    previous = _temperature_vector(build)

    while elapsed < duration_s - 1.0e-12:
        dt = min(dt_s, duration_s - elapsed)
        system.step(dt)
        elapsed += dt
        steps += 1
        current = _temperature_vector(build)
        last_rate = float(np.max(np.abs(current - previous)) / dt)
        previous = current
        if math.isfinite(last_rate) and last_rate <= steady_dtemp_k_s:
            satisfied += 1
        else:
            satisfied = 0
        if satisfied >= int(steady_window_steps):
            return {
                "steady_reached": True,
                "steps": steps,
                "final_time_s": float(system.global_time),
                "max_temperature_rate_k_s": last_rate,
                "steady_dtemp_k_s": float(steady_dtemp_k_s),
                "steady_window_steps": int(steady_window_steps),
            }

    return {
        "steady_reached": False,
        "steps": steps,
        "final_time_s": float(system.global_time),
        "max_temperature_rate_k_s": last_rate,
        "steady_dtemp_k_s": float(steady_dtemp_k_s),
        "steady_window_steps": int(steady_window_steps),
    }


def _temperature_vector(build: dict[str, Any]) -> np.ndarray:
    values = []
    for solid in build["tfe"].solids.values():
        if hasattr(solid, "T"):
            values.append(np.asarray(solid.T, dtype=float).ravel())
    if not values:
        return np.zeros(1, dtype=float)
    return np.concatenate(values)


def _build_summary(
    config: LuchauSingleTFEConfig,
    build: dict[str, Any],
    thermo_results: dict[str, Any] | None,
    transient_result: dict[str, Any] | None = None,
) -> dict[str, Any]:
    tfe = build["tfe"]
    pellet_power = float(np.sum(tfe.solids["pellet"].Q_source))
    emitter_t = np.asarray(tfe.solids["emitter"].T, dtype=float)
    collector_t = np.asarray(tfe.solids["collector"].T, dtype=float)
    coolant_t = np.array([float(vol.T) for vol in build["channel"].volumes], dtype=float)
    return {
        "case": "Luchau_V_single_TFE_fixed_u",
        "config": asdict(config),
        "full_length_m": float(tfe.geom.height),
        "axial_node_count": int(tfe.mesh.n_axial),
        "axial_length_allocation_m": [0.065, 0.377, 0.065],
        "axial_node_allocation": [6, 25, 6],
        "single_tfe_flow_kg_s": float(config.single_tfe_flow_kg_s),
        "coolant_inlet_temperature_k": float(config.inlet_temperature_k),
        "coolant_outlet_temperature_k": float(coolant_t[-1]),
        "heater_profile_sum": float(np.sum(build["heater_profile"])),
        "pellet_power_w": pellet_power,
        "cesium_reservoir_temperature_k": tcs_from_cesium_pressure(float(config.cesium_pressure_torr)),
        "emitter_temperature_mean_k": float(np.mean(emitter_t)),
        "emitter_temperature_max_k": float(np.max(emitter_t)),
        "collector_temperature_mean_k": float(np.mean(collector_t)),
        "collector_temperature_max_k": float(np.max(collector_t)),
        "thermocalc_global_results": thermo_results,
        "transient_result": transient_result,
        "final_time_s": float(build["system"].global_time),
        "finite": bool(np.isfinite([pellet_power, np.mean(emitter_t), np.mean(collector_t), coolant_t[-1]]).all()),
    }


def _write_summary(summary: dict[str, Any], label: str | None) -> Path:
    if label is None:
        label = datetime.now().strftime("%Y%m%d_%H%M%S_luchau_single_tfe")
    output_dir = CASE_DIR / "runs" / label
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=json_default)
    return output_dir


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    summary = run_case(args)
    print(json.dumps(summary, indent=2, ensure_ascii=False, default=json_default))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
