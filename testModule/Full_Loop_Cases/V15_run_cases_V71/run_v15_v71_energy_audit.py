"""Finite-difference global energy audit for the V15 V71 closed loop."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import numpy as np

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopPumpConfig,
    V15PipeFinRadiatorConfig,
    V15_V71_RADIATOR_TUBE_K_LOSS,
)
from testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_tuning import (
    _apply_core_power,
    _set_ode_method,
    build_v15_v71_case_a_system,
    calibrate_pump_head,
    configure_tec,
    final_metrics,
)


def _solid_multiplier(build: Dict[str, Any], name: str) -> float:
    for tfe_name, multiplier in build["ring_multipliers"].items():
        if name.startswith(f"{tfe_name}_"):
            return float(multiplier)
    return 1.0


def _volume_scale_map(build: Dict[str, Any]) -> Dict[int, float]:
    scale = {id(vol): 1.0 for vol in build["system"].fluid_solver.volumes_obj}
    for name, channel in build["fluid_channels"].items():
        for vol in channel.volumes:
            scale[id(vol)] = float(build["ring_multipliers"][name])
    return scale


def _capture_storage(build: Dict[str, Any]) -> Dict[str, Any]:
    solids = {}
    for name, solid in build["system"].solid_components.items():
        solid._update_properties()
        solids[name] = {
            "temperature_k": np.asarray(solid.T, dtype=float).copy(),
            "capacity_j_per_k": np.asarray(solid.thermal_capacitance, dtype=float).copy(),
            "multiplier": _solid_multiplier(build, name),
        }

    scale_map = _volume_scale_map(build)
    fluids = {}
    for index, vol in enumerate(build["system"].fluid_solver.volumes_obj):
        name = str(getattr(vol, "name", f"volume_{index}"))
        fluids[name] = {
            "enthalpy_j_per_kg": float(vol.h),
            "mass_kg": float(vol.rho * vol.vol),
            "multiplier": float(scale_map[id(vol)]),
            "fixed": bool(getattr(vol, "is_pressure_boundary", False)),
        }
    return {"solids": solids, "fluids": fluids}


def _storage_rates(before: Dict[str, Any], after: Dict[str, Any], dt: float) -> Dict[str, float]:
    solid = 0.0
    for name, old in before["solids"].items():
        new = after["solids"][name]
        solid += (
            float(
                np.sum(
                    0.5
                    * (old["capacity_j_per_k"] + new["capacity_j_per_k"])
                    * (new["temperature_k"] - old["temperature_k"])
                )
            )
            * old["multiplier"]
            / dt
        )

    fluid = 0.0
    fluid_exact = 0.0
    for name, old in before["fluids"].items():
        if old["fixed"]:
            continue
        new = after["fluids"][name]
        scale = old["multiplier"]
        fluid += (
            0.5
            * (old["mass_kg"] + new["mass_kg"])
            * (new["enthalpy_j_per_kg"] - old["enthalpy_j_per_kg"])
            * scale
            / dt
        )
        fluid_exact += (
            new["mass_kg"] * new["enthalpy_j_per_kg"]
            - old["mass_kg"] * old["enthalpy_j_per_kg"]
        ) * scale / dt
    return {
        "solid_storage_rate_W": solid,
        "fluid_storage_rate_W": fluid,
        "fluid_exact_delta_mh_rate_W": fluid_exact,
        "combined_storage_rate_W": solid + fluid,
    }


def _tec_heat_removed(build: Dict[str, Any]) -> float:
    total = 0.0
    for name, multiplier in build["ring_multipliers"].items():
        if int(build["tec_ring_multipliers"].get(name, 0)) <= 0:
            continue
        tfe = build["tfes"][name]
        area = np.asarray(tfe.solids["emitter"].boundaries["right"].area, dtype=float)
        emitter_removed = -float(np.sum(np.asarray(tfe.plasma_data.electron_cooling_flux) * area))
        collector_added = float(np.sum(np.asarray(tfe.plasma_data.electron_heating_flux) * area))
        joule = float(np.sum(tfe.electric_data.emitter_joule_heat)) + float(
            np.sum(tfe.electric_data.collector_joule_heat)
        )
        total += (emitter_removed - collector_added - joule) * float(multiplier)
    return total


def _power_snapshot(build: Dict[str, Any], flow_calibration: Dict[str, float]) -> Dict[str, float]:
    metrics = final_metrics(build, flow_calibration)
    return {
        "core_heat_W": float(build["core"].last_total_core_power),
        "external_heat_W": float(metrics["external_heat_input_W"]),
        "radiator_rejection_W": float(metrics["radiator_heat_rejection_W"]),
        "tec_thermal_heat_removed_W": _tec_heat_removed(build),
        "tec_terminal_power_W": float(metrics["tec_power_W"]),
        "coolant_core_enthalpy_gain_W": float(metrics["fluid_core_enthalpy_gain_W"]),
        "total_flow_kg_s": float(metrics["pump"]["total_flow_kg_s"]),
        "pump_total_head_Pa": float(metrics["pump"]["pump_total_head_pa"]),
    }


def _average(integral: Dict[str, float], dt: float) -> Dict[str, float]:
    return {key: value / dt for key, value in integral.items()}


def run(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v15_v71_case_a_system(
        core_config=FullLoopCoreConfig(inlet_temperature_k=723.0, main_tec_enabled=False),
        pump_config=FullLoopPumpConfig(pump_total_head_pa=float(args.pump_head)),
        radiator_config=V15PipeFinRadiatorConfig(
            t_space_k=4.0,
            tube_emissivity=float(args.emissivity),
            fin_emissivity=float(args.emissivity),
            external_heat_enabled=False,
            solid_ode_method="implicit_euler",
            radiator_tube_inlet_k_loss=V15_V71_RADIATOR_TUBE_K_LOSS,
            radiator_tube_outlet_k_loss=V15_V71_RADIATOR_TUBE_K_LOSS,
        ),
    )
    system = build["system"]
    system.load_global_state(str(args.restart))
    _set_ode_method(system)
    _apply_core_power(build, float(args.power))
    configure_tec(build, float(args.voltage), float(args.wire_scale))
    build["core"].post_step(0.0, float(system.global_time))
    build["core"].thermo_calc.calculate(verbose=False)
    flow_calibration = calibrate_pump_head(build, float(args.target_flow), float(args.pump_head))

    start_time = float(system.global_time)
    end_time = start_time + float(args.duration)
    storage_before = _capture_storage(build)
    power_before = _power_snapshot(build, flow_calibration)
    integral = {key: 0.0 for key in power_before}
    step_count = 0

    while system.global_time < end_time - 1.0e-12:
        dt = min(float(args.dt), end_time - float(system.global_time))
        _apply_core_power(build, float(args.power))
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=100)
        _apply_core_power(build, float(args.power))
        power_after = _power_snapshot(build, flow_calibration)
        for key in integral:
            integral[key] += 0.5 * (power_before[key] + power_after[key]) * dt
        power_before = power_after
        step_count += 1

    storage_after = _capture_storage(build)
    storage = _storage_rates(storage_before, storage_after, float(args.duration))
    average = _average(integral, float(args.duration))
    thermal_residual = (
        average["core_heat_W"]
        + average["external_heat_W"]
        - average["radiator_rejection_W"]
        - average["tec_thermal_heat_removed_W"]
        - storage["combined_storage_rate_W"]
    )
    terminal_residual = (
        average["core_heat_W"]
        + average["external_heat_W"]
        - average["radiator_rejection_W"]
        - average["tec_terminal_power_W"]
        - storage["combined_storage_rate_W"]
    )
    result = {
        "case": "V15_V71_global_energy_audit",
        "restart_in": str(args.restart),
        "absolute_start_time_s": start_time,
        "absolute_end_time_s": float(system.global_time),
        "audit_duration_s": float(args.duration),
        "dt_s": float(args.dt),
        "step_count": step_count,
        "average_powers": average,
        "storage_rates": storage,
        "residuals": {
            "thermal_model_residual_W": thermal_residual,
            "thermal_model_relative": thermal_residual / average["core_heat_W"],
            "terminal_power_residual_W": terminal_residual,
            "terminal_power_relative": terminal_residual / average["core_heat_W"],
            "fast_terminal_residual_without_storage_W": (
                average["core_heat_W"]
                + average["external_heat_W"]
                - average["radiator_rejection_W"]
                - average["tec_terminal_power_W"]
            ),
            "tec_thermal_minus_terminal_W": (
                average["tec_thermal_heat_removed_W"] - average["tec_terminal_power_W"]
            ),
        },
        "definitions": {
            "thermal_model_residual": "Qcore + Qexternal - Qradiator - Qtec_applied - dUsolid/dt - dUfluid/dt",
            "terminal_power_residual": "Qcore + Qexternal - Qradiator - Ptec_terminal - dUsolid/dt - dUfluid/dt",
            "solid_storage": "sum(average thermal capacitance * temperature change) / audit duration, with representative TFE multipliers",
            "fluid_storage": "sum(average mass * enthalpy change) / audit duration, with representative channel multipliers",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--power", type=float, default=98000.0)
    parser.add_argument("--emissivity", type=float, default=0.83)
    parser.add_argument("--voltage", type=float, default=28.0)
    parser.add_argument("--wire-scale", type=float, default=0.02)
    parser.add_argument("--pump-head", type=float, default=20003.891)
    parser.add_argument("--target-flow", type=float, default=1.18)
    return parser


if __name__ == "__main__":
    print(json.dumps(run(_parser().parse_args()), indent=2, sort_keys=True))
