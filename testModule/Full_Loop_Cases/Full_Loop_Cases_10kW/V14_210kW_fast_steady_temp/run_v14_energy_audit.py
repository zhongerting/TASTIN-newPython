"""Finite-difference global energy audit for the V14 210 kW closed loop."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import numpy as np

from testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_energy_audit import (
    _average,
    _storage_rates,
    _tec_heat_removed,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    _apply_fixed_core_power,
    _ring_rejection,
    _ring_wall_rejection,
    build_debug_case,
    collect_metrics,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    ReactivityControlRunConfig,
    load_baseline_debug_config,
)


def _solid_scale_map(build: Dict[str, Any]) -> Dict[int, float]:
    scale = {id(solid): 1.0 for solid in build["system"].solid_components.values()}
    for name, tfe in build["tfes"].items():
        for solid in tfe.solids.values():
            scale[id(solid)] = float(build["ring_multipliers"][name])
    for ring_hp in build["ring_hps"]:
        for _, hp, multiplier in ring_hp._iter_present_hp_units_with_multiplier():
            for solid in hp.get_solids():
                scale[id(solid)] = float(multiplier)
    return scale


def _fluid_scale_map(build: Dict[str, Any]) -> Dict[int, float]:
    scale = {id(vol): 1.0 for vol in build["system"].fluid_solver.volumes_obj}
    for name, channel in build["fluid_channels"].items():
        for vol in channel.volumes:
            scale[id(vol)] = float(build["ring_multipliers"][name])
    return scale


def _capture_storage(build: Dict[str, Any]) -> Dict[str, Any]:
    solid_scale = _solid_scale_map(build)
    solids = {}
    for name, solid in build["system"].solid_components.items():
        solid._update_properties()
        solids[name] = {
            "temperature_k": np.asarray(solid.T, dtype=float).copy(),
            "capacity_j_per_k": np.asarray(solid.thermal_capacitance, dtype=float).copy(),
            "multiplier": solid_scale[id(solid)],
        }

    fluid_scale = _fluid_scale_map(build)
    fluids = {}
    for index, vol in enumerate(build["system"].fluid_solver.volumes_obj):
        fluids[str(getattr(vol, "name", f"volume_{index}"))] = {
            "enthalpy_j_per_kg": float(vol.h),
            "mass_kg": float(vol.rho * vol.vol),
            "multiplier": fluid_scale[id(vol)],
            "fixed": bool(getattr(vol, "is_pressure_boundary", False)),
        }
    return {"solids": solids, "fluids": fluids}


def _power_snapshot(build: Dict[str, Any]) -> Dict[str, float]:
    metrics = collect_metrics(build, stage_index=1, dt_s=0.0)
    net = build["system"].fluid_solver
    pump = next(j for j in net.junctions_obj if getattr(j, "name", "") == "J_PumpA")
    coolant_gain = float(pump.W) * (
        float(build["core_outlet_connector"].h)
        - float(build["core_inlet_connector"].h)
    )
    return {
        "core_heat_W": float(build["core"].last_total_core_power),
        "radiator_rejection_W": (
            _ring_rejection(build["ring_hps"])
            + _ring_wall_rejection(build)
        ),
        "tec_thermal_heat_removed_W": _tec_heat_removed(build),
        "tec_terminal_power_W": float(metrics["tec_main_electric_power_W"]),
        "coolant_core_enthalpy_gain_W": coolant_gain,
        "total_flow_kg_s": float(pump.W),
    }


def run(
        restart: Path, output: Path, duration: float, dt: float,
        warmup: float = 1.0) -> Dict[str, Any]:
    runtime = ReactivityControlRunConfig(
        restart_in=restart,
        output_dir=output.parent,
        duration_s=duration,
        dt_s=dt,
    )
    debug, source = load_baseline_debug_config(runtime)
    build = build_debug_case(debug)
    system = build["system"]
    core = build["core"]
    if bool(source.get("external_heat_enabled", False)):
        raise ValueError("audit expects external heat to be disabled")

    current_time = float(system.global_time)
    for component in system.components:
        component.post_step(0.0, current_time)
    core.thermo_calc.calculate(verbose=False)
    core.pre_step(dt, current_time)

    warmup_end = current_time + warmup
    while float(system.global_time) < warmup_end - 1.0e-12:
        step_dt = min(dt, warmup_end - float(system.global_time))
        _apply_fixed_core_power(build, debug.power_w)
        system.step(
            step_dt,
            inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=True,
            fluid_max_iter=int(debug.fluid_max_iter),
        )
        _apply_fixed_core_power(build, debug.power_w)

    before = _capture_storage(build)
    power_before = _power_snapshot(build)
    integral = {key: 0.0 for key in power_before}
    audit_start_time = float(system.global_time)
    end_time = audit_start_time + duration
    steps = 0
    while float(system.global_time) < end_time - 1.0e-12:
        step_dt = min(dt, end_time - float(system.global_time))
        _apply_fixed_core_power(build, debug.power_w)
        system.step(
            step_dt,
            inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=True,
            fluid_max_iter=int(debug.fluid_max_iter),
        )
        _apply_fixed_core_power(build, debug.power_w)
        power_after = _power_snapshot(build)
        for key in integral:
            integral[key] += 0.5 * (power_before[key] + power_after[key]) * step_dt
        power_before = power_after
        steps += 1

    storage = _storage_rates(before, _capture_storage(build), duration)
    average = _average(integral, duration)
    thermal_residual = (
        average["core_heat_W"]
        - average["radiator_rejection_W"]
        - average["tec_thermal_heat_removed_W"]
        - storage["combined_storage_rate_W"]
    )
    terminal_residual = (
        average["core_heat_W"]
        - average["radiator_rejection_W"]
        - average["tec_terminal_power_W"]
        - storage["combined_storage_rate_W"]
    )
    result = {
        "case": "V14_210kW_global_energy_audit",
        "restart_in": str(restart),
        "restart_time_s": current_time,
        "warmup_duration_s": warmup,
        "absolute_start_time_s": audit_start_time,
        "absolute_end_time_s": float(system.global_time),
        "audit_duration_s": duration,
        "dt_s": dt,
        "step_count": steps,
        "average_powers": average,
        "storage_rates": storage,
        "residuals": {
            "thermal_model_residual_W": thermal_residual,
            "thermal_model_relative": thermal_residual / average["core_heat_W"],
            "terminal_power_residual_W": terminal_residual,
            "terminal_power_relative": terminal_residual / average["core_heat_W"],
            "tec_thermal_minus_terminal_W": (
                average["tec_thermal_heat_removed_W"]
                - average["tec_terminal_power_W"]
            ),
        },
        "definitions": {
            "thermal_model_residual": "Qcore - Qradiator - Qtec_applied - dUsolid/dt - dUfluid/dt",
            "terminal_power_residual": "Qcore - Qradiator - Ptec_terminal - dUsolid/dt - dUfluid/dt",
            "coolant_core_enthalpy_gain": "internal core-to-coolant transfer; excluded from the global residual",
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--warmup", type=float, default=1.0)
    args = parser.parse_args()
    print(json.dumps(
        run(args.restart, args.output, args.duration, args.dt, args.warmup),
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
