"""Run the V14 Stage 0 shield/radiator startup with restart support."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig
from .v14_case import build_v14_case_a_system
from .v14_heatpipe_radiator import V14HeatPipeRadiatorConfig
from .V14_210kW_fixed_power_LOCA_1.run_v14_210kw_fixed_power_loca_1 import (
    append_postprocessing_histories,
    build_snapshot_payload,
)


def _set_uniform_temperature(build, temperature_k: float) -> None:
    system = build["system"]
    net = system.fluid_solver
    for volume in net.volumes_obj:
        volume.T = float(temperature_k)
        volume.h = float(volume.material.enthalpy(volume.T, volume.P))
        volume.update_properties(volume.material)
    net._initialize_state_from_objects()
    net._update_fluid_properties()
    net._sync_vectors_to_objects()
    for solid in system.solid_components.values():
        solid.T[...] = float(temperature_k)
        if hasattr(solid, "dTdt"):
            solid.dTdt[...] = 0.0
        solid._update_properties()
        solid._update_boundaries_state(current_time=float(system.global_time))
        solid.set_ode_method("implicit_euler")
    for unit in build["radiator_units"]:
        unit.last_fin_temperature[...] = float(temperature_k)
        unit.last_fin_effective_temperature_distribution[...] = float(temperature_k)
        unit._has_valid_fin_temperature = True


def build_case(
    initial_temperature_k: float = 300.0,
    target_flow_kg_s: float = 0.615,
    orbit_phase_origin_s: float = 0.0,
):
    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(
            tec_gap_h_eq_w_m2_k=5678.0,
            tec_gap_gas="Helium",
            inlet_temperature_k=float(initial_temperature_k),
            main_tec_enabled=False,
        ),
        flow_config=FullLoopFlowConfig(total_flow_kg_s=float(target_flow_kg_s)),
        pump_config=FullLoopPumpConfig(
            pump_total_head_pa=1.0,
            pump_flow_control=True,
            target_flow_kg_s=float(target_flow_kg_s),
        ),
        radiator_config=V14HeatPipeRadiatorConfig(
            t_space_k=200.0,
            hp_initial_temp_k=float(initial_temperature_k),
            external_heat_enabled=True,
            external_heat_period_s=5668.14,
            external_heat_time_origin_s=float(orbit_phase_origin_s),
            external_heat_absorption_efficiency=0.992,
            thermal_shield_enabled=True,
            thermal_shield_initially_active=True,
        ),
    )
    build["core"].point_reactor = None
    build["core"].enable_tec_coupled = False
    build["core"].update_neutronic_power(p_total=0.0, p_fiss=0.0, p_decay=0.0, alpha=1.0)
    return build


def _metrics(build) -> dict:
    system = build["system"]
    shield = build["radiator_thermal_shield"]
    solids = [np.asarray(solid.T, dtype=float) for solid in system.solid_components.values()]
    direct_heat = sum(
        ring.get_total_external_heat_absorption_scaled(system.global_time)
        for ring in build["ring_hps"]
    )
    result = {
        "time_s": float(system.global_time),
        "shield_active": bool(shield.last_active),
        "orbit_phase_s": float(
            (system.global_time - shield.external_heat_source.time_origin_s) % 5668.14
        ),
        "fluid_min_k": float(np.min(system.fluid_solver.T_vec)),
        "fluid_max_k": float(np.max(system.fluid_solver.T_vec)),
        "solid_min_k": float(min(np.min(value) for value in solids)),
        "solid_max_k": float(max(np.max(value) for value in solids)),
        "shield_inner_mean_k": float(shield.last_inner_temperature_mean_k),
        "shield_outer_mean_k": float(shield.last_outer_temperature_mean_k),
        "radiator_background_mean_k": float(shield.last_effective_background_mean_k),
        "shield_qsss_max_w_m2": float(np.max(shield.qsss_w_m2)),
        "direct_radiator_external_heat_w": float(direct_heat),
        "shield_converged": bool(shield.last_shield2_converged),
    }
    if not all(math.isfinite(float(value)) for value in result.values()):
        raise FloatingPointError(f"Non-finite startup state: {result}")
    return result


def capture_v14_history_reference(build) -> dict:
    net = build["system"].fluid_solver
    return {
        "feedback_reference_total": float(build["core"].compute_reactivity_feedback().total),
        "reference_fluid": {
            "T": np.asarray(net.T_vec, dtype=float).copy(),
            "P": np.asarray(net.P_vec, dtype=float).copy(),
            "h": np.asarray(net.h_vec, dtype=float).copy(),
            "W": np.asarray(net.W_vec, dtype=float).copy(),
        },
        "tec_open_circuit_active": False,
        "tec_open_circuit_time_s": float("nan"),
    }


def append_v14_system_history(
    output_dir: Path,
    summary_row: dict,
    build,
    history_reference: dict,
    start_time_s: float,
    external_reactivity_dollars: float = 0.0,
) -> None:
    history_path = output_dir / "history.csv"
    with history_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_row))
        if handle.tell() == 0:
            writer.writeheader()
        writer.writerow(summary_row)
    core = build["core"]
    beta = core.point_reactor.beta_total if core.has_point_reactor else 0.0
    payload = build_snapshot_payload(
        build,
        history_reference,
        start_time=float(start_time_s),
        coolant_present=True,
        hydraulic_solve_enabled=True,
        external_reactivity=float(external_reactivity_dollars) * float(beta),
        external_reactivity_dollars=float(external_reactivity_dollars),
    )
    if not core.has_point_reactor:
        feedback = core.compute_reactivity_feedback()
        effective = float(feedback.total) - float(history_reference["feedback_reference_total"])
        payload.update({
            "fission_power_W": np.asarray([core.last_total_core_power]),
            "decay_power_W": np.asarray([0.0]),
            "external_reactivity": np.asarray([float(external_reactivity_dollars) * float(beta)]),
            "external_reactivity_dollars": np.asarray([float(external_reactivity_dollars)]),
            "effective_temperature_feedback": np.asarray([effective]),
            "total_reactivity": np.asarray([effective + core.get_control_drum_reactivity()]),
        })
    append_postprocessing_histories(output_dir, payload)


def run(
    output_dir: Path,
    duration_s: float,
    max_dt_s: float,
    restart_in: Path | None = None,
    withdraw_shield: bool = False,
    record_interval_s: float = 10.0,
    initial_temperature_k: float = 300.0,
    target_flow_kg_s: float = 0.615,
) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    output_dir.mkdir(parents=True)
    build = build_case(
        initial_temperature_k=initial_temperature_k,
        target_flow_kg_s=target_flow_kg_s,
    )
    system = build["system"]
    shield = build["radiator_thermal_shield"]
    if restart_in is None:
        _set_uniform_temperature(build, initial_temperature_k)
        system.initialize_system(dt_init=min(0.01, max_dt_s), tol=1.0e-5, max_iter=1000)
    else:
        system.load_global_state(str(restart_in))
        for solid in system.solid_components.values():
            solid.set_ode_method("implicit_euler")
    if withdraw_shield:
        shield.set_active(False)
    shield.pre_step(0.0, float(system.global_time))

    start_time = float(system.global_time)
    end_time = start_time + float(duration_s)
    next_record = start_time + float(record_interval_s)
    history_reference = capture_v14_history_reference(build)
    latest = _metrics(build)
    append_v14_system_history(
        output_dir, latest, build, history_reference, start_time,
    )
    print(json.dumps(latest, sort_keys=True), flush=True)
    while system.global_time < end_time - 1.0e-10:
        build["core"].update_neutronic_power(p_total=0.0, p_fiss=0.0, p_decay=0.0, alpha=1.0)
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=float(max_dt_s),
            safety_factor=0.8,
            respect_fluid_cfl=False,
        )
        dt = min(
            float(dt), float(max_dt_s), end_time - system.global_time,
            next_record - system.global_time,
        )
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=300)
        if system.global_time >= next_record - 1.0e-10:
            latest = _metrics(build)
            append_v14_system_history(
                output_dir, latest, build, history_reference, start_time,
            )
            print(json.dumps(latest, sort_keys=True), flush=True)
            next_record += float(record_interval_s)

    latest = _metrics(build)
    restart_path = output_dir / "final_restart.npz"
    system.save_global_state(str(restart_path))
    summary = {
        "initial_temperature_k": float(initial_temperature_k),
        "target_flow_kg_s": float(target_flow_kg_s),
        "duration_s": float(duration_s),
        "record_interval_s": float(record_interval_s),
        "orbit_period_s": 5668.14,
        "restart_in": None if restart_in is None else str(restart_in),
        "restart_out": str(restart_path),
        "withdraw_shield": bool(withdraw_shield),
        "latest": latest,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--duration", type=float, default=1800.0)
    parser.add_argument("--max-dt", type=float, default=0.1)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--restart-in", type=Path)
    parser.add_argument("--withdraw-shield", action="store_true")
    parser.add_argument("--initial-temperature", type=float, default=300.0)
    parser.add_argument("--target-flow", type=float, default=0.615)
    args = parser.parse_args()
    run(
        args.output_dir,
        args.duration,
        args.max_dt,
        args.restart_in,
        args.withdraw_shield,
        args.record_interval,
        args.initial_temperature,
        args.target_flow,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
