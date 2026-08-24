"""Ramp prescribed V14 thermal power from 10 kW to 70 kW at 600 W/s."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.run_v14_shield_radiator_startup import (
    _metrics,
    append_v14_system_history,
    build_case,
    capture_v14_history_reference,
)


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_RESTART = CASE_DIR / "stage_1_fixed_0p50_to_10kw" / "final_restart.npz"


def _apply_power(core, power_w: float) -> None:
    core.update_neutronic_power(
        p_total=float(power_w), p_fiss=float(power_w), p_decay=0.0, alpha=1.0,
    )


def _update_threshold_events(build, state: dict) -> None:
    system = build["system"]
    coolant_min_k = float(np.min(system.fluid_solver.T_vec))
    outlet_k = float(build["core_outlet_connector"].T)
    if not state["shield_jettisoned"] and coolant_min_k >= state["shield_threshold_k"]:
        build["radiator_thermal_shield"].set_active(False)
        build["radiator_thermal_shield"].pre_step(0.0, float(system.global_time))
        state["shield_jettisoned"] = True
        state["shield_jettison_time_s"] = float(system.global_time)
    if not state["flow_increased"] and outlet_k >= state["outlet_threshold_k"]:
        build["pump_a"].set_flow_rate(state["high_flow_kg_s"])
        state["flow_increased"] = True
        state["flow_increase_time_s"] = float(system.global_time)


def _collect(build, start_time_s: float, prescribed_power_w: float, state: dict) -> dict:
    row = _metrics(build)
    row.update({
        "elapsed_s": float(build["system"].global_time - start_time_s),
        "power_control_mode": "prescribed_linear_ramp",
        "prescribed_power_w": float(prescribed_power_w),
        "power_ramp_rate_w_s": float(state["ramp_rate_w_s"]),
        "coolant_min_k": float(np.min(build["system"].fluid_solver.T_vec)),
        "core_outlet_coolant_k": float(build["core_outlet_connector"].T),
        "target_flow_kg_s": float(build["pump_a"].target_W),
        "pump_a_flow_kg_s": float(build["pump_a"].W),
        "pump_b_flow_kg_s": float(build["pump_b"].W),
        "shield_jettisoned": bool(state["shield_jettisoned"]),
        "flow_increased": bool(state["flow_increased"]),
    })
    if not all(isinstance(value, (bool, str)) or math.isfinite(float(value)) for value in row.values()):
        raise FloatingPointError(f"Non-finite power-ramp state: {row}")
    return row


def run(
    *,
    restart_in: Path,
    output_dir: Path,
    target_power_w: float = 70_000.0,
    ramp_rate_w_s: float = 600.0,
    max_duration_s: float = 200.0,
    max_dt_s: float = 0.2,
    record_interval_s: float = 1.0,
    checkpoint_interval_s: float = 50.0,
    low_flow_kg_s: float = 0.615,
    high_flow_kg_s: float = 1.23,
    shield_threshold_k: float = 373.0,
    outlet_threshold_k: float = 500.0,
) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if ramp_rate_w_s <= 0.0 or record_interval_s <= 0.0 or checkpoint_interval_s <= 0.0:
        raise ValueError("Ramp and recording intervals must be positive.")

    build = build_case(initial_temperature_k=300.0, target_flow_kg_s=low_flow_kg_s)
    system = build["system"]
    core = build["core"]
    system.load_global_state(str(restart_in))
    start_power_w = float(core.point_reactor.total_power if core.has_point_reactor else core.last_total_core_power)
    core.point_reactor = None
    core.enable_tec_coupled = False
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    build["pump_a"].set_flow_rate(float(low_flow_kg_s))
    build["radiator_thermal_shield"].pre_step(0.0, float(system.global_time))
    _apply_power(core, start_power_w)

    output_dir.mkdir(parents=True)
    start_time_s = float(system.global_time)
    end_time_s = start_time_s + float(max_duration_s)
    next_record_s = start_time_s + float(record_interval_s)
    next_checkpoint_s = start_time_s + float(checkpoint_interval_s)
    state = {
        "ramp_rate_w_s": float(ramp_rate_w_s),
        "shield_threshold_k": float(shield_threshold_k),
        "outlet_threshold_k": float(outlet_threshold_k),
        "high_flow_kg_s": float(high_flow_kg_s),
        "shield_jettisoned": not bool(build["radiator_thermal_shield"].last_active),
        "shield_jettison_time_s": None,
        "flow_increased": False,
        "flow_increase_time_s": None,
    }
    _update_threshold_events(build, state)
    history_reference = capture_v14_history_reference(build)
    prescribed_power_w = start_power_w
    latest = _collect(build, start_time_s, prescribed_power_w, state)
    append_v14_system_history(output_dir, latest, build, history_reference, start_time_s)
    print(json.dumps(latest, sort_keys=True), flush=True)

    stop_reason = "completed"
    checkpoint_paths = []
    while system.global_time < end_time_s - 1.0e-10 and prescribed_power_w < target_power_w - 1.0e-9:
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4, max_dt=float(max_dt_s), safety_factor=0.8,
            respect_fluid_cfl=False,
        )
        dt = min(
            float(dt), end_time_s - float(system.global_time),
            next_record_s - float(system.global_time),
            next_checkpoint_s - float(system.global_time),
            (float(target_power_w) - prescribed_power_w) / float(ramp_rate_w_s),
        )
        prescribed_power_w = min(
            float(target_power_w),
            start_power_w + float(ramp_rate_w_s) * (float(system.global_time) + dt - start_time_s),
        )
        _apply_power(core, prescribed_power_w)
        system.step(
            dt, inner_iter=1, fail_on_fluid_nonconvergence=True, fluid_max_iter=300,
        )
        _update_threshold_events(build, state)

        final_step = prescribed_power_w >= target_power_w - 1.0e-9
        if system.global_time >= next_record_s - 1.0e-10 or final_step:
            latest = _collect(build, start_time_s, prescribed_power_w, state)
            append_v14_system_history(output_dir, latest, build, history_reference, start_time_s)
            print(json.dumps(latest, sort_keys=True), flush=True)
            while next_record_s <= system.global_time + 1.0e-10:
                next_record_s += float(record_interval_s)
        if system.global_time >= next_checkpoint_s - 1.0e-10 and not final_step:
            checkpoint = output_dir / f"checkpoint_tplus_{system.global_time - start_time_s:07.1f}s.npz"
            system.save_global_state(str(checkpoint))
            checkpoint_paths.append(str(checkpoint))
            while next_checkpoint_s <= system.global_time + 1.0e-10:
                next_checkpoint_s += float(checkpoint_interval_s)

    if prescribed_power_w >= target_power_w - 1.0e-9:
        stop_reason = "target_power"
    elif system.global_time >= end_time_s - 1.0e-10:
        stop_reason = "maximum_duration"

    restart_out = output_dir / "final_restart.npz"
    system.save_global_state(str(restart_out))
    summary = {
        "restart_in": str(restart_in),
        "restart_out": str(restart_out),
        "checkpoint_paths": checkpoint_paths,
        "start_time_s": start_time_s,
        "end_time_s": float(system.global_time),
        "elapsed_s": float(system.global_time - start_time_s),
        "start_power_w": start_power_w,
        "target_power_w": float(target_power_w),
        "ramp_rate_w_s": float(ramp_rate_w_s),
        "record_interval_s": float(record_interval_s),
        "checkpoint_interval_s": float(checkpoint_interval_s),
        "low_flow_kg_s": float(low_flow_kg_s),
        "high_flow_kg_s": float(high_flow_kg_s),
        "shield_threshold_k": float(shield_threshold_k),
        "outlet_threshold_k": float(outlet_threshold_k),
        "shield_jettison_time_s": state["shield_jettison_time_s"],
        "flow_increase_time_s": state["flow_increase_time_s"],
        "stop_reason": stop_reason,
        "latest": latest,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-power", type=float, default=70_000.0)
    parser.add_argument("--ramp-rate", type=float, default=600.0)
    parser.add_argument("--max-duration", type=float, default=200.0)
    parser.add_argument("--max-dt", type=float, default=0.2)
    parser.add_argument("--record-interval", type=float, default=1.0)
    parser.add_argument("--checkpoint-interval", type=float, default=50.0)
    parser.add_argument("--low-flow", type=float, default=0.615)
    parser.add_argument("--high-flow", type=float, default=1.23)
    parser.add_argument("--shield-threshold", type=float, default=373.0)
    parser.add_argument("--outlet-threshold", type=float, default=500.0)
    args = parser.parse_args()
    run(
        restart_in=args.restart_in, output_dir=args.output_dir,
        target_power_w=args.target_power, ramp_rate_w_s=args.ramp_rate,
        max_duration_s=args.max_duration, max_dt_s=args.max_dt,
        record_interval_s=args.record_interval,
        checkpoint_interval_s=args.checkpoint_interval,
        low_flow_kg_s=args.low_flow, high_flow_kg_s=args.high_flow,
        shield_threshold_k=args.shield_threshold,
        outlet_threshold_k=args.outlet_threshold,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
