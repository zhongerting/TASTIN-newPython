"""Hold +0.50 dollars from the Stage 0 restart until V14 reaches 10 kW."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.run_v14_shield_radiator_startup import (
    _metrics,
    append_v14_system_history,
    build_case,
    capture_v14_history_reference,
)


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_RESTART = CASE_DIR / "phase_0_initial_1800s" / "final_restart.npz"


def _collect(
    build,
    start_time_s: float,
    step_dollars: float,
    commanded_dollars: float,
    control_phase: str,
) -> dict:
    system = build["system"]
    core = build["core"]
    point = core.point_reactor
    feedback = core.compute_reactivity_feedback()
    effective_feedback = float(feedback.total - core.feedback_reference_result.total)
    row = _metrics(build)
    row.update({
        "elapsed_s": float(system.global_time - start_time_s),
        "core_power_w": float(point.total_power),
        "fission_power_w": float(point.fission_power),
        "decay_power_w": float(point.decay_power),
        "step_reactivity_dollars": float(step_dollars),
        "command_reactivity_dollars": float(commanded_dollars),
        "command_reactivity": float(commanded_dollars * point.beta_total),
        "control_phase": str(control_phase),
        "temperature_feedback": effective_feedback,
        "total_reactivity": float(commanded_dollars * point.beta_total + effective_feedback),
        "fuel_mean_k": float(feedback.temperatures.fuel),
        "moderator_mean_k": float(feedback.temperatures.moderator),
        "reflector_mean_k": float(feedback.temperatures.reflector),
        "pump_a_flow_kg_s": float(build["pump_a"].W),
        "pump_b_flow_kg_s": float(build["pump_b"].W),
    })
    if not all(
        isinstance(value, (bool, str)) or math.isfinite(float(value))
        for value in row.values()
    ):
        raise FloatingPointError(f"Non-finite reactivity-startup state: {row}")
    return row


def run(
    *,
    restart_in: Path,
    output_dir: Path,
    step_dollars: float,
    duration_s: float,
    max_dt_s: float,
    source_power_w: float,
    target_power_w: float,
    max_power_w: float,
    record_interval_s: float,
) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if source_power_w <= 0.0 or target_power_w <= source_power_w:
        raise ValueError("Require 0 < source_power_w < target_power_w.")
    if max_power_w <= target_power_w:
        raise ValueError("max_power_w must exceed target_power_w.")
    if record_interval_s <= 0.0:
        raise ValueError("record_interval_s must be positive.")

    build = build_case(initial_temperature_k=300.0, target_flow_kg_s=0.615)
    system = build["system"]
    core = build["core"]
    system.load_global_state(str(restart_in))
    core.enable_tec_coupled = False
    core.initialize_point_reactor(total_power_initial=float(source_power_w))
    build["radiator_thermal_shield"].pre_step(0.0, float(system.global_time))

    output_dir.mkdir(parents=True)
    start_time_s = float(system.global_time)
    end_time_s = start_time_s + float(duration_s)
    next_record_s = start_time_s + float(record_interval_s)
    target_crossing_time_s = None
    commanded_dollars = float(step_dollars)
    stop_reason = "completed"
    power_growth_w_s = 0.0
    history_reference = capture_v14_history_reference(build)
    latest = _collect(
        build,
        start_time_s,
        step_dollars,
        commanded_dollars,
        "fixed_step",
    )
    append_v14_system_history(
        output_dir, latest, build, history_reference, start_time_s,
        external_reactivity_dollars=commanded_dollars,
    )
    print(json.dumps(latest, sort_keys=True), flush=True)

    while system.global_time < end_time_s - 1.0e-10:
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=float(max_dt_s),
            safety_factor=0.8,
            respect_fluid_cfl=False,
        )
        dt = min(
            float(dt), end_time_s - float(system.global_time),
            next_record_s - float(system.global_time),
        )
        power_before_step_w = float(core.point_reactor.total_power)
        if power_growth_w_s > 0.0:
            dt = min(dt, (float(target_power_w) - power_before_step_w) / power_growth_w_s)
        system.step(
            dt,
            inner_iter=1,
            fail_on_fluid_nonconvergence=True,
            fluid_max_iter=300,
            reactivity_control=float(commanded_dollars * core.point_reactor.beta_total),
        )
        latest = _collect(
            build,
            start_time_s,
            step_dollars,
            commanded_dollars,
            "fixed_step",
        )
        power_w = float(latest["core_power_w"])
        power_growth_w_s = max(0.0, (power_w - power_before_step_w) / float(dt))
        if target_crossing_time_s is None and power_w >= target_power_w:
            target_crossing_time_s = float(system.global_time)
            stop_reason = "target_power"
        if power_w >= max_power_w:
            stop_reason = "maximum_power"

        if (
            system.global_time >= next_record_s - 1.0e-10
            or stop_reason != "completed"
            or system.global_time >= end_time_s - 1.0e-10
        ):
            append_v14_system_history(
                output_dir, latest, build, history_reference, start_time_s,
                external_reactivity_dollars=commanded_dollars,
            )
            print(json.dumps(latest, sort_keys=True), flush=True)
            while next_record_s <= system.global_time + 1.0e-10:
                next_record_s += float(record_interval_s)
        if stop_reason != "completed":
            break

    restart_out = output_dir / "final_restart.npz"
    system.save_global_state(str(restart_out))
    summary = {
        "restart_in": str(restart_in),
        "restart_out": str(restart_out),
        "start_time_s": start_time_s,
        "end_time_s": float(system.global_time),
        "elapsed_s": float(system.global_time - start_time_s),
        "source_power_w": float(source_power_w),
        "target_power_w": float(target_power_w),
        "max_power_w": float(max_power_w),
        "record_interval_s": float(record_interval_s),
        "step_reactivity_dollars": float(step_dollars),
        "target_crossing_time_s": target_crossing_time_s,
        "stop_reason": stop_reason,
        "latest": latest,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--step-dollars", type=float, default=0.50)
    parser.add_argument("--duration", type=float, default=300.0)
    parser.add_argument("--max-dt", type=float, default=0.2)
    parser.add_argument("--source-power", type=float, default=1.0)
    parser.add_argument("--target-power", type=float, default=10000.0)
    parser.add_argument("--max-power", type=float, default=12000.0)
    parser.add_argument("--record-interval", type=float, default=1.0)
    args = parser.parse_args()
    run(
        restart_in=args.restart_in,
        output_dir=args.output_dir,
        step_dollars=args.step_dollars,
        duration_s=args.duration,
        max_dt_s=args.max_dt,
        source_power_w=args.source_power,
        target_power_w=args.target_power,
        max_power_w=args.max_power,
        record_interval_s=args.record_interval,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
