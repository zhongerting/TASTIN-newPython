"""Whole-core TEC open-circuit accident runner shared by two V14 cases."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    _apply_fixed_core_power, build_debug_case, collect_metrics,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1.run_v14_210kw_fixed_power_loca_1 import (
    _neutronics_metrics, append_postprocessing_histories, build_snapshot_payload,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.run_v14_helium_depressurization import (
    collect_temperature_peaks, find_nonfinite_model_state,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    ReactivityControlRunConfig, load_baseline_debug_config, prepare_reactivity_control,
)


ORBIT_PERIOD_S = 5668.144369
HALF_ORBIT_S = ORBIT_PERIOD_S / 2.0
FIXED_POWER_W = 210000.0
SCRAM_DOLLARS = -2.0
DEFAULT_RESTART = (
    Path(__file__).resolve().parent
    / "V14_210kW_fixed_power_external_heat_2orbits"
    / "runs" / "two_orbits_from13864_20260720"
    / "checkpoint_t019865s.npz"
)
DEFAULT_LIMITS_K = {
    "channel_wall": 1058.0, "pellet": 2700.0, "collector": 1023.0,
    "moderator": 930.0, "reflector": 1000.0,
}
SUMMARY_FIELDS = [
    "time_s", "accident_elapsed_s", "dt_s", "case_mode",
    "external_heat_enabled", "external_heat_period_s", "external_heat_phase_s",
    "tec_open_circuit_active", "fixed_power_control_active",
    "point_kinetics_enabled", "scram_active", "scram_time_absolute_s",
    "scram_elapsed_s", "scram_trigger_component", "scram_trigger_representative",
    "scram_trigger_axial_position_m", "scram_trigger_actual_K",
    "scram_trigger_limit_K", "core_power_W", "fission_power_W", "decay_power_W",
    "external_reactivity", "external_reactivity_dollars",
    "effective_temperature_feedback", "total_reactivity",
    "feedback_fuel", "feedback_electrode", "feedback_moderator",
    "feedback_reflector", "feedback_total_absolute",
    "channel_wall_max_T_K", "fuel_max_T_K", "collector_max_T_K",
    "moderator_max_T_K", "reflector_max_T_K", "coolant_max_T_K",
    "core_inlet_T_K", "core_outlet_T_K", "radiator_heat_rejection_W",
    "total_external_heat_rejection_W", "pump_a_flow_kg_s", "pump_b_flow_kg_s",
    "fluid_converged", "stop_reason",
]


@dataclass(frozen=True)
class TecOpenCircuitRunConfig:
    mode: str
    case_name: str
    output_dir: Path
    restart_in: Path = DEFAULT_RESTART
    minimum_duration_s: float = ORBIT_PERIOD_S
    post_scram_duration_s: float = HALF_ORBIT_S
    dt_s: float = 0.05
    checkpoint_interval_s: float = 50.0
    wall_limit_k: float = 1058.0
    fuel_limit_k: float = 2700.0
    collector_limit_k: float = 1023.0
    moderator_limit_k: float = 930.0
    reflector_limit_k: float = 1000.0


def record_interval_s(phase_elapsed_s: float) -> float:
    phase_elapsed_s = round(float(phase_elapsed_s), 6)
    if phase_elapsed_s < 20.0:
        return 0.1
    if phase_elapsed_s < 100.0:
        return 1.0
    return 10.0


def required_end_elapsed_s(
    minimum_duration_s: float,
    scram_elapsed_s: Optional[float],
    post_scram_duration_s: float,
) -> float:
    if scram_elapsed_s is None:
        return float(minimum_duration_s)
    return max(float(minimum_duration_s), float(scram_elapsed_s) + float(post_scram_duration_s))


def disable_all_tec(build: Dict[str, Any]) -> None:
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    build["core"].enable_tec_coupled = False
    build["system"]._refresh_solid_boundary_cache(
        update_flux=True, current_time=float(build["system"].global_time)
    )


def _validate(config: TecOpenCircuitRunConfig) -> None:
    if config.mode not in {"fixed_power", "reactive_feedback"}:
        raise ValueError("mode must be fixed_power or reactive_feedback")
    if not Path(config.restart_in).is_file():
        raise FileNotFoundError(f"restart not found: {config.restart_in}")
    for name in (
        "minimum_duration_s", "post_scram_duration_s", "dt_s",
        "checkpoint_interval_s", "wall_limit_k", "fuel_limit_k",
        "collector_limit_k", "moderator_limit_k", "reflector_limit_k",
    ):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")


def _limits(config: TecOpenCircuitRunConfig) -> Dict[str, float]:
    return {
        "channel_wall": float(config.wall_limit_k),
        "pellet": float(config.fuel_limit_k),
        "collector": float(config.collector_limit_k),
        "moderator": float(config.moderator_limit_k),
        "reflector": float(config.reflector_limit_k),
    }


def _temperature_metrics(peaks: Sequence[Dict[str, Any]]) -> Dict[str, float]:
    result = {}
    for component, field in (
        ("channel_wall", "channel_wall_max_T_K"),
        ("pellet", "fuel_max_T_K"),
        ("collector", "collector_max_T_K"),
        ("moderator", "moderator_max_T_K"),
        ("reflector", "reflector_max_T_K"),
    ):
        result[field] = max(
            float(peak["actual_k"]) for peak in peaks
            if peak.get("source_component", peak["component"]) == component
        )
    return result


def _temperature_trip(peaks: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    violations = [
        peak for peak in peaks
        if peak["component"] != "nonfinite_temperature"
        and float(peak["actual_k"]) >= float(peak["limit_k"])
    ]
    return max(
        violations,
        key=lambda item: float(item["actual_k"]) / float(item["limit_k"]),
        default=None,
    )


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )


def _append_summary(path: Path, row: Dict[str, Any]) -> None:
    new_file = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        writer.writerow(row)


def _slope_summary(rows: Sequence[Dict[str, Any]], window_s: float = 300.0) -> Dict[str, Any]:
    end = float(rows[-1]["time_s"])
    selected = [row for row in rows if float(row["time_s"]) >= end - window_s]
    result: Dict[str, Any] = {"window_s": window_s, "sample_count": len(selected)}
    if len(selected) < 2:
        return result
    x = np.asarray([float(row["time_s"]) for row in selected])
    x -= x[0]
    for field in (
        "core_power_W", "channel_wall_max_T_K", "fuel_max_T_K",
        "collector_max_T_K", "moderator_max_T_K", "reflector_max_T_K",
        "coolant_max_T_K", "core_inlet_T_K", "core_outlet_T_K",
    ):
        y = np.asarray([float(row[field]) for row in selected])
        result[f"{field}_slope_per_s"] = float(np.polyfit(x, y, 1)[0])
    return result


def run_accident(config: TecOpenCircuitRunConfig) -> Dict[str, Any]:
    _validate(config)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / "history.csv"
    if history_path.exists():
        raise FileExistsError(f"output history already exists: {history_path}")

    source_path = Path(config.restart_in).parent / "run_config.json"
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if not bool(source.get("external_heat_enabled", False)):
        raise ValueError("source restart must retain orbital external heat")
    if bool(source.get("point_kinetics_enabled", False)):
        raise ValueError("source restart must be fixed-power without point kinetics")
    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in), output_dir=out_dir,
        duration_s=float(config.minimum_duration_s), dt_s=float(config.dt_s),
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=None, external_heat_enabled=True,
        external_heat_period_s=float(source["external_heat_period_s"]),
        external_heat_time_origin_s=float(source["external_heat_time_origin_s"]),
    )
    debug, source = load_baseline_debug_config(runtime)
    fixed_power_active = config.mode == "fixed_power"
    build = build_debug_case(debug, apply_fixed_power=fixed_power_active)
    system, core = build["system"], build["core"]
    if not bool(build.get("external_heat_enabled", False)):
        raise RuntimeError("external heat boundary was not attached")

    handoff_type = "fixed_power_control"
    if config.mode == "reactive_feedback":
        handoff_type = prepare_reactivity_control(
            core, source_point_kinetics_enabled=False, expected_power_w=FIXED_POWER_W,
        )

    start_time = float(system.global_time)
    disable_all_tec(build)
    feedback = core.compute_reactivity_feedback()
    accident = {
        "reference_fluid": {
            "T": np.asarray(system.fluid_solver.T_vec).copy(),
            "P": np.asarray(system.fluid_solver.P_vec).copy(),
            "h": np.asarray(system.fluid_solver.h_vec).copy(),
            "W": np.asarray(system.fluid_solver.W_vec).copy(),
        },
        "feedback_reference_total": float(feedback.total),
        "tec_open_circuit_active": True,
        "tec_open_circuit_time_s": 0.0,
        "tec_open_circuit_trigger_current_A": float(source["tec_current_guess_a"]),
    }
    limits = _limits(config)
    state: Dict[str, Any] = {
        "fixed_power_active": fixed_power_active, "scram_active": False,
        "scram_time_absolute_s": float("nan"), "scram_elapsed_s": float("nan"),
        "scram_trigger": None, "handoff_type": handoff_type,
    }
    event = {
        "case": config.case_name, "mode": config.mode,
        "accident_model": "instantaneous_whole_core_TEC_open_circuit",
        "accident_time_absolute_s": start_time,
        "physical_tfe_count": int(sum(build["ring_multipliers"].values())),
        "current_electric_power_electron_heat_and_joule_heat_cleared": True,
        "passive_tec_gap_heat_transfer_retained": True,
        "external_heat_enabled": True,
        "external_heat_period_s": float(source["external_heat_period_s"]),
        "external_heat_time_origin_s": float(source["external_heat_time_origin_s"]),
    }
    _write_json(out_dir / "accident_event.json", event)
    run_config = dict(source)
    run_config.update({
        **event, "restart_in": str(config.restart_in),
        "source_run_config": str(source_path),
        "point_kinetics_enabled": config.mode == "reactive_feedback",
        "fixed_power_control_active": fixed_power_active,
        "power_w": FIXED_POWER_W, "tec_electrical_calculation_enabled": False,
        "tec_open_circuit_accident_active": True,
        "scram_reactivity_dollars": SCRAM_DOLLARS,
        "minimum_duration_s": float(config.minimum_duration_s),
        "post_scram_duration_s": float(config.post_scram_duration_s),
        "dt_s": float(config.dt_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "temperature_limits_K": limits, "coolant_temperature_is_trip": False,
        "history_schedule_accident_s": [
            [0.0, 20.0, 0.1], [20.0, 100.0, 1.0], [100.0, None, 10.0]
        ],
        "history_schedule_resets_after_scram": True,
    })
    _write_json(out_dir / "run_config.json", run_config)

    rows: list[Dict[str, Any]] = []
    latest_checkpoint = ""

    def external_reactivity() -> tuple[float, float]:
        if not state["scram_active"]:
            return 0.0, 0.0
        return SCRAM_DOLLARS * float(core.point_reactor.beta_total), SCRAM_DOLLARS

    def collect_row(dt: float, stop_reason: str = ""):
        elapsed = float(system.global_time) - start_time
        rho, dollars = external_reactivity()
        base = collect_metrics(build, stage_index=1, dt_s=dt)
        peaks = collect_temperature_peaks(core, limits)
        current_feedback = core.compute_reactivity_feedback()
        trigger = state["scram_trigger"] or {}
        row = {
            **base, **_neutronics_metrics(core, rho, dollars),
            **_temperature_metrics(peaks),
            "time_s": float(system.global_time), "accident_elapsed_s": elapsed,
            "case_mode": config.mode, "external_heat_enabled": True,
            "external_heat_period_s": float(source["external_heat_period_s"]),
            "external_heat_phase_s": (
                float(system.global_time) - float(source["external_heat_time_origin_s"])
            ) % float(source["external_heat_period_s"]),
            "tec_open_circuit_active": True,
            "fixed_power_control_active": bool(state["fixed_power_active"]),
            "scram_active": bool(state["scram_active"]),
            "scram_time_absolute_s": float(state["scram_time_absolute_s"]),
            "scram_elapsed_s": float(state["scram_elapsed_s"]),
            "scram_trigger_component": trigger.get("component", ""),
            "scram_trigger_representative": trigger.get("representative", ""),
            "scram_trigger_axial_position_m": trigger.get("axial_position_m", float("nan")),
            "scram_trigger_actual_K": trigger.get("actual_k", float("nan")),
            "scram_trigger_limit_K": trigger.get("limit_k", float("nan")),
            "core_power_W": float(core.last_total_core_power),
            "feedback_fuel": float(current_feedback.fuel),
            "feedback_electrode": float(current_feedback.electrode),
            "feedback_moderator": float(current_feedback.moderator),
            "feedback_reflector": float(current_feedback.reflector),
            "feedback_total_absolute": float(current_feedback.total),
            "coolant_max_T_K": float(np.max(np.asarray(system.fluid_solver.T_vec))),
            "stop_reason": stop_reason,
        }
        return row, peaks, rho, dollars

    def record(row: Dict[str, Any], rho: float, dollars: float) -> None:
        payload = build_snapshot_payload(
            build, accident, start_time=start_time, coolant_present=True,
            hydraulic_solve_enabled=True, external_reactivity=rho,
            external_reactivity_dollars=dollars,
        )
        append_postprocessing_histories(out_dir, payload)
        _append_summary(history_path, row)
        rows.append(dict(row))
        _write_json(out_dir / "latest_state.json", {
            "case": config.case_name, "latest_checkpoint_path": latest_checkpoint,
            "latest_metrics": row,
        })
        print(
            f"[{config.mode} +{row['accident_elapsed_s']:.2f}s] "
            f"P={row['core_power_W']:.2f}W Twall={row['channel_wall_max_T_K']:.3f}K "
            f"Tfuel={row['fuel_max_T_K']:.3f}K "
            f"Tcollector={row['collector_max_T_K']:.3f}K "
            f"scram={row['scram_active']}", flush=True,
        )

    def trigger_scram(trigger: Dict[str, Any]) -> None:
        nonlocal latest_checkpoint
        elapsed = float(system.global_time) - start_time
        if config.mode == "fixed_power":
            core.initialize_point_reactor(total_power_initial=FIXED_POWER_W)
            state["handoff_type"] = "limit_state_fixed_power_handoff"
            state["fixed_power_active"] = False
            accident["feedback_reference_total"] = float(core.feedback_reference_result.total)
        state.update({
            "scram_active": True, "scram_time_absolute_s": float(system.global_time),
            "scram_elapsed_s": elapsed, "scram_trigger": dict(trigger),
        })
        run_config.update({
            "point_kinetics_enabled": True, "fixed_power_control_active": False,
            "scram_active": True, "scram_time_absolute_s": float(system.global_time),
            "scram_elapsed_s": elapsed, "scram_trigger": trigger,
            "handoff_type": state["handoff_type"],
        })
        _write_json(out_dir / "run_config.json", run_config)
        path = out_dir / "scram_restart.npz"
        system.save_global_state(str(path))
        latest_checkpoint = str(path)
        _write_json(out_dir / "scram_event.json", {
            **trigger, "time_s": float(system.global_time),
            "accident_elapsed_s": elapsed,
            "external_reactivity_dollars": SCRAM_DOLLARS,
        })

    next_record = start_time + 0.1
    next_checkpoint = start_time + float(config.checkpoint_interval_s)
    row, peaks, rho, dollars = collect_row(0.0)
    nonfinite = find_nonfinite_model_state(build)
    if nonfinite is not None:
        raise RuntimeError(f"initial model state is nonfinite: {nonfinite}")
    initial_trip = _temperature_trip(peaks)
    if initial_trip is not None:
        trigger_scram(initial_trip)
        next_record = float(system.global_time) + 0.1
        row, peaks, rho, dollars = collect_row(0.0)
    record(row, rho, dollars)
    last_record_time = float(system.global_time)
    stop_reason = "completed"

    try:
        while True:
            end_elapsed = required_end_elapsed_s(
                config.minimum_duration_s,
                float(state["scram_elapsed_s"]) if state["scram_active"] else None,
                config.post_scram_duration_s,
            )
            end_time = start_time + end_elapsed
            now = float(system.global_time)
            if now >= end_time - 1.0e-9:
                break
            dt = min(
                float(config.dt_s), end_time - now,
                max(1.0e-12, next_record - now),
                max(1.0e-12, next_checkpoint - now),
            )
            rho, _ = external_reactivity()
            if state["fixed_power_active"]:
                _apply_fixed_core_power(build, FIXED_POWER_W)
            system.step(
                dt, inner_iter=int(debug.inner_iter),
                fail_on_fluid_nonconvergence=True,
                fluid_max_iter=int(debug.fluid_max_iter),
                reactivity_control=rho,
            )
            if state["fixed_power_active"]:
                _apply_fixed_core_power(build, FIXED_POWER_W)

            nonfinite = find_nonfinite_model_state(build)
            if nonfinite is not None:
                stop_reason = "nonfinite_model_state"
                _write_json(out_dir / "failure.json", nonfinite)
            row, peaks, rho, dollars = collect_row(dt, stop_reason)
            if not bool(row["fluid_converged"]):
                stop_reason = "hydraulic_nonconvergence"
                row["stop_reason"] = stop_reason
            if not math.isfinite(float(row["core_power_W"])):
                stop_reason = "nonfinite_core_power"
                row["stop_reason"] = stop_reason

            trip_now = False
            if not state["scram_active"] and stop_reason == "completed":
                trigger = _temperature_trip(peaks)
                if trigger is not None:
                    trigger_scram(trigger)
                    trip_now = True
                    next_record = float(system.global_time) + 0.1
                    row, peaks, rho, dollars = collect_row(dt)

            now = float(system.global_time)
            scheduled = now >= next_record - 1.0e-9
            terminal = now >= end_time - 1.0e-9 or stop_reason != "completed"
            if trip_now or scheduled or terminal:
                record(row, rho, dollars)
                last_record_time = now
                if scheduled and not trip_now:
                    origin = (
                        float(state["scram_time_absolute_s"])
                        if state["scram_active"] else start_time
                    )
                    next_record = now + record_interval_s(now - origin)

            if now >= next_checkpoint - 1.0e-9:
                elapsed = now - start_time
                path = out_dir / f"checkpoint_tplus_{elapsed:09.3f}s.npz"
                system.save_global_state(str(path))
                latest_checkpoint = str(path)
                next_checkpoint += float(config.checkpoint_interval_s)

            if stop_reason != "completed":
                path = out_dir / "emergency_restart.npz"
                system.save_global_state(str(path))
                latest_checkpoint = str(path)
                break
    except Exception as exc:
        path = out_dir / "emergency_restart.npz"
        system.save_global_state(str(path))
        _write_json(out_dir / "failure.json", {
            "type": type(exc).__name__, "message": str(exc),
            "time_s": float(system.global_time),
            "accident_elapsed_s": float(system.global_time) - start_time,
        })
        raise

    if abs(float(system.global_time) - last_record_time) > 1.0e-9:
        row, _, rho, dollars = collect_row(
            float(system.global_time) - last_record_time, stop_reason
        )
        record(row, rho, dollars)

    final_restart = out_dir / "final_restart.npz"
    system.save_global_state(str(final_restart))
    result = {
        "case": config.case_name, "mode": config.mode,
        "output_dir": str(out_dir), "source_restart_path": str(config.restart_in),
        "start_time_s": start_time, "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "minimum_duration_s": float(config.minimum_duration_s),
        "post_scram_duration_s": float(config.post_scram_duration_s),
        "scram_active": bool(state["scram_active"]),
        "scram_time_absolute_s": float(state["scram_time_absolute_s"]),
        "scram_elapsed_s": float(state["scram_elapsed_s"]),
        "stop_reason": stop_reason, "final_restart_path": str(final_restart),
        "latest_metrics": rows[-1], "final_window_slopes": _slope_summary(rows),
    }
    _write_json(out_dir / "run_summary.json", result)
    run_config.update({
        "point_kinetics_enabled": bool(core.has_point_reactor),
        "fixed_power_control_active": bool(state["fixed_power_active"]),
        "scram_active": bool(state["scram_active"]),
        "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
    })
    _write_json(out_dir / "run_config.json", run_config)
    _write_json(out_dir / "latest_state.json", {
        "case": config.case_name, "latest_checkpoint_path": str(final_restart),
        "latest_metrics": rows[-1],
    })
    return result


def run_cli(
    mode: str,
    case_name: str,
    default_output_dir: Path,
    argv: Optional[Sequence[str]] = None,
) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=default_output_dir)
    parser.add_argument("--duration", type=float, default=ORBIT_PERIOD_S)
    parser.add_argument("--post-scram-duration", type=float, default=HALF_ORBIT_S)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--checkpoint-interval", type=float, default=50.0)
    parser.add_argument("--wall-limit", type=float, default=1058.0)
    parser.add_argument("--fuel-limit", type=float, default=2700.0)
    parser.add_argument("--collector-limit", type=float, default=1023.0)
    parser.add_argument("--moderator-limit", type=float, default=930.0)
    parser.add_argument("--reflector-limit", type=float, default=1000.0)
    args = parser.parse_args(argv)
    result = run_accident(TecOpenCircuitRunConfig(
        mode=mode, case_name=case_name, output_dir=args.output_dir,
        restart_in=args.restart_in, minimum_duration_s=float(args.duration),
        post_scram_duration_s=float(args.post_scram_duration), dt_s=float(args.dt),
        checkpoint_interval_s=float(args.checkpoint_interval),
        wall_limit_k=float(args.wall_limit), fuel_limit_k=float(args.fuel_limit),
        collector_limit_k=float(args.collector_limit),
        moderator_limit_k=float(args.moderator_limit),
        reflector_limit_k=float(args.reflector_limit),
    ))
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0
