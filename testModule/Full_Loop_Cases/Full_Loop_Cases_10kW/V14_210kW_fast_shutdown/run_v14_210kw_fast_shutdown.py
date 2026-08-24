"""Run V14 from its two-orbit steady state with an immediate -2 dollar scram."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Sequence

import numpy as np

REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    build_debug_case,
    collect_metrics,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1.run_v14_210kw_fixed_power_loca_1 import (
    _neutronics_metrics,
    append_postprocessing_histories,
    build_snapshot_payload,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.run_v14_helium_depressurization import (
    find_nonfinite_model_state,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    ReactivityControlRunConfig,
    load_baseline_debug_config,
    prepare_reactivity_control,
)


CASE_DIR = Path(__file__).resolve().parent
CASE_NAME = "V14_210kW_fast_shutdown"
ORBIT_PERIOD_S = 5668.144369
SCRAM_DOLLARS = -2.0
DEFAULT_RESTART = (
    CASE_DIR.parent / "V14_210kW_fixed_power_external_heat_2orbits" / "runs"
    / "two_orbits_from13864_20260720" / "stage_01_restart.npz"
)
DEFAULT_OUTPUT_DIR = CASE_DIR / "runs" / "one_orbit_minus2dollar"

SUMMARY_FIELDS = [
    "time_s", "shutdown_elapsed_s", "dt_s", "core_power_W",
    "fission_power_W", "decay_power_W", "external_reactivity",
    "external_reactivity_dollars", "effective_temperature_feedback",
    "total_reactivity", "feedback_fuel", "feedback_electrode",
    "feedback_moderator", "feedback_reflector", "feedback_total_absolute",
    "core_inlet_T_K", "core_outlet_T_K", "min_fluid_T_K", "max_fluid_T_K",
    "min_solid_T_K", "max_solid_T_K", "radiator_heat_rejection_W",
    "total_external_heat_rejection_W", "pump_a_flow_kg_s", "pump_b_flow_kg_s",
    "fluid_converged", "tec_main_current_A", "tec_main_voltage_V",
    "tec_main_electric_power_W", "tec_main_converged",
    "tec_open_circuit_active", "tec_open_circuit_time_s", "stop_reason",
]


@dataclass(frozen=True)
class FastShutdownRunConfig:
    restart_in: Path = DEFAULT_RESTART
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = ORBIT_PERIOD_S
    dt_s: float = 0.05
    history_interval_s: float = 1.0
    checkpoint_interval_s: float = 100.0
    scram_reactivity_dollars: float = SCRAM_DOLLARS
    tec_open_current_threshold_a: float = 0.01
    continuation_settle_duration_s: float = 2.0
    continuation_settle_dt_s: float = 0.01


def _validate(config: FastShutdownRunConfig) -> None:
    if not Path(config.restart_in).is_file():
        raise FileNotFoundError(f"restart not found: {config.restart_in}")
    for name in (
        "duration_s", "dt_s", "history_interval_s",
        "continuation_settle_dt_s",
    ):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(float(config.checkpoint_interval_s)) or float(config.checkpoint_interval_s) < 0.0:
        raise ValueError("checkpoint_interval_s must be finite and non-negative")
    if (
        not math.isfinite(float(config.continuation_settle_duration_s))
        or float(config.continuation_settle_duration_s) < 0.0
    ):
        raise ValueError("continuation_settle_duration_s must be finite and non-negative")
    if not math.isfinite(float(config.scram_reactivity_dollars)) or float(config.scram_reactivity_dollars) >= 0.0:
        raise ValueError("scram_reactivity_dollars must be finite and negative")
    if not math.isfinite(float(config.tec_open_current_threshold_a)) or float(config.tec_open_current_threshold_a) < 0.0:
        raise ValueError("tec_open_current_threshold_a must be finite and non-negative")


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def _append_summary(path: Path, row: Dict[str, Any]) -> None:
    new_file = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        writer.writerow(row)


def external_reactivity(core: Any, dollars: float) -> float:
    return float(dollars) * float(core.point_reactor.beta_total)


def maybe_open_tec(build: Dict[str, Any], state: Dict[str, Any], threshold_a: float, elapsed_s: float) -> bool:
    core = build["core"]
    if state["tec_open_circuit_active"] or not core.enable_tec_coupled:
        return False
    results = core.thermo_calc.get_global_results()
    current = float(results.get("Iout", float("nan")))
    if not math.isfinite(current) or current > float(threshold_a):
        return False
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    core.enable_tec_coupled = False
    state["tec_open_circuit_active"] = True
    state["tec_open_circuit_time_s"] = float(elapsed_s)
    state["tec_open_circuit_trigger_current_A"] = current
    return True


def _restore_open_tec_continuation(
    build: Dict[str, Any], source_path: Path, threshold_a: float,
) -> tuple[bool, float]:
    """Restore an already-zero TEC without reusing the fixed-U current seed."""
    history_path = source_path.parent / "history.csv"
    if not history_path.is_file():
        return False, float("nan")
    first_zero_elapsed = float("nan")
    latest_current = float("nan")
    with history_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            latest_current = float(row["tec_main_current_A"])
            if (
                not math.isfinite(first_zero_elapsed)
                and math.isfinite(latest_current)
                and latest_current <= float(threshold_a)
            ):
                first_zero_elapsed = float(row["shutdown_elapsed_s"])
    if not math.isfinite(latest_current) or latest_current > float(threshold_a):
        return False, float("nan")
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    build["core"].enable_tec_coupled = False
    return True, first_zero_elapsed


def run_fast_shutdown(config: FastShutdownRunConfig) -> Dict[str, Any]:
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

    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in), output_dir=out_dir,
        duration_s=float(config.duration_s), dt_s=float(config.dt_s),
        record_interval_s=float(config.history_interval_s),
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=None, external_heat_enabled=True,
        external_heat_period_s=float(source["external_heat_period_s"]),
        external_heat_time_origin_s=float(source["external_heat_time_origin_s"]),
    )
    debug, source = load_baseline_debug_config(runtime)
    build = build_debug_case(debug, apply_fixed_power=False)
    system, core = build["system"], build["core"]
    initial_power = float(core.last_total_core_power)
    source_point_kinetics_enabled = bool(source.get("point_kinetics_enabled", False))
    handoff_type = prepare_reactivity_control(
        core, source_point_kinetics_enabled=source_point_kinetics_enabled,
        expected_power_w=float(source["power_w"]),
    )
    segment_start_time = float(system.global_time)
    shutdown_origin_time = float(
        source.get("shutdown_time_absolute_s", segment_start_time)
    )
    tec_already_open, tec_open_elapsed = (
        _restore_open_tec_continuation(
            build, source_path, config.tec_open_current_threshold_a,
        )
        if source_point_kinetics_enabled else (False, float("nan"))
    )
    rho = external_reactivity(core, config.scram_reactivity_dollars)
    state: Dict[str, Any] = {
        "reference_fluid": {
            "T": np.asarray(system.fluid_solver.T_vec).copy(),
            "P": np.asarray(system.fluid_solver.P_vec).copy(),
            "h": np.asarray(system.fluid_solver.h_vec).copy(),
            "W": np.asarray(system.fluid_solver.W_vec).copy(),
        },
        "feedback_reference_total": float(core.feedback_reference_result.total),
        "tec_open_circuit_active": tec_already_open,
        "tec_open_circuit_time_s": tec_open_elapsed,
        "continuation_settle_duration_s": (
            float(config.continuation_settle_duration_s)
            if source_point_kinetics_enabled else 0.0
        ),
        "continuation_settle_dt_s": float(config.continuation_settle_dt_s),
        "tec_open_circuit_trigger_current_A": float("nan"),
    }

    event_restart = out_dir / (
        "continuation_start_restart.npz"
        if source_point_kinetics_enabled else "shutdown_event_restart.npz"
    )
    system.save_global_state(str(event_restart))
    run_config = dict(source)
    run_config.update({
        "case": CASE_NAME, "restart_in": str(config.restart_in),
        "source_run_config": str(source_path), "handoff_type": handoff_type,
        "shutdown_time_absolute_s": shutdown_origin_time,
        "segment_start_time_absolute_s": segment_start_time,
        "point_kinetics_enabled": True,
        "fixed_power_control_active": False,
        "scram_reactivity_dollars": float(config.scram_reactivity_dollars),
        "scram_reactivity": rho, "target_flow_kg_s": float(source["target_flow_kg_s"]),
        "external_heat_time_origin_s": float(source["external_heat_time_origin_s"]),
        "duration_s": float(config.duration_s), "dt_s": float(config.dt_s),
        "history_interval_s": float(config.history_interval_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "tec_open_circuit_current_threshold_A": float(config.tec_open_current_threshold_a),
        "tec_open_circuit_active_at_segment_start": tec_already_open,
        "tec_open_circuit_time_s": tec_open_elapsed,
        "segment_start_restart": str(event_restart),
    })
    _write_json(out_dir / "run_config.json", run_config)
    _write_json(out_dir / (
        "continuation_start_event.json"
        if source_point_kinetics_enabled else "shutdown_event.json"
    ), {
        "time_s": segment_start_time,
        "shutdown_elapsed_s": segment_start_time - shutdown_origin_time,
        "handoff_type": handoff_type,
        "external_reactivity_dollars": float(config.scram_reactivity_dollars),
        "external_reactivity": rho, "initial_total_power_W": initial_power,
        "initial_fission_power_W": float(core.point_reactor.fission_power),
        "initial_decay_power_W": float(core.point_reactor.decay_power),
        "restart_path": str(event_restart),
    })

    latest_checkpoint = str(event_restart)
    rows: list[Dict[str, Any]] = []

    def collect_row(dt_s: float, stop_reason: str = "") -> Dict[str, Any]:
        base = collect_metrics(build, stage_index=1, dt_s=dt_s)
        feedback = core.compute_reactivity_feedback()
        row = {
            **base, **_neutronics_metrics(core, rho, float(config.scram_reactivity_dollars)),
            "shutdown_elapsed_s": float(system.global_time) - shutdown_origin_time,
            "core_power_W": float(core.last_total_core_power),
            "feedback_fuel": float(feedback.fuel),
            "feedback_electrode": float(feedback.electrode),
            "feedback_moderator": float(feedback.moderator),
            "feedback_reflector": float(feedback.reflector),
            "feedback_total_absolute": float(feedback.total),
            "tec_open_circuit_active": bool(state["tec_open_circuit_active"]),
            "tec_open_circuit_time_s": float(state["tec_open_circuit_time_s"]),
            "stop_reason": stop_reason,
        }
        return row

    def record(row: Dict[str, Any]) -> None:
        payload = build_snapshot_payload(
            build, state, start_time=shutdown_origin_time, coolant_present=True,
            hydraulic_solve_enabled=True, external_reactivity=rho,
            external_reactivity_dollars=float(config.scram_reactivity_dollars),
        )
        payload["shutdown_elapsed_s"] = payload["accident_elapsed_s"].copy()
        append_postprocessing_histories(out_dir, payload)
        _append_summary(history_path, row)
        rows.append(dict(row))
        _write_json(out_dir / "latest_state.json", {
            "case": CASE_NAME, "latest_checkpoint_path": latest_checkpoint,
            "latest_metrics": row,
        })
        print(
            f"[shutdown +{row['shutdown_elapsed_s']:.2f}s] "
            f"P={row['core_power_W']:.3f} W "
            f"Pf={row['fission_power_W']:.3f} W Pd={row['decay_power_W']:.3f} W "
            f"I={row['tec_main_current_A']:.3f} A "
            f"Tfluid=[{row['min_fluid_T_K']:.3f}, {row['max_fluid_T_K']:.3f}] K",
            flush=True,
        )

    record(collect_row(0.0))
    end_time = segment_start_time + float(config.duration_s)
    next_record = segment_start_time + float(config.history_interval_s)
    next_checkpoint = (
        segment_start_time + float(config.checkpoint_interval_s)
        if config.checkpoint_interval_s > 0.0 else float("inf")
    )
    stop_reason = "completed"
    try:
        while float(system.global_time) < end_time - 1.0e-9:
            now = float(system.global_time)
            segment_elapsed = now - segment_start_time
            settling = (
                source_point_kinetics_enabled
                and segment_elapsed < float(config.continuation_settle_duration_s) - 1.0e-12
            )
            requested_dt = (
                float(config.continuation_settle_dt_s)
                if settling else float(config.dt_s)
            )
            dt = min(
                requested_dt, end_time - now,
                max(1.0e-12, next_record - now),
                max(1.0e-12, next_checkpoint - now),
            )
            if settling:
                dt = min(
                    dt,
                    segment_start_time
                    + float(config.continuation_settle_duration_s) - now,
                )
            system.step(
                dt, inner_iter=int(debug.inner_iter),
                fail_on_fluid_nonconvergence=True,
                fluid_max_iter=int(debug.fluid_max_iter),
                reactivity_control=rho,
            )
            elapsed = float(system.global_time) - shutdown_origin_time
            if maybe_open_tec(build, state, config.tec_open_current_threshold_a, elapsed):
                path = out_dir / "tec_open_event_restart.npz"
                system.save_global_state(str(path))
                latest_checkpoint = str(path)
                _write_json(out_dir / "tec_open_event.json", {
                    "time_s": float(system.global_time), "shutdown_elapsed_s": elapsed,
                    "trigger_current_A": state["tec_open_circuit_trigger_current_A"],
                    "restart_path": str(path),
                })
            nonfinite = find_nonfinite_model_state(build)
            if nonfinite is not None:
                stop_reason = "nonfinite_model_state"
                _write_json(out_dir / "failure.json", nonfinite)
            row = collect_row(dt, stop_reason)
            if not bool(row["fluid_converged"]):
                stop_reason = "hydraulic_nonconvergence"
                row["stop_reason"] = stop_reason
            now = float(system.global_time)
            if now >= next_record - 1.0e-9 or now >= end_time - 1.0e-9 or stop_reason != "completed":
                record(row)
                next_record += float(config.history_interval_s)
            if now >= next_checkpoint - 1.0e-9:
                path = out_dir / f"checkpoint_tplus_{now - shutdown_origin_time:09.3f}s.npz"
                system.save_global_state(str(path))
                latest_checkpoint = str(path)
                next_checkpoint += float(config.checkpoint_interval_s)
            if stop_reason != "completed":
                break
    except Exception as exc:
        path = out_dir / "emergency_restart.npz"
        system.save_global_state(str(path))
        _write_json(out_dir / "failure.json", {
            "type": type(exc).__name__, "message": str(exc),
            "time_s": float(system.global_time),
            "shutdown_elapsed_s": float(system.global_time) - shutdown_origin_time,
        })
        raise

    final_restart = out_dir / "final_restart.npz"
    system.save_global_state(str(final_restart))
    result = {
        "case": CASE_NAME, "source_restart_path": str(config.restart_in),
        "start_time_s": segment_start_time, "end_time_s": float(system.global_time),
        "shutdown_origin_time_s": shutdown_origin_time,
        "duration_s": float(system.global_time) - segment_start_time,
        "scram_reactivity_dollars": float(config.scram_reactivity_dollars),
        "stop_reason": stop_reason, "final_restart_path": str(final_restart),
        "latest_metrics": rows[-1],
    }
    _write_json(out_dir / "run_summary.json", result)
    return result


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--duration", type=float, default=ORBIT_PERIOD_S)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--history-interval", type=float, default=1.0)
    parser.add_argument("--checkpoint-interval", type=float, default=100.0)
    parser.add_argument("--scram-dollars", type=float, default=SCRAM_DOLLARS)
    parser.add_argument("--continuation-settle-duration", type=float, default=2.0)
    parser.add_argument("--continuation-settle-dt", type=float, default=0.01)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    result = run_fast_shutdown(FastShutdownRunConfig(
        restart_in=args.restart_in, output_dir=args.output_dir,
        duration_s=float(args.duration), dt_s=float(args.dt),
        history_interval_s=float(args.history_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        scram_reactivity_dollars=float(args.scram_dollars),
        continuation_settle_duration_s=float(args.continuation_settle_duration),
        continuation_settle_dt_s=float(args.continuation_settle_dt),
    ))
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
