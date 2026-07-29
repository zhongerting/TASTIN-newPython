"""Run V14 instantaneous helium-gap depressurization from the orbital checkpoint."""

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
    _tec_main_metrics, build_debug_case,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1.run_v14_210kw_fixed_power_loca_1 import (
    _failure_temperature_metrics, _neutronics_metrics, _radiator_rejection,
    _read_source_config, _solid_min_max_by_category,
    append_postprocessing_histories, evaluate_failure_reason,
    external_reactivity_for_elapsed, maybe_transition_tec_to_open_circuit,
    staged_record_interval, write_snapshot,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.run_v14_helium_depressurization import (
    HELIUM_H_FINAL_W_M2K, HELIUM_H_INITIAL_W_M2K,
    _refresh_gap_diagnostics, collect_helium_gaps, collect_helium_metrics,
    find_nonfinite_model_state, refresh_tec_now, set_helium_h_eq,
    set_tec_update_interval,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    ReactivityControlRunConfig, load_baseline_debug_config,
    prepare_reactivity_control,
)

CASE_DIR = Path(__file__).resolve().parent
CASE_NAME = "V14_10kW_210kW_helium_depressurization_1"
DEFAULT_RESTART = (
    CASE_DIR.parent / "V14_210kW_fixed_power_external_heat_2orbits" / "runs"
    / "two_orbits_from13864_20260720" / "checkpoint_t019265s.npz"
)
DEFAULT_OUTPUT_DIR = CASE_DIR / "runs" / "default"
SUMMARY_FIELDS = [
    "time_s", "accident_elapsed_s", "dt_s", "hydraulic_solve_enabled",
    "helium_h_eq_W_m2K", "helium_gap_heat_out_scaled_W",
    "core_power_W", "fission_power_W", "decay_power_W",
    "external_reactivity", "external_reactivity_dollars",
    "effective_temperature_feedback", "total_reactivity",
    "feedback_fuel", "feedback_electrode", "feedback_moderator",
    "feedback_reflector", "feedback_total_absolute",
    "feedback_total_change_from_accident",
    "tec_main_current_A", "tec_main_voltage_V", "tec_main_electric_power_W",
    "tec_main_converged", "tec_open_circuit_active",
    "tec_open_circuit_time_s", "tec_open_circuit_trigger_current_A",
    "coolant_max_T_K", "collector_max_T_K", "emitter_max_T_K",
    "moderator_max_T_K", "reflector_max_T_K",
    "core_structure_min_T_K", "core_structure_max_T_K",
    "pipe_wall_min_T_K", "pipe_wall_max_T_K",
    "heat_pipe_min_T_K", "heat_pipe_max_T_K",
    "radiator_net_rejection_W", "failure_reason", "snapshot_path",
]


@dataclass(frozen=True)
class HeliumRunConfig:
    restart_in: Path = DEFAULT_RESTART
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = 2000.0
    dt_s: float = 0.05
    checkpoint_interval_s: float = 100.0
    tec_update_interval_s: float = 0.05
    collector_failure_temperature_k: float = 1500.0
    emitter_failure_temperature_k: float = 3000.0
    coolant_failure_temperature_k: float = 1058.0
    moderator_failure_temperature_k: float = 930.0
    reflector_failure_temperature_k: float = 1000.0
    scram_time_s: Optional[float] = None
    scram_reactivity_dollars: float = -2.0
    tec_open_circuit_current_threshold_a: float = 0.01


def _validate_config(config: HeliumRunConfig) -> None:
    if not Path(config.restart_in).is_file():
        raise FileNotFoundError(f"restart not found: {config.restart_in}")
    for name in ("duration_s", "dt_s", "tec_update_interval_s"):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(float(config.checkpoint_interval_s)):
        raise ValueError("checkpoint_interval_s must be finite")
    for name in (
        "collector_failure_temperature_k", "emitter_failure_temperature_k",
        "coolant_failure_temperature_k", "moderator_failure_temperature_k",
        "reflector_failure_temperature_k",
    ):
        if not math.isfinite(float(getattr(config, name))):
            raise ValueError(f"{name} must be finite")
    if config.scram_time_s is not None and (
        not math.isfinite(float(config.scram_time_s))
        or float(config.scram_time_s) < 0.0
    ):
        raise ValueError("scram_time_s must be finite and non-negative")
    threshold = float(config.tec_open_circuit_current_threshold_a)
    if not math.isfinite(threshold) or threshold < 0.0:
        raise ValueError("TEC open-circuit threshold must be finite and non-negative")


def _append_summary(path: Path, row: Dict[str, Any]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _collect_summary(
    build: Dict[str, Any], gaps: Dict[str, tuple[Any, int]],
    accident: Dict[str, Any], *, start_time: float, dt_s: float,
    snapshot_path: str, failure_reason: str = "",
    external_reactivity: float = 0.0,
    external_reactivity_dollars: float = 0.0,
) -> Dict[str, Any]:
    system, core = build["system"], build["core"]
    feedback = core.compute_reactivity_feedback()
    row = {
        "time_s": float(system.global_time),
        "accident_elapsed_s": float(system.global_time) - start_time,
        "dt_s": float(dt_s), "hydraulic_solve_enabled": True,
        "core_power_W": float(core.last_total_core_power),
        "failure_reason": failure_reason,
        "feedback_fuel": float(feedback.fuel),
        "feedback_electrode": float(feedback.electrode),
        "feedback_moderator": float(feedback.moderator),
        "feedback_reflector": float(feedback.reflector),
        "feedback_total_absolute": float(feedback.total),
        "feedback_total_change_from_accident": (
            float(feedback.total) - accident["feedback_reference_total"]
        ),
        "tec_open_circuit_active": bool(accident["tec_open_circuit_active"]),
        "tec_open_circuit_time_s": float(accident["tec_open_circuit_time_s"]),
        "tec_open_circuit_trigger_current_A": float(
            accident["tec_open_circuit_trigger_current_A"]
        ),
        "radiator_net_rejection_W": _radiator_rejection(build),
        "snapshot_path": snapshot_path,
        **_tec_main_metrics(build),
        **_neutronics_metrics(core, external_reactivity, external_reactivity_dollars),
        **_solid_min_max_by_category(system),
        **_failure_temperature_metrics(system, coolant_present=True),
    }
    row.update(collect_helium_metrics(
        build, gaps, accident_time_s=start_time, active=True,
    ))
    return row


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )


def run_helium_accident(config: HeliumRunConfig) -> Dict[str, Any]:
    _validate_config(config)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / "history.csv"
    if history_path.exists():
        raise FileExistsError(f"output history already exists: {history_path}")

    source = _read_source_config(Path(config.restart_in))
    if not bool(source.get("external_heat_enabled", False)):
        raise ValueError("source checkpoint must retain orbital external heat")
    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in), output_dir=out_dir,
        duration_s=float(config.duration_s), dt_s=float(config.dt_s),
        record_interval_s=0.5,
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=None, external_heat_enabled=True,
        external_heat_period_s=float(source["external_heat_period_s"]),
        external_heat_time_origin_s=float(source["external_heat_time_origin_s"]),
    )
    debug, source = load_baseline_debug_config(runtime)
    build = build_debug_case(debug, apply_fixed_power=False)
    system, core = build["system"], build["core"]
    handoff_type = prepare_reactivity_control(
        core, source_point_kinetics_enabled=bool(source["point_kinetics_enabled"]),
        expected_power_w=float(source["power_w"]),
    )
    set_tec_update_interval(core, float(config.tec_update_interval_s))
    start_time = float(system.global_time)
    refresh_tec_now(core, start_time)
    gaps = collect_helium_gaps(build)
    accident = {
        "reference_fluid": {
            "T": np.asarray(system.fluid_solver.T_vec).copy(),
            "P": np.asarray(system.fluid_solver.P_vec).copy(),
            "h": np.asarray(system.fluid_solver.h_vec).copy(),
            "W": np.asarray(system.fluid_solver.W_vec).copy(),
        },
        "feedback_reference_total": float(core.compute_reactivity_feedback().total),
        "tec_open_circuit_active": False,
        "tec_open_circuit_time_s": float("nan"),
        "tec_open_circuit_trigger_current_A": float("nan"),
    }
    initial_reason, initial_metrics = evaluate_failure_reason(
        system, config, coolant_present=True,
    )
    if initial_reason:
        raise RuntimeError(f"initial state violates limit: {initial_reason}: {initial_metrics}")
    if find_nonfinite_model_state(build) is not None:
        raise RuntimeError("initial model state is nonfinite")

    write_snapshot(
        out_dir / "snapshot_pre_accident.npz", build, accident,
        start_time=start_time, coolant_present=True,
        hydraulic_solve_enabled=True,
    )
    set_helium_h_eq(gaps, HELIUM_H_FINAL_W_M2K)
    _refresh_gap_diagnostics(build, gaps)
    event = {
        "case": CASE_NAME, "accident_time_absolute_s": start_time,
        "accident_model": "instantaneous_total_loss_of_helium_gas_conduction",
        "helium_h_initial_W_m2K": HELIUM_H_INITIAL_W_M2K,
        "helium_h_final_W_m2K": HELIUM_H_FINAL_W_M2K,
        "affected_representative_tfes": list(gaps),
        "physical_tfe_count": sum(multiplier for _, multiplier in gaps.values()),
        "hydraulic_solve_enabled": True, "external_heat_enabled": True,
    }
    _write_json(out_dir / "accident_event.json", event)
    run_config = dict(source)
    run_config.update({
        **event, "restart_in": str(config.restart_in),
        "source_run_config": str(Path(config.restart_in).parent / "run_config.json"),
        "duration_s": float(config.duration_s), "dt_s": float(config.dt_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "tec_update_interval_s": float(config.tec_update_interval_s),
        "point_kinetics_enabled": True, "fixed_power_enabled": False,
        "reactivity_control_mode": (
            "temperature_feedback_plus_scram"
            if config.scram_time_s is not None else "temperature_feedback_only"
        ),
        "handoff_type": handoff_type, "scram_time_s": config.scram_time_s,
        "scram_reactivity_dollars": (
            float(config.scram_reactivity_dollars)
            if config.scram_time_s is not None else 0.0
        ),
        "tec_open_circuit_after_scram": config.scram_time_s is not None,
        "tec_open_circuit_current_threshold_A": float(
            config.tec_open_circuit_current_threshold_a
        ),
        "record_schedule_s": [
            [0.0, 20.0, 0.5], [20.0, 100.0, 2.0],
            [100.0, 400.0, 5.0], [400.0, 600.0, 10.0],
            [600.0, float(config.duration_s), 20.0],
        ],
        "failure_temperature_limits_K": {
            "collector": config.collector_failure_temperature_k,
            "emitter": config.emitter_failure_temperature_k,
            "coolant": config.coolant_failure_temperature_k,
            "moderator": config.moderator_failure_temperature_k,
            "reflector": config.reflector_failure_temperature_k,
        },
    })
    _write_json(out_dir / "run_config.json", run_config)

    post_path = out_dir / "snapshot_tplus_00000.000s.npz"
    payload = write_snapshot(
        post_path, build, accident, start_time=start_time,
        coolant_present=True, hydraulic_solve_enabled=True,
    )
    append_postprocessing_histories(out_dir, payload)
    latest = _collect_summary(
        build, gaps, accident, start_time=start_time, dt_s=0.0,
        snapshot_path=str(post_path),
    )
    _append_summary(history_path, latest)

    end_time = start_time + float(config.duration_s)
    next_record = staged_record_interval(0.0)
    last_checkpoint = 0.0
    stop_reason = "completed"
    trip_payload: Dict[str, Any] = {}
    while float(system.global_time) < end_time - 1.0e-9:
        dt = min(float(config.dt_s), end_time - float(system.global_time))
        elapsed_before = float(system.global_time) - start_time
        rho, dollars = external_reactivity_for_elapsed(core, config, elapsed_before)
        system.step(
            dt, inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=True,
            fluid_max_iter=int(debug.fluid_max_iter),
            reactivity_control=rho,
        )
        elapsed = float(system.global_time) - start_time
        if maybe_transition_tec_to_open_circuit(build, accident, config, elapsed):
            print(
                f"[He-loss +{elapsed:.3f}s] TEC open circuit: "
                f"I={accident['tec_open_circuit_trigger_current_A']:.6g} A",
                flush=True,
            )
        rho, dollars = external_reactivity_for_elapsed(core, config, elapsed)
        nonfinite = find_nonfinite_model_state(build)
        if nonfinite is not None:
            stop_reason, trip_payload = "nonfinite_model_state", nonfinite
        else:
            stop_reason, trip_payload = evaluate_failure_reason(
                system, config, coolant_present=True,
            )
        if (
            elapsed + 1.0e-9 >= next_record
            or float(system.global_time) >= end_time - 1.0e-9
            or stop_reason
        ):
            snapshot = out_dir / f"snapshot_tplus_{elapsed:09.3f}s.npz"
            payload = write_snapshot(
                snapshot, build, accident, start_time=start_time,
                coolant_present=True, hydraulic_solve_enabled=True,
                external_reactivity=rho, external_reactivity_dollars=dollars,
            )
            append_postprocessing_histories(out_dir, payload)
            latest = _collect_summary(
                build, gaps, accident, start_time=start_time, dt_s=dt,
                snapshot_path=str(snapshot), failure_reason=stop_reason,
                external_reactivity=rho, external_reactivity_dollars=dollars,
            )
            _append_summary(history_path, latest)
            print(
                f"[He-loss +{elapsed:.3f}s] Tc={latest['collector_max_T_K']:.3f} K "
                f"Te={latest['emitter_max_T_K']:.3f} K "
                f"Tcool={latest['coolant_max_T_K']:.3f} K "
                f"P={latest['core_power_W']:.3f} W "
                f"I={latest['tec_main_current_A']:.6g} A",
                flush=True,
            )
            while next_record <= elapsed + 1.0e-9:
                next_record += staged_record_interval(next_record)
        if stop_reason:
            system.save_global_state(str(out_dir / "emergency_restart.npz"))
            _write_json(out_dir / "limit_trip.json", {
                **trip_payload, "stop_reason": stop_reason,
                "time_s": float(system.global_time),
                "accident_elapsed_s": elapsed,
            })
            break
        if (
            float(config.checkpoint_interval_s) > 0.0
            and elapsed - last_checkpoint >= float(config.checkpoint_interval_s) - 1.0e-9
        ):
            system.save_global_state(
                str(out_dir / f"checkpoint_tplus_{elapsed:.3f}s.npz")
            )
            last_checkpoint = elapsed

    if not stop_reason:
        stop_reason = "completed"
    restart_path = out_dir / "stage_01_restart.npz"
    system.save_global_state(str(restart_path))
    result = {
        "case": CASE_NAME, "output_dir": str(out_dir),
        "source_restart_path": str(config.restart_in),
        "start_time_s": start_time, "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "history_path": str(history_path),
        "postprocessing_history_paths": {
            name: str(out_dir / f"history_{name}.csv")
            for name in ("coolant", "solids", "electrical", "reactivity")
        },
        "restart_path": str(restart_path), "latest_metrics": latest,
        "stop_reason": stop_reason,
    }
    _write_json(out_dir / "run_summary.json", result)
    _write_json(out_dir / "latest_state.json", result)
    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--duration", type=float, default=2000.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--checkpoint-interval", type=float, default=100.0)
    parser.add_argument("--tec-update-interval", type=float, default=0.05)
    parser.add_argument("--scram-time", type=float)
    parser.add_argument("--scram-reactivity-dollars", type=float, default=-2.0)
    parser.add_argument("--tec-open-circuit-current-threshold", type=float, default=0.01)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = run_helium_accident(HeliumRunConfig(
        restart_in=args.restart_in, output_dir=args.output_dir,
        duration_s=float(args.duration), dt_s=float(args.dt),
        checkpoint_interval_s=float(args.checkpoint_interval),
        tec_update_interval_s=float(args.tec_update_interval),
        scram_time_s=args.scram_time,
        scram_reactivity_dollars=float(args.scram_reactivity_dollars),
        tec_open_circuit_current_threshold_a=float(
            args.tec_open_circuit_current_threshold
        ),
    ))
    print(json.dumps(result["latest_metrics"], indent=2, sort_keys=True, ensure_ascii=False))
    print(f"Stop reason: {result['stop_reason']}")
    print(f"Saved outputs to: {result['output_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
