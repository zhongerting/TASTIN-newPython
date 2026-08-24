"""V14 heat-pipe transfer-failure accident runner.

The accident keeps the nominal RingHP objects and hydraulic network intact.  A
failure only scales the affected evaporator fluid-solid coupling multiplier,
so the nominal local-loss map and flow resistance remain active.
All three case directories call this module with a small, explicit failure
signature.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence
import sys

import numpy as np

REPO_ROOT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "testModule").is_dir()
)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (  # noqa: E402
    DebugRunConfig,
    _apply_fixed_core_power,
    build_debug_case,
    collect_metrics,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1.run_v14_210kw_fixed_power_loca_1 import (  # noqa: E402
    _neutronics_metrics,
    append_postprocessing_histories,
    build_snapshot_payload,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.run_v14_helium_depressurization import (  # noqa: E402
    collect_temperature_peaks,
    find_nonfinite_model_state,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (  # noqa: E402
    ReactivityControlRunConfig,
    load_baseline_debug_config,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.v14_heatpipe_radiator import (  # noqa: E402
    LOWER_HP_MULTIPLIERS,
    SEGMENT_PATH,
    UPPER_HP_MULTIPLIERS,
)


ORBIT_PERIOD_S = 5668.144369
HALF_ORBIT_S = ORBIT_PERIOD_S / 2.0
FIXED_POWER_W = 210000.0
SCRAM_DOLLARS = -2.0
DEFAULT_RESTART = (
    Path(__file__).resolve().parent
    / "V14_210kW_fixed_power_external_heat_2orbits"
    / "runs"
    / "two_orbits_from13864_20260720"
    / "checkpoint_t019865s.npz"
)
DEFAULT_LIMITS_K = {
    "channel_wall": 1058.0,
    "pellet": 2700.0,
    "collector": 1023.0,
    "moderator": 930.0,
    "reflector": 1000.0,
}

SUMMARY_FIELDS = [
    "time_s", "accident_elapsed_s", "dt_s", "case_mode",
    "failure_mode", "external_heat_enabled", "external_heat_period_s",
    "external_heat_phase_s", "tec_electrical_calculation_enabled",
    "fixed_power_control_active", "point_kinetics_enabled",
    "scram_active", "scram_time_absolute_s", "scram_elapsed_s",
    "scram_trigger_component", "scram_trigger_representative",
    "scram_trigger_axial_position_m", "scram_trigger_actual_K",
    "scram_trigger_limit_K", "core_power_W", "fission_power_W",
    "decay_power_W", "external_reactivity", "external_reactivity_dollars",
    "effective_temperature_feedback", "total_reactivity", "feedback_fuel",
    "feedback_electrode", "feedback_moderator", "feedback_reflector",
    "feedback_total_absolute", "channel_wall_max_T_K", "fuel_max_T_K",
    "collector_max_T_K", "moderator_max_T_K", "reflector_max_T_K",
    "coolant_max_T_K", "core_inlet_T_K", "core_outlet_T_K",
    "radiator_heat_rejection_W", "total_external_heat_rejection_W",
    "pump_a_flow_kg_s", "pump_b_flow_kg_s", "fluid_converged",
    "stop_reason",
]


@dataclass(frozen=True)
class HeatPipeFailureRunConfig:
    """One fixed-power V14 heat-pipe accident calculation."""

    case_name: str
    failure_mode: str
    output_dir: Path
    failure_map: Mapping[str, Mapping[int, Mapping[int, float]]]
    restart_in: Path = DEFAULT_RESTART
    minimum_duration_s: float = ORBIT_PERIOD_S
    post_scram_duration_s: float = HALF_ORBIT_S
    dt_s: float = 0.05
    history_interval_s: float = 1.0
    checkpoint_interval_s: float = 100.0
    wall_limit_k: float = 1058.0
    fuel_limit_k: float = 2700.0
    collector_limit_k: float = 1023.0
    moderator_limit_k: float = 930.0
    reflector_limit_k: float = 1000.0


def _finite(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def _limits(config: HeatPipeFailureRunConfig) -> Dict[str, float]:
    return {
        "channel_wall": float(config.wall_limit_k),
        "pellet": float(config.fuel_limit_k),
        "collector": float(config.collector_limit_k),
        "moderator": float(config.moderator_limit_k),
        "reflector": float(config.reflector_limit_k),
    }


def _temperature_metrics(peaks: Sequence[Dict[str, Any]]) -> Dict[str, float]:
    result: Dict[str, float] = {}
    for component, field in (
        ("channel_wall", "channel_wall_max_T_K"),
        ("pellet", "fuel_max_T_K"),
        ("collector", "collector_max_T_K"),
        ("moderator", "moderator_max_T_K"),
        ("reflector", "reflector_max_T_K"),
    ):
        values = [
            float(peak["actual_k"])
            for peak in peaks
            if peak.get("source_component", peak["component"]) == component
        ]
        result[field] = max(values) if values else float("nan")
    return result


def _temperature_trip(
        peaks: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
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
        writer = csv.DictWriter(
            stream, fieldnames=SUMMARY_FIELDS, extrasaction="ignore"
        )
        if new_file:
            writer.writeheader()
        writer.writerow(row)


def _slope_summary(rows: Sequence[Dict[str, Any]], window_s: float = 300.0) -> Dict[str, Any]:
    if not rows:
        return {"window_s": float(window_s), "sample_count": 0}
    end = float(rows[-1]["time_s"])
    selected = [row for row in rows if float(row["time_s"]) >= end - window_s]
    result: Dict[str, Any] = {
        "window_s": float(window_s), "sample_count": len(selected)
    }
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
        if np.all(np.isfinite(y)):
            result[f"{field}_slope_per_s"] = float(np.polyfit(x, y, 1)[0])
    return result


def _json_failure_map(
        failure_map: Mapping[str, Mapping[int, Mapping[int, float]]]) -> Dict[str, Dict[str, Dict[str, float]]]:
    return {
        str(ring): {
            str(sector): {
                str(node): float(factor)
                for node, factor in nodes.items()
            }
            for sector, nodes in sectors.items()
        }
        for ring, sectors in failure_map.items()
    }


def _thermal_transfer_maps(
        signature: Sequence[Mapping[str, Any]],
        ) -> tuple[Dict[str, list[list[float]]], Dict[str, list[list[float]]]]:
    """Build complete nominal/effective maps for restart provenance."""
    nominal = {
        "upper": [[1.0, 1.0, 1.0] for _ in UPPER_HP_MULTIPLIERS],
        "lower": [[1.0, 1.0, 1.0] for _ in LOWER_HP_MULTIPLIERS],
    }
    effective = {
        "upper": [list(row) for row in nominal["upper"]],
        "lower": [list(row) for row in nominal["lower"]],
    }
    for item in signature:
        effective[str(item["ring"])][int(item["sector_index"])][int(item["node_index"])] = float(
            item["effective_transfer_fraction"]
        )
    return nominal, effective


def failure_signature(
        failure_map: Mapping[str, Mapping[int, Mapping[int, float]]]) -> list[Dict[str, Any]]:
    """Return a JSON-safe, auditable map of affected nominal heat pipes."""
    nominal = {"upper": UPPER_HP_MULTIPLIERS, "lower": LOWER_HP_MULTIPLIERS}
    result: list[Dict[str, Any]] = []
    for ring, sectors in failure_map.items():
        ring_key = str(ring).lower()
        if ring_key not in nominal:
            raise ValueError(f"unknown heat-pipe ring {ring!r}")
        for sector, nodes in sectors.items():
            sector_index = int(sector)
            if not 0 <= sector_index < len(SEGMENT_PATH):
                raise ValueError(f"sector index out of range: {sector_index}")
            for node, factor in nodes.items():
                node_index = int(node)
                if not 0 <= node_index < 3:
                    raise ValueError(f"node index out of range: {node_index}")
                transfer_fraction = float(factor)
                if not math.isfinite(transfer_fraction) or not 0.0 <= transfer_fraction <= 1.0:
                    raise ValueError("heat-pipe transfer fraction must be in [0, 1]")
                result.append({
                    "ring": ring_key,
                    "sector_index": sector_index,
                    "sector_name": SEGMENT_PATH[sector_index][0],
                    "node_index": node_index,
                    "nominal_heatpipe_count": int(nominal[ring_key][sector_index][node_index]),
                    "effective_transfer_fraction": transfer_fraction,
                    "hydraulic_resistance_preserved": True,
                })
    if not result:
        raise ValueError("failure_map must contain at least one affected node")
    return result


def apply_heatpipe_failure(
        build: Dict[str, Any],
        failure_map: Mapping[str, Mapping[int, Mapping[int, float]]],
) -> list[Dict[str, Any]]:
    """Scale only evaporator coupling; leave correlation and hydraulic losses unchanged."""
    signature = failure_signature(failure_map)
    all_ring_hps = list(build.get("ring_hps", []))
    ring_sets = {
        "upper": list(build.get("upper_ring_hps", all_ring_hps[:6])),
        "lower": list(build.get("lower_ring_hps", all_ring_hps[6:])),
    }
    for item in signature:
        ring_hp = ring_sets[item["ring"]][item["sector_index"]]
        node_index = item["node_index"]
        couplers = list(getattr(ring_hp, "coupler_hps", []))
        if len(couplers) != 3:
            raise RuntimeError(
                f"{getattr(ring_hp, 'name', ring_hp)!r} does not expose three node couplers"
            )
        coupler = couplers[node_index]
        factor = float(item["effective_transfer_fraction"])
        if not hasattr(coupler, "coupling_multiplier"):
            raise RuntimeError(
                f"{coupler!r} does not expose coupling_multiplier; "
                "use the heat-pipe transfer-factor implementation first"
            )
        coupler.coupling_multiplier = factor
        coupler.heatpipe_failure_transfer_fraction = factor
        coupler.heatpipe_failure_signature = dict(item)
        item["coupler_name"] = str(getattr(coupler, "name", ""))
        item["ring_hp_name"] = str(getattr(ring_hp, "name", ""))
        item["coupling_multiplier"] = factor
        item["coupling_disabled"] = factor == 0.0
    return signature


def _validate(config: HeatPipeFailureRunConfig) -> None:
    if not str(config.failure_mode).strip():
        raise ValueError("failure_mode must not be empty")
    if not Path(config.restart_in).is_file():
        raise FileNotFoundError(f"restart not found: {config.restart_in}")
    for name in (
        "minimum_duration_s", "post_scram_duration_s", "dt_s",
        "history_interval_s", "checkpoint_interval_s", "wall_limit_k",
        "fuel_limit_k", "collector_limit_k", "moderator_limit_k",
        "reflector_limit_k",
    ):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    failure_signature(config.failure_map)


def _required_end_elapsed(
        minimum_duration_s: float,
        scram_elapsed_s: Optional[float],
        post_scram_duration_s: float) -> float:
    if scram_elapsed_s is None:
        return float(minimum_duration_s)
    return max(float(minimum_duration_s), float(scram_elapsed_s) + float(post_scram_duration_s))


def run_accident(config: HeatPipeFailureRunConfig) -> Dict[str, Any]:
    _validate(config)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / "history.csv"
    if history_path.exists():
        raise FileExistsError(f"output history already exists: {history_path}")

    source_path = Path(config.restart_in).parent / "run_config.json"
    if not source_path.is_file():
        raise FileNotFoundError(f"run_config.json not found beside restart: {source_path}")
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if not bool(source.get("external_heat_enabled", False)):
        raise ValueError("source restart must retain orbital external heat")
    if bool(source.get("point_kinetics_enabled", False)):
        raise ValueError("source restart must be fixed-power without point kinetics")
    if not bool(source.get("tec_electrical_calculation_enabled", False)):
        raise ValueError("heat-pipe accident requires the normal TEC calculation")

    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in), output_dir=out_dir,
        duration_s=float(config.minimum_duration_s), dt_s=float(config.dt_s),
        record_interval_s=float(config.history_interval_s),
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=None, external_heat_enabled=True,
        external_heat_period_s=float(source["external_heat_period_s"]),
        external_heat_time_origin_s=float(source["external_heat_time_origin_s"]),
    )
    debug, source = load_baseline_debug_config(runtime)
    build = build_debug_case(debug, apply_fixed_power=True)
    system, core = build["system"], build["core"]
    if not bool(build.get("external_heat_enabled", False)):
        raise RuntimeError("external heat boundary was not attached")
    if build.get("radiator_thermal_shield") is not None:
        raise RuntimeError("heat-pipe accident must run without a thermal shield")
    if bool(getattr(core, "has_point_reactor", False)):
        raise RuntimeError("source build unexpectedly has point kinetics enabled")

    signature = apply_heatpipe_failure(build, config.failure_map)
    nominal_transfer_map, effective_transfer_map = _thermal_transfer_maps(signature)
    start_time = float(system.global_time)
    refresh_boundary_cache = getattr(system, "_refresh_solid_boundary_cache", None)
    if callable(refresh_boundary_cache):
        refresh_boundary_cache(
            update_flux=True,
            current_time=float(system.global_time),
        )
    limits = _limits(config)
    accident: Dict[str, Any] = {
        "reference_fluid": {
            "T": np.asarray(system.fluid_solver.T_vec).copy(),
            "P": np.asarray(system.fluid_solver.P_vec).copy(),
            "h": np.asarray(system.fluid_solver.h_vec).copy(),
            "W": np.asarray(system.fluid_solver.W_vec).copy(),
        },
        "feedback_reference_total": float(core.compute_reactivity_feedback().total),
        "heatpipe_failure_active": True,
        "failure_mode": str(config.failure_mode),
        "failure_signature": signature,
        "nominal_thermal_transfer_map": nominal_transfer_map,
        "effective_thermal_transfer_map": effective_transfer_map,
    }
    state: Dict[str, Any] = {
        "fixed_power_active": True,
        "scram_active": False,
        "scram_time_absolute_s": float("nan"),
        "scram_elapsed_s": float("nan"),
        "scram_trigger": None,
    }
    event = {
        "case": config.case_name,
        "case_mode": "fixed_power_then_scram",
        "accident_model": "V14_heatpipe_transfer_failure",
        "accident_time_absolute_s": start_time,
        "failure_mode": str(config.failure_mode),
        "failure_signature": signature,
        "hydraulic_resistance_preserved": True,
        "nominal_upper_heatpipe_count": int(sum(sum(v) for v in UPPER_HP_MULTIPLIERS)),
        "nominal_lower_heatpipe_count": int(sum(sum(v) for v in LOWER_HP_MULTIPLIERS)),
        "thermal_shield_enabled": False,
        "external_heat_enabled": True,
        "external_heat_period_s": float(source["external_heat_period_s"]),
        "external_heat_time_origin_s": float(source["external_heat_time_origin_s"]),
        "temperature_limits_K": limits,
        "coolant_temperature_is_trip": False,
    }
    _write_json(out_dir / "accident_event.json", event)
    run_config = dict(source)
    run_config.update({
        **event,
        "restart_in": str(config.restart_in),
        "source_run_config": str(source_path),
        "point_kinetics_enabled": False,
        "fixed_power_control_active": True,
        "power_w": FIXED_POWER_W,
        "tec_electrical_calculation_enabled": True,
        "heatpipe_failure_active": True,
        "failure_signature": signature,
        "failure_map": _json_failure_map(config.failure_map),
        "nominal_thermal_transfer_map": nominal_transfer_map,
        "effective_thermal_transfer_map": effective_transfer_map,
        "history_interval_s": float(config.history_interval_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "history_schedule": [[0.0, None, float(config.history_interval_s)]],
        "checkpoint_schedule": {
            "periodic_interval_s": float(config.checkpoint_interval_s),
            "accident_start": "accident_start_restart.npz",
            "scram_event": "scram_restart.npz",
            "final": "final_restart.npz",
        },
        "temperature_limits_K": limits,
        "coolant_temperature_is_trip": False,
        "thermal_shield_enabled": False,
        "external_heat_enabled": True,
        "scram_reactivity_dollars": SCRAM_DOLLARS,
        "minimum_duration_s": float(config.minimum_duration_s),
        "post_scram_duration_s": float(config.post_scram_duration_s),
        "dt_s": float(config.dt_s),
    })
    _write_json(out_dir / "run_config.json", run_config)

    # This is deliberately separate from periodic checkpoints and the scram
    # event checkpoint so every run has an unambiguous accident start state.
    start_restart = out_dir / "accident_start_restart.npz"
    system.save_global_state(str(start_restart))
    latest_checkpoint = str(start_restart)

    rows: list[Dict[str, Any]] = []

    def external_reactivity() -> tuple[float, float]:
        if not state["scram_active"]:
            return 0.0, 0.0
        return SCRAM_DOLLARS * float(core.point_reactor.beta_total), SCRAM_DOLLARS

    def collect_row(dt_s: float, stop_reason: str = ""):
        elapsed = float(system.global_time) - start_time
        rho, dollars = external_reactivity()
        base = collect_metrics(build, stage_index=1, dt_s=dt_s)
        peaks = collect_temperature_peaks(core, limits)
        feedback = core.compute_reactivity_feedback()
        trigger = state["scram_trigger"] or {}
        neutronics = _neutronics_metrics(core, rho, dollars)
        row = {
            **base,
            **neutronics,
            **_temperature_metrics(peaks),
            "time_s": float(system.global_time),
            "accident_elapsed_s": elapsed,
            "case_mode": "fixed_power_then_scram",
            "failure_mode": str(config.failure_mode),
            "external_heat_enabled": True,
            "external_heat_period_s": float(source["external_heat_period_s"]),
            "external_heat_phase_s": (
                (float(system.global_time) - float(source["external_heat_time_origin_s"]))
                % float(source["external_heat_period_s"])
            ),
            "tec_electrical_calculation_enabled": True,
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
            "feedback_fuel": float(feedback.fuel),
            "feedback_electrode": float(feedback.electrode),
            "feedback_moderator": float(feedback.moderator),
            "feedback_reflector": float(feedback.reflector),
            "feedback_total_absolute": float(feedback.total),
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
            "case": config.case_name,
            "latest_checkpoint_path": latest_checkpoint,
            "latest_metrics": row,
        })
        print(
            f"[{config.failure_mode} +{row['accident_elapsed_s']:.2f}s] "
            f"P={row['core_power_W']:.2f}W Twall={row['channel_wall_max_T_K']:.3f}K "
            f"Tcollector={row['collector_max_T_K']:.3f}K scram={row['scram_active']}",
            flush=True,
        )

    def trigger_scram(trigger: Dict[str, Any]) -> None:
        nonlocal latest_checkpoint
        elapsed = float(system.global_time) - start_time
        core.initialize_point_reactor(total_power_initial=FIXED_POWER_W)
        state.update({
            "fixed_power_active": False,
            "scram_active": True,
            "scram_time_absolute_s": float(system.global_time),
            "scram_elapsed_s": elapsed,
            "scram_trigger": dict(trigger),
        })
        run_config.update({
            "point_kinetics_enabled": True,
            "fixed_power_control_active": False,
            "scram_active": True,
            "scram_time_absolute_s": float(system.global_time),
            "scram_elapsed_s": elapsed,
            "scram_trigger": dict(trigger),
        })
        _write_json(out_dir / "run_config.json", run_config)
        path = out_dir / "scram_restart.npz"
        system.save_global_state(str(path))
        latest_checkpoint = str(path)
        _write_json(out_dir / "scram_event.json", {
            **trigger,
            "time_s": float(system.global_time),
            "accident_elapsed_s": elapsed,
            "external_reactivity_dollars": SCRAM_DOLLARS,
            "restart_path": str(path),
        })

    next_record = start_time + float(config.history_interval_s)
    next_checkpoint = start_time + float(config.checkpoint_interval_s)
    row, peaks, rho, dollars = collect_row(0.0)
    nonfinite = find_nonfinite_model_state(build)
    if nonfinite is not None:
        raise RuntimeError(f"initial model state is nonfinite: {nonfinite}")
    initial_trip = _temperature_trip(peaks)
    if initial_trip is not None:
        trigger_scram(initial_trip)
        next_record = float(system.global_time) + float(config.history_interval_s)
        row, peaks, rho, dollars = collect_row(0.0)
    record(row, rho, dollars)
    last_record_time = float(system.global_time)
    stop_reason = "completed"

    try:
        while True:
            end_elapsed = _required_end_elapsed(
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
                dt,
                inner_iter=int(debug.inner_iter),
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
                    next_record = float(system.global_time) + float(config.history_interval_s)
                    row, peaks, rho, dollars = collect_row(dt)

            now = float(system.global_time)
            scheduled = now >= next_record - 1.0e-9
            terminal = now >= end_time - 1.0e-9 or stop_reason != "completed"
            if trip_now or scheduled or terminal:
                record(row, rho, dollars)
                last_record_time = now
                if scheduled and not trip_now:
                    next_record = now + float(config.history_interval_s)

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
            "type": type(exc).__name__,
            "message": str(exc),
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
        "case": config.case_name,
        "failure_mode": config.failure_mode,
        "output_dir": str(out_dir),
        "source_restart_path": str(config.restart_in),
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "minimum_duration_s": float(config.minimum_duration_s),
        "post_scram_duration_s": float(config.post_scram_duration_s),
        "scram_active": bool(state["scram_active"]),
        "scram_time_absolute_s": float(state["scram_time_absolute_s"]),
        "scram_elapsed_s": float(state["scram_elapsed_s"]),
        "stop_reason": stop_reason,
        "final_restart_path": str(final_restart),
        "latest_metrics": rows[-1] if rows else {},
        "final_window_slopes": _slope_summary(rows),
        "failure_signature": signature,
    }
    _write_json(out_dir / "run_summary.json", result)
    run_config.update({
        "point_kinetics_enabled": bool(core.has_point_reactor),
        "fixed_power_control_active": bool(state["fixed_power_active"]),
        "scram_active": bool(state["scram_active"]),
        "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "final_restart_path": str(final_restart),
    })
    _write_json(out_dir / "run_config.json", run_config)
    _write_json(out_dir / "latest_state.json", {
        "case": config.case_name,
        "latest_checkpoint_path": str(final_restart),
        "latest_metrics": rows[-1] if rows else {},
    })
    return result


def run_cli(
        *, case_name: str, failure_mode: str,
        failure_map: Mapping[str, Mapping[int, Mapping[int, float]]],
        default_output_dir: Path,
        argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=default_output_dir)
    parser.add_argument("--duration", type=float, default=ORBIT_PERIOD_S)
    parser.add_argument("--post-scram-duration", type=float, default=HALF_ORBIT_S)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--history-interval", type=float, default=1.0)
    parser.add_argument("--checkpoint-interval", type=float, default=100.0)
    parser.add_argument("--wall-limit", type=float, default=1058.0)
    parser.add_argument("--fuel-limit", type=float, default=2700.0)
    parser.add_argument("--collector-limit", type=float, default=1023.0)
    parser.add_argument("--moderator-limit", type=float, default=930.0)
    parser.add_argument("--reflector-limit", type=float, default=1000.0)
    args = parser.parse_args(argv)
    result = run_accident(HeatPipeFailureRunConfig(
        case_name=case_name,
        failure_mode=failure_mode,
        failure_map=failure_map,
        output_dir=args.output_dir,
        restart_in=args.restart_in,
        minimum_duration_s=float(args.duration),
        post_scram_duration_s=float(args.post_scram_duration),
        dt_s=float(args.dt),
        history_interval_s=float(args.history_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        wall_limit_k=float(args.wall_limit),
        fuel_limit_k=float(args.fuel_limit),
        collector_limit_k=float(args.collector_limit),
        moderator_limit_k=float(args.moderator_limit),
        reflector_limit_k=float(args.reflector_limit),
    ))
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0
