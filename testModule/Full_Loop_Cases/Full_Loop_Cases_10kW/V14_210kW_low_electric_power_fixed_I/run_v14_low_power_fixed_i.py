"""Run V14 from an orbital checkpoint with the main TEC at fixed current."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import traceback
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    DebugRunConfig,
    _apply_fixed_core_power,
    _apply_wire_resistance,
    build_debug_case,
    collect_metrics,
)


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_RESTART = (
    CASE_DIR.parent / "V14_210kW_fixed_power_external_heat_2orbits" / "runs"
    / "two_orbits_from13864_20260720" / "checkpoint_t019265s.npz"
)
DEFAULT_OUTPUT = CASE_DIR / "runs" / "stage_b_fixed_i_smoke_0p1s"


@dataclass(frozen=True)
class LowPowerRunConfig:
    restart_in: Path = DEFAULT_RESTART
    output_dir: Path = DEFAULT_OUTPUT
    duration_s: float = 0.1
    dt_s: float = 0.05
    tec_update_interval_s: float = 0.5
    record_interval_s: float = 0.05
    hold_before_ramp_s: float = 0.0
    ramp_duration_s: float = 0.0
    ramp_shape: str = "cubic"
    final_power_w: float = 210000.0
    final_flow_kg_s: float = 2.46
    fluid_max_iter: int = 100
    checkpoint_interval_s: float = 60.0
    staged_recording: bool = False
    resume_from: Optional[Path] = None
    fixed_current_a: Optional[float] = None
    external_heat_enabled: bool = True


def validate_run_config(config: LowPowerRunConfig) -> None:
    positive = {
        "duration_s": config.duration_s,
        "dt_s": config.dt_s,
        "tec_update_interval_s": config.tec_update_interval_s,
        "record_interval_s": config.record_interval_s,
        "final_power_w": config.final_power_w,
        "final_flow_kg_s": config.final_flow_kg_s,
    }
    for name, raw in positive.items():
        value = float(raw)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    for name in ("hold_before_ramp_s", "ramp_duration_s", "checkpoint_interval_s"):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name} must be finite and non-negative")
    if int(config.fluid_max_iter) <= 0:
        raise ValueError("fluid_max_iter must be positive")
    if config.fixed_current_a is not None:
        current = float(config.fixed_current_a)
        if not math.isfinite(current) or current <= 0.0:
            raise ValueError("fixed_current_a must be finite and positive")
    if config.ramp_shape not in ("cubic", "quintic"):
        raise ValueError("ramp_shape must be cubic or quintic")


def make_manifest(
    *, config: LowPowerRunConfig, source_config: Dict[str, Any],
    candidate_start_time_s: float, initial_current_a: float,
    initial_electric_power_w: float, initial_outlet_k: float,
    initial_power_w: float, initial_flow_kg_s: float,
) -> Dict[str, Any]:
    return {
        "source_build_config": dict(source_config),
        "candidate_start_time_s": float(candidate_start_time_s),
        "baseline": {
            "I0_A": float(initial_current_a),
            "Pe0_W": float(initial_electric_power_w),
            "Tout0_K": float(initial_outlet_k),
            "Q0_W": float(initial_power_w),
            "W0_kg_s": float(initial_flow_kg_s),
        },
        "trajectory": {
            "hold_before_ramp_s": float(config.hold_before_ramp_s),
            "ramp_duration_s": float(config.ramp_duration_s),
            "ramp_shape": config.ramp_shape,
            "final_power_w": float(config.final_power_w),
            "final_flow_kg_s": float(config.final_flow_kg_s),
            "duration_s": float(config.duration_s),
            "dt_s": float(config.dt_s),
            "tec_update_interval_s": float(config.tec_update_interval_s),
            "fixed_current_a": (
                None if config.fixed_current_a is None
                else float(config.fixed_current_a)
            ),
            "external_heat_enabled": bool(config.external_heat_enabled),
        },
    }

def restore_trajectory_config(
    config: LowPowerRunConfig, manifest: Dict[str, Any],
) -> LowPowerRunConfig:
    trajectory = manifest["trajectory"]
    restored = {
        name: trajectory[name] for name in (
            "hold_before_ramp_s", "ramp_duration_s", "final_power_w",
            "final_flow_kg_s", "duration_s", "dt_s",
        )
    }
    if "tec_update_interval_s" in trajectory:
        restored["tec_update_interval_s"] = trajectory["tec_update_interval_s"]
    if "ramp_shape" in trajectory:
        restored["ramp_shape"] = trajectory["ramp_shape"]
    for name in ("fixed_current_a", "external_heat_enabled"):
        if name in trajectory:
            restored[name] = trajectory[name]
    return replace(config, **restored)

def load_source_config(restart: Path) -> Dict[str, Any]:
    manifest_path = Path(restart).parent / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if "source_build_config" in manifest:
            return dict(manifest["source_build_config"])
    return json.loads(
        (Path(restart).parent / "run_config.json").read_text(encoding="utf-8"))

def smooth_ramp(time_s: float, start_s: float, duration_s: float,
                initial_value: float, final_value: float,
                *, shape: str = "cubic") -> float:
    if duration_s <= 0.0:
        return float(initial_value if time_s <= start_s else final_value)
    fraction = min(1.0, max(0.0, (float(time_s) - float(start_s)) / float(duration_s)))
    if shape == "quintic":
        weight = fraction**3 * (10.0 + fraction * (-15.0 + 6.0 * fraction))
    elif shape == "cubic":
        weight = fraction * fraction * (3.0 - 2.0 * fraction)
    else:
        raise ValueError("shape must be cubic or quintic")
    return float(initial_value) + weight * (float(final_value) - float(initial_value))


def control_setpoints_at_step_end(
    *, step_end_elapsed_s: float, hold_before_ramp_s: float, ramp_duration_s: float,
    ramp_shape: str = "cubic",
    initial_power_w: float, final_power_w: float,
    initial_flow_kg_s: float, final_flow_kg_s: float,
) -> tuple[float, float]:
    return (
        smooth_ramp(step_end_elapsed_s, hold_before_ramp_s, ramp_duration_s,
                    initial_power_w, final_power_w, shape=ramp_shape),
        smooth_ramp(step_end_elapsed_s, hold_before_ramp_s, ramp_duration_s,
                    initial_flow_kg_s, final_flow_kg_s, shape=ramp_shape),
    )


def record_interval_for_elapsed(
    elapsed_s: float, *, staged: bool, default_s: float = 0.5,
) -> float:
    if not staged:
        return float(default_s)
    elapsed = float(elapsed_s)
    if elapsed <= 20.0:
        return 0.5
    if elapsed <= 100.0:
        return 2.0
    if elapsed <= 400.0:
        return 5.0
    if elapsed <= 600.0:
        return 10.0
    return 20.0


def interval_due(elapsed_s: float, last_s: float, interval_s: float) -> bool:
    return float(elapsed_s) - float(last_s) >= float(interval_s) - 1.0e-9


def record_due(
    *, elapsed_s: float, last_record_s: float, interval_s: float,
    duration_s: float, stopped: bool,
) -> bool:
    tolerance = 1.0e-9
    return (
        interval_due(elapsed_s, last_record_s, interval_s)
        or float(elapsed_s) >= float(duration_s) - tolerance
        or bool(stopped)
    )


def resolve_history_fields(
    path: Path, proposed: Sequence[str], *, resuming: bool,
) -> tuple[str, ...]:
    fields = tuple(proposed)
    if not resuming:
        return fields
    if not path.is_file():
        raise RuntimeError(f"resume history is missing: {path}")
    with path.open(newline="", encoding="utf-8") as stream:
        existing = tuple(next(csv.reader(stream), ()))
    if existing != fields:
        raise RuntimeError(
            f"resume history schema mismatch: existing={existing}, current={fields}")
    return existing

def rewind_history_for_resume(path: Path, checkpoint_elapsed_s: float) -> None:
    history_path = Path(path)
    if not history_path.is_file():
        raise RuntimeError(f"resume history is missing: {history_path}")
    with history_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fields = tuple(reader.fieldnames or ())
        if "elapsed_s" not in fields:
            raise RuntimeError("resume history has no elapsed_s column")
        rows = list(reader)
    cutoff = float(checkpoint_elapsed_s) - 1.0e-9
    retained = []
    for row in rows:
        try:
            elapsed = float(row["elapsed_s"])
        except (TypeError, ValueError, OverflowError) as exc:
            raise RuntimeError("resume history has an invalid elapsed_s value") from exc
        if elapsed < cutoff:
            retained.append(row)
    temporary = history_path.with_name(f"{history_path.stem}.rewind.tmp{history_path.suffix}")
    try:
        with temporary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(retained)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, history_path)
    finally:
        if temporary.exists():
            temporary.unlink()


def append_history_row(
    path: Path, row: Dict[str, Any], fieldnames: Sequence[str],
) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


def save_checkpoint_atomic(system: Any, path: Path) -> None:
    target = Path(path)
    temporary = target.with_name(f"{target.stem}.tmp{target.suffix}")
    try:
        system.save_global_state(str(temporary))
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def outlet_limit_is_active(elapsed_s: float, grace_s: float = 1.0) -> bool:
    return float(elapsed_s) >= float(grace_s) - 1.0e-9


def evaluate_hard_trip(
    metrics: Dict[str, Any], *, initial_outlet_k: float,
    enforce_outlet_limit: bool = True,
) -> Optional[Dict[str, Any]]:
    for field, raw in metrics.items():
        if isinstance(raw, (str, bool, np.bool_)) or not np.isscalar(raw):
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError, OverflowError):
            continue
        if not math.isfinite(value):
            return {"stop_reason": "nonfinite_metric", "field": field, "actual": value}
    if not bool(metrics.get("fluid_converged", False)):
        return {"stop_reason": "hydraulic_nonconvergence"}
    if not bool(metrics.get("tec_main_converged", False)):
        return {"stop_reason": "tec_nonconvergence"}
    limits = {
        "fuel_max_T_K": ("fuel_temperature_limit", 2700.0),
        "collector_max_T_K": ("collector_temperature_limit", 1500.0),
        "emitter_max_T_K": ("emitter_temperature_limit", 3000.0),
        "coolant_max_T_K": ("coolant_temperature_limit", 1058.0),
        "moderator_max_T_K": ("moderator_temperature_limit", 930.0),
        "reflector_max_T_K": ("reflector_temperature_limit", 1000.0),
    }
    for field, (reason, limit) in limits.items():
        if float(metrics[field]) >= limit:
            return {"stop_reason": reason, "field": field,
                    "actual": float(metrics[field]), "limit": limit}
    outlet_limit = float(initial_outlet_k) + 0.5
    if enforce_outlet_limit and float(metrics["core_outlet_T_K"]) > outlet_limit:
        return {"stop_reason": "core_outlet_temperature_limit",
                "field": "core_outlet_T_K",
                "actual": float(metrics["core_outlet_T_K"]), "limit": outlet_limit}
    return None


def exit_code_for_result(result: Dict[str, Any]) -> int:
    return 0 if result.get("stop_reason") == "completed" else 1


def write_failure_artifacts(
        out_dir: Path, exc: BaseException, *, system: Optional[Any] = None) -> None:
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    restart_resumable = system is not None and (directory / "run_manifest.json").is_file()
    payload = {
        "exception_type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
        "restart_resumable": restart_resumable,
    }
    if restart_resumable:
        try:
            save_checkpoint_atomic(system, directory / "emergency_restart.npz")
        except Exception:
            payload["restart_resumable"] = False
            payload["emergency_restart_error"] = traceback.format_exc()
    _write_json(directory / "failure.json", payload)


def validate_final_setpoints(
        build: Dict[str, Any], core: Any, power_w: float, flow_kg_s: float) -> None:
    controlled_pumps = 0
    for key in ("pump_a", "pump_b"):
        pump = build.get(key)
        if pump is None or (
                not callable(getattr(pump, "set_flow_rate", None))
                and not hasattr(pump, "target_W")):
            continue
        controlled_pumps += 1
        actual = float(getattr(pump, "target_W", float("nan")))
        if not math.isclose(actual, float(flow_kg_s), rel_tol=0.0, abs_tol=1.0e-12):
            raise RuntimeError(f"{key} did not retain the exact final flow setpoint")
    if controlled_pumps == 0:
        raise RuntimeError("no pump supports prescribed-flow control")
    actual_power = float(getattr(core, "last_total_core_power", float("nan")))
    if not math.isclose(actual_power, float(power_w), rel_tol=1.0e-12, abs_tol=1.0e-6):
        raise RuntimeError("failed to apply exact final thermal-power setpoint")


def set_total_flow_target(build: Dict[str, Any], flow_kg_s: float) -> None:
    value = float(flow_kg_s)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("total flow target must be finite and positive")
    controlled_pumps = 0
    for key in ("pump_a", "pump_b"):
        pump = build.get(key)
        setter = getattr(pump, "set_flow_rate", None)
        if callable(setter):
            setter(value)
            controlled_pumps += 1
        elif pump is not None and hasattr(pump, "target_W"):
            pump.target_W = value
            controlled_pumps += 1
    if controlled_pumps == 0:
        raise RuntimeError("no pump supports prescribed-flow control")
    build["pump_target_flow_kg_s"] = value


def validate_fixed_i_result(metrics: Dict[str, Any], *, target_current_a: float,
                            relative_tolerance: float = 1.0e-6) -> None:
    current = float(metrics["tec_main_current_A"])
    voltage = float(metrics["tec_main_voltage_V"])
    target = float(target_current_a)
    if not bool(metrics["tec_main_converged"]):
        raise RuntimeError("TEC fixed-I solve did not converge")
    if not all(math.isfinite(value) for value in (current, voltage, target)):
        raise RuntimeError("TEC fixed-I result is non-finite")
    tolerance = max(1.0e-6, abs(target) * float(relative_tolerance))
    if abs(current - target) > tolerance:
        raise RuntimeError(
            f"TEC fixed-I current mismatch: target={target:.9g} A, actual={current:.9g} A"
        )
    if target > 0.0 and voltage <= 0.0:
        raise RuntimeError("TEC fixed-I state has no positive generated voltage")


def _refresh_tec(core: Any, current_time_s: float) -> None:
    for group in core.iter_tec_circuit_groups():
        if group.thermo_calc is None:
            continue
        core._sync_tec_group_temperatures(group)
        group.thermo_calc.calculate(verbose=False)
        core._apply_tec_group_results(group)
        group.last_update_time = float(current_time_s)
        if group.name == "main":
            core._last_thermo_update_time = float(current_time_s)


def build_case_for_run(debug: DebugRunConfig, *, resuming: bool) -> Dict[str, Any]:
    return build_debug_case(debug, apply_fixed_power=not resuming)

def prepare_candidate_state(
    build: Dict[str, Any], source: Dict[str, Any], config: LowPowerRunConfig,
    *, current_time_s: float, resume_manifest: Dict[str, Any],
) -> Dict[str, Any]:
    """Restore the saved fixed-I contract without recalculating a fixed-U baseline."""
    baseline = resume_manifest["baseline"]
    initial_current = float(baseline["I0_A"])
    candidate_start = float(resume_manifest["candidate_start_time_s"])
    elapsed = float(current_time_s) - candidate_start
    q_resume, w_resume = control_setpoints_at_step_end(
        step_end_elapsed_s=elapsed,
        hold_before_ramp_s=config.hold_before_ramp_s,
        ramp_duration_s=config.ramp_duration_s,
        ramp_shape=config.ramp_shape,
        initial_power_w=float(baseline["Q0_W"]),
        final_power_w=config.final_power_w,
        initial_flow_kg_s=float(baseline["W0_kg_s"]),
        final_flow_kg_s=config.final_flow_kg_s,
    )
    _apply_fixed_core_power(build, q_resume)
    set_total_flow_target(build, w_resume)
    core = build["core"]
    core.setup_tec_circuit(
        "fixed_i", initial_current, I_guess=initial_current, topology="series")
    _apply_wire_resistance(core, float(source["wire_resistance_scale"]))
    _refresh_tec(core, current_time_s)
    core.thermo_update_interval = float(config.tec_update_interval_s)
    fixed_i = collect_metrics(build, stage_index=0, dt_s=0.0)
    validate_fixed_i_result(fixed_i, target_current_a=initial_current)
    return {
        "fixed_i_metrics": fixed_i,
        "initial_current_A": initial_current,
        "initial_electric_power_W": float(baseline["Pe0_W"]),
        "initial_outlet_T_K": float(baseline["Tout0_K"]),
        "initial_power_W": float(baseline["Q0_W"]),
        "initial_flow_kg_s": float(baseline["W0_kg_s"]),
        "candidate_start_time_s": candidate_start,
        "elapsed_s": elapsed,
        "resume_power_setpoint_W": q_resume,
        "resume_flow_setpoint_kg_s": w_resume,
    }

def _load_debug_config(
    config: LowPowerRunConfig, *, restart_override: Optional[Path] = None,
    source_override: Optional[Dict[str, Any]] = None,
) -> tuple[DebugRunConfig, Dict[str, Any]]:
    restart = Path(restart_override if restart_override is not None else config.restart_in)
    if not restart.is_file():
        raise FileNotFoundError(restart)
    source = (
        dict(source_override) if source_override is not None
        else load_source_config(restart)
    )
    debug = DebugRunConfig(
        output_dir=Path(config.output_dir), stage_durations_s=(float(config.duration_s),),
        dt_s=float(config.dt_s), record_interval_s=float(config.record_interval_s),
        checkpoint_interval_s=0.0, min_fluid_temperature_stop_k=None,
        tec_electrical_enabled=True, tec_voltage_v=float(source["tec_voltage_v"]),
        tec_current_guess_a=float(source["tec_current_guess_a"]),
        tec_lookup_enabled=bool(source["tec_lookup_enabled"]),
        tec_lookup_db=source.get("tec_lookup_db"),
        tec_lookup_regions=tuple(source["tec_lookup_regions"]),
        wire_resistance_scale=float(source["wire_resistance_scale"]),
        radiator_emissivity=float(source["radiator_emissivity"]),
        hp_up_view_factor=float(source["hp_up_view_factor"]),
        upper_hp_down_view_factor=float(source["upper_hp_down_view_factor"]),
        lower_hp_down_view_factor=float(source["lower_hp_down_view_factor"]),
        inner_iter=int(source["inner_iter"]), fluid_max_iter=int(config.fluid_max_iter),
        power_w=float(source["power_w"]), target_flow_kg_s=float(source["target_flow_kg_s"]),
        initial_temperature_k=float(source["initial_temperature_k"]),
        space_temperature_k=float(source["space_temperature_k"]),
        external_heat_enabled=bool(config.external_heat_enabled),
        external_heat_period_s=float(source["external_heat_period_s"]),
        external_heat_time_origin_s=float(source["external_heat_time_origin_s"]),
        restart_in=restart, case_prefix="V14_210kW_low_electric_power_fixed_I",
    )
    return debug, source


def _collect_safety_metrics(build: Dict[str, Any]) -> Dict[str, float]:
    grouped = {
        "fuel": [], "collector": [], "emitter": [],
        "moderator": [], "reflector": [], "inner_clad": [], "outer_clad": [],
    }
    core = build["core"]
    for tfe in core.tfes.values():
        solids = tfe.solids
        grouped["fuel"].append(np.asarray(solids["pellet"].T, dtype=float))
        grouped["collector"].append(np.asarray(solids["collector"].T, dtype=float))
        grouped["emitter"].append(np.asarray(solids["emitter"].T, dtype=float))
        grouped["inner_clad"].append(
            np.asarray(solids["inner_clad"].T, dtype=float))
        grouped["outer_clad"].append(
            np.asarray(solids["outer_clad"].T, dtype=float))
        if "moderator" in solids:
            grouped["moderator"].append(np.asarray(solids["moderator"].T, dtype=float))
    grouped["moderator"].extend(
        np.asarray(solid.T, dtype=float) for solid in core.mod_rings)
    grouped["reflector"].append(np.asarray(core.reflector.T, dtype=float))
    result = {
        f"{name}_max_T_K": max(float(np.max(values)) for values in arrays)
        for name, arrays in grouped.items()
    }
    result["coolant_max_T_K"] = float(np.max(np.asarray(
        build["system"].fluid_solver.T_vec, dtype=float)))
    return result




def _weighted_stats(prefix: str, values: list[np.ndarray],
                    weights: list[np.ndarray]) -> Dict[str, float]:
    if not values:
        return {f"{prefix}_{stat}": float("nan") for stat in ("min_K", "mean_K", "max_K")}
    data = np.concatenate([np.asarray(value, dtype=float).ravel() for value in values])
    weight = np.concatenate([np.asarray(value, dtype=float).ravel() for value in weights])
    valid = np.isfinite(data) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return {f"{prefix}_{stat}": float("nan") for stat in ("min_K", "mean_K", "max_K")}
    data, weight = data[valid], weight[valid]
    return {
        f"{prefix}_min_K": float(np.min(data)),
        f"{prefix}_mean_K": float(np.average(data, weights=weight)),
        f"{prefix}_max_K": float(np.max(data)),
    }


def _refresh_ring_hp_diagnostic_state(build: Dict[str, Any], current_time_s: float) -> None:
    for ring_hp in build.get("ring_hps", []):
        ring_hp.pre_step(0.0, float(current_time_s))


def _collect_heat_rejection_diagnostics(
    build: Dict[str, Any], *, current_time_s: float,
    external_heat_period_s: float, external_heat_time_origin_s: float,
) -> Dict[str, float]:
    section_values = {name: [] for name in ("evaporator", "adiabatic", "condenser")}
    section_weights = {name: [] for name in section_values}
    wick_values: list[np.ndarray] = []
    wick_weights: list[np.ndarray] = []
    condenser_wall_values: list[np.ndarray] = []
    condenser_wall_weights: list[np.ndarray] = []
    fin_values: list[np.ndarray] = []
    fin_weights: list[np.ndarray] = []
    header_values: list[np.ndarray] = []
    header_weights: list[np.ndarray] = []
    evaporator_heat = gross_rejection = absorbed_heat = wall_rejection = 0.0

    for ring_hp in build.get("ring_hps", []):
        header = ring_hp.solid_header
        header_values.append(np.asarray(header.T, dtype=float))
        header_weights.append(np.asarray(header.mesh.geom_data.volumes, dtype=float))
        outer = getattr(header, "boundaries", {}).get("right")
        if outer is not None and hasattr(outer, "current_flux"):
            wall_rejection -= float(np.sum(np.asarray(outer.current_flux, dtype=float)))

        for _, hp_unit, multiplier in ring_hp._iter_present_hp_units_with_multiplier():
            hp = hp_unit.hp
            temperature = np.asarray(hp.T, dtype=float).reshape(hp.shape_nodes)
            volume = np.asarray(hp.mesh.geom_data.volumes, dtype=float).reshape(hp.shape_nodes)
            for name, section in (
                ("evaporator", hp._slice_eva),
                ("adiabatic", hp._slice_aba),
                ("condenser", hp._slice_con),
            ):
                section_values[name].append(temperature[:, section])
                section_weights[name].append(volume[:, section] * multiplier)

            wick_temperature = temperature[:hp.n_wick, :]
            wick_values.append(np.asarray(
                hp.wick_mat.conductivity_axial(wick_temperature), dtype=float))
            wick_weights.append(volume[:hp.n_wick, :] * multiplier)
            condenser_wall_values.append(temperature[-1, hp._slice_con])
            condenser_wall_weights.append(volume[-1, hp._slice_con] * multiplier)

            fin_temperature = np.asarray(hp_unit.last_fin_temperature, dtype=float)
            fin_cell_volume = (
                np.asarray(hp_unit.fin_width_array, dtype=float)[:, None]
                * float(hp_unit.fin_thickness) * float(hp_unit.fin_height)
                / int(hp_unit.n_fin_height)
            )
            fin_values.append(fin_temperature)
            fin_weights.append(np.broadcast_to(fin_cell_volume, fin_temperature.shape) * multiplier)

            evaporator_heat += multiplier * float(np.sum(
                np.asarray(hp.boundaries["outer_eva"].current_flux, dtype=float)))
            q_aba, _ = hp_unit.get_heat_rejection_distribution()
            breakdown = hp_unit.get_heat_exchange_breakdown()
            gross_rejection += multiplier * (
                float(np.sum(q_aba)) + float(np.sum(breakdown["gross_rejection"])))
            _, _, absorbed = hp_unit.get_external_heat_absorption_distribution(current_time_s)
            absorbed_heat += multiplier * float(np.sum(absorbed))

    result: Dict[str, float] = {}
    for name in section_values:
        result.update(_weighted_stats(
            f"hp_{name}_temperature", section_values[name], section_weights[name]))
    wick_stats = _weighted_stats("hp_wick_axial_conductivity", wick_values, wick_weights)
    result.update({key.replace("_K", "_W_mK"): value for key, value in wick_stats.items()})
    result.update(_weighted_stats(
        "hp_condenser_outer_wall_temperature", condenser_wall_values,
        condenser_wall_weights))
    result.update(_weighted_stats("radiator_fin_temperature", fin_values, fin_weights))
    result.update(_weighted_stats(
        "collector_ring_wall_temperature", header_values, header_weights))
    result.update({
        "hp_evaporator_minus_condenser_mean_K": (
            result["hp_evaporator_temperature_mean_K"]
            - result["hp_condenser_temperature_mean_K"]),
        "hp_evaporator_heat_input_W": evaporator_heat,
        "radiator_gross_heat_rejection_W": gross_rejection,
        "radiator_external_heat_absorption_W": absorbed_heat,
        "radiator_net_heat_rejection_W": gross_rejection - absorbed_heat,
        "collector_ring_wall_outward_rejection_W": wall_rejection,
    })
    period = float(external_heat_period_s)
    phase = (float(current_time_s) - float(external_heat_time_origin_s)) % period
    result["external_heat_phase_s"] = phase
    result["external_heat_phase_fraction"] = phase / period
    return result

def _find_nonfinite_model_state(build: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    system = build["system"]
    arrays = [
        (f"fluid.{field}", getattr(system.fluid_solver, field))
        for field in ("T_vec", "P_vec", "h_vec", "rho_vec", "W_vec")
        if hasattr(system.fluid_solver, field)
    ]
    arrays.extend(
        (f"solid:{name}.T", solid.T)
        for name, solid in system.solid_components.items()
        if hasattr(solid, "T")
    )
    for field, raw in arrays:
        values = np.asarray(raw, dtype=float).ravel()
        bad = np.flatnonzero(~np.isfinite(values))
        if bad.size:
            index = int(bad[0])
            return {"field": field, "flat_index": index,
                    "actual": float(values[index])}
    return None

def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def _run_low_power_case_impl(
        config: LowPowerRunConfig, *, failure_context: Dict[str, Any]) -> Dict[str, Any]:
    validate_run_config(config)
    resume_manifest: Optional[Dict[str, Any]] = None
    resume_path = Path(config.resume_from) if config.resume_from is not None else None
    if resume_path is not None:
        if not resume_path.is_file():
            raise FileNotFoundError(resume_path)
        manifest_path = resume_path.parent / "run_manifest.json"
        resume_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        config = restore_trajectory_config(config, resume_manifest)
        validate_run_config(config)
        source = dict(resume_manifest["source_build_config"])
        debug, _ = _load_debug_config(
            config, restart_override=resume_path, source_override=source)
        out_dir = resume_path.parent
    else:
        debug, source = _load_debug_config(config)
        out_dir = Path(config.output_dir)

    failure_context["out_dir"] = out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / "history_control.csv"
    if history_path.exists() and resume_manifest is None:
        raise FileExistsError(f"output history already exists: {history_path}")
    _write_json(out_dir / "run_config.json", config.__dict__)

    build = build_case_for_run(debug, resuming=resume_manifest is not None)
    system, core = build["system"], build["core"]
    failure_context["system"] = system
    current_time = float(system.global_time)
    _refresh_ring_hp_diagnostic_state(build, current_time)

    if resume_manifest is not None:
        state = prepare_candidate_state(
            build, source, config, current_time_s=current_time,
            resume_manifest=resume_manifest)
        fixed_u: Optional[Dict[str, Any]] = None
        fixed_i = dict(state["fixed_i_metrics"])
        fixed_i.update(_collect_safety_metrics(build))
        initial_current = float(state["initial_current_A"])
        initial_electric_power = float(state["initial_electric_power_W"])
        initial_outlet = float(state["initial_outlet_T_K"])
        initial_thermal_power = float(state["initial_power_W"])
        initial_flow = float(state["initial_flow_kg_s"])
        candidate_start_time = float(state["candidate_start_time_s"])
        initial_elapsed = float(state["elapsed_s"])
        current_q_set = float(state["resume_power_setpoint_W"])
        current_w_set = float(state["resume_flow_setpoint_kg_s"])
    else:
        if config.fixed_current_a is None:
            _refresh_tec(core, current_time)
            fixed_u = collect_metrics(build, stage_index=0, dt_s=0.0)
            if not fixed_u["tec_main_converged"]:
                raise RuntimeError("checkpoint fixed-U refresh did not converge")
            initial_current = float(fixed_u["tec_main_current_A"])
            if initial_current <= 0.0:
                raise RuntimeError("checkpoint does not provide a positive TEC current baseline")
            initial_thermal_power = float(source["power_w"])
            initial_flow = float(source["target_flow_kg_s"])
        else:
            fixed_u = None
            initial_current = float(config.fixed_current_a)
            initial_thermal_power = float(config.final_power_w)
            initial_flow = float(config.final_flow_kg_s)
            _apply_fixed_core_power(build, initial_thermal_power)
            set_total_flow_target(build, initial_flow)
        core.setup_tec_circuit(
            "fixed_i", initial_current, I_guess=initial_current, topology="series")
        _apply_wire_resistance(core, float(source["wire_resistance_scale"]))
        _refresh_tec(core, current_time)
        core.thermo_update_interval = float(config.tec_update_interval_s)
        fixed_i = collect_metrics(build, stage_index=0, dt_s=0.0)
        fixed_i.update(_collect_safety_metrics(build))
        validate_fixed_i_result(fixed_i, target_current_a=initial_current)
        initial_electric_power = float(fixed_i["tec_main_electric_power_W"])
        if initial_electric_power <= 0.0:
            raise RuntimeError("checkpoint does not provide a positive TEC power baseline")
        initial_outlet = float(fixed_i["core_outlet_T_K"])
        candidate_start_time = current_time
        initial_elapsed = 0.0
        current_q_set = initial_thermal_power
        current_w_set = initial_flow
        manifest = make_manifest(
            config=config, source_config=source,
            candidate_start_time_s=candidate_start_time,
            initial_current_a=initial_current,
            initial_electric_power_w=initial_electric_power,
            initial_outlet_k=initial_outlet,
            initial_power_w=initial_thermal_power,
            initial_flow_kg_s=initial_flow,
        )
        _write_json(out_dir / "run_manifest.json", manifest)
    if resume_manifest is not None:
        rewind_history_for_resume(history_path, initial_elapsed)
    history_fields: Optional[Sequence[str]] = None
    latest: Dict[str, Any] = {}

    def add_row(row: Dict[str, Any], elapsed: float, q_set: float, w_set: float) -> None:
        nonlocal history_fields, latest
        latest = dict(row)
        latest.update(
            elapsed_s=float(elapsed),
            thermal_power_setpoint_W=float(q_set),
            flow_setpoint_kg_s=float(w_set),
            electric_power_ratio=(
                float(row["tec_main_electric_power_W"]) / initial_electric_power),
        )
        latest.update(_collect_heat_rejection_diagnostics(
            build,
            current_time_s=float(row.get("time_s", system.global_time)),
            external_heat_period_s=float(debug.external_heat_period_s),
            external_heat_time_origin_s=float(debug.external_heat_time_origin_s),
        ))
        if history_fields is None:
            history_fields = resolve_history_fields(
                history_path, tuple(latest), resuming=resume_manifest is not None)
        append_history_row(history_path, latest, history_fields)

    add_row(fixed_i, initial_elapsed, current_q_set, current_w_set)
    initial_nonfinite = _find_nonfinite_model_state(build)
    initial_trip = (
        {"stop_reason": "nonfinite_model_state", **initial_nonfinite}
        if initial_nonfinite is not None
        else evaluate_hard_trip(
            fixed_i, initial_outlet_k=initial_outlet,
            enforce_outlet_limit=outlet_limit_is_active(initial_elapsed),
        )
    )
    if initial_trip is not None:
        restart_path = out_dir / "emergency_restart.npz"
        save_checkpoint_atomic(system, restart_path)
        _write_json(out_dir / "limit_trip.json", initial_trip)
        summary = {
            "candidate_start_time_s": candidate_start_time,
            "end_time_s": float(system.global_time),
            "initial_current_A": initial_current,
            "initial_electric_power_W": initial_electric_power,
            "fixed_u_metrics": fixed_u,
            "fixed_i_handoff_metrics": fixed_i,
            "latest_metrics": latest,
            "stop_reason": str(initial_trip["stop_reason"]),
            "trip": initial_trip,
            "restart_path": str(restart_path),
            "final_power_setpoint_W": current_q_set,
            "final_flow_setpoint_kg_s": current_w_set,
        }
        _write_json(out_dir / "run_summary.json", summary)
        return summary

    last_record = initial_elapsed
    last_checkpoint = initial_elapsed
    end_time = candidate_start_time + float(config.duration_s)
    stop_reason = "completed"
    trip_payload: Optional[Dict[str, Any]] = None
    q_set, w_set = current_q_set, current_w_set

    while float(system.global_time) < end_time - 1.0e-12:
        dt = min(float(config.dt_s), end_time - float(system.global_time))
        step_end_elapsed = float(system.global_time) + dt - candidate_start_time
        q_set, w_set = control_setpoints_at_step_end(
            step_end_elapsed_s=step_end_elapsed,
            hold_before_ramp_s=config.hold_before_ramp_s,
            ramp_duration_s=config.ramp_duration_s,
            ramp_shape=config.ramp_shape,
            initial_power_w=initial_thermal_power,
            final_power_w=config.final_power_w,
            initial_flow_kg_s=initial_flow,
            final_flow_kg_s=config.final_flow_kg_s,
        )
        set_total_flow_target(build, w_set)
        _apply_fixed_core_power(build, q_set)
        system.step(
            dt, inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=False,
            fluid_max_iter=int(config.fluid_max_iter),
        )
        set_total_flow_target(build, w_set)
        _apply_fixed_core_power(build, q_set)
        elapsed = float(system.global_time) - candidate_start_time

        row = collect_metrics(build, stage_index=1, dt_s=dt)
        row.update(_collect_safety_metrics(build))
        nonfinite = _find_nonfinite_model_state(build)
        if nonfinite is not None:
            trip_payload = {"stop_reason": "nonfinite_model_state", **nonfinite}
        else:
            try:
                validate_fixed_i_result(row, target_current_a=initial_current)
            except RuntimeError as exc:
                trip_payload = {"stop_reason": "tec_fixed_i_invalid", "detail": str(exc)}
            if trip_payload is None:
                trip_payload = evaluate_hard_trip(
                    row, initial_outlet_k=initial_outlet,
                    enforce_outlet_limit=outlet_limit_is_active(elapsed),
                )
        if trip_payload is not None:
            stop_reason = str(trip_payload["stop_reason"])

        interval = record_interval_for_elapsed(
            elapsed, staged=bool(config.staged_recording),
            default_s=float(config.record_interval_s),
        )
        if record_due(
            elapsed_s=elapsed,
            last_record_s=last_record,
            interval_s=interval,
            duration_s=float(config.duration_s),
            stopped=stop_reason != "completed",
        ):
            add_row(row, elapsed, q_set, w_set)
            last_record = elapsed
            print(
                f"[t+{elapsed:.3f}s] I={row['tec_main_current_A']:.6f} A "
                f"U={row['tec_main_voltage_V']:.6f} V "
                f"Pe/Pe0={latest['electric_power_ratio']:.6f} "
                f"Tout={row['core_outlet_T_K']:.3f} K",
                flush=True,
            )

        if stop_reason != "completed":
            save_checkpoint_atomic(system, out_dir / "emergency_restart.npz")
            _write_json(out_dir / "limit_trip.json", trip_payload or {})
            break

        if (
            float(config.checkpoint_interval_s) > 0.0
            and interval_due(
                elapsed, last_checkpoint, float(config.checkpoint_interval_s))
        ):
            save_checkpoint_atomic(system, out_dir / "latest_restart.npz")
            last_checkpoint = elapsed

    if stop_reason != "completed":
        restart_path = out_dir / "emergency_restart.npz"
        summary = {
            "candidate_start_time_s": candidate_start_time,
            "end_time_s": float(system.global_time),
            "initial_current_A": initial_current,
            "initial_electric_power_W": initial_electric_power,
            "fixed_u_metrics": fixed_u,
            "fixed_i_handoff_metrics": fixed_i,
            "latest_metrics": latest,
            "stop_reason": stop_reason,
            "trip": trip_payload,
            "restart_path": str(restart_path),
            "final_power_setpoint_W": q_set,
            "final_flow_setpoint_kg_s": w_set,
        }
        _write_json(out_dir / "run_summary.json", summary)
        return summary

    freeze_time = float(config.hold_before_ramp_s) + float(config.ramp_duration_s)
    elapsed = float(system.global_time) - candidate_start_time
    if elapsed >= freeze_time - 1.0e-12:
        q_set, w_set = float(config.final_power_w), float(config.final_flow_kg_s)
        set_total_flow_target(build, w_set)
        _apply_fixed_core_power(build, q_set)
        validate_final_setpoints(build, core, q_set, w_set)
        if (
            not math.isclose(float(latest.get("thermal_power_setpoint_W", float("nan"))), q_set)
            or not math.isclose(float(latest.get("flow_setpoint_kg_s", float("nan"))), w_set)
        ):
            terminal = collect_metrics(build, stage_index=2, dt_s=0.0)
            terminal.update(_collect_safety_metrics(build))
            add_row(terminal, elapsed, q_set, w_set)

    restart_path = out_dir / "final_restart.npz"
    save_checkpoint_atomic(system, restart_path)
    summary = {
        "candidate_start_time_s": candidate_start_time,
        "end_time_s": float(system.global_time),
        "initial_current_A": initial_current,
        "initial_electric_power_W": initial_electric_power,
        "fixed_u_metrics": fixed_u,
        "fixed_i_handoff_metrics": fixed_i,
        "latest_metrics": latest,
        "stop_reason": stop_reason,
        "trip": trip_payload,
        "restart_path": str(restart_path),
        "final_power_setpoint_W": q_set,
        "final_flow_setpoint_kg_s": w_set,
    }
    _write_json(out_dir / "run_summary.json", summary)
    return summary

def run_low_power_case(config: LowPowerRunConfig) -> Dict[str, Any]:
    failure_context: Dict[str, Any] = {"out_dir": Path(config.output_dir)}
    try:
        return _run_low_power_case_impl(config, failure_context=failure_context)
    except Exception as exc:
        try:
            write_failure_artifacts(
                Path(failure_context["out_dir"]), exc,
                system=failure_context.get("system"),
            )
        except Exception:
            pass
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--duration", type=float, default=0.1)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--tec-update-interval", type=float, default=0.5)
    parser.add_argument("--record-interval", type=float, default=0.05)
    parser.add_argument("--hold-before-ramp", type=float, default=0.0)
    parser.add_argument("--ramp-duration", type=float, default=0.0)
    parser.add_argument("--ramp-shape", choices=("cubic", "quintic"), default="cubic")
    parser.add_argument("--final-power", type=float, default=210000.0)
    parser.add_argument("--final-flow", type=float, default=2.46)
    parser.add_argument("--checkpoint-interval", type=float, default=60.0)
    parser.add_argument("--staged-recording", action="store_true")
    parser.add_argument("--resume-from", type=Path)
    parser.add_argument("--fixed-current", type=float)
    parser.add_argument("--disable-external-heat", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    result = run_low_power_case(LowPowerRunConfig(
        restart_in=args.restart_in, output_dir=args.output_dir, duration_s=args.duration,
        dt_s=args.dt, tec_update_interval_s=args.tec_update_interval,
        record_interval_s=args.record_interval,
        hold_before_ramp_s=args.hold_before_ramp, ramp_duration_s=args.ramp_duration,
        ramp_shape=args.ramp_shape,
        final_power_w=args.final_power, final_flow_kg_s=args.final_flow,
        checkpoint_interval_s=args.checkpoint_interval,
        staged_recording=args.staged_recording, resume_from=args.resume_from,
        fixed_current_a=args.fixed_current,
        external_heat_enabled=not args.disable_external_heat,
    ))
    return exit_code_for_result(result)


if __name__ == "__main__":
    raise SystemExit(main())
