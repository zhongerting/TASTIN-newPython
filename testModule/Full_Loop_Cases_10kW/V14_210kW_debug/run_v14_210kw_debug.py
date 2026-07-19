"""Staged powered debug run for V14_10kW heat-pipe full-loop case."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases_10kW import (  # noqa: E402
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    build_v14_case_a_system,
)

CASE_NAME = "V14_10kW_210kW_powered_debug"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "runs" / "default"
DEFAULT_LOOKUP_DB = str(REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr")
DEFAULT_LOOKUP_REGIONS = ("core", "startup", "high_power", "accident")
WIRE_RESISTANCE_OHM = (0.00155199999999970, 0.00102400000000000, 0.000336000000000000, 0.000608000000000000)
HISTORY_FIELDS = [
    "stage",
    "time_s",
    "dt_s",
    "core_inlet_T_K",
    "core_outlet_T_K",
    "core_delta_T_K",
    "radiator_heat_rejection_W",
    "upper_ring_heatpipe_rejection_W",
    "lower_ring_heatpipe_rejection_W",
    "ring_wall_heat_rejection_W",
    "total_external_heat_rejection_W",
    "net_power_estimate_W",
    "tec_main_current_A",
    "tec_main_voltage_V",
    "tec_main_electric_power_W",
    "tec_main_converged",
    "pump_required_head_total_Pa",
    "pump_a_required_head_Pa",
    "pump_b_required_head_Pa",
    "pump_a_flow_kg_s",
    "pump_b_flow_kg_s",
    "min_fluid_T_K",
    "max_fluid_T_K",
    "min_solid_T_K",
    "max_solid_T_K",
    "fluid_converged",
]


@dataclass(frozen=True)
class DebugRunConfig:
    output_dir: Path = DEFAULT_OUTPUT_DIR
    stage_durations_s: Sequence[float] = (1.0, 10.0, 100.0)
    dt_s: float = 0.05
    record_interval_s: float = 2.0
    checkpoint_interval_s: float = 10.0
    min_fluid_temperature_stop_k: Optional[float] = 500.0
    tec_electrical_enabled: bool = False
    tec_voltage_v: float = 50.0
    tec_current_guess_a: float = 150.0
    tec_lookup_enabled: bool = False
    tec_lookup_db: Optional[str] = DEFAULT_LOOKUP_DB
    tec_lookup_regions: Sequence[str] = DEFAULT_LOOKUP_REGIONS
    wire_resistance_scale: float = 1.0
    radiator_emissivity: float = 0.75
    hp_up_view_factor: float = 0.0
    upper_hp_down_view_factor: float = 0.3
    lower_hp_down_view_factor: float = 0.3
    inner_iter: int = 1
    fluid_max_iter: int = 100
    init_dt_s: float = 0.01
    init_tol_kg_s: float = 1.0e-4
    init_max_iter: int = 1000
    power_w: float = 210000.0
    target_flow_kg_s: float = 2.46
    initial_temperature_k: float = 754.15
    space_temperature_k: float = 4.0
    restart_in: Optional[Path] = None
    case_prefix: str = CASE_NAME


def _parse_stage_durations(text: str) -> List[float]:
    values = [float(item.strip()) for item in str(text).split(",") if item.strip()]
    if not values or any(value <= 0.0 for value in values):
        raise ValueError("stage durations must be positive comma-separated seconds")
    return values


def _finite_float(value: Any) -> float:
    try:
        out = float(value)
    except Exception:
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def _sync_fluid_temperature(net: Any, temperature_k: float) -> None:
    t_value = float(temperature_k)
    for vol in net.volumes_obj:
        vol.T = t_value
        if hasattr(vol.material, "enthalpy"):
            vol.h = float(vol.material.enthalpy(t_value, vol.P))
        if hasattr(vol, "update_properties"):
            vol.update_properties(vol.material)
    if hasattr(net, "_initialize_state_from_objects"):
        net._initialize_state_from_objects()
    if hasattr(net, "_update_fluid_properties"):
        net._update_fluid_properties()
    if hasattr(net, "_sync_vectors_to_objects"):
        net._sync_vectors_to_objects()


def _sync_solid_temperature(system: Any, temperature_k: float) -> None:
    t_value = float(temperature_k)
    for solid in system.solid_components.values():
        if hasattr(solid, "T"):
            solid.T[:] = t_value
        if hasattr(solid, "current_time"):
            solid.current_time = float(system.global_time)
        if hasattr(solid, "_update_properties"):
            solid._update_properties()
        if hasattr(solid, "_update_boundaries_state"):
            solid._update_boundaries_state(current_time=float(system.global_time))


def _set_implicit_euler(system: Any) -> None:
    for solid in system.solid_components.values():
        if hasattr(solid, "set_ode_method"):
            solid.set_ode_method("implicit_euler")


def _apply_wire_resistance(core: Any, scale: float) -> List[float]:
    thermo_calc = getattr(core, "thermo_calc", None)
    wire_res = np.asarray(WIRE_RESISTANCE_OHM, dtype=float) * float(scale)
    if thermo_calc is None:
        return [float(value) for value in wire_res]
    thermo_calc._input_data.resistanceWire = np.tile(wire_res, (thermo_calc.N_elem, 1))
    thermo_calc.build()
    return [float(value) for value in wire_res]


def _apply_fixed_core_power(build: Dict[str, Any], power_w: float) -> None:
    build["core"].update_neutronic_power(
        p_total=float(power_w),
        p_fiss=float(power_w),
        p_decay=0.0,
        alpha=1.0,
    )


def _junction_by_name(net: Any, name: str) -> Optional[Any]:
    for junction in net.junctions_obj:
        if getattr(junction, "name", "") == name:
            return junction
    return None


def _pump_head(junction: Optional[Any]) -> float:
    if junction is None:
        return float("nan")
    return _finite_float(getattr(junction.to_vol, "P", np.nan)) - _finite_float(getattr(junction.from_vol, "P", np.nan))


def _ring_rejection(ring_hps: Iterable[Any]) -> float:
    total = 0.0
    for ring_hp in ring_hps:
        getter = getattr(ring_hp, "get_total_heat_rejection_scaled", None)
        if callable(getter):
            total += _finite_float(getter())
    return total


def _ring_wall_rejection(build: Dict[str, Any]) -> float:
    total = 0.0
    for solid in build.get("ring_solids", []):
        boundary = getattr(solid, "boundaries", {}).get("right")
        if boundary is None or not hasattr(boundary, "current_flux"):
            continue
        outward = -float(np.sum(np.asarray(boundary.current_flux, dtype=float)))
        total += max(0.0, outward)
    return total



def _tec_main_metrics(build: Dict[str, Any]) -> Dict[str, Any]:
    core = build.get("core")
    if core is None or not bool(getattr(core, "enable_tec_coupled", False)):
        return {
            "tec_main_current_A": 0.0,
            "tec_main_voltage_V": 0.0,
            "tec_main_electric_power_W": 0.0,
            "tec_main_converged": True,
        }
    try:
        results = core.get_tec_circuit_global_results().get("main")
    except Exception:
        results = None
    if not results:
        return {
            "tec_main_current_A": float("nan"),
            "tec_main_voltage_V": float("nan"),
            "tec_main_electric_power_W": float("nan"),
            "tec_main_converged": False,
        }
    current = _finite_float(results.get("Iout", np.nan))
    voltage = _finite_float(results.get("Uout", np.nan))
    return {
        "tec_main_current_A": current,
        "tec_main_voltage_V": voltage,
        "tec_main_electric_power_W": current * voltage,
        "tec_main_converged": bool(results.get("converged", True)),
    }

def _solid_min_max(system: Any) -> tuple[float, float]:
    values = []
    for solid in system.solid_components.values():
        if hasattr(solid, "T"):
            arr = np.asarray(solid.T, dtype=float)
            if arr.size:
                values.append(arr.reshape(-1))
    if not values:
        return float("nan"), float("nan")
    all_values = np.concatenate(values)
    return float(np.nanmin(all_values)), float(np.nanmax(all_values))


def collect_metrics(build: Dict[str, Any], stage_index: int, dt_s: float) -> Dict[str, Any]:
    system = build["system"]
    net = system.fluid_solver
    pump_a = _junction_by_name(net, "J_PumpA")
    pump_b = _junction_by_name(net, "J_PumpB")
    pump_a_head = _pump_head(pump_a)
    pump_b_head = _pump_head(pump_b)
    upper_q = _ring_rejection(build.get("ring_hps", [])[:6])
    lower_q = _ring_rejection(build.get("ring_hps", [])[6:])
    heatpipe_q = upper_q + lower_q
    ring_wall_q = _ring_wall_rejection(build)
    total_external_q = heatpipe_q + ring_wall_q
    core_power = _finite_float(getattr(build["core"], "last_total_core_power", np.nan))
    solid_min, solid_max = _solid_min_max(system)
    tec = _tec_main_metrics(build)
    diag = getattr(system, "last_step_diagnostics", None) or {}
    fluid_flags = diag.get("fluid_converged_by_iteration", [])
    return {
        "stage": int(stage_index),
        "time_s": float(system.global_time),
        "dt_s": float(dt_s),
        "core_inlet_T_K": _finite_float(build["core_inlet_connector"].T),
        "core_outlet_T_K": _finite_float(build["core_outlet_connector"].T),
        "core_delta_T_K": _finite_float(build["core_outlet_connector"].T) - _finite_float(build["core_inlet_connector"].T),
        "radiator_heat_rejection_W": heatpipe_q,
        "upper_ring_heatpipe_rejection_W": upper_q,
        "lower_ring_heatpipe_rejection_W": lower_q,
        "ring_wall_heat_rejection_W": ring_wall_q,
        "total_external_heat_rejection_W": total_external_q,
        "net_power_estimate_W": core_power - total_external_q,
        **tec,
        "pump_required_head_total_Pa": pump_a_head + pump_b_head,
        "pump_a_required_head_Pa": pump_a_head,
        "pump_b_required_head_Pa": pump_b_head,
        "pump_a_flow_kg_s": _finite_float(getattr(pump_a, "W", np.nan)),
        "pump_b_flow_kg_s": _finite_float(getattr(pump_b, "W", np.nan)),
        "min_fluid_T_K": float(np.nanmin(net.T_vec)),
        "max_fluid_T_K": float(np.nanmax(net.T_vec)),
        "min_solid_T_K": solid_min,
        "max_solid_T_K": solid_max,
        "fluid_converged": bool(all(fluid_flags)) if fluid_flags else True,
    }


def _format_progress(row: Dict[str, Any]) -> str:
    return (
        f"[t={row['time_s']:.2f}s stage={row['stage']}] "
        f"Tin={row['core_inlet_T_K']:.3f}K Tout={row['core_outlet_T_K']:.3f}K "
        f"Qhp={row['radiator_heat_rejection_W']:.2f}W "
        f"Qwall={row['ring_wall_heat_rejection_W']:.2f}W "
        f"Qext={row['total_external_heat_rejection_W']:.2f}W "
        f"Pnet={row['net_power_estimate_W']:.2f}W "
        f"Ptec={row['tec_main_electric_power_W']:.2f}W "
        f"Itec={row['tec_main_current_A']:.3f}A "
        f"TminF={row['min_fluid_T_K']:.3f}K "
        f"dPpump={row['pump_required_head_total_Pa']:.2f}Pa "
        f"Wpump=({row['pump_a_flow_kg_s']:.4f},{row['pump_b_flow_kg_s']:.4f})kg/s "
        f"fluid_ok={row['fluid_converged']}"
    )


def _load_restart_with_coupling_dt(
        system: Any,
        restart_path: Any,
        *,
        fallback_dt: float) -> None:
    had_instance_override = '_run_couplers' in vars(system)
    previous_instance_override = vars(system).get('_run_couplers')
    original_run_couplers = system._run_couplers

    def run_couplers_with_dt(*args: Any, **kwargs: Any) -> Any:
        if kwargs.get('dt') is None:
            sync_dt = float(getattr(system, '_last_dt', fallback_dt))
            if not math.isfinite(sync_dt) or sync_dt <= 0.0:
                sync_dt = float(fallback_dt)
            kwargs['dt'] = sync_dt
        return original_run_couplers(*args, **kwargs)

    system._run_couplers = run_couplers_with_dt
    try:
        system.load_global_state(str(restart_path))
    finally:
        if had_instance_override:
            system._run_couplers = previous_instance_override
        else:
            del system._run_couplers


def build_debug_case(config: DebugRunConfig, *, apply_fixed_power: bool = True) -> Dict[str, Any]:
    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(
            inlet_temperature_k=float(config.initial_temperature_k),
            main_tec_enabled=bool(config.tec_electrical_enabled),
            main_tec_target_value=float(config.tec_voltage_v),
            main_tec_current_guess_a=float(config.tec_current_guess_a),
            tec_lookup_enabled=bool(config.tec_lookup_enabled),
            tec_lookup_db=config.tec_lookup_db,
            tec_lookup_regions=tuple(config.tec_lookup_regions),
        ),
        flow_config=FullLoopFlowConfig(total_flow_kg_s=float(config.target_flow_kg_s)),
        pump_config=FullLoopPumpConfig(
            pump_total_head_pa=1.0,
            pump_flow_control=True,
            target_flow_kg_s=float(config.target_flow_kg_s),
        ),
        radiator_config=V14HeatPipeRadiatorConfig(
            t_space_k=float(config.space_temperature_k),
            hp_initial_temp_k=float(config.initial_temperature_k),
            hp_emissivity=float(config.radiator_emissivity),
            fin_emissivity=float(config.radiator_emissivity),
            hp_up_view_factor=float(config.hp_up_view_factor),
            upper_hp_down_view_factor=float(config.upper_hp_down_view_factor),
            lower_hp_down_view_factor=float(config.lower_hp_down_view_factor),
        ),
    )
    system = build["system"]
    _set_implicit_euler(system)
    if config.restart_in is not None:
        _load_restart_with_coupling_dt(
            system,
            config.restart_in,
            fallback_dt=float(config.init_dt_s),
        )
    else:
        _sync_fluid_temperature(system.fluid_solver, config.initial_temperature_k)
        _sync_solid_temperature(system, config.initial_temperature_k)
        system.global_time = 0.0
    _set_implicit_euler(system)
    if bool(config.tec_electrical_enabled):
        build["core"].enable_tec_coupled = True
        build["core"].setup_tec_circuit(
            mode_str="fixed_u",
            target_value=float(config.tec_voltage_v),
            I_guess=float(config.tec_current_guess_a),
            topology="series",
        )
        build["wire_resistance_ohm"] = _apply_wire_resistance(build["core"], config.wire_resistance_scale)
    else:
        build["wire_resistance_ohm"] = [float(value) * float(config.wire_resistance_scale) for value in WIRE_RESISTANCE_OHM]
    build["wire_resistance_scale"] = float(config.wire_resistance_scale)
    if apply_fixed_power:
        _apply_fixed_core_power(build, config.power_w)
    system.initialize_system(
        dt_init=float(config.init_dt_s),
        tol=float(config.init_tol_kg_s),
        max_iter=int(config.init_max_iter),
    )
    if apply_fixed_power:
        _apply_fixed_core_power(build, config.power_w)
    return build


def _write_history_header(path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        csv.DictWriter(f, fieldnames=HISTORY_FIELDS).writeheader()


def _append_history(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=HISTORY_FIELDS, extrasaction="ignore")
        for row in rows:
            writer.writerow(row)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def _write_latest_state(out_dir: Path, stage: int, restart_path: Path, summary_path: Optional[Path], metrics: Dict[str, Any]) -> None:
    _write_json(out_dir / "latest_state.json", {
        "case": CASE_NAME,
        "latest_stage": int(stage),
        "latest_restart_path": str(restart_path),
        "latest_summary_path": str(summary_path) if summary_path is not None else None,
        "latest_metrics": metrics,
    })


def _record_row(history_path: Path, rows: List[Dict[str, Any]], row: Dict[str, Any]) -> None:
    rows.append(row)
    _append_history(history_path, [row])
    print(_format_progress(row), flush=True)


def _low_temperature_triggered(config: DebugRunConfig, row: Dict[str, Any]) -> bool:
    threshold = config.min_fluid_temperature_stop_k
    return threshold is not None and float(threshold) > 0.0 and row["min_fluid_T_K"] < float(threshold)


def run_debug_case(config: DebugRunConfig) -> Dict[str, Any]:
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / "history.csv"
    _write_history_header(history_path)
    _write_json(out_dir / "run_config.json", {
        "case": CASE_NAME,
        "stage_durations_s": [float(v) for v in config.stage_durations_s],
        "dt_s": float(config.dt_s),
        "record_interval_s": float(config.record_interval_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "min_fluid_temperature_stop_k": None if config.min_fluid_temperature_stop_k is None else float(config.min_fluid_temperature_stop_k),
        "tec_electrical_calculation_enabled": bool(config.tec_electrical_enabled),
        "tec_voltage_v": float(config.tec_voltage_v),
        "tec_current_guess_a": float(config.tec_current_guess_a),
        "tec_lookup_enabled": bool(config.tec_lookup_enabled),
        "tec_lookup_db": config.tec_lookup_db,
        "tec_lookup_regions": list(config.tec_lookup_regions),
        "wire_resistance_scale": float(config.wire_resistance_scale),
        "wire_resistance_ohm": [float(value) * float(config.wire_resistance_scale) for value in WIRE_RESISTANCE_OHM],
        "radiator_emissivity": float(config.radiator_emissivity),
        "hp_up_view_factor": float(config.hp_up_view_factor),
        "upper_hp_down_view_factor": float(config.upper_hp_down_view_factor),
        "lower_hp_down_view_factor": float(config.lower_hp_down_view_factor),
        "inner_iter": int(config.inner_iter),
        "fluid_max_iter": int(config.fluid_max_iter),
        "power_w": float(config.power_w),
        "target_flow_kg_s": float(config.target_flow_kg_s),
        "initial_temperature_k": float(config.initial_temperature_k),
        "space_temperature_k": float(config.space_temperature_k),
        "restart_in": str(config.restart_in) if config.restart_in is not None else None,
        "solid_ode_method": "implicit_euler",
        "external_heat_enabled": False,
        "point_kinetics_enabled": False,
        "pump_mode": "fixed_target_flow",
    })

    build = build_debug_case(config)
    system = build["system"]
    latest_rows: List[Dict[str, Any]] = []
    stage_summaries: List[Dict[str, Any]] = []
    last_record_time = -float("inf")
    last_checkpoint_time = float(system.global_time)
    stop_reason = "completed"
    stopped = False

    initial = collect_metrics(build, stage_index=0, dt_s=0.0)
    _record_row(history_path, latest_rows, initial)

    for stage_index, duration_s in enumerate(config.stage_durations_s, start=1):
        stage_start = float(system.global_time)
        stage_end = stage_start + float(duration_s)
        stage_rows: List[Dict[str, Any]] = []
        while float(system.global_time) < stage_end - 1.0e-12:
            dt = min(float(config.dt_s), stage_end - float(system.global_time))
            _apply_fixed_core_power(build, config.power_w)
            system.step(
                dt,
                inner_iter=int(config.inner_iter),
                fail_on_fluid_nonconvergence=False,
                fluid_max_iter=int(config.fluid_max_iter),
            )
            _apply_fixed_core_power(build, config.power_w)

            row_for_checks: Optional[Dict[str, Any]] = None
            should_record = (
                float(system.global_time) - last_record_time >= float(config.record_interval_s) - 1.0e-12
                or float(system.global_time) >= stage_end - 1.0e-12
            )
            if should_record:
                row_for_checks = collect_metrics(build, stage_index=stage_index, dt_s=dt)
                stage_rows.append(row_for_checks)
                _record_row(history_path, latest_rows, row_for_checks)
                last_record_time = float(system.global_time)

            if row_for_checks is None:
                row_for_checks = collect_metrics(build, stage_index=stage_index, dt_s=dt)

            if _low_temperature_triggered(config, row_for_checks):
                if not latest_rows or abs(latest_rows[-1]["time_s"] - row_for_checks["time_s"]) > 1.0e-12:
                    stage_rows.append(row_for_checks)
                    _record_row(history_path, latest_rows, row_for_checks)
                emergency_path = out_dir / "emergency_low_temp_restart.npz"
                system.save_global_state(str(emergency_path))
                stop_reason = "low_fluid_temperature"
                stopped = True
                _write_json(out_dir / "low_temperature_stop.json", {
                    "case": CASE_NAME,
                    "time_s": float(system.global_time),
                    "threshold_k": float(config.min_fluid_temperature_stop_k),
                    "restart_path": str(emergency_path),
                    "latest_metrics": row_for_checks,
                })
                _write_latest_state(out_dir, stage_index, emergency_path, None, row_for_checks)
                print(f"[stop] min_fluid_T_K below {float(config.min_fluid_temperature_stop_k):.3f} K; saved {emergency_path}", flush=True)
                break

            if float(config.checkpoint_interval_s) > 0.0 and float(system.global_time) - last_checkpoint_time >= float(config.checkpoint_interval_s) - 1.0e-12:
                checkpoint_path = out_dir / f"checkpoint_t{int(round(system.global_time)):06d}s.npz"
                system.save_global_state(str(checkpoint_path))
                last_checkpoint_time = float(system.global_time)
                _write_json(out_dir / "latest_checkpoint.json", {
                    "case": CASE_NAME,
                    "time_s": float(system.global_time),
                    "checkpoint_path": str(checkpoint_path),
                    "latest_metrics": row_for_checks,
                })
                _write_latest_state(out_dir, stage_index, checkpoint_path, None, row_for_checks)
                print(f"[checkpoint] saved {checkpoint_path}", flush=True)

        if stopped:
            break

        if not stage_rows or stage_rows[-1]["time_s"] < float(system.global_time) - 1.0e-12:
            row = collect_metrics(build, stage_index=stage_index, dt_s=0.0)
            stage_rows.append(row)
            _record_row(history_path, latest_rows, row)

        restart_path = out_dir / f"stage_{stage_index:02d}_restart.npz"
        system.save_global_state(str(restart_path))
        summary = {
            "case": CASE_NAME,
            "stage": int(stage_index),
            "stage_duration_s": float(duration_s),
            "stage_start_time_s": stage_start,
            "stage_end_time_s": float(system.global_time),
            "restart_path": str(restart_path),
            "latest_metrics": stage_rows[-1],
        }
        summary_path = out_dir / f"stage_{stage_index:02d}_summary.json"
        _write_json(summary_path, summary)
        stage_summaries.append(summary)
        _write_latest_state(out_dir, stage_index, restart_path, summary_path, stage_rows[-1])

    result = {
        "case": CASE_NAME,
        "output_dir": str(out_dir),
        "history_path": str(history_path),
        "stage_summaries": stage_summaries,
        "latest_metrics": latest_rows[-1],
        "stop_reason": stop_reason,
    }
    _write_json(out_dir / "run_summary.json", result)
    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--stage-durations", default="1,10,100", help="Comma-separated stage durations in seconds.")
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--record-interval", type=float, default=2.0)
    parser.add_argument("--checkpoint-interval", type=float, default=10.0)
    parser.add_argument("--min-fluid-temperature-stop", type=float, default=500.0, help="Emergency stop threshold in K; <=0 disables it.")
    parser.add_argument("--enable-tec", action="store_true", help="Enable main series TEC electrical calculation.")
    parser.add_argument("--tec-voltage", type=float, default=50.0)
    parser.add_argument("--tec-current-guess", type=float, default=150.0)
    parser.add_argument("--enable-tec-lookup", action="store_true", help="Enable ThermoCalc emission lookup tables.")
    parser.add_argument("--tec-lookup-db", default=DEFAULT_LOOKUP_DB)
    parser.add_argument("--tec-lookup-regions", default=",".join(DEFAULT_LOOKUP_REGIONS))
    parser.add_argument("--wire-resistance-scale", type=float, default=1.0)
    parser.add_argument("--radiator-emissivity", type=float, default=0.75)
    parser.add_argument("--hp-up-view-factor", type=float, default=0.0)
    parser.add_argument("--upper-hp-down-view-factor", type=float, default=0.3)
    parser.add_argument("--lower-hp-down-view-factor", type=float, default=0.3)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--fluid-max-iter", type=int, default=100)
    parser.add_argument("--init-dt", type=float, default=0.01)
    parser.add_argument("--init-tol", type=float, default=1.0e-4)
    parser.add_argument("--init-max-iter", type=int, default=1000)
    parser.add_argument("--power-w", type=float, default=210000.0)
    parser.add_argument("--target-flow", type=float, default=2.46)
    parser.add_argument("--initial-temperature-k", type=float, default=754.15)
    parser.add_argument("--space-temperature-k", type=float, default=4.0)
    parser.add_argument("--restart-in", type=Path, default=None)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    stop_threshold = None if float(args.min_fluid_temperature_stop) <= 0.0 else float(args.min_fluid_temperature_stop)
    config = DebugRunConfig(
        output_dir=args.output_dir,
        stage_durations_s=_parse_stage_durations(args.stage_durations),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        min_fluid_temperature_stop_k=stop_threshold,
        tec_electrical_enabled=bool(args.enable_tec),
        tec_voltage_v=float(args.tec_voltage),
        tec_current_guess_a=float(args.tec_current_guess),
        tec_lookup_enabled=bool(args.enable_tec_lookup),
        tec_lookup_db=str(args.tec_lookup_db) if args.tec_lookup_db else None,
        tec_lookup_regions=tuple(part.strip() for part in str(args.tec_lookup_regions).split(",") if part.strip()),
        wire_resistance_scale=float(args.wire_resistance_scale),
        radiator_emissivity=float(args.radiator_emissivity),
        hp_up_view_factor=float(args.hp_up_view_factor),
        upper_hp_down_view_factor=float(args.upper_hp_down_view_factor),
        lower_hp_down_view_factor=float(args.lower_hp_down_view_factor),
        inner_iter=int(args.inner_iter),
        fluid_max_iter=int(args.fluid_max_iter),
        init_dt_s=float(args.init_dt),
        init_tol_kg_s=float(args.init_tol),
        init_max_iter=int(args.init_max_iter),
        power_w=float(args.power_w),
        target_flow_kg_s=float(args.target_flow),
        initial_temperature_k=float(args.initial_temperature_k),
        space_temperature_k=float(args.space_temperature_k),
        restart_in=args.restart_in,
    )
    result = run_debug_case(config)
    print(json.dumps(result["latest_metrics"], indent=2, sort_keys=True, ensure_ascii=False))
    print(f"Stop reason: {result['stop_reason']}")
    print(f"Saved outputs to: {result['output_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
