import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from test_core_assemble_v7_caseA import (
    build_v7_case_a_system,
    _case_a_electric_diagnostics,
    _case_a_flow_diagnostics,
    _case_a_reset_design_flows_after_restart,
)
from test_core_assemble_v7_caseA_faststeady import compute_faststeady_energy_audit


DEFAULT_RESTART_IN = "test_core_assemble_v7_caseA_faststeady_restart_t8800.npz"
DEFAULT_OUTPUT_DIR = "testModule/v7_caseA_electric_dt08_outerCp001_overnight"
CASE_PREFIX = "test_core_assemble_v7_caseA_electric_dt08_outerCp001"
TOTAL_POWER_W = 115000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run V7 CaseA electric-coupled overnight continuation with dt <= 0.8 s."
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART_IN)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--case-prefix", default=CASE_PREFIX)
    parser.add_argument("--initial-duration", type=float, default=2000.0)
    parser.add_argument("--followup-duration", type=float, default=200.0)
    parser.add_argument("--max-followup-segments", type=int, default=None)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--outer-cp-scale", type=float, default=0.01)
    parser.add_argument("--thermo-update-interval", type=float, default=0.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    return parser.parse_args()


def _json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return str(value)


def _safe_time_label(value: float) -> str:
    if abs(value - round(value)) < 1.0e-6:
        return str(int(round(value)))
    return f"{value:.3f}".replace(".", "p")


def _remove_if_present(path: Optional[Path]) -> None:
    if path is not None and path.exists():
        path.unlink()


def build_loaded_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v7_case_a_system(
        pipe_n_nodes=args.pipe_n_nodes,
        solid_heat_capacity_scale=args.outer_cp_scale,
        solid_heat_capacity_scale_scope="global_outer",
    )
    system = build["system"]
    core = build["core"]

    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = float(args.thermo_update_interval)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )

    system.initialize_system()
    if not os.path.exists(args.restart_in):
        raise FileNotFoundError(f"Restart file not found: {args.restart_in}")
    system.load_global_state(args.restart_in)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = float(args.thermo_update_interval)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)
    return build


def advance_to(build: Dict[str, Any], stop_time: float, args: argparse.Namespace) -> None:
    system = build["system"]
    core = build["core"]
    current_time = float(system.global_time)
    while current_time < stop_time - 1.0e-10:
        core.update_neutronic_power(
            p_total=TOTAL_POWER_W,
            p_fiss=TOTAL_POWER_W,
            p_decay=0.0,
            alpha=1.0,
        )
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=float(args.max_dt),
            safety_factor=float(args.safety_factor),
        )
        dt = min(float(dt), float(args.max_dt), stop_time - current_time)
        if dt <= 0.0:
            raise RuntimeError(f"Non-positive dt={dt} at t={current_time}")
        system.step(dt, inner_iter=1)
        current_time = float(system.global_time)


def collect_state(
    build: Dict[str, Any],
    relative_elapsed_s: float,
    phase: str,
    segment_index: int,
    restart_path: Path,
    state_path: Path,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    system = build["system"]
    core = build["core"]
    tfes = build["tfes"]

    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is not None:
        core.thermo_calc.calculate(verbose=False)

    flow = _case_a_flow_diagnostics(build)
    electric = _case_a_electric_diagnostics(core)
    audit = compute_faststeady_energy_audit(
        {
            "build": build,
            "system": system,
            "core": core,
        }
    )

    tfe_summary: Dict[str, Dict[str, float]] = {}
    for name, tfe in tfes.items():
        solids = tfe.solids
        coolant_outlet = float(tfe.coolant.volumes[-1].T)
        component_max = {
            "pellet_max_k": float(np.max(solids["pellet"].T)),
            "emitter_max_k": float(np.max(solids["emitter"].T)),
            "collector_max_k": float(np.max(solids["collector"].T)),
            "inner_clad_max_k": float(np.max(solids["inner_clad"].T)),
            "outer_clad_max_k": float(np.max(solids["outer_clad"].T)),
            "moderator_max_k": float(np.max(solids["moderator"].T)),
        }
        tfe_summary[name] = {
            **component_max,
            "solid_max_k": max(component_max.values()),
            "coolant_outlet_k": coolant_outlet,
        }

    state = {
        "case": args.case_prefix,
        "phase": phase,
        "segment_index": int(segment_index),
        "restart_in": args.restart_in,
        "restart_file": str(restart_path),
        "state_file": str(state_path),
        "relative_elapsed_s": float(relative_elapsed_s),
        "absolute_time_s": float(system.global_time),
        "last_dt_s": float(getattr(system, "_last_dt", np.nan)),
        "max_dt_setting_s": float(args.max_dt),
        "safety_factor": float(args.safety_factor),
        "thermo_update_interval_s": float(core.thermo_update_interval),
        "outer_cp_scale": float(args.outer_cp_scale),
        "tfe_cp_scale": 1.0,
        "total_power_w": TOTAL_POWER_W,
        "inlet_plenum_temperature_k": float(build["inlet_plenum"].T),
        "outlet_plenum_temperature_k": float(build["outlet_plenum"].T),
        "center_coolant_outlet_k": float(tfes["Center"].coolant.volumes[-1].T),
        "core_delta_p_pa": float(flow["core_plenum_delta_p_pa"]),
        "total_macro_flow_kg_s": float(flow["inlet_total_macro_flow_kg_s"]),
        "tfe_total_macro_flow_kg_s": float(flow["tfe_total_macro_flow_kg_s"]),
        **electric,
        "energy": audit,
        "flow": flow,
        "tfes": tfe_summary,
    }
    return state


def flatten_for_csv(state: Dict[str, Any]) -> Dict[str, Any]:
    row = {
        "phase": state["phase"],
        "segment_index": state["segment_index"],
        "relative_elapsed_s": state["relative_elapsed_s"],
        "absolute_time_s": state["absolute_time_s"],
        "last_dt_s": state["last_dt_s"],
        "max_dt_setting_s": state["max_dt_setting_s"],
        "safety_factor": state["safety_factor"],
        "thermo_update_interval_s": state["thermo_update_interval_s"],
        "outer_cp_scale": state["outer_cp_scale"],
        "tfe_cp_scale": state["tfe_cp_scale"],
        "restart_file": state["restart_file"],
        "state_file": state["state_file"],
        "total_power_w": state["total_power_w"],
        "coolant_heat_pickup_w": state["energy"]["coolant_heat_pickup_w"],
        "electric_power_w": state["energy"]["electric_power_w"],
        "outer_wall_radiation_w": state["energy"]["outer_wall_radiation_w"],
        "balance_residual_w": state["energy"]["balance_residual_w"],
        "balance_residual_percent": state["energy"]["balance_residual_percent"],
        "tec_total_voltage_v": state["tec_total_voltage_v"],
        "tec_total_current_a": state["tec_total_current_a"],
        "tec_total_electric_power_w": state["tec_total_electric_power_w"],
        "inlet_plenum_temperature_k": state["inlet_plenum_temperature_k"],
        "outlet_plenum_temperature_k": state["outlet_plenum_temperature_k"],
        "center_coolant_outlet_k": state["center_coolant_outlet_k"],
        "core_delta_p_pa": state["core_delta_p_pa"],
        "total_macro_flow_kg_s": state["total_macro_flow_kg_s"],
        "tfe_total_macro_flow_kg_s": state["tfe_total_macro_flow_kg_s"],
    }
    for tfe_name, summary in state["tfes"].items():
        for key, value in summary.items():
            row[f"{tfe_name}_{key}"] = value
    return row


def append_history(history_path: Path, row: Dict[str, Any]) -> None:
    write_header = not history_path.exists()
    fieldnames = list(row.keys())
    with history_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def save_checkpoint(
    build: Dict[str, Any],
    output_dir: Path,
    history_path: Path,
    relative_elapsed_s: float,
    phase: str,
    segment_index: int,
    restart_path: Path,
    state_path: Path,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    build["system"].save_global_state(str(restart_path))
    state = collect_state(
        build=build,
        relative_elapsed_s=relative_elapsed_s,
        phase=phase,
        segment_index=segment_index,
        restart_path=restart_path,
        state_path=state_path,
        args=args,
    )
    with state_path.open("w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2, default=_json_default)
    append_history(history_path, flatten_for_csv(state))
    print(
        f"[Checkpoint] {phase} rel={relative_elapsed_s:.3f}s "
        f"abs={state['absolute_time_s']:.3f}s restart={restart_path}"
    )
    return state


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_history.csv"

    build = build_loaded_case(args)
    system = build["system"]
    start_time = float(system.global_time)

    initial_stop_time = start_time + float(args.initial_duration)
    advance_to(build, initial_stop_time, args)
    rel = float(system.global_time - start_time)
    rel_label = _safe_time_label(rel)
    abs_label = _safe_time_label(float(system.global_time))
    fixed_restart = output_dir / f"{args.case_prefix}_restart_rel{rel_label}_abs{abs_label}.npz"
    fixed_state = output_dir / f"{args.case_prefix}_state_rel{rel_label}_abs{abs_label}.json"
    save_checkpoint(
        build,
        output_dir,
        history_path,
        rel,
        "initial_2000s",
        0,
        fixed_restart,
        fixed_state,
        args,
    )

    latest_followup_restart: Optional[Path] = None
    segment_index = 0
    while args.max_followup_segments is None or segment_index < args.max_followup_segments:
        segment_index += 1
        stop_time = float(system.global_time) + float(args.followup_duration)
        advance_to(build, stop_time, args)
        rel = float(system.global_time - start_time)
        rel_label = _safe_time_label(rel)
        abs_label = _safe_time_label(float(system.global_time))
        restart_path = output_dir / f"{args.case_prefix}_restart_rel{rel_label}_abs{abs_label}.npz"
        state_path = output_dir / f"{args.case_prefix}_state_rel{rel_label}_abs{abs_label}.json"
        save_checkpoint(
            build,
            output_dir,
            history_path,
            rel,
            "followup_200s",
            segment_index,
            restart_path,
            state_path,
            args,
        )
        if latest_followup_restart is not None and latest_followup_restart != restart_path:
            _remove_if_present(latest_followup_restart)
        latest_followup_restart = restart_path


if __name__ == "__main__":
    main()
