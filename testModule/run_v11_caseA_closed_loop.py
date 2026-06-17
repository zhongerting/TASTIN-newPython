import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from run_v8_caseA_common import (
    DEFAULT_COOLANT_MATERIAL,
    DEFAULT_SOLID_ODE_METHOD,
    TOTAL_POWER_W,
    apply_solid_ode_method,
    apply_wire_resistance,
    get_solid_ode_methods,
    get_wire_resistance,
    json_default,
    parse_solid_ode_method,
    parse_v8_multipliers,
    passive_tec_source_totals,
)
from test_core_assemble_v10_caseA import build_v10_case_a_system
from test_core_assemble_v11_caseA import (
    V11_CASE_VERSION,
    V11_DEFAULT_PUMP_TOTAL_HEAD_PA,
    build_v11_case_a_system,
    reset_v11_design_flows,
    set_v11_pump_total_head,
    v11_basic_diagnostics,
)


DEFAULT_V10_RESTART = (
    "testModule/v10_caseA_tune727_ring020_hpfin075_steady_plus1000s/"
    "v10_caseA_tune727_ring020_hpfin075_steady_plus1000s_latest_restart.npz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V11 CaseA closed core + collector-ring pumped loop.")
    parser.add_argument("--restart-in", default=None)
    parser.add_argument("--init-from-v10", default=DEFAULT_V10_RESTART)
    parser.add_argument("--create-init-only", action="store_true")
    parser.add_argument("--output-dir", default="testModule/v11_caseA_closed_loop")
    parser.add_argument("--case-prefix", default="v11_caseA_closed_loop")
    parser.add_argument("--duration", type=float, default=200.0)
    parser.add_argument("--record-interval", type=float, default=50.0)
    parser.add_argument("--restart-interval", type=float, default=50.0)
    parser.add_argument("--max-dt", type=float, default=0.1)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--external-pipe-n-nodes", type=int, default=5)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--thermo-update-interval", type=float, default=0.8)
    parser.add_argument("--disable-tec-coupled", action="store_true")
    parser.add_argument("--inlet-temperature-k", type=float, default=727.0)
    parser.add_argument("--wire-resistance-scale", type=float, default=0.5)
    parser.add_argument("--ring-emissivity", type=float, default=0.2)
    parser.add_argument("--hp-emissivity", type=float, default=0.75)
    parser.add_argument("--fin-emissivity", type=float, default=0.75)
    parser.add_argument("--outer-header-emissivity", type=float, default=0.2)
    parser.add_argument("--outer-header-t-space-k", type=float, default=None)
    parser.add_argument("--outer-header-area-scale", type=float, default=1.0)
    parser.add_argument("--total-inlet-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--reference-pressure-pa", type=float, default=None)
    parser.add_argument("--pump-total-head-pa", type=float, default=V11_DEFAULT_PUMP_TOTAL_HEAD_PA)
    parser.add_argument("--target-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--enable-pump-head-control", action="store_true")
    parser.add_argument("--pump-control-interval", type=float, default=50.0)
    parser.add_argument("--pump-control-max-fraction", type=float, default=0.10)
    parser.add_argument("--connector-volume-m3", type=float, default=1.0e-5)
    parser.add_argument("--connector-length-m", type=float, default=0.02)
    parser.add_argument("--coolant-material", default=DEFAULT_COOLANT_MATERIAL)
    parser.add_argument(
        "--fluid-solid-coupling-scheme",
        choices=("current", "local_implicit"),
        default="current",
    )
    parser.add_argument(
        "--solid-ode-method",
        type=parse_solid_ode_method,
        default=DEFAULT_SOLID_ODE_METHOD,
    )
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=False),
        default=parse_v8_multipliers("1,6,12,15,3"),
    )
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=True),
        default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True),
    )
    args = parser.parse_args()
    if args.reference_pressure_pa is None:
        args.reference_pressure_pa = 166471.52
    if args.solid_ode_method == DEFAULT_SOLID_ODE_METHOD:
        args.solid_ode_method = parse_solid_ode_method("RK45")
    if args.fluid_solid_coupling_scheme == "current":
        args.fluid_solid_coupling_scheme = "local_implicit"
    return args


def flatten_for_csv(record: Dict[str, Any]) -> Dict[str, Any]:
    flat: Dict[str, Any] = {}
    for key, value in record.items():
        if isinstance(value, dict):
            for sub_key, sub_value in value.items():
                flat[f"{key}_{sub_key}"] = sub_value
        elif isinstance(value, (list, tuple)):
            for idx, sub_value in enumerate(value, start=1):
                flat[f"{key}_{idx}"] = sub_value
        else:
            flat[key] = value
    return flat


def append_row(path: Path, fieldnames: List[str], row: Dict[str, Any], *, write_header: bool) -> None:
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def apply_fluid_solid_coupling_scheme(system: Any, scheme: str) -> int:
    count = 0
    missing_capacitance: List[str] = []
    for coupler in getattr(system, "couplers", []):
        setter = getattr(coupler, "set_coupling_time_scheme", None)
        if not callable(setter):
            continue
        if scheme == "local_implicit" and getattr(coupler, "solid_node_capacitance", None) is None:
            missing_capacitance.append(getattr(coupler, "name", type(coupler).__name__))
            continue
        setter(scheme)
        count += 1
    if missing_capacitance:
        names = ", ".join(missing_capacitance)
        raise ValueError(
            "local_implicit fluid-solid coupling requires solid_node_capacitance for: "
            f"{names}"
        )
    return count


def _load_npz_dict(path: str) -> Dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def _copy_solid_states(build: Dict[str, Any], state: Dict[str, np.ndarray]) -> int:
    copied = 0
    solids = build["system"].solid_components
    solid_iter = solids.values() if isinstance(solids, dict) else solids
    for solid in solid_iter:
        prefix = f"Solid_{solid.name}"
        if f"{prefix}/T" in state and hasattr(solid, "load_state_dict"):
            solid.load_state_dict(state, prefix=prefix)
            copied += 1
    return copied


def _state_volume_map(system) -> Dict[str, Any]:
    return {getattr(vol, "name", ""): vol for vol in system.fluid_solver.volumes_obj}


def _state_junction_map(system) -> Dict[str, Any]:
    return {getattr(junc, "name", ""): junc for junc in system.fluid_solver.junctions_obj}


def _sync_objects_to_vectors(build: Dict[str, Any]) -> None:
    net = build["system"].fluid_solver
    for idx, vol in enumerate(net.volumes_obj):
        net.P_vec[idx] = float(vol.P)
        net.T_vec[idx] = float(vol.T)
        net.h_vec[idx] = float(vol.h)
    for idx, junc in enumerate(net.junctions_obj):
        net.W_vec[idx] = float(junc.W)
        if hasattr(junc, "update_velocity"):
            junc.update_velocity()
    if hasattr(net, "W_old"):
        net.W_old[:] = net.W_vec
    if hasattr(net, "W_iterate"):
        net.W_iterate[:] = net.W_vec
    net._update_fluid_properties()
    net._refresh_cached_pressure_targets()
    net._refresh_cached_boundary_targets()


def _copy_volume_state(dst, src) -> None:
    dst.P = float(src.P)
    dst.T = float(src.T)
    dst.h = float(src.h)
    if hasattr(dst, "rho") and hasattr(src, "rho"):
        dst.rho = float(src.rho)
    if hasattr(dst, "mu") and hasattr(src, "mu"):
        dst.mu = float(src.mu)


def _build_loaded_v10(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v10_case_a_system(
        inlet_temperature_k=args.inlet_temperature_k,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        outlet_pressure_pa=160000.0,
        pipe_n_nodes=args.pipe_n_nodes,
        external_pipe_n_nodes=args.external_pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        ring_emissivity=args.ring_emissivity,
        hp_emissivity=args.hp_emissivity,
        fin_emissivity=args.fin_emissivity,
        outer_header_emissivity=float(args.outer_header_emissivity),
        outer_header_t_space_k=args.outer_header_t_space_k,
        outer_header_area_scale=float(args.outer_header_area_scale),
    )
    build["system"].initialize_system()
    build["system"].load_global_state(args.init_from_v10)
    return build


def _average_volume_state(dst, a, b) -> None:
    dst.P = 0.5 * (float(a.P) + float(b.P))
    dst.T = 0.5 * (float(a.T) + float(b.T))
    if getattr(dst, "material", None) is not None:
        dst.h = dst.material.enthalpy(dst.T, dst.P)
        dst.update_properties(dst.material)


def inject_v10_restart(build: Dict[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    state = _load_npz_dict(args.init_from_v10)
    v10_loaded = _build_loaded_v10(args)
    copied_solids = _copy_solid_states(build, state)
    build["core"].load_state_dict(state, prefix="Macro_TASTIN_Core_V8_CaseA")

    dst_vols = _state_volume_map(build["system"])
    src_vols = _state_volume_map(v10_loaded["system"])
    copied_fluid = 0
    renamed_volume_map = {
        "V11_PumpOutletDistributor_51": "V10_HeaderToColdReturnSplit",
    }
    for name, dst in dst_vols.items():
        src_name = renamed_volume_map.get(name, name)
        if src_name in src_vols:
            _copy_volume_state(dst, src_vols[src_name])
            copied_fluid += 1

    h52_out = v10_loaded["radiator_outer_header_52"].volumes[-1]
    pump_outlet_source = src_vols.get("V10_HeaderToColdReturnSplit")
    if pump_outlet_source is not None:
        _average_volume_state(build["pump_mid_node"], h52_out, pump_outlet_source)

    if "CoreInletConnector" in dst_vols and "CoreInletConnector" in src_vols:
        reference_pressure = float(src_vols["CoreInletConnector"].P)
        dst_vols["CoreInletConnector"].P = reference_pressure
        dst_vols["CoreInletConnector"].target_P = reference_pressure
        args.reference_pressure_pa = reference_pressure

    dst_juncs = _state_junction_map(build["system"])
    src_juncs = _state_junction_map(v10_loaded["system"])
    renamed_junction_map = {
        "J_PumpOutletDistributor_51_to_ColdReturnBranch_1": "J_ColdReturnSplit_to_ColdReturnBranch_1",
        "J_PumpOutletDistributor_51_to_ColdReturnBranch_2_3_Rep": "J_ColdReturnSplit_to_ColdReturnBranch_2_3_Rep",
    }
    copied_junctions = 0
    for name, dst in dst_juncs.items():
        src_name = renamed_junction_map.get(name, name)
        if src_name in src_juncs:
            dst.W = float(src_juncs[src_name].W)
            copied_junctions += 1
    old_h52 = src_juncs.get("J_RadiatorOuterHeader_52_to_ColdReturnSplit")
    old_outlet = src_juncs.get("J_ColdReturnOutletMerge_to_OutletBoundary")
    if old_h52 is not None:
        build["pump_a"].W = float(old_h52.W)
        build["pump_b"].W = float(old_h52.W)
    if old_outlet is not None:
        build["j_cold_merge_to_core_inlet"].W = float(old_outlet.W)

    old_split = src_vols.get("V10_HeaderToColdReturnSplit")
    h52_dst = build["radiator_outer_header_52"].volumes[-1]
    if old_split is not None:
        downstream_offset = (
            float(h52_dst.P)
            + float(args.pump_total_head_pa)
            - float(old_split.P)
        )
        for name, dst in dst_vols.items():
            if (
                name == "V11_PumpOutletDistributor_51"
                or name == "V10_ColdReturnOutletMerge"
                or name.startswith("ColdReturnBranch_1")
                or name.startswith("ColdReturnBranch_2_3_Rep")
            ):
                src_name = renamed_volume_map.get(name, name)
                src = src_vols.get(src_name)
                if src is not None:
                    dst.P = float(src.P) + downstream_offset
        build["pump_mid_node"].P = float(h52_dst.P) + 0.5 * float(args.pump_total_head_pa)

    build["system"].global_time = float(state["System/global_time"][0])
    if "System/last_dt" in state:
        build["system"]._last_dt = float(state["System/last_dt"][0])

    reset_v11_design_flows(build, preserve_ring_restart_flows=True)
    set_v11_pump_total_head(build, float(args.pump_total_head_pa))
    _sync_objects_to_vectors(build)
    return {
        "copied_solids": copied_solids,
        "copied_fluid_volumes": copied_fluid,
        "copied_junctions": copied_junctions,
        "v10_restart": args.init_from_v10,
        "reference_pressure_pa": float(args.reference_pressure_pa),
    }


def build_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v11_case_a_system(
        inlet_temperature_k=args.inlet_temperature_k,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        reference_pressure_pa=args.reference_pressure_pa,
        pump_total_head_pa=args.pump_total_head_pa,
        pipe_n_nodes=args.pipe_n_nodes,
        external_pipe_n_nodes=args.external_pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        ring_emissivity=args.ring_emissivity,
        hp_emissivity=args.hp_emissivity,
        fin_emissivity=args.fin_emissivity,
        outer_header_emissivity=float(args.outer_header_emissivity),
        outer_header_t_space_k=args.outer_header_t_space_k,
        outer_header_area_scale=float(args.outer_header_area_scale),
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = not bool(args.disable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    apply_solid_ode_method(build, args.solid_ode_method)

    migration: Optional[Dict[str, Any]] = None
    if args.restart_in:
        system.load_global_state(args.restart_in)
        apply_solid_ode_method(build, args.solid_ode_method)
        set_v11_pump_total_head(build, float(args.pump_total_head_pa))
    else:
        migration = inject_v10_restart(build, args)

    core.point_reactor = None
    core.enable_tec_coupled = not bool(args.disable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    if core.enable_tec_coupled:
        core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=150.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    preserve_ring_restart_flows = bool(args.restart_in or migration is not None)
    reset_v11_design_flows(build, preserve_ring_restart_flows=preserve_ring_restart_flows)
    set_v11_pump_total_head(build, float(args.pump_total_head_pa))
    core.post_step(0.0, float(system.global_time))
    if core.enable_tec_coupled and core.thermo_calc is not None:
        apply_wire_resistance(core, scale=float(args.wire_resistance_scale))
        core._last_thermo_update_time = float(system.global_time)
    core.pre_step(0.0, float(system.global_time))
    build["fluid_solid_coupling_scheme"] = args.fluid_solid_coupling_scheme
    build["fluid_solid_coupler_count"] = apply_fluid_solid_coupling_scheme(
        system,
        args.fluid_solid_coupling_scheme,
    )
    build["solid_ode_method"] = args.solid_ode_method
    build["solid_ode_methods"] = get_solid_ode_methods(build)
    build["wire_resistance_scale"] = float(args.wire_resistance_scale)
    build["wire_resistance_ohm"] = get_wire_resistance(core)
    if core.thermo_calc is None:
        build["wire_resistance_ohm"] = [
            float(value) * float(args.wire_resistance_scale)
            for value in build["wire_resistance_ohm"]
        ]
    build["tec_coupled_enabled"] = bool(core.enable_tec_coupled)
    build["migration_summary"] = migration
    return build


def write_latest_state(path: Path, build: Dict[str, Any], args: argparse.Namespace, latest_restart: Path, history_path: Path, start_time: float, target_time: float, latest_record: Dict[str, Any]) -> None:
    latest = {
        "case_version": V11_CASE_VERSION,
        "restart_in": args.restart_in,
        "init_from_v10": args.init_from_v10,
        "restart_out": str(latest_restart),
        "history_csv": str(history_path),
        "start_time_s": start_time,
        "end_time_s": float(build["system"].global_time),
        "target_time_s": target_time,
        "record_interval_s": float(args.record_interval),
        "restart_interval_s": float(args.restart_interval),
        "max_dt_s": float(args.max_dt),
        "inner_iter": int(args.inner_iter),
        "coolant_material": build["coolant_material"],
        "inlet_temperature_k": float(args.inlet_temperature_k),
        "reference_pressure_pa": float(args.reference_pressure_pa),
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "pump_single_head_pa": float(build["pump_single_head_pa"]),
        "target_flow_kg_s": float(args.target_flow_kg_s),
        "enable_pump_head_control": bool(args.enable_pump_head_control),
        "wire_resistance_scale": float(args.wire_resistance_scale),
        "ring_emissivity": build["ring_emissivity"],
        "hp_emissivity": build["hp_emissivity"],
        "fin_emissivity": build["fin_emissivity"],
        "outer_header_emissivity": build["outer_header_emissivity"],
        "outer_header_t_space_k": build["outer_header_t_space_k"],
        "outer_header_area_scale": build["outer_header_area_scale"],
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "fluid_solid_coupler_count": build["fluid_solid_coupler_count"],
        "solid_ode_method": build["solid_ode_method"],
        "solid_ode_methods": get_solid_ode_methods(build),
        "wire_resistance_ohm": build["wire_resistance_ohm"],
        "tec_coupled_enabled": build["tec_coupled_enabled"],
        "ring_multipliers": build["ring_multipliers"],
        "tec_ring_multipliers": build["tec_ring_multipliers"],
        "migration_summary": build.get("migration_summary"),
        "latest": latest_record,
    }
    with path.open("w", encoding="utf-8") as f:
        json.dump(latest, f, indent=2, ensure_ascii=False, default=json_default)


def maybe_adjust_pump_head(build: Dict[str, Any], args: argparse.Namespace) -> Optional[Dict[str, float]]:
    if not bool(args.enable_pump_head_control):
        return None
    measured = float(build["pump_a"].W + build["pump_b"].W) * 0.5
    target = float(args.target_flow_kg_s)
    if measured <= 1.0e-9 or target <= 0.0:
        return None
    old_head = float(build["pump_total_head_pa"])
    raw_ratio = (target / measured) ** 2
    max_fraction = max(0.0, float(args.pump_control_max_fraction))
    ratio = float(np.clip(raw_ratio, 1.0 - max_fraction, 1.0 + max_fraction))
    new_head = max(1.0, old_head * ratio)
    set_v11_pump_total_head(build, new_head)
    return {
        "pump_control_old_head_pa": old_head,
        "pump_control_new_head_pa": new_head,
        "pump_control_measured_flow_kg_s": measured,
        "pump_control_target_flow_kg_s": target,
    }


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0 and not args.create_init_only:
        raise ValueError("duration must be positive unless --create-init-only is used.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_history.csv"
    latest_state_path = output_dir / f"{args.case_prefix}_latest_state.json"
    latest_restart_path = output_dir / f"{args.case_prefix}_latest_restart.npz"

    build = build_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)

    print("=== V11 CaseA closed core + collector ring pumped loop ===", flush=True)
    print(f"case_version={V11_CASE_VERSION}", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"init_from_v10={args.init_from_v10}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)
    print(f"latest_restart={latest_restart_path}", flush=True)
    print(f"fluid_solid_coupling_scheme={build['fluid_solid_coupling_scheme']}", flush=True)
    print(
        "tuning="
        f"inlet={args.inlet_temperature_k:.3f}K "
        f"wire_scale={args.wire_resistance_scale:.3f} "
        f"ring_eps={build['ring_emissivity']:.3f} "
        f"hp_eps={build['hp_emissivity']:.3f} "
        f"fin_eps={build['fin_emissivity']:.3f} "
        f"outer_header_eps={build['outer_header_emissivity']:.3f}",
        flush=True,
    )
    print(
        f"pump_total_head={build['pump_total_head_pa']:.3f}Pa "
        f"target_flow={args.target_flow_kg_s:.6f}kg/s "
        f"pump_control={bool(args.enable_pump_head_control)}",
        flush=True,
    )
    print(f"migration_summary={build.get('migration_summary')}", flush=True)

    initial_record = {
        **v11_basic_diagnostics(build),
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "wire_resistance_scale": build["wire_resistance_scale"],
    }
    system.save_global_state(str(latest_restart_path))
    write_latest_state(
        latest_state_path,
        build,
        args,
        latest_restart_path,
        history_path,
        start_time,
        target_time,
        initial_record,
    )
    if args.create_init_only:
        print("V11 initialized and saved without time advancement.", flush=True)
        return

    next_record_time = min(start_time + float(args.record_interval), target_time)
    next_restart_time = min(start_time + float(args.restart_interval), target_time)
    next_control_time = min(start_time + float(args.pump_control_interval), target_time)
    fieldnames: Optional[List[str]] = None
    last_record: Dict[str, Any] = initial_record

    while float(system.global_time) < target_time - 1.0e-10:
        build["core"].update_neutronic_power(
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
        dt = min(
            float(dt),
            float(args.max_dt),
            next_record_time - float(system.global_time),
            next_control_time - float(system.global_time),
            target_time - float(system.global_time),
        )
        system.step(dt, inner_iter=int(args.inner_iter))

        if float(system.global_time) >= next_control_time - 1.0e-10:
            control = maybe_adjust_pump_head(build, args)
            if control is not None:
                print(
                    "pump_control "
                    f"flow={control['pump_control_measured_flow_kg_s']:.6f}kg/s "
                    f"head={control['pump_control_old_head_pa']:.3f}->"
                    f"{control['pump_control_new_head_pa']:.3f}Pa",
                    flush=True,
                )
            next_control_time = min(next_control_time + float(args.pump_control_interval), target_time)

        if float(system.global_time) >= next_record_time - 1.0e-10:
            passive = passive_tec_source_totals(build)
            if any(value != 0.0 for value in passive.values()):
                raise RuntimeError(f"Ring3_Open TEC sources are not zero: {passive}")
            record = {
                **v11_basic_diagnostics(build),
                "relative_time_s": float(system.global_time) - start_time,
                "max_dt_s": float(args.max_dt),
                "inner_iter": int(args.inner_iter),
                "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
                "solid_ode_method": build["solid_ode_method"],
                "wire_resistance_scale": build["wire_resistance_scale"],
                "wire_resistance_ohm": build["wire_resistance_ohm"],
                "tec_coupled_enabled": build["tec_coupled_enabled"],
                "target_flow_kg_s": float(args.target_flow_kg_s),
                "flow_error_kg_s": float(build["pump_a"].W + build["pump_b"].W) * 0.5 - float(args.target_flow_kg_s),
            }
            flat = flatten_for_csv(record)
            if fieldnames is None:
                fieldnames = list(flat.keys())
            append_row(history_path, fieldnames, flat, write_header=not history_path.exists())
            last_record = record
            print(
                f"t_rel={record['relative_time_s']:.1f}s "
                f"Tin={record['core_inlet_connector_t_k']:.3f}K "
                f"Tcore_out={record['core_outlet_connector_t_k']:.3f}K "
                f"Touter_out={record['radiator_outer_header_52_t_out_k']:.3f}K "
                f"Wpump={record['pump_mean_flow_kg_s']:.6f}kg/s "
                f"H={record['pump_total_head_pa']:.2f}Pa "
                f"Pel={float(record['tec_total_electric_power_w'] or 0.0):.3f}W",
                flush=True,
            )
            next_record_time = min(next_record_time + float(args.record_interval), target_time)

        if float(system.global_time) >= next_restart_time - 1.0e-10:
            system.save_global_state(str(latest_restart_path))
            write_latest_state(
                latest_state_path,
                build,
                args,
                latest_restart_path,
                history_path,
                start_time,
                target_time,
                last_record,
            )
            next_restart_time = min(next_restart_time + float(args.restart_interval), target_time)

    system.save_global_state(str(latest_restart_path))
    write_latest_state(
        latest_state_path,
        build,
        args,
        latest_restart_path,
        history_path,
        start_time,
        target_time,
        last_record,
    )
    print("V11 closed-loop run completed.", flush=True)


if __name__ == "__main__":
    main()
