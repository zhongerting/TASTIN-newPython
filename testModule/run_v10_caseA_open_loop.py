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

import CoolantLoop.model_collector_ring_6segment_v9_interface as ring_cfg
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
from test_core_assemble_v9_caseA import build_v9_case_a_system
from test_core_assemble_v10_caseA import (
    V10_CASE_VERSION,
    build_v10_case_a_system,
    reset_v10_design_flows,
    v10_basic_diagnostics,
)


DEFAULT_V9_RESTART = (
    "testModule/v9_caseA_open_loop_tec_3000s/"
    "v9_caseA_open_loop_tec_3000s_latest_restart.npz"
)
DEFAULT_RING_RESTART = (
    "CoolantLoop/collector_ring_6segment_v9_interface_500s_from200s_restart.npz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V10 CaseA open core + collector-ring model.")
    parser.add_argument("--restart-in", default=None)
    parser.add_argument("--init-from-v9", default=DEFAULT_V9_RESTART)
    parser.add_argument("--init-from-ring", default=DEFAULT_RING_RESTART)
    parser.add_argument("--create-init-only", action="store_true")
    parser.add_argument("--output-dir", default="testModule/v10_caseA_open_loop")
    parser.add_argument("--case-prefix", default="v10_caseA_open_loop")
    parser.add_argument("--duration", type=float, default=200.0)
    parser.add_argument("--record-interval", type=float, default=20.0)
    parser.add_argument("--restart-interval", type=float, default=20.0)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--external-pipe-n-nodes", type=int, default=5)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--thermo-update-interval", type=float, default=0.8)
    parser.add_argument("--disable-tec-coupled", action="store_true")
    parser.add_argument("--inlet-temperature-k", type=float, default=753.330663091)
    parser.add_argument("--total-inlet-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--outlet-pressure-pa", type=float, default=160000.0)
    parser.add_argument("--connector-volume-m3", type=float, default=1.0e-5)
    parser.add_argument("--connector-length-m", type=float, default=0.02)
    parser.add_argument("--coolant-material", default=DEFAULT_COOLANT_MATERIAL)
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
    return parser.parse_args()


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


def _build_loaded_v9(args: argparse.Namespace):
    build = build_v9_case_a_system(
        inlet_temperature_k=743.0,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        outlet_pressure_pa=args.outlet_pressure_pa,
        pipe_n_nodes=args.pipe_n_nodes,
        external_pipe_n_nodes=args.external_pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    build["system"].initialize_system()
    build["system"].load_global_state(args.init_from_v9)
    return build


def _build_loaded_ring(args: argparse.Namespace):
    model = ring_cfg.build_model()
    model["sys_mgr"].initialize_system()
    model["sys_mgr"].load_global_state(args.init_from_ring)
    return model


def inject_v9_and_ring_restart(build: Dict[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    v9_state = _load_npz_dict(args.init_from_v9)
    ring_state = _load_npz_dict(args.init_from_ring)
    v9_loaded = _build_loaded_v9(args)
    ring_loaded = _build_loaded_ring(args)

    copied_core_solids = _copy_solid_states(build, v9_state)
    copied_ring_solids = _copy_solid_states(build, ring_state)
    build["core"].load_state_dict(v9_state, prefix="Macro_TASTIN_Core_V8_CaseA")

    dst_vols = _state_volume_map(build["system"])
    v9_vols = _state_volume_map(v9_loaded["system"])
    ring_vols = _state_volume_map(ring_loaded["sys_mgr"])
    copied_fluid_v9 = 0
    copied_fluid_ring = 0

    v9_excluded_after_ring = {
        "RadiatorInnerHeader_53",
        "RadiatorOuterHeader_52",
        "ColdReturnBranch_1",
        "ColdReturnBranch_2_3_Rep",
    }
    for name, src in v9_vols.items():
        if any(name == item or name.startswith(f"{item}_") for item in v9_excluded_after_ring):
            continue
        if name in dst_vols:
            _copy_volume_state(dst_vols[name], src)
            copied_fluid_v9 += 1
    for name, src in ring_vols.items():
        if name in dst_vols and (
            name.startswith("InletMix_")
            or name.startswith("OutletMix_")
            or name.startswith("A")
            or name.startswith("Manifold_")
        ):
            _copy_volume_state(dst_vols[name], src)
            copied_fluid_ring += 1

    if "CoreInletConnector" in dst_vols and "V10_InletBoundary_FixedFlow" in dst_vols:
        _copy_volume_state(dst_vols["V10_InletBoundary_FixedFlow"], dst_vols["CoreInletConnector"])

    ring_outlet_source = ring_vols.get("OutletBuffer_Vol_05")
    if ring_outlet_source is not None:
        downstream_prefixes = (
            "V10_RadiatorManifoldMerge",
            "RadiatorInnerHeader_53",
            "RadiatorOuterHeader_52",
            "V10_HeaderToColdReturnSplit",
            "ColdReturnBranch_1",
            "ColdReturnBranch_2_3_Rep",
            "V10_ColdReturnOutletMerge",
            "V10_OutletBoundary_FixedPressure",
        )
        for name, dst in dst_vols.items():
            if any(name == prefix or name.startswith(f"{prefix}_") for prefix in downstream_prefixes):
                _copy_volume_state(dst, ring_outlet_source)
        dst_vols["V10_OutletBoundary_FixedPressure"].P = float(args.outlet_pressure_pa)
        if hasattr(dst_vols["V10_OutletBoundary_FixedPressure"], "target_P"):
            dst_vols["V10_OutletBoundary_FixedPressure"].target_P = float(args.outlet_pressure_pa)

    dst_juncs = _state_junction_map(build["system"])
    for name, src in _state_junction_map(v9_loaded["system"]).items():
        if name in dst_juncs:
            dst_juncs[name].W = float(src.W)
    for name, src in _state_junction_map(ring_loaded["sys_mgr"]).items():
        if name in dst_juncs and (
            name.startswith("A")
            or name.startswith("Manifold_")
            or name.startswith("J_InletMix_")
            or name.startswith("J_A")
            or name.startswith("J_I")
            or name.startswith("J_OutletMix_")
            or name.startswith("J_Manifold_")
        ):
            dst_juncs[name].W = float(src.W)

    build["system"].global_time = float(v9_state["System/global_time"][0])
    if "System/last_dt" in v9_state:
        build["system"]._last_dt = float(v9_state["System/last_dt"][0])

    reset_v10_design_flows(build, preserve_ring_restart_flows=True)
    _sync_objects_to_vectors(build)
    return {
        "copied_core_solids": copied_core_solids,
        "copied_ring_solids": copied_ring_solids,
        "copied_fluid_v9": copied_fluid_v9,
        "copied_fluid_ring": copied_fluid_ring,
        "v9_restart": args.init_from_v9,
        "ring_restart": args.init_from_ring,
    }


def build_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v10_case_a_system(
        inlet_temperature_k=args.inlet_temperature_k,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        outlet_pressure_pa=args.outlet_pressure_pa,
        pipe_n_nodes=args.pipe_n_nodes,
        external_pipe_n_nodes=args.external_pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = not bool(args.disable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    apply_solid_ode_method(build, args.solid_ode_method)
    system.initialize_system()

    migration: Optional[Dict[str, Any]] = None
    if args.restart_in:
        system.load_global_state(args.restart_in)
        apply_solid_ode_method(build, args.solid_ode_method)
    else:
        migration = inject_v9_and_ring_restart(build, args)

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
    reset_v10_design_flows(build, preserve_ring_restart_flows=preserve_ring_restart_flows)
    core.post_step(0.0, float(system.global_time))
    if core.enable_tec_coupled and core.thermo_calc is not None:
        apply_wire_resistance(core)
        core._last_thermo_update_time = float(system.global_time)
    core.pre_step(0.0, float(system.global_time))
    build["solid_ode_method"] = args.solid_ode_method
    build["solid_ode_methods"] = get_solid_ode_methods(build)
    build["wire_resistance_ohm"] = get_wire_resistance(core)
    build["tec_coupled_enabled"] = bool(core.enable_tec_coupled)
    build["migration_summary"] = migration
    return build


def write_latest_state(path: Path, build: Dict[str, Any], args: argparse.Namespace, latest_restart: Path, history_path: Path, start_time: float, target_time: float, latest_record: Dict[str, Any]) -> None:
    latest = {
        "case_version": V10_CASE_VERSION,
        "restart_in": args.restart_in,
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

    print("=== V10 CaseA open core + collector ring ===", flush=True)
    print(f"case_version={V10_CASE_VERSION}", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"init_from_v9={args.init_from_v9}", flush=True)
    print(f"init_from_ring={args.init_from_ring}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)
    print(f"latest_restart={latest_restart_path}", flush=True)
    print(f"migration_summary={build.get('migration_summary')}", flush=True)

    initial_record = v10_basic_diagnostics(build)
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
        print("V10 initialized and saved without time advancement.", flush=True)
        return

    next_record_time = min(start_time + float(args.record_interval), target_time)
    next_restart_time = min(start_time + float(args.restart_interval), target_time)
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
            target_time - float(system.global_time),
        )
        system.step(dt, inner_iter=int(args.inner_iter))

        if float(system.global_time) >= next_record_time - 1.0e-10:
            passive = passive_tec_source_totals(build)
            if any(value != 0.0 for value in passive.values()):
                raise RuntimeError(f"Ring3_Open TEC sources are not zero: {passive}")
            record = {
                **v10_basic_diagnostics(build),
                "relative_time_s": float(system.global_time) - start_time,
                "max_dt_s": float(args.max_dt),
                "inner_iter": int(args.inner_iter),
                "solid_ode_method": build["solid_ode_method"],
                "wire_resistance_ohm": build["wire_resistance_ohm"],
                "tec_coupled_enabled": build["tec_coupled_enabled"],
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
                f"Tring_out={record['manifold_1_t_out_k']:.3f}K "
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
    print("V10 open-loop run completed.", flush=True)


if __name__ == "__main__":
    main()
