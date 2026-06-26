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
    load_tec_load_curve,
    parse_solid_ode_method,
    parse_v8_multipliers,
    passive_tec_source_totals,
)
from test_core_assemble_v12_caseA import build_v12_case_a_system
from test_core_assemble_v13_caseA import (
    V13_CASE_VERSION,
    V13_DEFAULT_INLET_TEMPERATURE_K,
    V13_DEFAULT_PUMP_TOTAL_HEAD_PA,
    V13_DEFAULT_REFERENCE_PRESSURE_PA,
    attach_radiator_thermal_shield,
    build_v13_case_a_system,
    reset_v13_design_flows,
    set_v13_pump_total_head,
    v13_basic_diagnostics,
)


DEFAULT_V12_RESTART = (
    "testModule/v12_caseA_open_loop_tec_wire05_200s_from1500s/"
    "v12_caseA_open_loop_tec_wire05_200s_from1500s_latest_restart.npz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V13 CaseA closed core + TOPAZ-II pipe-fin radiator loop.")
    parser.add_argument("--restart-in", default=None)
    parser.add_argument("--init-from-v12", default=DEFAULT_V12_RESTART)
    parser.add_argument("--create-init-only", action="store_true")
    parser.add_argument("--output-dir", default="testModule/v13_caseA_closed_loop")
    parser.add_argument("--case-prefix", default="v13_caseA_closed_loop")
    parser.add_argument("--duration", type=float, default=100.0)
    parser.add_argument("--record-interval", type=float, default=50.0)
    parser.add_argument("--restart-interval", type=float, default=50.0)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--init-dt", type=float, default=0.05)
    parser.add_argument("--hydraulic-tol", type=float, default=1.0e-6)
    parser.add_argument("--hydraulic-max-iter", type=int, default=800)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--inlet-temperature-k", type=float, default=V13_DEFAULT_INLET_TEMPERATURE_K)
    parser.add_argument("--total-inlet-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--reference-pressure-pa", type=float, default=None)
    parser.add_argument("--pump-total-head-pa", type=float, default=V13_DEFAULT_PUMP_TOTAL_HEAD_PA)
    parser.add_argument("--target-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--enable-pump-head-control", action="store_true")
    parser.add_argument("--disable-pump-flow-control", action="store_true")
    parser.add_argument("--pump-control-interval", type=float, default=50.0)
    parser.add_argument("--pump-control-max-fraction", type=float, default=0.10)
    parser.add_argument("--connector-volume-m3", type=float, default=1.0e-5)
    parser.add_argument("--connector-length-m", type=float, default=0.02)
    parser.add_argument("--coolant-material", default=DEFAULT_COOLANT_MATERIAL)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--enable-reserved-parallel-tec", action="store_true")
    parser.add_argument("--reserved-parallel-mode", choices=("fixed_u", "fixed_i", "load_curve"), default="fixed_u")
    parser.add_argument("--reserved-parallel-voltage", type=float, default=0.8)
    parser.add_argument("--reserved-parallel-current", type=float, default=6000.0)
    parser.add_argument("--reserved-parallel-load-curve", default=None)
    parser.add_argument("--thermo-update-interval", type=float, default=1.0)
    parser.add_argument("--disable-tec-coupled", action="store_true")
    parser.add_argument("--wire-resistance-scale", type=float, default=0.5)
    parser.add_argument("--total-power-w", type=float, default=TOTAL_POWER_W)
    parser.add_argument("--n-tubes", type=int, default=78)
    parser.add_argument("--n-axial", type=int, default=8)
    parser.add_argument("--n-radial-wall", type=int, default=1)
    parser.add_argument("--n-fin-width", type=int, default=12)
    parser.add_argument("--tube-length-m", type=float, default=1.85)
    parser.add_argument("--tube-inner-diameter-m", type=float, default=0.007)
    parser.add_argument("--tube-outer-diameter-m", type=float, default=0.008)
    parser.add_argument("--upper-header-centerline-diameter-m", type=float, default=0.824)
    parser.add_argument("--lower-header-centerline-diameter-m", type=float, default=1.346)
    parser.add_argument("--header-inner-diameter-m", type=float, default=0.020)
    parser.add_argument("--fin-thickness-m", type=float, default=0.0004)
    parser.add_argument("--fin-width-upper-m", type=float, default=0.03319)
    parser.add_argument("--fin-width-lower-m", type=float, default=0.05421)
    parser.add_argument("--tube-emissivity", type=float, default=0.80)
    parser.add_argument("--fin-emissivity", type=float, default=0.80)
    parser.add_argument("--tube-area-scale", type=float, default=1.0)
    parser.add_argument("--fin-area-scale", type=float, default=0.35)
    parser.add_argument("--t-space-k", type=float, default=3.0)
    parser.add_argument("--fin-conductivity-w-m-k", type=float, default=348.9)
    parser.add_argument("--fin-view-factor", type=float, default=1.0)
    parser.add_argument("--fin-contact-resistance-m2k-w", type=float, default=0.0)
    parser.add_argument("--enable-radiation-shield", "--enable-radiator-shield", dest="enable_radiation_shield", action="store_true")
    parser.add_argument("--shield-active-until-s", type=float, default=None)
    parser.add_argument("--shield-inner-emissivity", type=float, default=0.8)
    parser.add_argument("--shield-outer-emissivity", type=float, default=0.8)
    parser.add_argument("--shield-conductivity-w-m-k", type=float, default=1.0)
    parser.add_argument("--shield-thickness-m", type=float, default=0.002)
    parser.add_argument("--shield-view-factor", type=float, default=0.8)
    parser.add_argument("--shield-solar-heat-flux-w-m2", type=float, default=0.0)
    parser.add_argument("--shield-background-temperature-k", type=float, default=3.0)
    parser.add_argument("--shield-relaxation", type=float, default=1.0)
    parser.add_argument(
        "--shield-model",
        choices=("segment_balance", "fortran_shield2"),
        default="segment_balance",
    )
    parser.add_argument("--radiator-header-k-loss", type=float, default=1.0)
    parser.add_argument("--radiator-tube-inlet-k-loss", type=float, default=100.0)
    parser.add_argument("--radiator-tube-outlet-k-loss", type=float, default=100.0)
    parser.add_argument("--connector-k-loss", type=float, default=0.0)
    parser.add_argument(
        "--fluid-solid-coupling-scheme",
        choices=("current", "local_implicit"),
        default="local_implicit",
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
        args.reference_pressure_pa = V13_DEFAULT_REFERENCE_PRESSURE_PA
    if args.solid_ode_method == DEFAULT_SOLID_ODE_METHOD:
        args.solid_ode_method = parse_solid_ode_method("RK45")
    return args


def configure_core_tec_circuits(core, args: argparse.Namespace) -> None:
    """Configure the legacy main series TEC circuit and optional reserved parallel circuit."""
    core.setup_tec_circuit("fixed_u", float(args.target_voltage), I_guess=150.0, topology="series")
    if not bool(args.enable_reserved_parallel_tec):
        if hasattr(core, "disable_reserved_parallel_tec_circuit"):
            core.disable_reserved_parallel_tec_circuit()
        return

    mode = str(args.reserved_parallel_mode)
    load_curve = load_tec_load_curve(args.reserved_parallel_load_curve)
    if mode == "fixed_u":
        target_value = float(args.reserved_parallel_voltage)
    elif mode == "fixed_i":
        target_value = float(args.reserved_parallel_current)
    elif mode == "load_curve":
        if load_curve is None:
            raise ValueError("--reserved-parallel-load-curve is required when --reserved-parallel-mode=load_curve.")
        target_value = float(args.reserved_parallel_voltage)
    else:
        raise ValueError(f"Unsupported reserved parallel TEC mode: {mode}")

    core.setup_reserved_parallel_tec_circuit(
        mode_str=mode,
        target_value=target_value,
        I_guess=float(args.reserved_parallel_current),
        multipliers={"Ring3_Open": int(core.tfe_multipliers.get("Ring3_Open", 0))},
        load_curve=load_curve,
    )


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
        raise ValueError(f"local_implicit fluid-solid coupling requires solid_node_capacitance for: {names}")
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


def _state_volume_map(system: Any) -> Dict[str, Any]:
    return {getattr(vol, "name", ""): vol for vol in system.fluid_solver.volumes_obj}


def _state_junction_map(system: Any) -> Dict[str, Any]:
    return {getattr(junc, "name", ""): junc for junc in system.fluid_solver.junctions_obj}


def _copy_volume_state(dst: Any, src: Any) -> None:
    dst.P = float(src.P)
    dst.T = float(src.T)
    dst.h = float(src.h)
    if hasattr(dst, "rho") and hasattr(src, "rho"):
        dst.rho = float(src.rho)
    if hasattr(dst, "mu") and hasattr(src, "mu"):
        dst.mu = float(src.mu)


def _set_volume_thermo(dst: Any, *, pressure_pa: float, temperature_k: float) -> None:
    dst.P = float(pressure_pa)
    dst.T = float(temperature_k)
    if getattr(dst, "material", None) is not None:
        dst.h = dst.material.enthalpy(dst.T, dst.P)
        dst.update_properties(dst.material)


def _set_volume_pressure_preserve_temperature(vol: Any, pressure_pa: float) -> None:
    _set_volume_thermo(vol, pressure_pa=float(pressure_pa), temperature_k=float(vol.T))


def _linearize_channel_pressure(channel: Any, p_start: float, p_end: float) -> None:
    volumes = list(getattr(channel, "volumes", []))
    if not volumes:
        return
    if len(volumes) == 1:
        _set_volume_pressure_preserve_temperature(volumes[0], 0.5 * (float(p_start) + float(p_end)))
        return
    for idx, vol in enumerate(volumes):
        frac = idx / float(len(volumes) - 1)
        _set_volume_pressure_preserve_temperature(vol, (1.0 - frac) * float(p_start) + frac * float(p_end))


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


def _build_loaded_v12(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v12_case_a_system(
        inlet_temperature_k=args.inlet_temperature_k,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        pipe_n_nodes=args.pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        enable_tec_coupled=not bool(args.disable_tec_coupled),
        n_tubes=args.n_tubes,
        n_axial=args.n_axial,
        n_radial_wall=args.n_radial_wall,
        n_fin_width=args.n_fin_width,
        tube_length_m=args.tube_length_m,
        tube_inner_diameter_m=args.tube_inner_diameter_m,
        tube_outer_diameter_m=args.tube_outer_diameter_m,
        upper_header_centerline_diameter_m=args.upper_header_centerline_diameter_m,
        lower_header_centerline_diameter_m=args.lower_header_centerline_diameter_m,
        header_inner_diameter_m=args.header_inner_diameter_m,
        fin_thickness_m=args.fin_thickness_m,
        fin_width_upper_m=args.fin_width_upper_m,
        fin_width_lower_m=args.fin_width_lower_m,
        tube_emissivity=args.tube_emissivity,
        fin_emissivity=args.fin_emissivity,
        tube_area_scale=args.tube_area_scale,
        fin_area_scale=args.fin_area_scale,
        t_space_k=args.t_space_k,
        fin_conductivity_w_m_k=args.fin_conductivity_w_m_k,
        fin_view_factor=args.fin_view_factor,
        fin_contact_resistance_m2k_w=args.fin_contact_resistance_m2k_w,
        radiator_header_k_loss=args.radiator_header_k_loss,
        radiator_tube_inlet_k_loss=args.radiator_tube_inlet_k_loss,
        radiator_tube_outlet_k_loss=args.radiator_tube_outlet_k_loss,
        connector_k_loss=args.connector_k_loss,
        fluid_solid_coupling_scheme="current",
        solid_ode_method=args.solid_ode_method,
    )
    build["system"].initialize_system()
    build["system"].load_global_state(args.init_from_v12)
    return build


def inject_v12_restart(build: Dict[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    state = _load_npz_dict(args.init_from_v12)
    v12_loaded = _build_loaded_v12(args)
    copied_solids = _copy_solid_states(build, state)
    build["core"].load_state_dict(state, prefix="Macro_TASTIN_Core_V8_CaseA")

    dst_vols = _state_volume_map(build["system"])
    src_vols = _state_volume_map(v12_loaded["system"])
    copied_fluid = 0
    for name, dst in dst_vols.items():
        if name in src_vols:
            _copy_volume_state(dst, src_vols[name])
            copied_fluid += 1

    core_inlet = dst_vols.get("V12_CoreInletConnector")
    src_core_inlet = src_vols.get("V12_CoreInletConnector")
    if core_inlet is not None and src_core_inlet is not None:
        args.reference_pressure_pa = float(src_core_inlet.P)
        core_inlet.P = float(src_core_inlet.P)
        core_inlet.target_P = float(src_core_inlet.P)

    pipe09_last = build["flow_network_cold_pipes"][-1].volumes[-1]
    pipe11_first = build["pipe11_core_inlet_header"].volumes[0]
    p_ref = float(args.reference_pressure_pa)
    pump_head = float(args.pump_total_head_pa)
    p_feed = p_ref + 100.0
    p_core_out = p_ref - 0.45 * pump_head
    p_rad_in = p_ref - 0.55 * pump_head
    p_rad_out = p_ref - 0.88 * pump_head
    p_suction = p_feed - pump_head
    for channel in build["fluid_channels"].values():
        _linearize_channel_pressure(channel, p_ref - 100.0, p_core_out + 100.0)
    _set_volume_pressure_preserve_temperature(build["core_outlet_connector"], p_core_out)
    _linearize_channel_pressure(build["pipe05_core_outlet_to_radiator"], p_core_out - 100.0, p_rad_in + 100.0)
    _set_volume_pressure_preserve_temperature(build["radiator_inlet_split"], p_rad_in)
    for channel in build["radiator_upper_nodes"]:
        _linearize_channel_pressure(channel, p_rad_in + 50.0, p_rad_in - 50.0)
    for channel in build["radiator_lower_nodes"]:
        _linearize_channel_pressure(channel, p_rad_out + 50.0, p_rad_out - 50.0)
    for channel in build["radiator_tube_channels"]:
        _linearize_channel_pressure(channel, p_rad_in - 100.0, p_rad_out + 100.0)
    _set_volume_pressure_preserve_temperature(build["radiator_outlet_mix"], p_rad_out)
    cold_channels = list(build["flow_network_cold_pipes"])
    cold_bounds = np.linspace(p_rad_out - 100.0, p_suction, len(cold_channels) + 1)
    for idx, channel in enumerate(cold_channels):
        _linearize_channel_pressure(channel, float(cold_bounds[idx]), float(cold_bounds[idx + 1]))
    _set_volume_thermo(
        build["pump_mid_node"],
        pressure_pa=p_suction + 0.5 * pump_head,
        temperature_k=0.5 * (float(pipe09_last.T) + float(pipe11_first.T)),
    )
    _set_volume_thermo(
        build["pump_outlet_node"],
        pressure_pa=p_feed,
        temperature_k=float(pipe11_first.T),
    )
    _linearize_channel_pressure(build["pipe11_core_inlet_header"], p_feed, p_ref + 20.0)
    _set_volume_pressure_preserve_temperature(build["core_inlet_distribution"], p_ref + 15.0)
    _linearize_channel_pressure(build["inlet_pipe_1"], p_ref + 10.0, p_ref + 1.0)
    _linearize_channel_pressure(build["inlet_pipe_23"], p_ref + 10.0, p_ref + 1.0)
    if core_inlet is not None:
        _set_volume_pressure_preserve_temperature(core_inlet, p_ref)
        core_inlet.target_P = p_ref

    dst_juncs = _state_junction_map(build["system"])
    src_juncs = _state_junction_map(v12_loaded["system"])
    copied_junctions = 0
    for name, dst in dst_juncs.items():
        if name in src_juncs:
            dst.W = float(src_juncs[name].W)
            copied_junctions += 1

    old_return = src_juncs.get("J_Pipe09_to_OutletBoundary")
    old_inlet = src_juncs.get("J_V12_InletBoundary_to_Pipe11")
    pump_seed = float(args.total_inlet_flow_kg_s)
    if old_return is not None:
        pump_seed = float(old_return.W)
    elif old_inlet is not None:
        pump_seed = float(old_inlet.W)
    build["pump_a"].W = pump_seed
    build["pump_b"].W = pump_seed
    build["j_pump_to_pipe11"].W = pump_seed

    build["system"].global_time = float(state["System/global_time"][0])
    if "System/last_dt" in state:
        build["system"]._last_dt = float(state["System/last_dt"][0])

    reset_v13_design_flows(build)
    set_v13_pump_total_head(build, float(args.pump_total_head_pa))
    _sync_objects_to_vectors(build)
    return {
        "copied_solids": copied_solids,
        "copied_fluid_volumes": copied_fluid,
        "copied_junctions": copied_junctions,
        "v12_restart": args.init_from_v12,
        "reference_pressure_pa": float(args.reference_pressure_pa),
    }


def build_case(args: argparse.Namespace) -> Dict[str, Any]:
    build = build_v13_case_a_system(
        inlet_temperature_k=args.inlet_temperature_k,
        total_inlet_flow_kg_s=args.total_inlet_flow_kg_s,
        reference_pressure_pa=args.reference_pressure_pa,
        pump_total_head_pa=args.pump_total_head_pa,
        pipe_n_nodes=args.pipe_n_nodes,
        connector_volume_m3=args.connector_volume_m3,
        connector_length_m=args.connector_length_m,
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        enable_tec_coupled=not bool(args.disable_tec_coupled),
        n_tubes=args.n_tubes,
        n_axial=args.n_axial,
        n_radial_wall=args.n_radial_wall,
        n_fin_width=args.n_fin_width,
        tube_length_m=args.tube_length_m,
        tube_inner_diameter_m=args.tube_inner_diameter_m,
        tube_outer_diameter_m=args.tube_outer_diameter_m,
        upper_header_centerline_diameter_m=args.upper_header_centerline_diameter_m,
        lower_header_centerline_diameter_m=args.lower_header_centerline_diameter_m,
        header_inner_diameter_m=args.header_inner_diameter_m,
        fin_thickness_m=args.fin_thickness_m,
        fin_width_upper_m=args.fin_width_upper_m,
        fin_width_lower_m=args.fin_width_lower_m,
        tube_emissivity=args.tube_emissivity,
        fin_emissivity=args.fin_emissivity,
        tube_area_scale=args.tube_area_scale,
        fin_area_scale=args.fin_area_scale,
        t_space_k=args.t_space_k,
        fin_conductivity_w_m_k=args.fin_conductivity_w_m_k,
        fin_view_factor=args.fin_view_factor,
        fin_contact_resistance_m2k_w=args.fin_contact_resistance_m2k_w,
        radiator_header_k_loss=args.radiator_header_k_loss,
        radiator_tube_inlet_k_loss=args.radiator_tube_inlet_k_loss,
        radiator_tube_outlet_k_loss=args.radiator_tube_outlet_k_loss,
        connector_k_loss=args.connector_k_loss,
        pump_flow_control=not bool(args.disable_pump_flow_control),
        fluid_solid_coupling_scheme="current",
        solid_ode_method=args.solid_ode_method,
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
        set_v13_pump_total_head(build, float(args.pump_total_head_pa))
    else:
        migration = inject_v12_restart(build, args)

    if bool(args.enable_radiation_shield):
        active_until = None
        if args.shield_active_until_s is not None:
            active_until = float(system.global_time) + float(args.shield_active_until_s)
        shield = attach_radiator_thermal_shield(
            build,
            active_until_s=active_until,
            background_temperature_k=float(args.shield_background_temperature_k),
            shield_view_factor=float(args.shield_view_factor),
            inner_emissivity=float(args.shield_inner_emissivity),
            outer_emissivity=float(args.shield_outer_emissivity),
            conductivity_w_m_k=float(args.shield_conductivity_w_m_k),
            thickness_m=float(args.shield_thickness_m),
            solar_heat_flux_w_m2=float(args.shield_solar_heat_flux_w_m2),
            relaxation=float(args.shield_relaxation),
            model=args.shield_model,
        )
        shield.pre_step(0.0, float(system.global_time))

    core.point_reactor = None
    core.enable_tec_coupled = not bool(args.disable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    if core.enable_tec_coupled:
        configure_core_tec_circuits(core, args)
    core.update_neutronic_power(
        p_total=float(args.total_power_w),
        p_fiss=float(args.total_power_w),
        p_decay=0.0,
        alpha=1.0,
    )
    reset_v13_design_flows(build)
    set_v13_pump_total_head(build, float(args.pump_total_head_pa))
    core.post_step(0.0, float(system.global_time))
    if core.enable_tec_coupled and core.thermo_calc is not None:
        apply_wire_resistance(core, scale=float(args.wire_resistance_scale))
        core.set_thermo_update_time(float(system.global_time))
    core.pre_step(0.0, float(system.global_time))
    build["fluid_solid_coupling_scheme"] = args.fluid_solid_coupling_scheme
    build["fluid_solid_coupler_count"] = apply_fluid_solid_coupling_scheme(system, args.fluid_solid_coupling_scheme)
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
    build["tec_topology"] = str(getattr(core, "tec_topology", "series"))
    build["tec_circuit_mode"] = str(getattr(core, "tec_circuit_mode", "fixed_u"))
    build["reserved_parallel_tec_enabled"] = bool(getattr(core, "reserved_parallel_tec_enabled", False))
    build["reserved_parallel_tec_mode"] = str(args.reserved_parallel_mode)
    build["radiation_shield_enabled"] = bool(args.enable_radiation_shield)
    build["radiation_shield_model"] = args.shield_model
    build["radiation_shield_active_until_abs_s"] = (
        None
        if not bool(args.enable_radiation_shield) or args.shield_active_until_s is None
        else float(system.global_time) + float(args.shield_active_until_s)
    )
    build["migration_summary"] = migration
    return build


def write_latest_state(
    path: Path,
    build: Dict[str, Any],
    args: argparse.Namespace,
    latest_restart: Path,
    history_path: Path,
    start_time: float,
    target_time: float,
    latest_record: Dict[str, Any],
) -> None:
    latest = {
        "case_version": V13_CASE_VERSION,
        "restart_in": args.restart_in,
        "init_from_v12": args.init_from_v12,
        "restart_out": str(latest_restart),
        "history_csv": str(history_path),
        "start_time_s": float(start_time),
        "target_time_s": float(target_time),
        "absolute_time_s": float(build["system"].global_time),
        "duration_s": float(args.duration),
        "record_interval_s": float(args.record_interval),
        "restart_interval_s": float(args.restart_interval),
        "max_dt_s": float(args.max_dt),
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "pump_single_head_pa": float(build["pump_single_head_pa"]),
        "pump_flow_control": bool(build.get("pump_flow_control", False)),
        "target_flow_kg_s": float(args.target_flow_kg_s),
        "enable_pump_head_control": bool(args.enable_pump_head_control),
        "wire_resistance_scale": float(args.wire_resistance_scale),
        "wire_resistance_ohm": build["wire_resistance_ohm"],
        "tec_coupled_enabled": build["tec_coupled_enabled"],
        "tec_topology": build.get("tec_topology"),
        "tec_circuit_mode": build.get("tec_circuit_mode"),
        "reserved_parallel_tec_enabled": build.get("reserved_parallel_tec_enabled"),
        "reserved_parallel_tec_mode": build.get("reserved_parallel_tec_mode"),
        "radiation_shield_enabled": build.get("radiation_shield_enabled", False),
        "radiation_shield_model": build.get("radiation_shield_model"),
        "radiation_shield_active_until_abs_s": build.get("radiation_shield_active_until_abs_s"),
        "solid_ode_method": str(build["solid_ode_method"]),
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "latest_record": latest_record,
        "migration_summary": build.get("migration_summary"),
    }
    path.write_text(json.dumps(latest, indent=2, default=json_default), encoding="utf-8")


def maybe_adjust_pump_head(build: Dict[str, Any], args: argparse.Namespace) -> Optional[Dict[str, float]]:
    if not bool(args.enable_pump_head_control):
        return None
    if bool(build.get("pump_flow_control", False)):
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
    set_v13_pump_total_head(build, new_head)
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

    print("=== V13 CaseA closed core + TOPAZ-II pipe-fin radiator pumped loop ===", flush=True)
    print(f"case_version={V13_CASE_VERSION}", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"init_from_v12={args.init_from_v12}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)
    print(f"latest_restart={latest_restart_path}", flush=True)
    print(f"fluid_solid_coupling_scheme={build['fluid_solid_coupling_scheme']}", flush=True)
    print(
        f"radiator=n_tubes={args.n_tubes} n_axial={args.n_axial} "
        f"tube_eps={args.tube_emissivity:.3f} fin_eps={args.fin_emissivity:.3f}",
        flush=True,
    )
    print(
        f"pump_total_head={build['pump_total_head_pa']:.3f}Pa "
        f"target_flow={args.target_flow_kg_s:.6f}kg/s "
        f"pump_flow_control={bool(build.get('pump_flow_control', False))} "
        f"pump_control={bool(args.enable_pump_head_control)}",
        flush=True,
    )
    print(f"migration_summary={build.get('migration_summary')}", flush=True)

    initial_record = {
        **v13_basic_diagnostics(build),
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "wire_resistance_scale": build["wire_resistance_scale"],
        "tec_topology": build.get("tec_topology"),
        "tec_circuit_mode": build.get("tec_circuit_mode"),
        "reserved_parallel_tec_enabled": build.get("reserved_parallel_tec_enabled"),
        "reserved_parallel_tec_mode": build.get("reserved_parallel_tec_mode"),
        "radiation_shield_enabled": build.get("radiation_shield_enabled", False),
        "radiation_shield_model": build.get("radiation_shield_model"),
        "radiation_shield_active_until_abs_s": build.get("radiation_shield_active_until_abs_s"),
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
        print("V13 initialized and saved without time advancement.", flush=True)
        return

    next_record_time = min(start_time + float(args.record_interval), target_time)
    next_restart_time = min(start_time + float(args.restart_interval), target_time)
    next_control_time = min(start_time + float(args.pump_control_interval), target_time)
    fieldnames: Optional[List[str]] = None
    last_record: Dict[str, Any] = initial_record

    while float(system.global_time) < target_time - 1.0e-10:
        build["core"].update_neutronic_power(
            p_total=float(args.total_power_w),
            p_fiss=float(args.total_power_w),
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
            if not bool(build.get("reserved_parallel_tec_enabled", False)):
                passive = passive_tec_source_totals(build)
                if any(value != 0.0 for value in passive.values()):
                    raise RuntimeError(f"Ring3_Open TEC sources are not zero: {passive}")
            record = {
                **v13_basic_diagnostics(build),
                "relative_time_s": float(system.global_time) - start_time,
                "max_dt_s": float(args.max_dt),
                "inner_iter": int(args.inner_iter),
                "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
                "solid_ode_method": build["solid_ode_method"],
                "wire_resistance_scale": build["wire_resistance_scale"],
                "wire_resistance_ohm": build["wire_resistance_ohm"],
                "tec_coupled_enabled": build["tec_coupled_enabled"],
                "tec_topology": build.get("tec_topology"),
                "tec_circuit_mode": build.get("tec_circuit_mode"),
                "reserved_parallel_tec_enabled": build.get("reserved_parallel_tec_enabled"),
                "reserved_parallel_tec_mode": build.get("reserved_parallel_tec_mode"),
                "radiation_shield_enabled": build.get("radiation_shield_enabled", False),
                "radiation_shield_model": build.get("radiation_shield_model"),
                "radiation_shield_active_until_abs_s": build.get("radiation_shield_active_until_abs_s"),
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
                f"Trad_out={record['radiator_outlet_mix_t_k']:.3f}K "
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
    print("V13 closed-loop run completed.", flush=True)


if __name__ == "__main__":
    main()
