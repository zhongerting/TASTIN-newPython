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
    TOTAL_POWER_W,
    apply_wire_resistance,
    apply_solid_ode_method,
    get_wire_resistance,
    json_default,
    passive_tec_source_totals,
    parse_solid_ode_method,
    parse_v8_multipliers,
)
from test_core_assemble_v12_caseA import (
    V12_CASE_VERSION,
    V12_DEFAULT_INLET_TEMPERATURE_K,
    V12_DEFAULT_OUTLET_PRESSURE_PA,
    build_v12_case_a_system,
    reset_v12_design_flows,
    v12_basic_diagnostics,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V12 CaseA open core + TOPAZ-II pipe-fin radiator model.")
    parser.add_argument("--restart-in", default=None)
    parser.add_argument("--output-dir", default="testModule/v12_caseA_open_loop_no_tec")
    parser.add_argument("--case-prefix", default="v12_caseA_open_loop_no_tec")
    parser.add_argument("--duration", type=float, default=10.0)
    parser.add_argument("--record-interval", type=float, default=1.0)
    parser.add_argument("--restart-interval", type=float, default=1.0)
    parser.add_argument("--max-dt", type=float, default=0.05)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--init-dt", type=float, default=0.05)
    parser.add_argument("--hydraulic-tol", type=float, default=1.0e-6)
    parser.add_argument("--hydraulic-max-iter", type=int, default=800)
    parser.add_argument("--convergence-tol", type=float, default=1.0e-3)
    parser.add_argument("--create-init-only", action="store_true")

    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--inlet-temperature-k", type=float, default=V12_DEFAULT_INLET_TEMPERATURE_K)
    parser.add_argument("--total-inlet-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--outlet-pressure-pa", type=float, default=V12_DEFAULT_OUTLET_PRESSURE_PA)
    parser.add_argument("--connector-volume-m3", type=float, default=1.0e-5)
    parser.add_argument("--connector-length-m", type=float, default=0.02)
    parser.add_argument("--coolant-material", default=DEFAULT_COOLANT_MATERIAL)
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=False),
        default=parse_v8_multipliers("1,6,12,15,3"),
    )
    parser.add_argument("--total-power-w", type=float, default=TOTAL_POWER_W)
    parser.add_argument("--enable-tec-coupled", action="store_true")
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--thermo-update-interval", type=float, default=0.5)
    parser.add_argument("--wire-resistance-scale", type=float, default=1.0)
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=True),
        default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True),
    )

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
    parser.add_argument("--radiator-header-k-loss", type=float, default=1.0)
    parser.add_argument("--radiator-tube-inlet-k-loss", type=float, default=100.0)
    parser.add_argument("--radiator-tube-outlet-k-loss", type=float, default=100.0)
    parser.add_argument("--connector-k-loss", type=float, default=0.0)
    parser.add_argument(
        "--fluid-solid-coupling-scheme",
        choices=("current", "local_implicit"),
        default="current",
    )
    parser.add_argument(
        "--solid-ode-method",
        type=parse_solid_ode_method,
        default=parse_solid_ode_method("RK45"),
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


def enforce_inlet_boundary_temperature(build: Dict[str, Any], temperature_k: float) -> None:
    inlet = build["inlet_boundary"]
    temperature = float(temperature_k)
    set_boundary_state = getattr(inlet, "set_boundary_state", None)
    if callable(set_boundary_state):
        set_boundary_state(T=temperature)
    else:
        inlet.T = temperature
        if getattr(inlet, "material", None) is not None:
            inlet.h = inlet.material.enthalpy(inlet.T, inlet.P)
            inlet.update_properties(inlet.material)
    net = build["system"].fluid_solver
    for idx, vol in enumerate(getattr(net, "volumes_obj", [])):
        if vol is inlet:
            net.T_vec[idx] = float(inlet.T)
            net.h_vec[idx] = float(inlet.h)
            net.rho_vec[idx] = float(inlet.rho)
            net.mu_vec[idx] = float(getattr(inlet, "mu", net.mu_vec[idx]))
            if hasattr(net, "T_backup"):
                net.T_backup[idx] = net.T_vec[idx]
            if hasattr(net, "h_backup"):
                net.h_backup[idx] = net.h_vec[idx]
            if hasattr(net, "rho_backup"):
                net.rho_backup[idx] = net.rho_vec[idx]
            if hasattr(net, "mu_backup"):
                net.mu_backup[idx] = net.mu_vec[idx]
            break
    build["target_core_inlet_t_k"] = temperature


def build_case(args: argparse.Namespace) -> Dict[str, Any]:
    init_coupling_scheme = "current" if args.fluid_solid_coupling_scheme == "local_implicit" else args.fluid_solid_coupling_scheme
    build = build_v12_case_a_system(
        inlet_temperature_k=float(args.inlet_temperature_k),
        total_inlet_flow_kg_s=float(args.total_inlet_flow_kg_s),
        outlet_pressure_pa=float(args.outlet_pressure_pa),
        pipe_n_nodes=int(args.pipe_n_nodes),
        connector_volume_m3=float(args.connector_volume_m3),
        connector_length_m=float(args.connector_length_m),
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        enable_tec_coupled=bool(args.enable_tec_coupled),
        n_tubes=int(args.n_tubes),
        n_axial=int(args.n_axial),
        n_radial_wall=int(args.n_radial_wall),
        n_fin_width=int(args.n_fin_width),
        tube_length_m=float(args.tube_length_m),
        tube_inner_diameter_m=float(args.tube_inner_diameter_m),
        tube_outer_diameter_m=float(args.tube_outer_diameter_m),
        upper_header_centerline_diameter_m=float(args.upper_header_centerline_diameter_m),
        lower_header_centerline_diameter_m=float(args.lower_header_centerline_diameter_m),
        header_inner_diameter_m=float(args.header_inner_diameter_m),
        fin_thickness_m=float(args.fin_thickness_m),
        fin_width_upper_m=float(args.fin_width_upper_m),
        fin_width_lower_m=float(args.fin_width_lower_m),
        tube_emissivity=float(args.tube_emissivity),
        fin_emissivity=float(args.fin_emissivity),
        tube_area_scale=float(args.tube_area_scale),
        fin_area_scale=float(args.fin_area_scale),
        t_space_k=float(args.t_space_k),
        fin_conductivity_w_m_k=float(args.fin_conductivity_w_m_k),
        fin_view_factor=float(args.fin_view_factor),
        fin_contact_resistance_m2k_w=float(args.fin_contact_resistance_m2k_w),
        radiator_header_k_loss=float(args.radiator_header_k_loss),
        radiator_tube_inlet_k_loss=float(args.radiator_tube_inlet_k_loss),
        radiator_tube_outlet_k_loss=float(args.radiator_tube_outlet_k_loss),
        connector_k_loss=float(args.connector_k_loss),
        fluid_solid_coupling_scheme=init_coupling_scheme,
        solid_ode_method=args.solid_ode_method,
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = bool(args.enable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    if not args.enable_tec_coupled:
        core.thermo_calc = None
    apply_solid_ode_method(build, args.solid_ode_method)
    reset_v12_design_flows(build)
    system.initialize_system(
        dt_init=float(args.init_dt),
        tol=float(args.hydraulic_tol),
        max_iter=int(args.hydraulic_max_iter),
    )
    enforce_inlet_boundary_temperature(build, float(args.inlet_temperature_k))
    fluid_solid_coupler_count = apply_fluid_solid_coupling_scheme(system, args.fluid_solid_coupling_scheme)
    if args.restart_in:
        system.load_global_state(args.restart_in)
        enforce_inlet_boundary_temperature(build, float(args.inlet_temperature_k))
        apply_solid_ode_method(build, args.solid_ode_method)
    core.point_reactor = None
    core.enable_tec_coupled = bool(args.enable_tec_coupled)
    core.thermo_update_interval = float(args.thermo_update_interval)
    if not args.enable_tec_coupled:
        core.thermo_calc = None
    elif core.thermo_calc is not None:
        core.setup_tec_circuit("fixed_u", float(args.target_voltage), I_guess=150.0)
    core.update_neutronic_power(
        p_total=float(args.total_power_w),
        p_fiss=float(args.total_power_w),
        p_decay=0.0,
        alpha=1.0,
    )
    reset_v12_design_flows(build)
    core.post_step(0.0, float(system.global_time))
    if core.enable_tec_coupled and core.thermo_calc is not None:
        apply_wire_resistance(core, scale=float(args.wire_resistance_scale))
        core._last_thermo_update_time = float(system.global_time)
    core.pre_step(0.0, float(system.global_time))
    build["solid_ode_method"] = args.solid_ode_method
    build["fluid_solid_coupling_scheme"] = args.fluid_solid_coupling_scheme
    build["fluid_solid_coupler_count"] = fluid_solid_coupler_count
    build["tec_coupled_enabled"] = bool(core.enable_tec_coupled)
    build["target_voltage_v"] = float(args.target_voltage)
    build["thermo_update_interval_s"] = float(args.thermo_update_interval)
    build["wire_resistance_scale"] = float(args.wire_resistance_scale)
    build["wire_resistance_ohm"] = get_wire_resistance(core)
    return build


def write_latest_state(
    *,
    path: Path,
    args: argparse.Namespace,
    restart_path: Path,
    history_path: Path,
    start_time: float,
    target_time: float,
    build: Dict[str, Any],
    latest_record: Dict[str, Any],
) -> None:
    latest = {
        "case_version": V12_CASE_VERSION,
        "restart_in": args.restart_in,
        "restart_out": str(restart_path),
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
        "solid_ode_method": build["solid_ode_method"],
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "fluid_solid_coupler_count": build["fluid_solid_coupler_count"],
        "tec_coupled_enabled": build["tec_coupled_enabled"],
        "target_voltage_v": build["target_voltage_v"],
        "thermo_update_interval_s": build["thermo_update_interval_s"],
        "wire_resistance_scale": build["wire_resistance_scale"],
        "wire_resistance_ohm": build["wire_resistance_ohm"],
        "ring_multipliers": build["ring_multipliers"],
        "tec_ring_multipliers": build["tec_ring_multipliers"],
        "radiator_geometry": build["radiator_geometry"],
        "flow_network_pipe_specs": build["flow_network_pipe_specs"],
        "latest": latest_record,
    }
    with path.open("w", encoding="utf-8") as f:
        json.dump(latest, f, indent=2, ensure_ascii=False, default=json_default)


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0 and not args.create_init_only:
        raise ValueError("duration must be positive unless --create-init-only is used.")
    if args.record_interval <= 0.0:
        raise ValueError("record-interval must be positive.")
    if args.restart_interval <= 0.0:
        raise ValueError("restart-interval must be positive.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_history.csv"
    latest_state_path = output_dir / f"{args.case_prefix}_latest_state.json"
    latest_restart_path = output_dir / f"{args.case_prefix}_latest_restart.npz"

    build = build_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time if args.create_init_only else start_time + float(args.duration)
    fieldnames: Optional[List[str]] = None
    last_record = {
        **v12_basic_diagnostics(build),
        "relative_time_s": 0.0,
        "max_dt_s": float(args.max_dt),
        "inner_iter": int(args.inner_iter),
        "solid_ode_method": build["solid_ode_method"],
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "fluid_solid_coupler_count": build["fluid_solid_coupler_count"],
    }

    print("=== V12 CaseA open core + TOPAZ-II pipe-fin radiator ===", flush=True)
    print(f"case_version={V12_CASE_VERSION}", flush=True)
    print(f"restart_in={args.restart_in}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"history_csv={history_path}", flush=True)
    print(f"latest_restart={latest_restart_path}", flush=True)
    print(
        f"tec_coupled_enabled={build['tec_coupled_enabled']}, "
        f"thermo_update_interval={build['thermo_update_interval_s']:.3f}s",
        flush=True,
    )

    if args.create_init_only:
        system.save_global_state(str(latest_restart_path))
        write_latest_state(
            path=latest_state_path,
            args=args,
            restart_path=latest_restart_path,
            history_path=history_path,
            start_time=start_time,
            target_time=target_time,
            build=build,
            latest_record=last_record,
        )
        print(json.dumps(last_record, indent=2, ensure_ascii=False, default=json_default), flush=True)
        print("V12 init-only run completed.", flush=True)
        return

    next_record_time = min(start_time + float(args.record_interval), target_time)
    next_restart_time = min(start_time + float(args.restart_interval), target_time)
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
            target_time - float(system.global_time),
        )
        system.step(
            dt,
            inner_iter=int(args.inner_iter),
            convergence_tol=float(args.convergence_tol),
        )
        passive = passive_tec_source_totals(build)
        if any(value != 0.0 for value in passive.values()):
            raise RuntimeError(f"Passive TFE TEC sources are not zero: {passive}")

        if float(system.global_time) >= next_record_time - 1.0e-10:
            record = {
                **v12_basic_diagnostics(build),
                "relative_time_s": float(system.global_time) - start_time,
                "max_dt_s": float(args.max_dt),
                "inner_iter": int(args.inner_iter),
                "solid_ode_method": build["solid_ode_method"],
                "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
                "fluid_solid_coupler_count": build["fluid_solid_coupler_count"],
            }
            flat = flatten_for_csv(record)
            if fieldnames is None:
                fieldnames = list(flat.keys())
            append_row(history_path, fieldnames, flat, write_header=not history_path.exists())
            last_record = record
            pel = float(record["tec_total_electric_power_w"] or 0.0)
            print(
                f"t_rel={record['relative_time_s']:.3f}s "
                f"Tin={record['core_inlet_connector_t_k']:.3f}K "
                f"Tcore_out={record['core_outlet_connector_t_k']:.3f}K "
                f"Trad_out={record['radiator_outlet_mix_t_k']:.3f}K "
                f"tubeWmean={record['radiator_tube_mean_flow_kg_s']:.6f}kg/s "
                f"Qrad={record['q_radiator_total_w']:.3f}W "
                f"Pel={pel:.3f}W",
                flush=True,
            )
            next_record_time = min(next_record_time + float(args.record_interval), target_time)

        if float(system.global_time) >= next_restart_time - 1.0e-10:
            system.save_global_state(str(latest_restart_path))
            write_latest_state(
                path=latest_state_path,
                args=args,
                restart_path=latest_restart_path,
                history_path=history_path,
                start_time=start_time,
                target_time=target_time,
                build=build,
                latest_record=last_record,
            )
            next_restart_time = min(next_restart_time + float(args.restart_interval), target_time)

    if not latest_restart_path.exists():
        system.save_global_state(str(latest_restart_path))
    write_latest_state(
        path=latest_state_path,
        args=args,
        restart_path=latest_restart_path,
        history_path=history_path,
        start_time=start_time,
        target_time=target_time,
        build=build,
        latest_record=last_record,
    )
    print("V12 open-loop run completed.", flush=True)


if __name__ == "__main__":
    main()
