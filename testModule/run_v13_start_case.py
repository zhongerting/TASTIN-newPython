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

from run_v13_caseA_closed_loop import (  # noqa: E402
    apply_fluid_solid_coupling_scheme,
    flatten_for_csv,
    json_default,
    maybe_adjust_pump_head,
)
from run_v8_caseA_common import (  # noqa: E402
    DEFAULT_COOLANT_MATERIAL,
    DEFAULT_SOLID_ODE_METHOD,
    apply_solid_ode_method,
    apply_wire_resistance,
    get_solid_ode_methods,
    WIRE_RESISTANCE_OHM,
    get_wire_resistance,
    parse_solid_ode_method,
    parse_v8_multipliers,
)
from test_core_assemble_v13_caseA import (  # noqa: E402
    V13_DEFAULT_PUMP_TOTAL_HEAD_PA,
    V13_DEFAULT_REFERENCE_PRESSURE_PA,
    attach_radiator_thermal_shield,
    build_v13_case_a_system,
    reset_v13_design_flows,
    set_v13_pump_total_head,
    v13_basic_diagnostics,
)
from v13_startup_control import (  # noqa: E402
    V13StartupControlConfig,
    V13StartupController,
    apply_tec_gap_h_eq,
    attach_radiator_tube_external_heat,
    reset_solid_temperatures,
    radiator_external_heat_power_w,
    reset_fluid_temperatures,
    shield_qsss_from_matrix,
)


V13_START_CASE_VERSION = "v13_start_cold_startup_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run V13-start TOPAZ-II cold startup case.", allow_abbrev=False)
    parser.add_argument("--output-dir", default="testModule/v13_start_cold_start_smoke")
    parser.add_argument("--case-prefix", default="v13_start_cold_start_smoke")
    parser.add_argument("--restart-in", default=None)
    parser.add_argument("--duration", type=float, default=10.0)
    parser.add_argument("--record-interval", type=float, default=1.0)
    parser.add_argument("--restart-interval", type=float, default=10.0)
    parser.add_argument("--max-dt", type=float, default=0.1)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--init-dt", type=float, default=0.05)
    parser.add_argument("--hydraulic-tol", type=float, default=1.0e-6)
    parser.add_argument("--hydraulic-max-iter", type=int, default=800)
    parser.add_argument("--step-hydraulic-max-iter", type=int, default=300)
    parser.add_argument("--initial-temperature-k", type=float, default=373.0)
    parser.add_argument("--solid-initial-temperature-k", type=float, default=None)
    parser.add_argument("--startup-profile", choices=("engineering", "titam"), default="engineering")
    parser.add_argument("--reference-pressure-pa", type=float, default=V13_DEFAULT_REFERENCE_PRESSURE_PA)
    parser.add_argument("--pump-total-head-pa", type=float, default=V13_DEFAULT_PUMP_TOTAL_HEAD_PA)
    parser.add_argument("--target-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--total-inlet-flow-kg-s", type=float, default=1.3)
    parser.add_argument("--enable-pump-head-control", action="store_true")
    parser.add_argument("--disable-pump-flow-control", action="store_true")
    parser.add_argument("--pump-control-interval", type=float, default=50.0)
    parser.add_argument("--pump-control-max-fraction", type=float, default=0.10)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--connector-volume-m3", type=float, default=1.0e-5)
    parser.add_argument("--connector-length-m", type=float, default=0.02)
    parser.add_argument("--coolant-material", default=DEFAULT_COOLANT_MATERIAL)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--startup-main-tec-initial-mode", choices=("fixed_r", "fixed_u"), default="fixed_r")
    parser.add_argument("--startup-main-tec-load-resistance-ohm", type=float, default=0.0044)
    parser.add_argument("--startup-main-tec-load-resistance-scope", choices=("total", "per_tec"), default="total")
    parser.add_argument("--startup-main-tec-switch-voltage-v", type=float, default=27.2)
    parser.add_argument("--startup-main-tec-i-guess-a", type=float, default=150.0)
    parser.add_argument("--enable-reserved-parallel-tec", action="store_true")
    parser.add_argument("--reserved-parallel-mode", choices=("fixed_u", "fixed_i", "load_curve"), default="fixed_u")
    parser.add_argument("--reserved-parallel-voltage", type=float, default=0.8)
    parser.add_argument("--reserved-parallel-current", type=float, default=6000.0)
    parser.add_argument("--reserved-parallel-load-curve", default=None)
    parser.add_argument("--thermo-update-interval", type=float, default=0.5)
    parser.add_argument("--wire-resistance-scale", type=float, default=0.5)
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
    parser.add_argument("--shield-model", choices=("segment_balance", "fortran_shield2"), default="fortran_shield2")
    parser.add_argument("--shield-inner-emissivity", type=float, default=0.8)
    parser.add_argument("--shield-outer-emissivity", type=float, default=0.1)
    parser.add_argument("--shield-conductivity-w-m-k", type=float, default=0.0008)
    parser.add_argument("--shield-thickness-m", type=float, default=0.01)
    parser.add_argument("--shield-view-factor", type=float, default=0.8)
    parser.add_argument("--shield-background-temperature-k", type=float, default=3.0)
    parser.add_argument("--shield-relaxation", type=float, default=1.0)
    parser.add_argument("--tube-external-heat-matrix-key", default="is58p5_w0_8p12_N78_sum")
    parser.add_argument("--shield-external-heat-matrix-key", default="is58p5_w0_8p12_N6_sum")
    parser.add_argument("--external-heat-scale", type=float, default=1.0)
    parser.add_argument("--tube-external-heat-area-fraction", type=float, default=1.0)
    parser.add_argument("--tube-external-heat-units", choices=("W/m", "W/m2"), default="W/m2")
    parser.add_argument("--radiator-header-k-loss", type=float, default=1.0)
    parser.add_argument("--radiator-tube-inlet-k-loss", type=float, default=100.0)
    parser.add_argument("--radiator-tube-outlet-k-loss", type=float, default=100.0)
    parser.add_argument("--connector-k-loss", type=float, default=0.0)
    parser.add_argument("--fluid-solid-coupling-scheme", choices=("current", "local_implicit"), default="local_implicit")
    parser.add_argument("--solid-ode-method", type=parse_solid_ode_method, default=parse_solid_ode_method("implicit_euler"))
    parser.add_argument("--ring-multipliers", type=lambda text: parse_v8_multipliers(text, allow_zero=False), default=parse_v8_multipliers("1,6,12,15,3"))
    parser.add_argument("--tec-ring-multipliers", type=lambda text: parse_v8_multipliers(text, allow_zero=True), default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True))
    parser.add_argument("--source-power-w", type=float, default=1.0)
    parser.add_argument("--fixed-startup-power-w", type=float, default=None)
    parser.add_argument("--shield-jettison-temperature-k", type=float, default=400.0)
    parser.add_argument("--tfe-start-after-critical-s", type=float, default=1500.0)
    parser.add_argument("--tfe-start-emitter-temperature-k", type=float, default=1050.0)
    parser.add_argument("--tec-electrical-start-after-cesium-s", type=float, default=0.0)
    parser.add_argument("--tec-electrical-start-cs-fraction", type=float, default=0.0)
    parser.add_argument("--tec-electrical-start-emitter-temperature-k", type=float, default=0.0)
    parser.add_argument("--helium-gap-h-eq-w-m2-k", type=float, default=1200.0)
    parser.add_argument("--cesium-gap-h-eq-w-m2-k", type=float, default=29.0)
    parser.add_argument("--cs-transition-tau-s", type=float, default=120.0)
    parser.add_argument("--initial-cs-fraction", type=float, default=0.0)
    args = parser.parse_args()
    if args.startup_profile == "titam":
        args.initial_temperature_k = 300.0
        args.total_inlet_flow_kg_s = 1.5
        args.target_flow_kg_s = 1.5
    return args


def append_row(path: Path, fieldnames: List[str], row: Dict[str, Any], *, write_header: bool) -> None:
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def mean_emitter_temperature_k(core: Any) -> float:
    values = []
    for tfe in getattr(core, "tfes", {}).values():
        emitter = getattr(tfe, "solids", {}).get("emitter")
        if emitter is None:
            continue
        boundary = emitter.boundaries.get("right")
        if boundary is not None and hasattr(boundary, "T_surface"):
            values.extend(np.ravel(np.asarray(boundary.T_surface, dtype=float)).tolist())
        elif hasattr(emitter, "T"):
            values.extend(np.ravel(np.asarray(emitter.T, dtype=float)).tolist())
    if not values:
        return float("nan")
    return float(np.mean(values))





def apply_startup_wire_resistance_without_calculate(core: Any, scale: float = 1.0) -> None:
    """Set main TEC wire resistance for startup without forcing an unsynchronized first calculate."""
    if getattr(core, "thermo_calc", None) is None:
        return
    wire_res = np.asarray(WIRE_RESISTANCE_OHM, dtype=float) * float(scale)
    thermo_calc = core.thermo_calc
    thermo_calc._input_data.resistanceWire = np.tile(wire_res, (thermo_calc.N_elem, 1))
    thermo_calc.build()


def _main_tec_series_element_count(core: Any) -> int:
    thermo_calc = getattr(core, "thermo_calc", None)
    if thermo_calc is not None and hasattr(thermo_calc, "N_elem"):
        return max(1, int(getattr(thermo_calc, "N_elem")))
    groups = getattr(core, "tec_circuit_groups", {})
    main_group = groups.get("main") if isinstance(groups, dict) else None
    if main_group is not None and hasattr(main_group, "total_virtual_elements"):
        return max(1, int(getattr(main_group, "total_virtual_elements")))
    return 1


def effective_startup_main_tec_load_resistance(core: Any, args: argparse.Namespace) -> float:
    resistance = float(getattr(args, "startup_main_tec_load_resistance_ohm", 0.0044))
    scope = str(getattr(args, "startup_main_tec_load_resistance_scope", "total")).lower()
    if scope == "total":
        return resistance
    if scope == "per_tec":
        return resistance * float(_main_tec_series_element_count(core))
    raise ValueError(f"Unsupported startup main TEC load resistance scope: {scope!r}")


def configure_startup_main_tec_circuit(core: Any, args: argparse.Namespace) -> None:
    """Configure the startup main series TEC circuit before first coupling."""
    mode = str(getattr(args, "startup_main_tec_initial_mode", "fixed_r")).lower()
    if mode == "fixed_r":
        target = effective_startup_main_tec_load_resistance(core, args)
    elif mode == "fixed_u":
        target = float(getattr(args, "target_voltage", 27.2))
    else:
        raise ValueError(f"Unsupported startup main TEC initial mode: {mode!r}")
    i_guess = float(getattr(args, "startup_main_tec_i_guess_a", 150.0))
    core.setup_tec_circuit(mode, target, I_guess=i_guess, topology="series")


def _finite_record_value(record: Dict[str, Any], *keys: str) -> Optional[float]:
    for key in keys:
        value = record.get(key)
        if value is None:
            continue
        try:
            value_f = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value_f):
            return value_f
    return None


def maybe_switch_startup_tec_to_fixed_voltage(
        build: Dict[str, Any],
        args: argparse.Namespace,
        record: Dict[str, Any]) -> bool:
    """Switch startup TEC from load resistance to fixed voltage once bus voltage reaches target."""
    if bool(build.get("startup_main_tec_switched_to_fixed_voltage", False)):
        return False
    if str(getattr(args, "startup_main_tec_initial_mode", "fixed_r")).lower() != "fixed_r":
        return False
    voltage = _finite_record_value(record, "tec_main_voltage_v", "tec_total_voltage_v")
    threshold = float(getattr(args, "startup_main_tec_switch_voltage_v", getattr(args, "target_voltage", 27.2)))
    if voltage is None or voltage < threshold:
        return False
    current = _finite_record_value(record, "tec_main_current_a", "tec_total_current_a")
    if current is None or current <= 0.0:
        current = float(getattr(args, "startup_main_tec_i_guess_a", 150.0))
    core = build["core"]
    system = build.get("system")
    time_s = float(getattr(system, "global_time", float("nan")))
    target_voltage = float(getattr(args, "target_voltage", threshold))
    time_text = "unknown" if not np.isfinite(time_s) else f"{time_s:.3f}"
    print(
        "switching startup main TEC to fixed voltage "
        f"at t={time_text} s: "
        f"U={voltage:.6g} V -> target={target_voltage:.6g} V, I_guess={current:.6g} A",
        flush=True,
    )
    core.setup_tec_circuit("fixed_u", target_voltage, I_guess=current, topology="series")
    if system is not None and hasattr(core, "set_thermo_update_time"):
        core.set_thermo_update_time(time_s - float(getattr(args, "thermo_update_interval", 0.5)))
    switch_time = None if not np.isfinite(time_s) else time_s
    build["startup_main_tec_switched_to_fixed_voltage"] = True
    build["startup_main_tec_switch_time_s"] = switch_time
    record["startup_main_tec_switched_to_fixed_voltage"] = True
    record["startup_main_tec_switch_time_s"] = switch_time
    return True


def startup_config_from_args(args: argparse.Namespace) -> V13StartupControlConfig:
    return V13StartupControlConfig(
        source_power_w=float(args.source_power_w),
        fixed_power_w=None if args.fixed_startup_power_w is None else float(args.fixed_startup_power_w),
        initial_temperature_k=float(args.initial_temperature_k),
        shield_jettison_temperature_k=float(args.shield_jettison_temperature_k),
        tfe_start_after_critical_s=float(args.tfe_start_after_critical_s),
        tfe_start_emitter_temperature_k=float(args.tfe_start_emitter_temperature_k),
        tec_electrical_start_after_cesium_s=float(args.tec_electrical_start_after_cesium_s),
        tec_electrical_start_cs_fraction=float(args.tec_electrical_start_cs_fraction),
        tec_electrical_start_emitter_temperature_k=float(args.tec_electrical_start_emitter_temperature_k),
        helium_gap_h_eq_w_m2_k=float(args.helium_gap_h_eq_w_m2_k),
        cesium_gap_h_eq_w_m2_k=float(args.cesium_gap_h_eq_w_m2_k),
        cs_transition_tau_s=float(args.cs_transition_tau_s),
    )


def build_case(args: argparse.Namespace) -> Dict[str, Any]:
    print("building cold V13-start system...", flush=True)
    build = build_v13_case_a_system(
        inlet_temperature_k=float(args.initial_temperature_k),
        total_inlet_flow_kg_s=float(args.total_inlet_flow_kg_s),
        reference_pressure_pa=float(args.reference_pressure_pa),
        pump_total_head_pa=float(args.pump_total_head_pa),
        pipe_n_nodes=int(args.pipe_n_nodes),
        connector_volume_m3=float(args.connector_volume_m3),
        connector_length_m=float(args.connector_length_m),
        coolant_material=args.coolant_material,
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
        enable_tec_coupled=True,
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
        pump_flow_control=not bool(args.disable_pump_flow_control),
        fluid_solid_coupling_scheme="current",
        solid_ode_method=args.solid_ode_method,
    )
    system = build["system"]
    core = build["core"]
    apply_solid_ode_method(build, args.solid_ode_method)
    core.enable_tec_coupled = False
    core.thermo_update_interval = float(args.thermo_update_interval)
    apply_tec_gap_h_eq(core, float(args.helium_gap_h_eq_w_m2_k))

    shield = attach_radiator_thermal_shield(
        build,
        active_until_s=None,
        background_temperature_k=float(args.shield_background_temperature_k),
        shield_view_factor=float(args.shield_view_factor),
        inner_emissivity=float(args.shield_inner_emissivity),
        outer_emissivity=float(args.shield_outer_emissivity),
        conductivity_w_m_k=float(args.shield_conductivity_w_m_k),
        thickness_m=float(args.shield_thickness_m),
        solar_heat_flux_w_m2=0.0,
        relaxation=float(args.shield_relaxation),
        model=args.shield_model,
    )
    build["radiator_thermal_shield"] = shield
    build["radiator_tube_external_heat_count"] = attach_radiator_tube_external_heat(
        build["radiator_units"],
        matrix_key=str(args.tube_external_heat_matrix_key),
        scale_factor=float(args.external_heat_scale),
        periodic=True,
        area_fraction=float(args.tube_external_heat_area_fraction),
        input_units=str(args.tube_external_heat_units),
    )

    solid_initial_temperature = (
        float(args.initial_temperature_k)
        if args.solid_initial_temperature_k is None
        else float(args.solid_initial_temperature_k)
    )
    build["solid_initial_temperature_k"] = solid_initial_temperature
    if args.restart_in:
        build["cold_initialized_fluid_count"] = 0
        build["cold_initialized_solid_count"] = 0
        print(f"loading V13-start restart: {args.restart_in}", flush=True)
        system.load_global_state(args.restart_in)
        print(f"V13-start restart loaded at t={system.global_time:.6f} s.", flush=True)
        reset_v13_design_flows(build)
        set_v13_pump_total_head(build, float(args.pump_total_head_pa))
    else:
        build["cold_initialized_fluid_count"] = reset_fluid_temperatures(system, solid_initial_temperature)
        build["cold_initialized_solid_count"] = reset_solid_temperatures(system, solid_initial_temperature)
        for unit in build.get("radiator_units", []):
            if hasattr(unit, "last_fin_temperature"):
                unit.last_fin_temperature[...] = solid_initial_temperature
            if hasattr(unit, "last_fin_effective_temperature_distribution"):
                unit.last_fin_effective_temperature_distribution[...] = solid_initial_temperature
            if hasattr(unit, "_has_valid_fin_temperature"):
                unit._has_valid_fin_temperature = False
        reset_v13_design_flows(build)
        set_v13_pump_total_head(build, float(args.pump_total_head_pa))
        print("initializing cold V13-start hydraulic/thermal system...", flush=True)
        system.initialize_system(dt_init=float(args.init_dt), tol=float(args.hydraulic_tol), max_iter=int(args.hydraulic_max_iter))
        print("cold V13-start initialization complete.", flush=True)
    core.update_neutronic_power(p_total=float(args.source_power_w), p_fiss=float(args.source_power_w), p_decay=0.0, alpha=1.0)
    build["fluid_solid_coupling_scheme"] = args.fluid_solid_coupling_scheme
    build["fluid_solid_coupler_count"] = apply_fluid_solid_coupling_scheme(system, args.fluid_solid_coupling_scheme)
    build["solid_ode_method"] = args.solid_ode_method
    build["solid_ode_methods"] = get_solid_ode_methods(build)
    build["wire_resistance_scale"] = float(args.wire_resistance_scale)
    build["wire_resistance_ohm"] = []
    build["startup_controller"] = V13StartupController(startup_config_from_args(args))
    if float(args.initial_cs_fraction) > 0.0:
        build["startup_controller"].seed_cesium_conditioning(float(system.global_time), float(args.initial_cs_fraction))
    build["tec_has_been_enabled"] = False
    return build


def prepare_startup_step(build: Dict[str, Any], args: argparse.Namespace):
    system = build["system"]
    core = build["core"]
    controller: V13StartupController = build["startup_controller"]
    inlet_t = float(build["core_inlet_connector"].T)
    emitter_t = mean_emitter_temperature_k(core)
    command = controller.evaluate(
        float(system.global_time),
        core_inlet_temperature_k=inlet_t,
        emitter_temperature_k=emitter_t,
    )
    shield = build.get("radiator_thermal_shield")
    if shield is not None:
        shield.qsss_w_m2 = shield_qsss_from_matrix(
            str(args.shield_external_heat_matrix_key),
            float(system.global_time),
            scale_factor=float(args.external_heat_scale),
        )
        if command.radiation_shield_active:
            shield.active_until_s = None
        else:
            shield.active_until_s = min(float(system.global_time), float(system.global_time) - 1.0e-9)
    apply_tec_gap_h_eq(core, command.tec_gap_h_eq_w_m2_k)
    if command.tec_enabled and not bool(build.get("tec_has_been_enabled", False)):
        print(f"enabling TEC coupling at t={system.global_time:.3f} s...", flush=True)
        configure_startup_main_tec_circuit(core, args)
        apply_startup_wire_resistance_without_calculate(core, scale=float(args.wire_resistance_scale))
        build["wire_resistance_ohm"] = get_wire_resistance(core)
        core.enable_tec_coupled = True
        core.set_thermo_update_time(float(system.global_time) - float(args.thermo_update_interval))
        build["tec_has_been_enabled"] = True
        build["startup_main_tec_switched_to_fixed_voltage"] = str(getattr(args, "startup_main_tec_initial_mode", "fixed_r")).lower() == "fixed_u"
        build["startup_main_tec_switch_time_s"] = None
    elif not command.tec_enabled:
        core.enable_tec_coupled = False
    core.update_neutronic_power(
        p_total=command.thermal_power_w,
        p_fiss=command.fission_power_w,
        p_decay=command.decay_power_w,
        alpha=1.0,
    )
    return command


def startup_tec_global_diagnostics(core: Any) -> Dict[str, Any]:
    thermo_calc = getattr(core, "thermo_calc", None)
    if thermo_calc is None or not hasattr(thermo_calc, "get_global_results"):
        return {
            "tec_solver_mode": None,
            "tec_solver_converged": None,
            "tec_solver_iteration_count": None,
            "tec_solver_zero_emission_skipped": None,
            "tec_solver_zero_emission_reason": None,
            "tec_solver_output_finite": None,
        }
    results = thermo_calc.get_global_results()
    if results is None:
        return {
            "tec_solver_mode": None,
            "tec_solver_converged": None,
            "tec_solver_iteration_count": None,
            "tec_solver_zero_emission_skipped": None,
            "tec_solver_zero_emission_reason": None,
            "tec_solver_output_finite": None,
        }
    finite_values = []
    for key in ("Iout", "Uout", "Rload"):
        try:
            finite_values.append(np.isfinite(float(results.get(key))))
        except (TypeError, ValueError):
            finite_values.append(False)
    return {
        "tec_solver_mode": results.get("mode"),
        "tec_solver_converged": bool(results.get("converged")) if results.get("converged") is not None else None,
        "tec_solver_iteration_count": int(results.get("iteration_count")) if results.get("iteration_count") is not None else None,
        "tec_solver_zero_emission_skipped": bool(results.get("zero_emission_skipped")) if results.get("zero_emission_skipped") is not None else None,
        "tec_solver_zero_emission_reason": results.get("zero_emission_reason"),
        "tec_solver_output_finite": bool(all(finite_values)),
    }



def startup_energy_residual_diagnostics(record: Dict[str, Any]) -> Dict[str, float]:
    """Derived V13-start energy diagnostics for transient storage tracking."""
    core_heat = _finite_record_value(record, "core_heat_power_w", "startup_thermal_power_w")
    coolant_dh = _finite_record_value(record, "coolant_enthalpy_rise_w")
    electric = _finite_record_value(record, "tec_total_electric_power_w", "tec_main_electric_power_w")
    q_radiator = _finite_record_value(record, "q_radiator_total_w")
    if electric is None:
        electric = 0.0

    core_minus_coolant_minus_electric = float("nan")
    core_minus_radiator_minus_electric = float("nan")
    radiator_minus_coolant = float("nan")
    core_residual_rel = float("nan")
    radiator_balance_rel = float("nan")

    if core_heat is not None and coolant_dh is not None:
        core_minus_coolant_minus_electric = core_heat - coolant_dh - electric
        denominator = max(abs(core_heat), 1.0)
        core_residual_rel = core_minus_coolant_minus_electric / denominator
    if q_radiator is not None and coolant_dh is not None:
        radiator_minus_coolant = q_radiator - coolant_dh
        denominator = max(abs(q_radiator), abs(coolant_dh), 1.0)
        radiator_balance_rel = radiator_minus_coolant / denominator
    if core_heat is not None and q_radiator is not None:
        core_minus_radiator_minus_electric = core_heat - q_radiator - electric

    return {
        "core_heat_minus_coolant_enthalpy_minus_electric_w": float(core_minus_coolant_minus_electric),
        "core_heat_minus_radiator_minus_electric_w": float(core_minus_radiator_minus_electric),
        "radiator_minus_coolant_enthalpy_w": float(radiator_minus_coolant),
        "core_energy_storage_residual_rel": float(core_residual_rel),
        "radiator_coolant_balance_rel": float(radiator_balance_rel),
    }

def collect_record(build: Dict[str, Any], args: argparse.Namespace, command, start_time: float) -> Dict[str, Any]:
    core = build["core"]
    record = {
        **v13_basic_diagnostics(build),
        **command.diagnostics(),
        **startup_tec_global_diagnostics(build["core"]),
        "relative_time_s": float(build["system"].global_time) - float(start_time),
        "max_dt_s": float(args.max_dt),
        "solid_ode_method": str(build["solid_ode_method"]),
        "fluid_solid_coupling_scheme": build["fluid_solid_coupling_scheme"],
        "wire_resistance_scale": build["wire_resistance_scale"],
        "wire_resistance_ohm": build["wire_resistance_ohm"],
        "solid_initial_temperature_k": float(build.get("solid_initial_temperature_k", float("nan"))),
        "cold_initialized_solid_count": int(build.get("cold_initialized_solid_count", 0)),
        "cold_initialized_fluid_count": int(build.get("cold_initialized_fluid_count", 0)),
        "tec_coupled_enabled": bool(getattr(core, "enable_tec_coupled", False)),
        "tec_has_been_enabled": bool(build.get("tec_has_been_enabled", False)),
        "startup_main_tec_initial_mode": str(args.startup_main_tec_initial_mode),
        "startup_main_tec_load_resistance_ohm": float(args.startup_main_tec_load_resistance_ohm),
        "startup_main_tec_load_resistance_scope": str(args.startup_main_tec_load_resistance_scope),
        "startup_main_tec_effective_load_resistance_ohm": effective_startup_main_tec_load_resistance(core, args),
        "startup_main_tec_switch_voltage_v": float(args.startup_main_tec_switch_voltage_v),
        "startup_main_tec_switched_to_fixed_voltage": bool(build.get("startup_main_tec_switched_to_fixed_voltage", False)),
        "startup_main_tec_switch_time_s": build.get("startup_main_tec_switch_time_s"),
        "mean_emitter_temperature_k": mean_emitter_temperature_k(core),
        "target_flow_kg_s": float(args.target_flow_kg_s),
        "radiator_tube_external_heat_count": int(build.get("radiator_tube_external_heat_count", 0)),
        "radiator_tube_external_heat_w": radiator_external_heat_power_w(build.get("radiator_units", [])),
        "shield_external_qsss_mean_w_m2": float(np.mean(np.asarray(getattr(build.get("radiator_thermal_shield"), "qsss_w_m2", np.zeros(8)), dtype=float))),
        "shield_external_qsss_max_w_m2": float(np.max(np.asarray(getattr(build.get("radiator_thermal_shield"), "qsss_w_m2", np.zeros(8)), dtype=float))),
    }
    record.update(startup_energy_residual_diagnostics(record))
    return record


def write_latest_state(path: Path, build: Dict[str, Any], args: argparse.Namespace, latest_restart: Path, history_path: Path, start_time: float, target_time: float, latest_record: Dict[str, Any]) -> None:
    latest = {
        "case_version": V13_START_CASE_VERSION,
        "restart_out": str(latest_restart),
        "history_csv": str(history_path),
        "start_time_s": float(start_time),
        "target_time_s": float(target_time),
        "absolute_time_s": float(build["system"].global_time),
        "duration_s": float(args.duration),
        "initial_temperature_k": float(args.initial_temperature_k),
        "startup_profile": str(args.startup_profile),
        "max_dt_s": float(args.max_dt),
        "pump_total_head_pa": float(build["pump_total_head_pa"]),
        "target_flow_kg_s": float(args.target_flow_kg_s),
        "tube_external_heat_matrix_key": str(args.tube_external_heat_matrix_key),
        "shield_external_heat_matrix_key": str(args.shield_external_heat_matrix_key),
        "external_heat_scale": float(args.external_heat_scale),
        "restart_in": None if args.restart_in is None else str(args.restart_in),
        "tube_external_heat_units": str(args.tube_external_heat_units),
        "initial_cs_fraction": float(args.initial_cs_fraction),
        "latest_record": latest_record,
    }
    path.write_text(json.dumps(latest, indent=2, default=json_default), encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.duration <= 0.0:
        raise ValueError("duration must be positive.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    history_path = output_dir / f"{args.case_prefix}_history.csv"
    latest_state_path = output_dir / f"{args.case_prefix}_latest_state.json"
    latest_restart_path = output_dir / f"{args.case_prefix}_latest_restart.npz"

    build = build_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)
    command = prepare_startup_step(build, args)
    last_record = collect_record(build, args, command, start_time)
    system.save_global_state(str(latest_restart_path))
    write_latest_state(latest_state_path, build, args, latest_restart_path, history_path, start_time, target_time, last_record)

    print("=== V13-start TOPAZ-II cold startup ===", flush=True)
    print(f"case_version={V13_START_CASE_VERSION}", flush=True)
    print(f"start_time={start_time:.6f}, target_time={target_time:.6f}", flush=True)
    print(f"initial_temperature={args.initial_temperature_k:.3f}K profile={args.startup_profile}", flush=True)
    print(f"history_csv={history_path}", flush=True)

    next_record_time = min(start_time + float(args.record_interval), target_time)
    next_restart_time = min(start_time + float(args.restart_interval), target_time)
    next_control_time = min(start_time + float(args.pump_control_interval), target_time)
    fieldnames: Optional[List[str]] = None

    while float(system.global_time) < target_time - 1.0e-10:
        command = prepare_startup_step(build, args)
        dt = system.compute_adaptive_dt(min_dt=1.0e-4, max_dt=float(args.max_dt), safety_factor=float(args.safety_factor))
        dt = min(
            float(dt),
            float(args.max_dt),
            next_record_time - float(system.global_time),
            next_control_time - float(system.global_time),
            target_time - float(system.global_time),
        )
        system.step(dt, inner_iter=int(args.inner_iter), fluid_max_iter=int(args.step_hydraulic_max_iter))

        if float(system.global_time) >= next_control_time - 1.0e-10:
            control = maybe_adjust_pump_head(build, args)
            if control is not None:
                print(
                    "pump_control "
                    f"flow={control['pump_control_measured_flow_kg_s']:.6f}kg/s "
                    f"head={control['pump_control_old_head_pa']:.3f}->{control['pump_control_new_head_pa']:.3f}Pa",
                    flush=True,
                )
            next_control_time = min(next_control_time + float(args.pump_control_interval), target_time)

        if float(system.global_time) >= next_record_time - 1.0e-10:
            command = prepare_startup_step(build, args)
            record = collect_record(build, args, command, start_time)
            flat = flatten_for_csv(record)
            if fieldnames is None:
                fieldnames = list(flat.keys())
            append_row(history_path, fieldnames, flat, write_header=not history_path.exists())
            maybe_switch_startup_tec_to_fixed_voltage(build, args, record)
            last_record = record
            print(
                f"t_rel={record['relative_time_s']:.1f}s "
                f"phase={record['startup_phase']} "
                f"Q={record['startup_thermal_power_w']:.1f}W "
                f"Tin={record['core_inlet_connector_t_k']:.3f}K "
                f"Temit={record['mean_emitter_temperature_k']:.3f}K "
                f"shield={record['radiation_shield_active']} "
                f"tec={record['tec_coupled_enabled']}",
                flush=True,
            )
            next_record_time = min(next_record_time + float(args.record_interval), target_time)

        if float(system.global_time) >= next_restart_time - 1.0e-10:
            system.save_global_state(str(latest_restart_path))
            write_latest_state(latest_state_path, build, args, latest_restart_path, history_path, start_time, target_time, last_record)
            next_restart_time = min(next_restart_time + float(args.restart_interval), target_time)

    system.save_global_state(str(latest_restart_path))
    write_latest_state(latest_state_path, build, args, latest_restart_path, history_path, start_time, target_time, last_record)
    print("V13-start cold startup run completed.", flush=True)


if __name__ == "__main__":
    main()







