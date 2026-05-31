import os
import sys
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - plotting is optional
    plt = None

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from profiler import TEASAProfiler
from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.BoundaryVolume import (
    IncompressibleBoundaryVolume,
    InletJunction,
)
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidVolume,
    MacroFlowJunction,
    NonUniformIncompressibleFluidChannel,
)
from Materials.Solids.UO2 import UO2
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.ZrH import ZirconiumHydride
from Materials.Solids.BerylliumOxide import BerylliumOxide
from Materials.Solids.GasGaps import Xenon, Cesium, CarbonDioxide, Helium
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Components.TFEUnit import TFEUnit, TFEGeometry, TFEMeshParams, GapConfig
from Components.ReactorCore import (
    ReactorCore,
    GlobalAnnulusStructureConfig,
    GlobalGapStructureConfig,
)
from test_core_assemble_v5 import (
    build_axial_power_profile,
    build_global_moderator_meshes,
    build_ring_power_factors,
)
from test_core_assemble_v6 import get_time_dependent_dt_cap_v6


def _validate_power_table(
    times_s: Optional[Sequence[float]],
    values_w: Optional[Sequence[float]],
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    if times_s is None and values_w is None:
        return None, None
    if times_s is None or values_w is None:
        raise ValueError("power_table_times_s and power_table_values_w must be provided together.")

    times_arr = np.asarray(times_s, dtype=float)
    values_arr = np.asarray(values_w, dtype=float)

    if times_arr.ndim != 1 or values_arr.ndim != 1:
        raise ValueError("Power table times and values must be one-dimensional arrays.")
    if times_arr.size == 0:
        raise ValueError("Power table cannot be empty.")
    if times_arr.size != values_arr.size:
        raise ValueError("Power table times and values must have the same length.")
    if not np.all(np.isfinite(times_arr)) or not np.all(np.isfinite(values_arr)):
        raise ValueError("Power table times and values must be finite.")
    if np.any(np.diff(times_arr) <= 0.0):
        raise ValueError("Power table times must be strictly increasing.")
    if np.any(values_arr < 0.0):
        raise ValueError("Power table values must be non-negative.")

    return times_arr.copy(), values_arr.copy()


def _power_from_table(
    relative_time_s: float,
    total_power_w: float,
    table_times_s: Optional[np.ndarray],
    table_values_w: Optional[np.ndarray],
) -> float:
    if table_times_s is None:
        return float(total_power_w)
    return float(np.interp(float(relative_time_s), table_times_s, table_values_w))


def build_v7_system(
    inlet_temperature_k: float,
    channel_inlet_flow_kg_s: float,
    enable_tec_coupled: bool,
    tec_update_interval_s: float,
    tec_target_voltage_v: float,
    tec_initial_current_a: float,
    inlet_plenum_volume_m3: float = 5.0e-3,
    outlet_plenum_volume_m3: float = 5.0e-3,
    plenum_length_m: float = 0.05,
    outlet_boundary_k_loss: float = 0.884,
) -> Dict[str, Any]:
    if inlet_plenum_volume_m3 <= 0.0 or outlet_plenum_volume_m3 <= 0.0:
        raise ValueError("Plenum volumes must be positive.")
    if plenum_length_m <= 0.0:
        raise ValueError("plenum_length_m must be positive.")

    l_lower = 0.065
    l_active = 0.377
    l_upper = 0.065
    n_lower = 6
    n_active = 25
    n_upper = 6
    n_total = n_lower + n_active + n_upper
    total_height = l_lower + l_active + l_upper

    node_lengths = np.array(
        [l_lower / n_lower] * n_lower
        + [l_active / n_active] * n_active
        + [l_upper / n_upper] * n_upper,
        dtype=float,
    )
    common_y_faces = np.insert(np.cumsum(node_lengths), 0, 0.0)

    geom_params = TFEGeometry(
        r_pellet_inner=4.0e-3,
        r_pellet_outer=8.5e-3,
        r_fission_gas_outer=8.65e-3,
        r_emitter_outer=9.8e-3,
        r_collector_inner=10.3e-3,
        r_collector_outer=11.85e-3,
        r_inner_clad_inner=11.90e-3,
        r_inner_clad_outer=12.25e-3,
        r_coolant_inner=12.25e-3,
        r_coolant_outer=12.95e-3,
        r_outer_clad_outer=13.30e-3,
        r_moderator_inner=13.52e-3,
        r_moderator_outer=16.27e-3,
        height=total_height,
    )

    mesh_params = TFEMeshParams(
        n_axial=n_total,
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=1,
        n_r_outer_clad=1,
        n_r_moderator=3,
    )

    axial_power_profile = build_axial_power_profile(
        n_lower=n_lower,
        n_active=n_active,
        n_upper=n_upper,
    )

    materials_dict = {
        "UO2": UO2(),
        "MoNb": MoNb(),
        "Molybdenum": Molybdenum(),
        "StainlessSteel": AusteniticStainlessSteel(),
        "ZrH": ZirconiumHydride(),
        "BerylliumOxide": BerylliumOxide(),
        "Sodium": SodiumPotassium78(),
    }

    cfg_fg = GapConfig(
        mode="simplified",
        h_eq=5678.0,
        material=Xenon(),
        emissivity_inner=0.15,
        emissivity_outer=0.15,
    )
    cfg_tec = GapConfig(
        mode="simplified",
        h_eq=29.0,
        material=Cesium(),
        emissivity_inner=0.15,
        emissivity_outer=0.60,
    )
    cfg_he = GapConfig(
        mode="simplified",
        h_eq=5678.0,
        material=Helium(),
        emissivity_inner=0.60,
        emissivity_outer=0.80,
    )
    cfg_co2 = GapConfig(
        mode="simplified",
        h_eq=53.6,
        material=CarbonDioxide(),
        emissivity_inner=0.80,
        emissivity_outer=0.80,
    )

    sodium = materials_dict["Sodium"]
    p_inlet_sys = 165370.0
    p_outlet_sys = 161270.0
    w_single_design = float(channel_inlet_flow_kg_s)

    a_flow = np.pi * (geom_params.r_coolant_outer**2 - geom_params.r_coolant_inner**2)
    d_h = 2.0 * (geom_params.r_coolant_outer - geom_params.r_coolant_inner)

    ring_names = ["Center", "Ring1", "Ring2", "Ring3"]
    multipliers = [1, 6, 12, 18]
    total_multiplier = float(sum(multipliers))
    total_flow_design = w_single_design * total_multiplier
    total_flow_area = a_flow * total_multiplier

    all_fluid_vols = []
    all_fluid_juncs = []
    fluid_channels = {}

    inlet_boundary = IncompressibleBoundaryVolume(
        name="Global_Inlet_Boundary",
        material=sodium,
        P=p_inlet_sys,
        T=inlet_temperature_k,
        flow_area=max(total_flow_area, 1.0),
        hydraulic_diam=d_h,
    )
    outlet_boundary = IncompressibleBoundaryVolume(
        name="Global_Outlet_Boundary",
        material=sodium,
        P=p_outlet_sys,
        T=inlet_temperature_k,
        flow_area=max(total_flow_area, 1.0),
        hydraulic_diam=d_h,
    )
    inlet_boundary.is_pressure_boundary = True
    outlet_boundary.is_pressure_boundary = True

    inlet_plenum = IncompressibleFluidVolume(
        name="Global_Inlet_Plenum",
        volume=float(inlet_plenum_volume_m3),
        length=float(plenum_length_m),
        flow_area=max(float(inlet_plenum_volume_m3) / float(plenum_length_m), total_flow_area),
        hydraulic_diam=d_h,
        material=sodium,
        initial_P=p_inlet_sys,
        initial_T=inlet_temperature_k,
    )
    outlet_plenum = IncompressibleFluidVolume(
        name="Global_Outlet_Plenum",
        volume=float(outlet_plenum_volume_m3),
        length=float(plenum_length_m),
        flow_area=max(float(outlet_plenum_volume_m3) / float(plenum_length_m), total_flow_area),
        hydraulic_diam=d_h,
        material=sodium,
        initial_P=p_outlet_sys,
        initial_T=inlet_temperature_k,
    )

    all_fluid_vols.extend([inlet_boundary, inlet_plenum, outlet_plenum, outlet_boundary])

    j_inlet_boundary = InletJunction(
        name="J_Inlet_Boundary_to_Plenum",
        from_vol=inlet_boundary,
        to_vol=inlet_plenum,
        W_initial=total_flow_design,
    )
    all_fluid_juncs.append(j_inlet_boundary)

    for name, mult in zip(ring_names, multipliers):
        chan = NonUniformIncompressibleFluidChannel(
            name=f"Chan_{name}",
            node_lengths=node_lengths,
            flow_area=a_flow,
            hydraulic_diam=d_h,
            initial_P=p_inlet_sys,
            initial_T=inlet_temperature_k,
            material=sodium,
        )
        fluid_channels[name] = chan
        all_fluid_vols.extend(chan.volumes)
        all_fluid_juncs.extend(chan.internal_junctions)

        j_macro_in = MacroFlowJunction(
            name=f"J_PlenumIn_{name}",
            from_vol=inlet_plenum,
            to_vol=chan.volumes[0],
            macro_vol=inlet_plenum,
            multiplier=mult,
            flow_area=a_flow,
        )
        j_macro_out = MacroFlowJunction(
            name=f"J_PlenumOut_{name}",
            from_vol=chan.volumes[-1],
            to_vol=outlet_plenum,
            macro_vol=outlet_plenum,
            multiplier=mult,
            flow_area=a_flow,
        )
        j_macro_in.W = w_single_design
        j_macro_out.W = w_single_design
        all_fluid_juncs.extend([j_macro_in, j_macro_out])

    j_outlet_boundary = FlowJunction(
        name="J_Outlet_Plenum_to_Boundary",
        from_vol=outlet_plenum,
        to_vol=outlet_boundary,
        flow_area=total_flow_area,
        k_loss=float(outlet_boundary_k_loss),
    )
    j_outlet_boundary.W = total_flow_design
    all_fluid_juncs.append(j_outlet_boundary)

    tfes = {}
    for name in ring_names:
        tfes[name] = TFEUnit(
            name=name,
            geometry=geom_params,
            mesh_params=mesh_params,
            materials=materials_dict,
            coolant_channel=fluid_channels[name],
            fission_gas_config=cfg_fg,
            tec_gap_config=cfg_tec,
            he_gap_config=cfg_he,
            co2_gap_config=cfg_co2,
            power_fraction=1.0,
            axial_power_profile=axial_power_profile,
            axial_length_allocation=[l_lower, l_active, l_upper],
            axial_node_allocation=[n_lower, n_active, n_upper],
            axial_contact_resistance=0.0,
        )

    tfe_multipliers = {name: mult for name, mult in zip(ring_names, multipliers)}
    tfe_power_factors = build_ring_power_factors(ring_names, multipliers)
    ring_mapping = {name: i for i, name in enumerate(ring_names)}

    mod_meshes = build_global_moderator_meshes(
        inner_radius=geom_params.r_moderator_outer,
        outer_radius=60.0e-3,
        n_rings=len(ring_names),
        y_faces=common_y_faces,
        height=geom_params.height,
        n_axial=n_total,
    )

    moderator_barrel_gap_cfg = GlobalGapStructureConfig(
        mode="simplified",
        width=5.0e-3,
        h_eq=0.0,
        emissivity_inner=0.8,
        emissivity_outer=0.8,
    )
    barrel_cfg = GlobalAnnulusStructureConfig(
        material=materials_dict["StainlessSteel"],
        inner_radius=65.0e-3,
        thickness=3.0e-3,
        n_radial=3,
        initial_temp=inlet_temperature_k,
        outer_surface_emissivity=0.05,
    )
    barrel_reflector_gap_cfg = GlobalGapStructureConfig(
        mode="simplified",
        width=2.0e-3,
        h_eq=0.0,
        emissivity_inner=0.8,
        emissivity_outer=0.8,
    )
    reflector_cfg = GlobalAnnulusStructureConfig(
        material=materials_dict["BerylliumOxide"],
        outer_radius=102.0e-3,
        n_radial=8,
        initial_temp=inlet_temperature_k,
        outer_surface_emissivity=0.6,
    )

    core = ReactorCore(
        name="TASTIN_Core_V7",
        tfe_dict=tfes,
        tfe_multipliers=tfe_multipliers,
        tfe_power_factors=tfe_power_factors,
        mod_meshes=mod_meshes,
        mod_material=materials_dict["ZrH"],
        ring_mapping=ring_mapping,
        barrel_config=barrel_cfg,
        reflector_config=reflector_cfg,
        moderator_barrel_gap_config=moderator_barrel_gap_cfg,
        barrel_reflector_gap_config=barrel_reflector_gap_cfg,
        T_space=250.0,
        alpha_tec=0.5,
        enable_tec_coupled=enable_tec_coupled,
    )

    if enable_tec_coupled:
        if core.thermo_calc is None:
            raise RuntimeError(
                "TEC coupling is enabled, but ThermoCalc failed to initialize."
            )
        core.thermo_update_interval = float(tec_update_interval_s)
        core.setup_tec_circuit(
            mode_str="fixed_u",
            target_value=float(tec_target_voltage_v),
            I_guess=float(tec_initial_current_a),
        )

    hydraulic_net = HydraulicNetwork(
        volumes=all_fluid_vols,
        junctions=all_fluid_juncs,
        gravity_vector=0.0,
    )
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)

    return {
        "system": system,
        "core": core,
        "tfes": tfes,
        "ring_names": ring_names,
        "axial_power_profile": axial_power_profile,
        "fluid_channels": fluid_channels,
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "inlet_plenum": inlet_plenum,
        "outlet_plenum": outlet_plenum,
        "inlet_boundary_junction": j_inlet_boundary,
        "outlet_boundary_junction": j_outlet_boundary,
        "total_flow_design_kg_s": total_flow_design,
    }


def run_test_v7(
    run_duration_s: float = 20.0,
    total_power_w: float = 115000.0,
    power_table_times_s: Optional[Sequence[float]] = None,
    power_table_values_w: Optional[Sequence[float]] = None,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: float = 0.0351,
    restart_file: Optional[str] = None,
    save_interval: float = 0.0,
    enable_plot: bool = False,
    max_dt: float = 1.0,
    safety_factor: float = 20.0,
    enable_tec_coupled: bool = True,
    tec_update_interval_s: float = 1.0,
    tec_target_voltage_v: float = 27.2,
    tec_initial_current_a: float = 220.0,
    inlet_plenum_volume_m3: float = 5.0e-3,
    outlet_plenum_volume_m3: float = 5.0e-3,
    plenum_length_m: float = 0.05,
) -> Dict[str, Any]:
    print("=== TASTIN System Test V7: table-power core with inlet/outlet plena ===")

    table_times_s, table_values_w = _validate_power_table(
        power_table_times_s,
        power_table_values_w,
    )

    build = build_v7_system(
        inlet_temperature_k=inlet_temperature_k,
        channel_inlet_flow_kg_s=channel_inlet_flow_kg_s,
        enable_tec_coupled=enable_tec_coupled,
        tec_update_interval_s=tec_update_interval_s,
        tec_target_voltage_v=tec_target_voltage_v,
        tec_initial_current_a=tec_initial_current_a,
        inlet_plenum_volume_m3=inlet_plenum_volume_m3,
        outlet_plenum_volume_m3=outlet_plenum_volume_m3,
        plenum_length_m=plenum_length_m,
    )
    system = build["system"]
    core = build["core"]
    tfes = build["tfes"]
    ring_names = build["ring_names"]
    inlet_plenum = build["inlet_plenum"]
    outlet_plenum = build["outlet_plenum"]
    inlet_boundary_junction = build["inlet_boundary_junction"]
    outlet_boundary_junction = build["outlet_boundary_junction"]

    system.initialize_system()
    if restart_file:
        if not os.path.exists(restart_file):
            raise FileNotFoundError(f"Restart file not found: {restart_file}")
        try:
            system.load_global_state(restart_file)
        except ValueError as exc:
            raise ValueError(
                "Case7 restart_file must have the same V7 hydraulic topology. "
                "Old v5/v6 restart files are not compatible with the added plena."
            ) from exc
        current_time = float(system.global_time)
    else:
        current_time = 0.0

    # Case7 intentionally bypasses point kinetics.  This also protects against
    # accidental legacy neutronics state in a restart file.
    core.point_reactor = None
    current_power = _power_from_table(0.0, total_power_w, table_times_s, table_values_w)
    core.update_neutronic_power(
        p_total=current_power,
        p_fiss=current_power,
        p_decay=0.0,
        alpha=1.0,
    )

    simulation_start_time = current_time
    stop_time = simulation_start_time + float(run_duration_s)

    if enable_tec_coupled:
        core.enable_tec_coupled = True
        core.thermo_update_interval = float(tec_update_interval_s)
        core.setup_tec_circuit(
            mode_str="fixed_u",
            target_value=float(tec_target_voltage_v),
            I_guess=float(tec_initial_current_a),
        )
        core.post_step(0.0, current_time)
        core._last_thermo_update_time = current_time - core.thermo_update_interval
    else:
        core.enable_tec_coupled = False

    next_save_time = current_time + float(save_interval) if save_interval > 0.0 else None

    history_time = []
    history_relative_time = []
    history_max_fuel = []
    history_outlet_temp = {name: [] for name in ring_names}
    history_flow = {name: [] for name in ring_names}
    history_core_power = []
    history_power_table_target = []
    history_power_eval_time = []
    history_inlet_plenum_pressure = []
    history_outlet_plenum_pressure = []
    history_core_delta_p = []
    history_inlet_plenum_temp = []
    history_outlet_plenum_temp = []
    history_total_inlet_flow = []
    history_total_outlet_flow = []
    history_tec_current = []
    history_tec_voltage = []

    print(
        f"Running Case7 from t={simulation_start_time:.3f} s "
        f"to {stop_time:.3f} s..."
    )
    while current_time < stop_time:
        relative_time = current_time - simulation_start_time
        current_power = _power_from_table(
            relative_time,
            total_power_w,
            table_times_s,
            table_values_w,
        )
        core.update_neutronic_power(
            p_total=current_power,
            p_fiss=current_power,
            p_decay=0.0,
            alpha=1.0,
        )

        dt_cap = get_time_dependent_dt_cap_v6(relative_time, max_dt)
        if enable_tec_coupled:
            dt_cap = min(dt_cap, float(tec_update_interval_s))

        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=dt_cap,
            safety_factor=safety_factor,
        )
        dt = min(dt, stop_time - current_time)

        system.step(dt, inner_iter=1)
        current_time = float(system.global_time)

        global_max_fuel = max(float(np.max(tfe.solids["pellet"].T)) for tfe in tfes.values())
        tec_res = None
        if enable_tec_coupled and core.thermo_calc is not None:
            tec_res = core.thermo_calc.get_global_results()

        delta_p = float(inlet_plenum.P - outlet_plenum.P)
        history_time.append(current_time)
        history_relative_time.append(current_time - simulation_start_time)
        history_max_fuel.append(global_max_fuel)
        history_core_power.append(float(core.last_total_core_power))
        history_power_table_target.append(float(current_power))
        history_power_eval_time.append(float(relative_time))
        history_inlet_plenum_pressure.append(float(inlet_plenum.P))
        history_outlet_plenum_pressure.append(float(outlet_plenum.P))
        history_core_delta_p.append(delta_p)
        history_inlet_plenum_temp.append(float(inlet_plenum.T))
        history_outlet_plenum_temp.append(float(outlet_plenum.T))
        history_total_inlet_flow.append(float(inlet_boundary_junction.W))
        history_total_outlet_flow.append(float(outlet_boundary_junction.W))
        history_tec_current.append(0.0 if tec_res is None else float(tec_res.get("Iout", 0.0)))
        history_tec_voltage.append(0.0 if tec_res is None else float(tec_res.get("Uout", 0.0)))

        for name, tfe in tfes.items():
            history_outlet_temp[name].append(float(tfe.coolant.volumes[-1].T))
            history_flow[name].append(float(tfe.coolant.internal_junctions[0].W))

        if len(history_time) % 50 == 0:
            status_msg = (
                f"t_abs={current_time:9.3f} s | "
                f"t_rel={history_relative_time[-1]:8.3f} s | "
                f"P={history_core_power[-1] / 1000.0:8.3f} kW | "
                f"dP_core={delta_p:9.2f} Pa | "
                f"Tfuel_max={global_max_fuel:8.2f} K | "
                f"Tin_plenum={history_inlet_plenum_temp[-1]:7.2f} K | "
                f"Tout_plenum={history_outlet_plenum_temp[-1]:7.2f} K"
            )
            if enable_tec_coupled:
                status_msg += (
                    f" | TEC U={history_tec_voltage[-1]:6.2f} V"
                    f" | TEC I={history_tec_current[-1]:7.2f} A"
                )
            print(status_msg)

        if next_save_time is not None and current_time >= next_save_time:
            save_path = f"test_core_assemble_v7_restart_t{int(next_save_time)}.npz"
            print(f"[Checkpoint] Saving restart file: {save_path}")
            system.save_global_state(save_path)
            next_save_time += float(save_interval)

    if history_core_power:
        power_series = np.asarray(history_core_power, dtype=float)
        target_series = np.asarray(history_power_table_target, dtype=float)
        max_abs_power_error_w = float(np.max(np.abs(power_series - target_series)))
        final_table_power_w = _power_from_table(
            current_time - simulation_start_time,
            total_power_w,
            table_times_s,
            table_values_w,
        )
        final_summary = {
            "case_name": "V7",
            "restart_file": restart_file,
            "simulation_start_time_s": simulation_start_time,
            "final_time_s": current_time,
            "run_duration_s": current_time - simulation_start_time,
            "final_power_w": history_core_power[-1],
            "final_table_power_w_at_final_time": final_table_power_w,
            "max_abs_power_table_error_w": max_abs_power_error_w,
            "final_max_fuel_k": history_max_fuel[-1],
            "final_inlet_plenum_k": history_inlet_plenum_temp[-1],
            "final_outlet_plenum_k": history_outlet_plenum_temp[-1],
            "final_core_delta_p_pa": history_core_delta_p[-1],
            "final_inlet_plenum_pressure_pa": history_inlet_plenum_pressure[-1],
            "final_outlet_plenum_pressure_pa": history_outlet_plenum_pressure[-1],
            "final_total_inlet_flow_kg_s": history_total_inlet_flow[-1],
            "final_total_outlet_flow_kg_s": history_total_outlet_flow[-1],
            "final_ring_outlet_k": {
                name: history_outlet_temp[name][-1] for name in ring_names
            },
            "final_ring_flow_kg_s": {
                name: history_flow[name][-1] for name in ring_names
            },
        }
    else:
        final_summary = {
            "case_name": "V7",
            "restart_file": restart_file,
            "simulation_start_time_s": simulation_start_time,
            "final_time_s": current_time,
            "run_duration_s": 0.0,
            "final_power_w": float(core.last_total_core_power),
            "final_core_delta_p_pa": float(inlet_plenum.P - outlet_plenum.P),
        }

    print("Case7 run completed.")
    print(final_summary)

    if enable_plot and plt is not None and history_time:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        axes[0, 0].plot(history_relative_time, np.asarray(history_core_power) / 1000.0, color="tab:red")
        axes[0, 0].plot(
            history_relative_time,
            np.asarray(history_power_table_target) / 1000.0,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label="table target",
        )
        axes[0, 0].set_title("Core Total Power")
        axes[0, 0].set_xlabel("Relative Time [s]")
        axes[0, 0].set_ylabel("Power [kW]")
        axes[0, 0].grid(True)
        axes[0, 0].legend()

        axes[0, 1].plot(history_relative_time, history_core_delta_p, color="tab:blue")
        axes[0, 1].set_title("Plenum-to-Plenum Pressure Difference")
        axes[0, 1].set_xlabel("Relative Time [s]")
        axes[0, 1].set_ylabel("Delta P [Pa]")
        axes[0, 1].grid(True)

        axes[1, 0].plot(history_relative_time, history_max_fuel, color="tab:orange")
        axes[1, 0].set_title("Global Max Fuel Temperature")
        axes[1, 0].set_xlabel("Relative Time [s]")
        axes[1, 0].set_ylabel("Temperature [K]")
        axes[1, 0].grid(True)

        for name in ring_names:
            axes[1, 1].plot(history_relative_time, history_outlet_temp[name], label=name)
        axes[1, 1].plot(
            history_relative_time,
            history_outlet_plenum_temp,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label="Outlet Plenum",
        )
        axes[1, 1].set_title("Coolant Outlet Temperature")
        axes[1, 1].set_xlabel("Relative Time [s]")
        axes[1, 1].set_ylabel("Temperature [K]")
        axes[1, 1].grid(True)
        axes[1, 1].legend()

        fig.tight_layout()
        plt.show()
    elif enable_plot:
        print("matplotlib is not available or no history was recorded; skipping plots.")

    return {
        "history_time": history_time,
        "history_relative_time": history_relative_time,
        "history_max_fuel": history_max_fuel,
        "history_outlet_temp": history_outlet_temp,
        "history_flow": history_flow,
        "history_core_power": history_core_power,
        "history_power_table_target": history_power_table_target,
        "history_power_eval_time": history_power_eval_time,
        "history_inlet_plenum_pressure": history_inlet_plenum_pressure,
        "history_outlet_plenum_pressure": history_outlet_plenum_pressure,
        "history_core_delta_p": history_core_delta_p,
        "history_inlet_plenum_temp": history_inlet_plenum_temp,
        "history_outlet_plenum_temp": history_outlet_plenum_temp,
        "history_total_inlet_flow": history_total_inlet_flow,
        "history_total_outlet_flow": history_total_outlet_flow,
        "history_tec_current": history_tec_current,
        "history_tec_voltage": history_tec_voltage,
        "final_summary": final_summary,
        "core": core,
        "system": system,
        "tfes": tfes,
        "inlet_plenum": inlet_plenum,
        "outlet_plenum": outlet_plenum,
    }


if __name__ == "__main__":
    run_test_v7()
    TEASAProfiler.report()
