import os
import sys
from typing import Any, Dict, Optional

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - plotting is optional
    plt = None

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from profiler import TEASAProfiler
from Solvers.SystemManager import SystemManager
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidVolume,
    NonUniformIncompressibleFluidChannel,
    MacroFlowJunction,
)
from Materials.Solids.UO2 import UO2
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.ZrH import ZirconiumHydride
from Materials.Solids.BerylliumOxide import BerylliumOxide
from Materials.Solids.GasGaps import Xenon, Cesium, CarbonDioxide, Helium
from Materials.Fluids.Sodium import Sodium
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


def get_time_dependent_dt_cap_v6(current_time: float, default_max_dt: float) -> float:
    if current_time < 1.0:
        return min(default_max_dt, 0.01)
    if current_time < 10.0:
        return min(default_max_dt, 0.05)
    if current_time < 100.0:
        return min(default_max_dt, 0.2)
    if current_time < 1000.0:
        return min(default_max_dt, 0.5)
    return min(default_max_dt, 1.0)


def get_case_reactivity_control(
    case_name: str,
    relative_time: float,
    rho_step: float,
    step_time_s: float,
) -> float:
    case_key = case_name.strip().upper()
    if case_key in {"A", "D"}:
        return 0.0
    if case_key in {"B", "C"}:
        return float(rho_step) if relative_time >= step_time_s else 0.0
    raise ValueError(f"Unsupported case name: {case_name}")


def build_v6_system(
    inlet_temperature_k: float,
    channel_inlet_flow_kg_s: float,
    enable_tec_coupled: bool,
    tec_update_interval_s: float,
    tec_target_voltage_v: float,
    tec_initial_current_a: float,
):
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
        "Sodium": Sodium(),
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
    k_outlet = 0.884

    a_flow = np.pi * (geom_params.r_coolant_outer**2 - geom_params.r_coolant_inner**2)
    d_h = 2.0 * (geom_params.r_coolant_outer - geom_params.r_coolant_inner)

    ring_names = ["Center", "Ring1", "Ring2", "Ring3"]
    multipliers = [1, 6, 12, 18]

    all_fluid_vols = []
    all_fluid_juncs = []
    fluid_channels = {}

    inlet_plenum = IncompressibleBoundaryVolume(
        name="Global_Inlet",
        material=sodium,
        P=p_inlet_sys,
        T=inlet_temperature_k,
    )
    outlet_plenum = IncompressibleBoundaryVolume(
        name="Global_Outlet",
        material=sodium,
        P=p_outlet_sys,
        T=inlet_temperature_k,
    )
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum.is_pressure_boundary = True
    all_fluid_vols.extend([inlet_plenum, outlet_plenum])

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

        l_buffer = 0.01
        inter_in = IncompressibleFluidVolume(
            name=f"InterIn_{name}",
            volume=a_flow * mult * l_buffer,
            length=l_buffer,
            flow_area=a_flow * mult,
            hydraulic_diam=d_h,
            material=sodium,
            initial_P=p_inlet_sys,
            initial_T=inlet_temperature_k,
        )
        inter_out = IncompressibleFluidVolume(
            name=f"InterOut_{name}",
            volume=a_flow * mult * l_buffer,
            length=l_buffer,
            flow_area=a_flow * mult,
            hydraulic_diam=d_h,
            material=sodium,
            initial_P=p_outlet_sys,
            initial_T=inlet_temperature_k,
        )
        all_fluid_vols.extend([inter_in, inter_out])

        j_in = InletJunction(
            name=f"J_In_{name}",
            from_vol=inlet_plenum,
            to_vol=inter_in,
            W_initial=w_single_design * mult,
        )
        j_macro_in = MacroFlowJunction(
            name=f"J_MacroIn_{name}",
            from_vol=inter_in,
            to_vol=chan.volumes[0],
            macro_vol=inter_in,
            multiplier=mult,
            flow_area=a_flow,
        )
        j_macro_out = MacroFlowJunction(
            name=f"J_MacroOut_{name}",
            from_vol=chan.volumes[-1],
            to_vol=inter_out,
            macro_vol=inter_out,
            multiplier=mult,
            flow_area=a_flow,
        )
        j_out = FlowJunction(
            name=f"J_Out_{name}",
            from_vol=inter_out,
            to_vol=outlet_plenum,
            flow_area=a_flow * mult,
            k_loss=k_outlet,
        )
        all_fluid_juncs.extend([j_in, j_macro_in, j_macro_out, j_out])

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

    # Keep the same core name as v5 so the v5 steady-state restart file can be
    # loaded without triggering shape-fingerprint mismatches on global solids.
    core = ReactorCore(
        name="TASTIN_Core_V5",
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
    }


def run_test_v6_case_a(
    run_duration_s: float = 200.0,
    total_power_w: float = 115000.0,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: float = 0.0351,
    restart_file: Optional[str] = "test_core_assemble_v5_restart_t5000.npz",
    save_interval: float = 0.0,
    enable_plot: bool = False,
    max_dt: float = 1.0,
    safety_factor: float = 20.0,
    enable_tec_coupled: bool = False,
    tec_update_interval_s: float = 1.0,
    tec_target_voltage_v: float = 27.2,
    tec_initial_current_a: float = 220.0,
    case_name: str = "A",
    rho_step: float = 0.0,
    step_time_s: float = 5.0,
) -> Dict[str, Any]:
    print(
        f"=== TASTIN System Test V6 Case {case_name.upper()}: "
        "Point-kinetics restart verification ==="
    )

    build = build_v6_system(
        inlet_temperature_k=inlet_temperature_k,
        channel_inlet_flow_kg_s=channel_inlet_flow_kg_s,
        enable_tec_coupled=enable_tec_coupled,
        tec_update_interval_s=tec_update_interval_s,
        tec_target_voltage_v=tec_target_voltage_v,
        tec_initial_current_a=tec_initial_current_a,
    )
    system = build["system"]
    core = build["core"]
    tfes = build["tfes"]
    ring_names = build["ring_names"]

    system.initialize_system()
    if restart_file:
        if not os.path.exists(restart_file):
            raise FileNotFoundError(f"Restart file not found: {restart_file}")
        system.load_global_state(restart_file)
        current_time = float(system.global_time)
    else:
        current_time = 0.0
        core.update_neutronic_power(p_total=total_power_w, alpha=1.0)

    simulation_start_time = current_time
    stop_time = simulation_start_time + float(run_duration_s)

    core.initialize_point_reactor(total_power_initial=total_power_w)

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
    history_power_fission = []
    history_power_decay = []
    history_reactivity_control = []
    history_reactivity_feedback_total = []
    history_reactivity_feedback_fuel = []
    history_reactivity_feedback_electrode = []
    history_reactivity_feedback_moderator = []
    history_reactivity_feedback_reflector = []
    history_effective_reactivity_feedback = []
    history_tec_current = []
    history_tec_voltage = []

    print(
        f"Running Case {case_name.upper()} from restart time {simulation_start_time:.3f} s "
        f"to {stop_time:.3f} s..."
    )
    while current_time < stop_time:
        relative_time = current_time - simulation_start_time
        rho_control = get_case_reactivity_control(
            case_name=case_name,
            relative_time=relative_time,
            rho_step=rho_step,
            step_time_s=step_time_s,
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

        system.step(dt, inner_iter=1, reactivity_control=rho_control)
        current_time = float(system.global_time)

        global_max_fuel = max(float(np.max(tfe.solids["pellet"].T)) for tfe in tfes.values())
        feedback = core.last_feedback_result
        tec_res = None
        if enable_tec_coupled and core.thermo_calc is not None:
            tec_res = core.thermo_calc.get_global_results()

        history_time.append(current_time)
        history_relative_time.append(current_time - simulation_start_time)
        history_max_fuel.append(global_max_fuel)
        history_core_power.append(float(core.last_total_core_power))
        history_power_fission.append(float(core.point_reactor.fission_power))
        history_power_decay.append(float(core.point_reactor.decay_power))
        history_reactivity_control.append(float(rho_control))
        history_reactivity_feedback_total.append(float(feedback.total))
        history_reactivity_feedback_fuel.append(float(feedback.fuel))
        history_reactivity_feedback_electrode.append(float(feedback.electrode))
        history_reactivity_feedback_moderator.append(float(feedback.moderator))
        history_reactivity_feedback_reflector.append(float(feedback.reflector))
        history_effective_reactivity_feedback.append(float(core.last_effective_reactivity_feedback))
        history_tec_current.append(0.0 if tec_res is None else float(tec_res.get("Iout", 0.0)))
        history_tec_voltage.append(0.0 if tec_res is None else float(tec_res.get("Uout", 0.0)))

        for name, tfe in tfes.items():
            history_outlet_temp[name].append(float(tfe.coolant.volumes[-1].T))
            history_flow[name].append(float(tfe.coolant.internal_junctions[0].W))

        if len(history_time) % 50 == 0:
            status_msg = (
                f"t_abs={current_time:9.3f} s | "
                f"t_rel={history_relative_time[-1]:8.3f} s | "
                f"rho_ctl={rho_control:+.6e} | "
                f"rho_fb_eff={history_effective_reactivity_feedback[-1]:+.6e} | "
                f"Ptot={history_core_power[-1] / 1000.0:8.3f} kW | "
                f"Tfuel_max={global_max_fuel:8.2f} K"
            )
            if enable_tec_coupled:
                status_msg += (
                    f" | TEC U={history_tec_voltage[-1]:6.2f} V"
                    f" | TEC I={history_tec_current[-1]:7.2f} A"
                )
            print(status_msg)

        if next_save_time is not None and current_time >= next_save_time:
            save_path = f"test_core_assemble_v6_case_{case_name.lower()}_restart_t{int(next_save_time)}.npz"
            print(f"[Checkpoint] Saving restart file: {save_path}")
            system.save_global_state(save_path)
            next_save_time += float(save_interval)

    power_series = np.asarray(history_core_power, dtype=float)
    power_deviation_w = power_series - float(total_power_w)
    final_summary = {
        "case_name": case_name.upper(),
        "restart_file": restart_file,
        "simulation_start_time_s": simulation_start_time,
        "final_time_s": current_time,
        "run_duration_s": current_time - simulation_start_time,
        "initial_power_w": float(total_power_w),
        "final_power_w": history_core_power[-1],
        "max_abs_power_deviation_w": float(np.max(np.abs(power_deviation_w))),
        "max_abs_power_deviation_pct": float(np.max(np.abs(power_deviation_w)) / total_power_w * 100.0),
        "final_reactivity_control": history_reactivity_control[-1],
        "final_effective_feedback": history_effective_reactivity_feedback[-1],
        "final_total_feedback": history_reactivity_feedback_total[-1],
        "final_max_fuel_k": history_max_fuel[-1],
        "final_center_outlet_k": history_outlet_temp["Center"][-1],
        "final_ring_outlet_k": {
            name: history_outlet_temp[name][-1] for name in ring_names
        },
        "final_ring_flow_kg_s": {
            name: history_flow[name][-1] for name in ring_names
        },
    }

    print("Case run completed.")
    print(final_summary)

    if enable_plot and plt is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        axes[0, 0].plot(history_relative_time, np.asarray(history_core_power) / 1000.0, color="tab:red")
        axes[0, 0].axhline(total_power_w / 1000.0, color="black", linestyle="--", linewidth=1.0)
        axes[0, 0].set_title("Core Total Power")
        axes[0, 0].set_xlabel("Relative Time [s]")
        axes[0, 0].set_ylabel("Power [kW]")
        axes[0, 0].grid(True)

        axes[0, 1].plot(history_relative_time, history_effective_reactivity_feedback, label="effective feedback")
        axes[0, 1].plot(history_relative_time, history_reactivity_control, label="rho control")
        axes[0, 1].set_title("Reactivity")
        axes[0, 1].set_xlabel("Relative Time [s]")
        axes[0, 1].set_ylabel("rho [-]")
        axes[0, 1].grid(True)
        axes[0, 1].legend()

        axes[1, 0].plot(history_relative_time, history_max_fuel, color="tab:orange")
        axes[1, 0].set_title("Global Max Fuel Temperature")
        axes[1, 0].set_xlabel("Relative Time [s]")
        axes[1, 0].set_ylabel("Temperature [K]")
        axes[1, 0].grid(True)

        for name in ring_names:
            axes[1, 1].plot(history_relative_time, history_outlet_temp[name], label=name)
        axes[1, 1].set_title("Coolant Outlet Temperature by Ring")
        axes[1, 1].set_xlabel("Relative Time [s]")
        axes[1, 1].set_ylabel("Temperature [K]")
        axes[1, 1].grid(True)
        axes[1, 1].legend()

        fig.tight_layout()
        plt.show()
    elif enable_plot:
        print("matplotlib is not available; skipping plots.")

    return {
        "history_time": history_time,
        "history_relative_time": history_relative_time,
        "history_max_fuel": history_max_fuel,
        "history_outlet_temp": history_outlet_temp,
        "history_flow": history_flow,
        "history_core_power": history_core_power,
        "history_power_fission": history_power_fission,
        "history_power_decay": history_power_decay,
        "history_reactivity_control": history_reactivity_control,
        "history_reactivity_feedback_total": history_reactivity_feedback_total,
        "history_reactivity_feedback_fuel": history_reactivity_feedback_fuel,
        "history_reactivity_feedback_electrode": history_reactivity_feedback_electrode,
        "history_reactivity_feedback_moderator": history_reactivity_feedback_moderator,
        "history_reactivity_feedback_reflector": history_reactivity_feedback_reflector,
        "history_effective_reactivity_feedback": history_effective_reactivity_feedback,
        "history_tec_current": history_tec_current,
        "history_tec_voltage": history_tec_voltage,
        "final_summary": final_summary,
        "core": core,
        "system": system,
        "tfes": tfes,
    }


if __name__ == "__main__":
    run_test_v6_case_a()
    TEASAProfiler.report()
