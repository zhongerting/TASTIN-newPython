import os
import sys
from typing import Dict, Any, Optional

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


COREINPUT_RING_SHARES_RAW = np.array(
    [0.030269999, 0.166800001, 0.335999999, 0.4644],
    dtype=float,
)

# Fitted from the normalized 25-point active-region power-shape table extracted
# from TASTIN-normal/CoreInput.txt (lines 593-596 after removing 6+6 reflector zeros).
# shape(z) = a0 + a2*z^2 + a4*z^4 + a6*z^6 + a8*z^8, z in [-1, 1]
AXIAL_SHAPE_COEFFS = np.array(
    [
        6.27905178e-02,
        -7.13913811e-02,
        1.35276842e-02,
        -1.02326367e-03,
        3.90936491e-05,
    ],
    dtype=float,
)


def coreinput_axial_shape(z: np.ndarray) -> np.ndarray:
    z2 = np.asarray(z, dtype=float) ** 2
    shape = (
        AXIAL_SHAPE_COEFFS[0]
        + AXIAL_SHAPE_COEFFS[1] * z2
        + AXIAL_SHAPE_COEFFS[2] * z2**2
        + AXIAL_SHAPE_COEFFS[3] * z2**3
        + AXIAL_SHAPE_COEFFS[4] * z2**4
    )
    return np.maximum(shape, 0.0)


def build_axial_power_profile(n_lower: int, n_active: int, n_upper: int) -> np.ndarray:
    z_centers = np.linspace(-1.0, 1.0, n_active)
    active_profile = coreinput_axial_shape(z_centers)
    active_profile /= np.sum(active_profile)
    return np.concatenate((np.zeros(n_lower), active_profile, np.zeros(n_upper)))


def build_ring_power_factors(ring_names, multipliers) -> Dict[str, float]:
    normalized_shares = COREINPUT_RING_SHARES_RAW / np.sum(COREINPUT_RING_SHARES_RAW)
    return {
        name: float(share / mult)
        for name, share, mult in zip(ring_names, normalized_shares, multipliers)
    }


def build_global_moderator_meshes(
    inner_radius: float,
    outer_radius: float,
    n_rings: int,
    y_faces: np.ndarray,
    height: float,
    n_axial: int,
) -> list[Mesh2D]:
    meshes = []
    radial_edges = np.linspace(inner_radius, outer_radius, n_rings + 1)
    for r_in, r_out in zip(radial_edges[:-1], radial_edges[1:]):
        meshes.append(
            Mesh2D(
                x_dim=r_out - r_in,
                n_x=3,
                y_dim=height,
                n_y=n_axial,
                y_faces=y_faces,
                geometry_type="cylindrical",
                inner_radius=r_in,
            )
        )
    return meshes


def get_time_dependent_dt_cap(current_time: float, default_max_dt: float) -> float:
    if current_time < 20.0:
        return min(default_max_dt, 0.2)
    if current_time < 200.0:
        return min(default_max_dt, 1.0)
    if current_time < 1000.0:
        return min(default_max_dt, 2.0)
    return min(default_max_dt, 5.0)


def run_test_v5(
    t_end: float = 10000.0,
    total_power_w: float = 115000.0,
    inlet_temperature_k: float = 743.0,
    channel_inlet_flow_kg_s: float = 0.0351,
    save_interval: float = 2500,
    # restart_file: Optional[str] = None,
    restart_file: Optional[str] = "test_core_assemble_v5_restart_t5000.npz",
    enable_plot: bool = False,
    max_dt: float = 1.0,
    safety_factor: float = 200.0,
    enable_tec_coupled: bool = True,
    tec_update_interval_s: float = 1.0,
    tec_target_voltage_v: float = 27.2,
    tec_initial_current_a: float = 220.0,
) -> Dict[str, Any]:
    print("=== TASTIN System Test V5: CoreInput-based thermal steady run ===")

    # ---------------------------------------------------------------------
    # 1. CoreInput-based geometry and mesh
    # ---------------------------------------------------------------------
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

    # ---------------------------------------------------------------------
    # 2. Hydraulic network
    # ---------------------------------------------------------------------
    sodium = materials_dict["Sodium"]
    p_inlet_sys = 165370.0
    p_outlet_sys = 161270.0
    # CoreInput.txt contains 0.003519 on the flow line, but in the current
    # Python hydraulic model that value drives the long transient unstable.
    # The v5 baseline therefore uses the stable design flow carried in the
    # earlier Python assembly scripts and keeps the CoreInput geometry/power map.
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

    # ---------------------------------------------------------------------
    # 3. TFE units and ReactorCore
    # ---------------------------------------------------------------------
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

    # ---------------------------------------------------------------------
    # 4. Assemble system
    # ---------------------------------------------------------------------
    hydraulic_net = HydraulicNetwork(
        volumes=all_fluid_vols,
        junctions=all_fluid_juncs,
        gravity_vector=0.0,
    )
    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)

    if restart_file and os.path.exists(restart_file):
        system.initialize_system()
        system.load_global_state(restart_file)
        current_time = system.global_time
    else:
        system.initialize_system()
        current_time = 0.0
        core.update_neutronic_power(p_total=total_power_w, alpha=1.0)

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

    next_save_time = current_time + float(save_interval) if save_interval > 0.0 else None

    history_time = []
    history_max_fuel = []
    history_outlet_temp = {name: [] for name in ring_names}
    history_flow = {name: [] for name in ring_names}
    history_mod_avg = {f"ModRing_{i}": [] for i in range(len(core.mod_rings))}
    history_core_power = []
    history_tec_current = []
    history_tec_voltage = []

    # ---------------------------------------------------------------------
    # 5. Transient to thermal steady state
    # ---------------------------------------------------------------------
    if enable_tec_coupled:
        print("Running thermo-hydraulic + TEC transient...")
    else:
        print("Running thermal-only transient...")
    while current_time < t_end:
        core.update_neutronic_power(p_total=total_power_w, alpha=1.0)

        dt_cap = get_time_dependent_dt_cap(current_time, max_dt)
        if enable_tec_coupled:
            dt_cap = min(dt_cap, float(tec_update_interval_s))
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-3,
            max_dt=dt_cap,
            safety_factor=safety_factor,
        )
        dt = min(dt, t_end - current_time)
        system.step(dt)
        current_time = system.global_time

        global_max_fuel = max(float(np.max(tfe.solids["pellet"].T)) for tfe in tfes.values())
        history_time.append(current_time)
        history_max_fuel.append(global_max_fuel)
        history_core_power.append(core.last_total_core_power)

        for name, tfe in tfes.items():
            history_outlet_temp[name].append(float(tfe.coolant.volumes[-1].T))
            history_flow[name].append(float(tfe.coolant.internal_junctions[0].W))

        for i, ring in enumerate(core.mod_rings):
            history_mod_avg[f"ModRing_{i}"].append(float(np.mean(ring.T)))

        tec_res = None
        if enable_tec_coupled and core.thermo_calc is not None:
            tec_res = core.thermo_calc.get_global_results()
        history_tec_current.append(0.0 if tec_res is None else float(tec_res.get("Iout", 0.0)))
        history_tec_voltage.append(0.0 if tec_res is None else float(tec_res.get("Uout", 0.0)))

        if len(history_time) % 50 == 0:
            status_msg = (
                f"Time: {current_time:8.2f} s | "
                f"Power: {core.last_total_core_power / 1000.0:7.2f} kW | "
                f"Max Fuel T: {global_max_fuel:8.2f} K | "
                f"Outlet[Center]: {history_outlet_temp['Center'][-1]:8.2f} K"
            )
            if enable_tec_coupled:
                status_msg += (
                    f" | TEC U: {history_tec_voltage[-1]:6.2f} V"
                    f" | TEC I: {history_tec_current[-1]:7.2f} A"
                )
            print(status_msg)

        if next_save_time is not None and current_time >= next_save_time:
            save_path = f"test_core_assemble_v5_restart_t{int(next_save_time)}.npz"
            print(f"[Checkpoint] Saving restart file: {save_path}")
            system.save_global_state(save_path)
            next_save_time += float(save_interval)

    # ---------------------------------------------------------------------
    # 6. Post-process
    # ---------------------------------------------------------------------
    power_summary = core.get_power_distribution_summary()
    final_summary = {
        "final_time_s": current_time,
        "final_max_fuel_k": history_max_fuel[-1],
        "final_center_outlet_k": history_outlet_temp["Center"][-1],
        "final_ring_outlet_k": {
            name: history_outlet_temp[name][-1] for name in ring_names
        },
        "final_ring_flow_kg_s": {
            name: history_flow[name][-1] for name in ring_names
        },
        "ring_power_factors": power_summary["tfe_power_factors"],
        "ring_group_power_shares": power_summary["tfe_group_power_shares"],
        "reconstructed_total_power_w": power_summary["reconstructed_total_real_power"],
    }

    if enable_tec_coupled:
        print("Thermo-hydraulic + TEC run completed.")
    else:
        print("Thermal-only run completed.")
    print(final_summary)

    if enable_plot and plt is not None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        axes[0, 0].plot(history_time, history_max_fuel, color="red")
        axes[0, 0].set_title("Global Max Fuel Temperature")
        axes[0, 0].set_xlabel("Time [s]")
        axes[0, 0].set_ylabel("Temperature [K]")
        axes[0, 0].grid(True)

        for name in ring_names:
            axes[0, 1].plot(history_time, history_outlet_temp[name], label=name)
        axes[0, 1].set_title("Coolant Outlet Temperature by Ring")
        axes[0, 1].set_xlabel("Time [s]")
        axes[0, 1].set_ylabel("Temperature [K]")
        axes[0, 1].grid(True)
        axes[0, 1].legend()

        for name in ring_names:
            axes[1, 0].plot(history_time, history_flow[name], label=name)
        axes[1, 0].set_title("Representative Channel Flow Rate")
        axes[1, 0].set_xlabel("Time [s]")
        axes[1, 0].set_ylabel("Mass Flow [kg/s]")
        axes[1, 0].grid(True)
        axes[1, 0].legend()

        for key, values in history_mod_avg.items():
            axes[1, 1].plot(history_time, values, label=key)
        axes[1, 1].set_title("Global Moderator Ring Average Temperature")
        axes[1, 1].set_xlabel("Time [s]")
        axes[1, 1].set_ylabel("Temperature [K]")
        axes[1, 1].grid(True)
        axes[1, 1].legend()

        fig.tight_layout()
        plt.show()
    elif enable_plot:
        print("matplotlib is not available; skipping plots.")

    return {
        "history_time": history_time,
        "history_max_fuel": history_max_fuel,
        "history_outlet_temp": history_outlet_temp,
        "history_flow": history_flow,
        "history_mod_avg": history_mod_avg,
        "history_core_power": history_core_power,
        "history_tec_current": history_tec_current,
        "history_tec_voltage": history_tec_voltage,
        "final_summary": final_summary,
        "core": core,
        "system": system,
        "tfes": tfes,
        "axial_power_profile": axial_power_profile,
    }


if __name__ == "__main__":
    run_test_v5()
    TEASAProfiler.report()
