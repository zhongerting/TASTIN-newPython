import logging
import os
import sys
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import (
    IncompressibleFluidVolume,
    IncompressibleFluidChannel,
    FlowJunction,
    MacroFlowJunction,
)
from Solvers.Hydrodynamics.BoundaryVolume import (
    IncompressibleBoundaryVolume,
    InletJunction,
)
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Solids.WallMaterial import SS316

import CoolantLoop.model_collector_ring_full_ringhp_geometry100hp_potassium_mixed as cfg


logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------
# Geometry and topology
# ---------------------------------------------------------
T_SPACE = cfg.T_SPACE
T_INLET = cfg.T_INLET
P_OUTLET = cfg.P_OUTLET
T_INIT = cfg.T_INIT
W_TOTAL = cfg.W_TOTAL

L_SECTOR = cfg.L_SECTOR
N_SECTOR = 3
AREA_RING = cfg.AREA_RING
DH_RING = cfg.DH_RING
PERIM_HEADER = cfg.PERIM_HEADER
R_IN_RING = cfg.R_IN_RING
R_OUT_RING = cfg.R_OUT_RING

AREA_HOT_LEG = cfg.AREA_HOT_LEG
DH_HOT_LEG = cfg.DH_HOT_LEG
L_HOT_LEG = cfg.L_HOT_LEG
N_HOT_LEG = cfg.N_HOT_LEG

AREA_MANIFOLD = cfg.AREA_MANIFOLD
DH_MANIFOLD = cfg.DH_MANIFOLD
L_MANIFOLD = cfg.L_MANIFOLD
N_MANIFOLD = cfg.N_MANIFOLD

L_INLET_BUFFER = cfg.L_INLET_BUFFER
N_INLET_BUFFER = cfg.N_INLET_BUFFER
AREA_INLET_BUFFER = cfg.AREA_INLET_BUFFER
DH_INLET_BUFFER = cfg.DH_INLET_BUFFER

L_OUTLET_BUFFER = cfg.L_OUTLET_BUFFER
N_OUTLET_BUFFER = cfg.N_OUTLET_BUFFER
AREA_OUTLET_BUFFER = cfg.AREA_OUTLET_BUFFER
DH_OUTLET_BUFFER = cfg.DH_OUTLET_BUFFER

RING_EMISSIVITY = cfg.RING_EMISSIVITY
HP_EMISSIVITY = cfg.HP_EMISSIVITY
FIN_EMISSIVITY = cfg.FIN_EMISSIVITY
HP_TOTAL_COUNT = 100

SEGMENT_SPECS = [
    ("A1_I1_to_O1", "I1", "O1", [5, 6, 6]),
    ("A2_O1_to_I2", "O1", "I2", [5, 5, 6]),
    ("A3_I2_to_O2", "I2", "O2", [5, 6, 6]),
    ("A4_O2_to_I3", "O2", "I3", [5, 5, 6]),
    ("A5_I3_to_O3", "I3", "O3", [5, 6, 6]),
    ("A6_O3_to_I1", "O3", "I1", [5, 6, 6]),
]
INLET_MIX_KEYS = ["I1", "I2", "I3"]
OUTLET_MIX_KEYS = ["O1", "O2", "O3"]

DEFAULT_T_END = 50.0
DEFAULT_PRINT_EVERY_TIME = 1.0
DEFAULT_RESTART_SAVE_EVERY = 10.0
RING_WALL_ODE_METHOD = "LSODA"
HP_WICK_ANISOTROPIC = True
HP_WICK_CONDUCTIVITY_CAP = None
HP_TIME_INTEGRATOR = "theta_implicit"
HP_THETA_IMPLICIT_VALUE = 0.6
HP_IMPLICIT_BOUNDARY_LINEARIZATION = True

if sum(sum(spec[3]) for spec in SEGMENT_SPECS) != HP_TOTAL_COUNT:
    raise ValueError("Segment heat-pipe multipliers must sum to HP_TOTAL_COUNT.")


nak = cfg.nak


def build_sector_solid(name):
    mesh = Mesh2D(
        x_dim=R_OUT_RING - R_IN_RING,
        n_x=1,
        y_dim=L_SECTOR,
        n_y=N_SECTOR,
        geometry_type="cylindrical",
        inner_radius=R_IN_RING,
    )
    solid = HeatConduction2D(
        mesh=mesh,
        material=SS316(),
        name=name,
        initial_temp=T_INIT,
    )
    solid.set_ode_method(RING_WALL_ODE_METHOD)
    bare_area_array = solid.boundaries["right"].area / 2.0
    solid.boundaries["right"].add_dynamic_radiation_condition(
        emissivity=RING_EMISSIVITY,
        bare_area_array=bare_area_array,
        T_env=T_SPACE,
    )
    return solid


def build_mix_node(name, kind):
    if kind == "inlet":
        area = AREA_HOT_LEG
        dh = DH_HOT_LEG
    elif kind == "outlet":
        area = AREA_MANIFOLD
        dh = DH_MANIFOLD
    else:
        raise ValueError(f"Unknown mix node kind: {kind}")

    length = 2.0 * dh
    return IncompressibleFluidVolume(
        name=name,
        volume=area * length,
        length=length,
        flow_area=area,
        hydraulic_diam=dh,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=nak,
    )


def configure_ring_hp_heat_pipe_solver(ring_hp):
    for hp_unit in ring_hp.hp_units:
        hp = hp_unit.hp
        hp.set_wick_conductivity_mode(HP_WICK_ANISOTROPIC)
        hp.wick_mat.set_conductivity_cap(HP_WICK_CONDUCTIVITY_CAP)
        hp.set_time_integrator(HP_TIME_INTEGRATOR)
        if HP_TIME_INTEGRATOR == "theta_implicit":
            hp.set_theta_implicit_value(HP_THETA_IMPLICIT_VALUE)
        hp.set_implicit_boundary_linearization(HP_IMPLICIT_BOUNDARY_LINEARIZATION)


def build_model():
    inlet_boundary = IncompressibleBoundaryVolume(
        name="InletBoundary",
        material=nak,
        P=P_OUTLET + 5000.0,
        T=T_INLET,
    )
    outlet_boundary = IncompressibleBoundaryVolume(
        name="OutletBoundary",
        material=nak,
        P=P_OUTLET,
        T=T_INIT,
    )
    outlet_boundary.is_pressure_boundary = True

    inlet_buffer_channel = IncompressibleFluidChannel(
        name="InletBuffer",
        n_nodes=N_INLET_BUFFER,
        total_length=L_INLET_BUFFER,
        flow_area=AREA_INLET_BUFFER,
        hydraulic_diam=DH_INLET_BUFFER,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=nak,
    )
    outlet_buffer_channel = IncompressibleFluidChannel(
        name="OutletBuffer",
        n_nodes=N_OUTLET_BUFFER,
        total_length=L_OUTLET_BUFFER,
        flow_area=AREA_OUTLET_BUFFER,
        hydraulic_diam=DH_OUTLET_BUFFER,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=nak,
    )

    hot_legs = [
        IncompressibleFluidChannel(
            name=f"HotLeg_{idx}",
            n_nodes=N_HOT_LEG,
            total_length=L_HOT_LEG,
            flow_area=AREA_HOT_LEG,
            hydraulic_diam=DH_HOT_LEG,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        for idx in range(1, 4)
    ]
    manifolds = [
        IncompressibleFluidChannel(
            name=f"Manifold_{idx}",
            n_nodes=N_MANIFOLD,
            total_length=L_MANIFOLD,
            flow_area=AREA_MANIFOLD,
            hydraulic_diam=DH_MANIFOLD,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        for idx in range(1, 4)
    ]

    inlet_mix_nodes = {
        key: build_mix_node(f"InletMix_{key}", "inlet")
        for key in INLET_MIX_KEYS
    }
    outlet_mix_nodes = {
        key: build_mix_node(f"OutletMix_{key}", "outlet")
        for key in OUTLET_MIX_KEYS
    }
    mix_nodes = {**inlet_mix_nodes, **outlet_mix_nodes}

    sectors = []
    solids = []
    ring_hps = []
    segment_links = []
    segment_entry_links = []
    segment_exit_links = []

    for sector_name, start_key, end_key, multipliers in SEGMENT_SPECS:
        channel = IncompressibleFluidChannel(
            name=f"{sector_name}_Channel",
            n_nodes=N_SECTOR,
            total_length=L_SECTOR,
            flow_area=AREA_RING,
            hydraulic_diam=DH_RING,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        solid = build_sector_solid(f"{sector_name}_Solid")
        ring_hp = cfg.build_ring_hp(
            name=f"{sector_name}_RingHP",
            fluid_channel=channel,
            solid_header=solid,
            hp_multipliers=multipliers,
        )
        configure_ring_hp_heat_pipe_solver(ring_hp)
        sectors.append(channel)
        solids.append(solid)
        ring_hps.append(ring_hp)

        entry_link = FlowJunction(
            name=f"J_{start_key}_to_{sector_name}",
            from_vol=mix_nodes[start_key],
            to_vol=channel.volumes[0],
            flow_area=AREA_RING,
            k_loss=0.0,
        )
        exit_link = FlowJunction(
            name=f"J_{sector_name}_to_{end_key}",
            from_vol=channel.volumes[-1],
            to_vol=mix_nodes[end_key],
            flow_area=AREA_RING,
            k_loss=ring_hp.outlet_k_loss,
            dynamic_loss_params=ring_hp.outlet_dynamic_loss_params,
        )
        segment_entry_links.append(entry_link)
        segment_exit_links.append(exit_link)
        segment_links.extend([entry_link, exit_link])

    inlet_junction = InletJunction(
        name="J_InletBoundary_InletBuffer",
        from_vol=inlet_boundary,
        to_vol=inlet_buffer_channel.volumes[0],
        W_initial=W_TOTAL,
    )
    outlet_junction = FlowJunction(
        name="J_OutletBuffer_OutletBoundary",
        from_vol=outlet_buffer_channel.volumes[-1],
        to_vol=outlet_boundary,
        flow_area=AREA_OUTLET_BUFFER,
        k_loss=0.0,
    )

    inlet_buffer_to_hot_leg = []
    hot_leg_to_inlet_mix = []
    outlet_mix_to_manifold = []
    manifold_to_outlet_buffer = []

    for idx, key in enumerate(INLET_MIX_KEYS):
        inlet_buffer_to_hot_leg.append(
            MacroFlowJunction(
                name=f"J_InletBuffer_HotLeg_{idx + 1}",
                from_vol=inlet_buffer_channel.volumes[-1],
                to_vol=hot_legs[idx].volumes[0],
                macro_vol=inlet_buffer_channel.volumes[-1],
                multiplier=2,
                flow_area=AREA_HOT_LEG,
                k_loss=0.0,
            )
        )
        hot_leg_to_inlet_mix.append(
            FlowJunction(
                name=f"J_HotLeg_{idx + 1}_InletMix_{key}",
                from_vol=hot_legs[idx].volumes[-1],
                to_vol=inlet_mix_nodes[key],
                flow_area=AREA_HOT_LEG,
                k_loss=0.0,
            )
        )

    for idx, key in enumerate(OUTLET_MIX_KEYS):
        outlet_mix_to_manifold.append(
            FlowJunction(
                name=f"J_OutletMix_{key}_Manifold_{idx + 1}",
                from_vol=outlet_mix_nodes[key],
                to_vol=manifolds[idx].volumes[0],
                flow_area=AREA_MANIFOLD,
                k_loss=0.0,
            )
        )
        manifold_to_outlet_buffer.append(
            MacroFlowJunction(
                name=f"J_Manifold_{idx + 1}_OutletBuffer",
                from_vol=manifolds[idx].volumes[-1],
                to_vol=outlet_buffer_channel.volumes[0],
                macro_vol=outlet_buffer_channel.volumes[0],
                multiplier=2,
                flow_area=AREA_MANIFOLD,
                k_loss=0.0,
            )
        )

    all_vols = [inlet_boundary, outlet_boundary]
    all_vols.extend(inlet_buffer_channel.volumes)
    for channel in hot_legs:
        all_vols.extend(channel.volumes)
    all_vols.extend(inlet_mix_nodes.values())
    all_vols.extend(outlet_mix_nodes.values())
    for channel in sectors:
        all_vols.extend(channel.volumes)
    for channel in manifolds:
        all_vols.extend(channel.volumes)
    all_vols.extend(outlet_buffer_channel.volumes)

    all_juncs = [inlet_junction]
    all_juncs.extend(inlet_buffer_to_hot_leg)
    all_juncs.extend(hot_leg_to_inlet_mix)
    all_juncs.extend(segment_links)
    all_juncs.extend(outlet_mix_to_manifold)
    all_juncs.extend(manifold_to_outlet_buffer)
    all_juncs.append(outlet_junction)
    all_juncs.extend(inlet_buffer_channel.internal_junctions)
    for channel in hot_legs:
        all_juncs.extend(channel.internal_junctions)
    for channel in sectors:
        all_juncs.extend(channel.internal_junctions)
    for channel in manifolds:
        all_juncs.extend(channel.internal_junctions)
    all_juncs.extend(outlet_buffer_channel.internal_junctions)

    network = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    sys_mgr = SystemManager(fluid_network=network)
    for ring_hp in ring_hps:
        sys_mgr.add_component(ring_hp)

    return {
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "inlet_buffer_channel": inlet_buffer_channel,
        "outlet_buffer_channel": outlet_buffer_channel,
        "hot_legs": hot_legs,
        "manifolds": manifolds,
        "inlet_mix_nodes": inlet_mix_nodes,
        "outlet_mix_nodes": outlet_mix_nodes,
        "sectors": sectors,
        "solids": solids,
        "ring_hps": ring_hps,
        "inlet_junction": inlet_junction,
        "outlet_junction": outlet_junction,
        "inlet_buffer_to_hot_leg": inlet_buffer_to_hot_leg,
        "hot_leg_to_inlet_mix": hot_leg_to_inlet_mix,
        "segment_links": segment_links,
        "segment_entry_links": segment_entry_links,
        "segment_exit_links": segment_exit_links,
        "outlet_mix_to_manifold": outlet_mix_to_manifold,
        "manifold_to_outlet_buffer": manifold_to_outlet_buffer,
        "all_vols": all_vols,
        "all_juncs": all_juncs,
        "network": network,
        "sys_mgr": sys_mgr,
    }


def get_model_statistics(model):
    solid_nodes = int(sum(np.size(solid.T) for solid in model["solids"]))
    fluid_nodes = int(len(model["all_vols"]) - 2)
    flow_junctions = int(len(model["all_juncs"]))
    return solid_nodes, fluid_nodes, flow_junctions


def print_model_summary(model):
    print("Model: 6 RingHP segments + independent inlet/outlet mixing volumes")
    print(f"Volumes: {len(model['all_vols'])}")
    print(f"Junctions: {len(model['all_juncs'])}")
    print(f"RingHP segment nodes : {N_SECTOR} per segment, {6 * N_SECTOR} total")
    print(f"HP segment totals    : {[sum(spec[3]) for spec in SEGMENT_SPECS]}")
    print(f"HP multipliers       : {[spec[3] for spec in SEGMENT_SPECS]}")
    print(f"HPs in whole ring    : {sum(sum(spec[3]) for spec in SEGMENT_SPECS)}")
    print(f"HP/fin emissivity    : {HP_EMISSIVITY:.3f} / {FIN_EMISSIVITY:.3f}")
    print(f"Ring emissivity      : {RING_EMISSIVITY:.3f}, T_space = {T_SPACE:.1f} K")
    print(f"Ring wall ODE method : {RING_WALL_ODE_METHOD}")
    print(
        "HP solver            : "
        f"wick_anisotropic={HP_WICK_ANISOTROPIC}, "
        f"k_cap={HP_WICK_CONDUCTIVITY_CAP}, "
        f"integrator={HP_TIME_INTEGRATOR}, "
        f"theta={HP_THETA_IMPLICIT_VALUE:.3f}, "
        f"implicit_boundary={HP_IMPLICIT_BOUNDARY_LINEARIZATION}"
    )
    print("  Ring topology: I1 -> A1 -> O1 -> A2 -> I2 -> A3 -> O2")
    print("                 -> A4 -> I3 -> A5 -> O3 -> A6 -> I1")
    print("  MacroFlow is only between macro buffers and single-ring branches.")


def print_pre_run_summary(model, case_name):
    solid_nodes, fluid_nodes, flow_junctions = get_model_statistics(model)
    print("=" * 78)
    print(f"Case: {case_name}")
    print("=" * 78)
    print_model_summary(model)
    print("-" * 78)
    print("Pre-run Summary")
    print(f"  Solid nodes    : {solid_nodes}")
    print(f"  Fluid nodes    : {fluid_nodes}")
    print(f"  Flow junctions : {flow_junctions}")
    print("=" * 78)


def run_case(
    case_name="collector_ring_6segment_geometry100hp_potassium_mixed_50s",
    t_end=DEFAULT_T_END,
    min_dt=1.0e-3,
    max_dt=0.5,
    safety_factor=1.0,
    inner_iter=2,
    print_every_time=DEFAULT_PRINT_EVERY_TIME,
    csv_path=None,
    restart_from=None,
    restart_save_path=None,
    restart_save_every=DEFAULT_RESTART_SAVE_EVERY,
):
    model = build_model()
    sys_mgr = model["sys_mgr"]
    network = model["network"]
    inlet_boundary = model["inlet_boundary"]
    outlet_boundary = model["outlet_boundary"]
    inlet_junction = model["inlet_junction"]
    outlet_junction = model["outlet_junction"]
    inlet_buffer_channel = model["inlet_buffer_channel"]
    outlet_buffer_channel = model["outlet_buffer_channel"]
    hot_legs = model["hot_legs"]
    manifolds = model["manifolds"]
    inlet_mix_nodes = model["inlet_mix_nodes"]
    outlet_mix_nodes = model["outlet_mix_nodes"]
    hot_leg_to_inlet_mix = model["hot_leg_to_inlet_mix"]
    outlet_mix_to_manifold = model["outlet_mix_to_manifold"]
    inlet_buffer_to_hot_leg = model["inlet_buffer_to_hot_leg"]
    manifold_to_outlet_buffer = model["manifold_to_outlet_buffer"]
    sectors = model["sectors"]
    segment_entry_links = model["segment_entry_links"]
    segment_exit_links = model["segment_exit_links"]

    print_pre_run_summary(model, case_name)

    if restart_from is not None:
        sys_mgr.load_global_state(restart_from)
        inlet_boundary.set_boundary_state(P=P_OUTLET + 5000.0, T=T_INLET)
        cfg.sync_boundary_to_network(network, inlet_boundary)
        outlet_boundary.set_boundary_state(P=P_OUTLET)
        cfg.sync_boundary_to_network(network, outlet_boundary)
        print(f"Restart loaded from: {restart_from}")
        print(f"Restart time: {sys_mgr.global_time:.6f} s")
    else:
        sys_mgr.initialize_system()
        print("System initialized from initial condition.")

    history = []
    next_print_time = cfg.next_event_time(sys_mgr.global_time, print_every_time)
    next_restart_save_time = cfg.next_event_time(sys_mgr.global_time, restart_save_every)

    while sys_mgr.global_time < t_end:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=min_dt,
            max_dt=max_dt,
            safety_factor=safety_factor,
        )
        if next_restart_save_time is not None:
            if sys_mgr.global_time < next_restart_save_time < sys_mgr.global_time + dt:
                dt = next_restart_save_time - sys_mgr.global_time
        if next_print_time is not None:
            if sys_mgr.global_time < next_print_time < sys_mgr.global_time + dt:
                dt = min(dt, next_print_time - sys_mgr.global_time)
        dt = min(dt, t_end - sys_mgr.global_time)

        sys_mgr.step(dt=dt, inner_iter=inner_iter)
        current_t = sys_mgr.global_time

        w_ring_in_total = float(sum(j.W for j in hot_leg_to_inlet_mix))
        w_ring_out_total = float(sum(j.W for j in outlet_mix_to_manifold))
        t_out_avg = float(np.mean([channel.volumes[-1].T for channel in manifolds]))

        row = {
            "time": current_t,
            "dt": dt,
            "W_in_total": float(inlet_junction.W),
            "W_out_total": float(outlet_junction.W),
            "W_ring_in_total": w_ring_in_total,
            "W_ring_out_total": w_ring_out_total,
            "T_out_avg": t_out_avg,
            "T_inlet_buffer_out": float(inlet_buffer_channel.volumes[-1].T),
            "T_outlet_buffer_out": float(outlet_buffer_channel.volumes[-1].T),
        }
        for idx, key in enumerate(INLET_MIX_KEYS, start=1):
            node = inlet_mix_nodes[key]
            row[f"T_inlet_mix_{key}"] = float(node.T)
            row[f"P_inlet_mix_{key}"] = float(node.P)
            row[f"W_macro_inlet_to_hotleg_{idx}"] = float(
                inlet_buffer_to_hot_leg[idx - 1].get_mass_flow_for(inlet_buffer_channel.volumes[-1])
            )
            row[f"W_hotleg_to_inlet_mix_{key}"] = float(hot_leg_to_inlet_mix[idx - 1].W)
        for idx, key in enumerate(OUTLET_MIX_KEYS, start=1):
            node = outlet_mix_nodes[key]
            row[f"T_outlet_mix_{key}"] = float(node.T)
            row[f"P_outlet_mix_{key}"] = float(node.P)
            row[f"W_outlet_mix_to_manifold_{key}"] = float(outlet_mix_to_manifold[idx - 1].W)
            row[f"W_macro_manifold_to_outlet_{idx}"] = float(
                manifold_to_outlet_buffer[idx - 1].get_mass_flow_for(outlet_buffer_channel.volumes[0])
            )
        for idx, channel in enumerate(hot_legs, start=1):
            row[f"T_hotleg_{idx}_out"] = float(channel.volumes[-1].T)
        for idx, channel in enumerate(manifolds, start=1):
            row[f"T_manifold_{idx}_out"] = float(channel.volumes[-1].T)
        for idx, channel in enumerate(sectors, start=1):
            row[f"T_A{idx}_in"] = float(channel.volumes[0].T)
            row[f"T_A{idx}_out"] = float(channel.volumes[-1].T)
            row[f"W_A{idx}_entry"] = float(segment_entry_links[idx - 1].W)
            row[f"W_A{idx}_exit"] = float(segment_exit_links[idx - 1].W)

        history.append(row)

        should_print = current_t >= t_end
        if next_print_time is not None and current_t >= next_print_time - 1.0e-12:
            should_print = True
        if should_print:
            print(
                f"t = {current_t:8.3f} s | "
                f"T_out_avg = {t_out_avg:.3f} K | "
                f"W_in_total = {inlet_junction.W:.4f} kg/s | "
                f"W_ring_in_total = {w_ring_in_total:.4f} kg/s"
            )
            while next_print_time is not None and current_t >= next_print_time - 1.0e-12:
                next_print_time += print_every_time

        if (
            restart_save_path is not None
            and restart_save_every > 0.0
            and next_restart_save_time is not None
            and current_t >= next_restart_save_time - 1.0e-12
        ):
            checkpoint_path = cfg.restart_checkpoint_path(restart_save_path, next_restart_save_time)
            sys_mgr.save_global_state(checkpoint_path)
            print(f"Restart saved at t={current_t:.3f} s: {checkpoint_path}")
            while next_restart_save_time is not None and current_t >= next_restart_save_time - 1.0e-12:
                next_restart_save_time += restart_save_every

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}_history.csv")
    cfg.write_history_csv(csv_path, history)

    if restart_save_path is not None:
        sys_mgr.save_global_state(restart_save_path)
        print(f"Final restart saved: {restart_save_path}")

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print("=" * 70)
    return model, history


if __name__ == "__main__":
    run_case(
        case_name="collector_ring_6segment_geometry100hp_potassium_mixed_50s",
        t_end=DEFAULT_T_END,
        print_every_time=DEFAULT_PRINT_EVERY_TIME,
        csv_path=os.path.join(current_dir, "collector_ring_6segment_geometry100hp_potassium_mixed_50s_history.csv"),
        restart_save_path=os.path.join(current_dir, "collector_ring_6segment_geometry100hp_potassium_mixed_50s_restart.npz"),
        restart_save_every=DEFAULT_RESTART_SAVE_EVERY,
        restart_from=None,
    )
