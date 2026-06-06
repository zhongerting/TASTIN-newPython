import csv
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
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.KHP import PotassiumHP
from Materials.Solids.WickMaterial import WickMaterial
from Components.RingHP import RingHP


logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------
# 0. Geometry and operating parameters
# ---------------------------------------------------------
T_SPACE = 200.0
T_INLET = 843.0
P_OUTLET = 160000.0
T_INIT = T_INLET
HP_INITIAL_TEMP = 800.0

W_TOTAL = 1.3
W_BRANCH_TOTAL = W_TOTAL / 3.0
W_INLET_LEG_INIT = W_BRANCH_TOTAL

L_INLET_BUFFER = 0.20
N_INLET_BUFFER = 5

L_SECTOR = 0.793
N_SECTOR = 4
L_RING = 6.0 * L_SECTOR
N_RING = 6 * N_SECTOR
RECT_LENGTH = 0.110
RECT_WIDTH = 0.040
AREA_RING = RECT_LENGTH * RECT_WIDTH
PERIM_HEADER = 2.0 * (RECT_LENGTH + RECT_WIDTH)
DH_RING = 4.0 * AREA_RING / PERIM_HEADER
HEADER_WALL_THICKNESS = 0.002
R_IN_RING = PERIM_HEADER / (2.0 * np.pi)
R_OUT_RING = R_IN_RING + HEADER_WALL_THICKNESS
HP_TOTAL_COUNT = 100
HP_COUNT_PER_SECTOR = HP_TOTAL_COUNT / 6.0
HP_MULTIPLIERS_RING = [5 if idx in (0, 6, 12, 18) else 4 for idx in range(N_RING)]
HP_MULTIPLIERS_SECTOR = HP_MULTIPLIERS_RING[:N_SECTOR]

L_HOT_LEG = 2.19632
R_IN_HOT_LEG = 0.0138
DH_HOT_LEG = 2.0 * R_IN_HOT_LEG
AREA_HOT_LEG = np.pi * R_IN_HOT_LEG**2
N_HOT_LEG = 28

AREA_INLET_BUFFER = 3.0 * AREA_HOT_LEG
DH_INLET_BUFFER = 2.0 * np.sqrt(AREA_INLET_BUFFER / np.pi)

L_MANIFOLD = 0.40911
R_IN_MANIFOLD = 0.009
DH_MANIFOLD = 2.0 * R_IN_MANIFOLD
AREA_MANIFOLD = np.pi * R_IN_MANIFOLD**2
N_MANIFOLD = 5

L_OUTLET_BUFFER = 0.20
AREA_OUTLET_BUFFER = 3.0 * AREA_MANIFOLD
DH_OUTLET_BUFFER = 2.0 * np.sqrt(AREA_OUTLET_BUFFER / np.pi)
N_OUTLET_BUFFER = 5

R_VAPOR_HP = 0.0075
WICK_THICKNESS = 0.0005
HP_WALL_THICKNESS = 0.0010
R_IN_HP = R_VAPOR_HP + WICK_THICKNESS
R_OUT_HP = R_IN_HP + HP_WALL_THICKNESS
L_EVA = 0.100
L_ABA = 0.0
L_CON = 0.500
HP_N_EVA = 1
HP_N_ABA = 0
HP_N_CON = 12
HP_N_WICK = 1
HP_N_WALL = 2
POROSITY = 0.6

THIN_FIN = 0.0004
FIN_HEIGHT = 0.020
N_FIN_HEIGHT = 15
HP_EMISSIVITY = 0.6
FIN_EMISSIVITY = 0.6
RING_EMISSIVITY = 0.05
HP_CROSSFLOW_C = 0.65
HP_CROSSFLOW_K_CAL = 1.0
HP_CROSSFLOW_WAKE_FACTOR = 1.0
MIX_NODE_VOLUME_FACTOR = 0.10

DEFAULT_T_END = 50.0
DEFAULT_PRINT_EVERY_TIME = 1.0
DEFAULT_RESTART_SAVE_EVERY = 10.0

if len(HP_MULTIPLIERS_SECTOR) != N_SECTOR:
    raise ValueError("HP_MULTIPLIERS_SECTOR length must equal N_SECTOR.")
if len(HP_MULTIPLIERS_RING) != N_RING:
    raise ValueError("HP_MULTIPLIERS_RING length must equal N_RING.")
if sum(HP_MULTIPLIERS_RING) != HP_TOTAL_COUNT:
    raise ValueError("HP_MULTIPLIERS_RING sum must equal HP_TOTAL_COUNT.")

INLET_NODE_INDICES = [0, 2 * N_SECTOR, 4 * N_SECTOR]
OUTLET_NODE_INDICES = [1 * N_SECTOR, 3 * N_SECTOR, 5 * N_SECTOR]


# ---------------------------------------------------------
# 0.1 Materials
# ---------------------------------------------------------
nak = SodiumPotassium78()
mat_wall = SS316(name="SS316_wall")
mat_hp_fluid = PotassiumHP(name="HP_Fluid_K")
mat_wick = WickMaterial(
    name="WickMaterial_K",
    solid_mat=mat_wall,
    fluid_mat=mat_hp_fluid,
    porosity=POROSITY,
    r_vapor=R_VAPOR_HP,
    r_in_wall=R_IN_HP,
)


def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    pe = np.maximum(Re * Pr, 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


def build_mix_node(name):
    ring_node_length = L_RING / N_RING
    mix_length = MIX_NODE_VOLUME_FACTOR * ring_node_length
    return IncompressibleFluidVolume(
        name=name,
        volume=AREA_RING * mix_length,
        length=mix_length,
        flow_area=AREA_RING,
        hydraulic_diam=DH_RING,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=nak,
    )


def build_ring_hp(name, fluid_channel, solid_header, hp_multipliers):
    return RingHP(
        name=name,
        fluid_channel=fluid_channel,
        solid_header=solid_header,
        hp_multipliers=hp_multipliers,
        header_flow_area=AREA_RING,
        header_dh=DH_RING,
        header_heated_perimeter=PERIM_HEADER,
        hp_r_out=R_OUT_HP,
        hp_r_in=R_IN_HP,
        hp_r_vapor=R_VAPOR_HP,
        hp_L_eva=L_EVA,
        hp_L_con=L_CON,
        hp_L_aba=L_ABA,
        hp_n_eva=HP_N_EVA,
        hp_n_con=HP_N_CON,
        hp_n_aba=HP_N_ABA,
        hp_n_wick=HP_N_WICK,
        hp_n_wall=HP_N_WALL,
        porosity_hp=POROSITY,
        HP_initial_temp=HP_INITIAL_TEMP,
        hp_wall_mat=mat_wall,
        hp_fluid_mat=mat_hp_fluid,
        hp_wick_mat=mat_wick,
        fin_thickness=THIN_FIN,
        fin_height=FIN_HEIGHT,
        n_fin_height=N_FIN_HEIGHT,
        fin_wrap_ratio=(2.0 * THIN_FIN) / (2.0 * np.pi * R_OUT_HP),
        emissivity=HP_EMISSIVITY,
        fin_emissivity=FIN_EMISSIVITY,
        up_view_factor=0.0,
        down_view_factor=0.3,
        T_space=T_SPACE,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=None,
        hp_crossflow_c=HP_CROSSFLOW_C,
        hp_crossflow_k_cal=HP_CROSSFLOW_K_CAL,
        hp_crossflow_wake_factor=HP_CROSSFLOW_WAKE_FACTOR,
    )


def build_model():
    # -----------------------------------------------------
    # 1. External boundaries
    # -----------------------------------------------------
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

    # -----------------------------------------------------
    # 2. Adiabatic inlet/outlet buffers + three inlet legs/manifolds
    # -----------------------------------------------------
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

    hot_legs = [
        IncompressibleFluidChannel(
            name=f"HotLeg_{i}",
            n_nodes=N_HOT_LEG,
            total_length=L_HOT_LEG,
            flow_area=AREA_HOT_LEG,
            hydraulic_diam=DH_HOT_LEG,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        for i in range(1, 4)
    ]

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

    manifolds = [
        IncompressibleFluidChannel(
            name=f"Manifold_{i}",
            n_nodes=N_MANIFOLD,
            total_length=L_MANIFOLD,
            flow_area=AREA_MANIFOLD,
            hydraulic_diam=DH_MANIFOLD,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        for i in range(1, 4)
    ]

    # -----------------------------------------------------
    # 3. One full 360 deg ring channel + one RingHP
    # -----------------------------------------------------
    ring_channel = IncompressibleFluidChannel(
        name="CollectorRing360_Channel",
        n_nodes=N_RING,
        total_length=L_RING,
        flow_area=AREA_RING,
        hydraulic_diam=DH_RING,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=nak,
    )

    ring_mesh = Mesh2D(
        x_dim=R_OUT_RING - R_IN_RING,
        n_x=1,
        y_dim=L_RING,
        n_y=N_RING,
        geometry_type="cylindrical",
        inner_radius=R_IN_RING,
    )
    ring_solid = HeatConduction2D(
        mesh=ring_mesh,
        material=SS316(),
        name="CollectorRing360_Solid",
        initial_temp=T_INIT,
    )
    bare_area_array = ring_solid.boundaries["right"].area / 2.0
    ring_solid.boundaries["right"].add_dynamic_radiation_condition(
        emissivity=RING_EMISSIVITY,
        bare_area_array=bare_area_array,
        T_env=T_SPACE,
    )

    ring_hp = build_ring_hp(
        name="CollectorRing360_RingHP",
        fluid_channel=ring_channel,
        solid_header=ring_solid,
        hp_multipliers=HP_MULTIPLIERS_RING,
    )

    # -----------------------------------------------------
    # 4. Connect macro headers, single-ring legs, mixing nodes, and ring closure
    # -----------------------------------------------------
    inlet_mix_nodes = [
        build_mix_node(name=f"InletMix_{idx + 1}")
        for idx in range(len(INLET_NODE_INDICES))
    ]
    outlet_mix_nodes = [
        build_mix_node(name=f"OutletMix_{idx + 1}")
        for idx in range(len(OUTLET_NODE_INDICES))
    ]

    inlet_junction = InletJunction(
        name="J_InletBoundary_InletBuffer",
        from_vol=inlet_boundary,
        to_vol=inlet_buffer_channel.volumes[0],
        W_initial=W_TOTAL,
    )

    inlet_buffer_to_hot_leg = []
    hot_leg_to_inlet_mix = []
    inlet_mix_to_ring = []
    ring_to_outlet_mix = []
    outlet_mix_to_manifold = []
    manifold_to_outlet_buffer = []

    outlet_junction = FlowJunction(
        name="J_OutletBuffer_OutletBoundary",
        from_vol=outlet_buffer_channel.volumes[-1],
        to_vol=outlet_boundary,
        flow_area=AREA_OUTLET_BUFFER,
        k_loss=0.0,
    )

    for idx, ring_node_idx in enumerate(INLET_NODE_INDICES):
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
                name=f"J_HotLeg_{idx + 1}_InletMix_{idx + 1}",
                from_vol=hot_legs[idx].volumes[-1],
                to_vol=inlet_mix_nodes[idx],
                flow_area=AREA_HOT_LEG,
                k_loss=0.0,
            )
        )
        inlet_mix_to_ring.append(
            FlowJunction(
                name=f"J_InletMix_{idx + 1}_RingNode_{ring_node_idx}",
                from_vol=inlet_mix_nodes[idx],
                to_vol=ring_channel.volumes[ring_node_idx],
                flow_area=AREA_RING,
                k_loss=0.0,
            )
        )

    for idx, ring_node_idx in enumerate(OUTLET_NODE_INDICES):
        ring_to_outlet_mix.append(
            FlowJunction(
                name=f"J_RingNode_{ring_node_idx}_OutletMix_{idx + 1}",
                from_vol=ring_channel.volumes[ring_node_idx],
                to_vol=outlet_mix_nodes[idx],
                flow_area=AREA_RING,
                k_loss=0.0,
            )
        )
        outlet_mix_to_manifold.append(
            FlowJunction(
                name=f"J_OutletMix_{idx + 1}_Manifold_{idx + 1}",
                from_vol=outlet_mix_nodes[idx],
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

    ring_closure = FlowJunction(
        name="J_Ring_Closure",
        from_vol=ring_channel.volumes[-1],
        to_vol=ring_channel.volumes[0],
        flow_area=AREA_RING,
        k_loss=ring_hp.outlet_k_loss,
        dynamic_loss_params=ring_hp.outlet_dynamic_loss_params,
    )

    # -----------------------------------------------------
    # 5. Assemble network and system manager
    # -----------------------------------------------------
    all_vols = []
    all_vols.append(inlet_boundary)
    all_vols.append(outlet_boundary)
    all_vols.extend(inlet_buffer_channel.volumes)
    for channel in hot_legs:
        all_vols.extend(channel.volumes)
    all_vols.extend(inlet_mix_nodes)
    all_vols.extend(outlet_mix_nodes)
    for channel in manifolds:
        all_vols.extend(channel.volumes)
    all_vols.extend(outlet_buffer_channel.volumes)
    all_vols.extend(ring_channel.volumes)

    all_juncs = []
    all_juncs.append(inlet_junction)
    all_juncs.extend(inlet_buffer_to_hot_leg)
    all_juncs.extend(hot_leg_to_inlet_mix)
    all_juncs.extend(inlet_mix_to_ring)
    all_juncs.extend(ring_to_outlet_mix)
    all_juncs.extend(outlet_mix_to_manifold)
    all_juncs.extend(manifold_to_outlet_buffer)
    all_juncs.append(outlet_junction)
    all_juncs.append(ring_closure)
    all_juncs.extend(inlet_buffer_channel.internal_junctions)
    for channel in hot_legs:
        all_juncs.extend(channel.internal_junctions)
    for channel in manifolds:
        all_juncs.extend(channel.internal_junctions)
    all_juncs.extend(outlet_buffer_channel.internal_junctions)
    all_juncs.extend(ring_channel.internal_junctions)

    network = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    sys_mgr = SystemManager(fluid_network=network)
    sys_mgr.add_component(ring_hp)

    return {
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "inlet_buffer_channel": inlet_buffer_channel,
        "hot_legs": hot_legs,
        "inlet_mix_nodes": inlet_mix_nodes,
        "outlet_mix_nodes": outlet_mix_nodes,
        "manifolds": manifolds,
        "outlet_buffer_channel": outlet_buffer_channel,
        "ring_channel": ring_channel,
        "ring_solid": ring_solid,
        "ring_hp": ring_hp,
        "inlet_junction": inlet_junction,
        "inlet_buffer_to_hot_leg": inlet_buffer_to_hot_leg,
        "hot_leg_to_inlet_mix": hot_leg_to_inlet_mix,
        "inlet_mix_to_ring": inlet_mix_to_ring,
        "ring_to_outlet_mix": ring_to_outlet_mix,
        "outlet_mix_to_manifold": outlet_mix_to_manifold,
        "manifold_to_outlet_buffer": manifold_to_outlet_buffer,
        "outlet_junction": outlet_junction,
        "ring_closure": ring_closure,
        "all_vols": all_vols,
        "all_juncs": all_juncs,
        "network": network,
        "sys_mgr": sys_mgr,
    }


def write_history_csv(csv_path, history):
    if not history:
        print("[WARN] No history rows to write.")
        return

    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)

    print(f"CSV saved: {csv_path}")


def restart_checkpoint_path(base_path, checkpoint_time):
    stem, ext = os.path.splitext(base_path)
    if not ext:
        ext = ".npz"
    return f"{stem}_t{int(round(checkpoint_time)):04d}s{ext}"


def sync_boundary_to_network(network, boundary):
    vol_idx = network.vol_to_idx[boundary]
    network.P_vec[vol_idx] = boundary.P
    network.T_vec[vol_idx] = boundary.T
    network.h_vec[vol_idx] = boundary.h


def next_event_time(current_time, interval):
    if interval is None or interval <= 0.0:
        return None
    return (np.floor((current_time + 1.0e-12) / interval) + 1.0) * interval


def get_model_statistics(model):
    solid_nodes = int(np.size(model["ring_solid"].T))
    fluid_nodes = int(len(model["all_vols"]) - 2)
    flow_junctions = int(len(model["all_juncs"]))
    return solid_nodes, fluid_nodes, flow_junctions


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
    case_name="collector_ring_full_ringhp_buffered_50s",
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
    inlet_buffer_channel = model["inlet_buffer_channel"]
    hot_legs = model["hot_legs"]
    inlet_mix_nodes = model["inlet_mix_nodes"]
    outlet_mix_nodes = model["outlet_mix_nodes"]
    manifolds = model["manifolds"]
    outlet_buffer_channel = model["outlet_buffer_channel"]
    ring_channel = model["ring_channel"]
    inlet_junction = model["inlet_junction"]
    inlet_buffer_to_hot_leg = model["inlet_buffer_to_hot_leg"]
    hot_leg_to_inlet_mix = model["hot_leg_to_inlet_mix"]
    inlet_mix_to_ring = model["inlet_mix_to_ring"]
    ring_to_outlet_mix = model["ring_to_outlet_mix"]
    outlet_mix_to_manifold = model["outlet_mix_to_manifold"]
    manifold_to_outlet_buffer = model["manifold_to_outlet_buffer"]
    outlet_junction = model["outlet_junction"]
    ring_closure = model["ring_closure"]

    print_pre_run_summary(model, case_name)

    if restart_from is not None:
        sys_mgr.load_global_state(restart_from)
        inlet_boundary.set_boundary_state(P=P_OUTLET + 5000.0, T=T_INLET)
        sync_boundary_to_network(network, inlet_boundary)
        outlet_boundary.set_boundary_state(P=P_OUTLET)
        sync_boundary_to_network(network, outlet_boundary)
        print(f"Restart loaded from: {restart_from}")
        print(f"Restart time: {sys_mgr.global_time:.6f} s")
    else:
        sys_mgr.initialize_system()
        print("System initialized from initial condition.")

    history = []
    next_print_time = next_event_time(sys_mgr.global_time, print_every_time)
    next_restart_save_time = next_event_time(sys_mgr.global_time, restart_save_every)

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

        w_in_total = float(inlet_junction.W)
        w_out_total = float(outlet_junction.W)
        w_ring_in_total = float(sum(j.W for j in inlet_mix_to_ring))
        w_ring_out_total = float(sum(j.W for j in ring_to_outlet_mix))
        t_outlet_list = [channel.volumes[-1].T for channel in manifolds]
        t_out_avg = float(np.mean(t_outlet_list))

        row = {
            "time": current_t,
            "dt": dt,
            "W_in_total": w_in_total,
            "W_out_total": w_out_total,
            "W_ring_in_total": w_ring_in_total,
            "W_ring_out_total": w_ring_out_total,
            "W_ring_closure": float(ring_closure.W),
            "T_out_avg": t_out_avg,
            "T_inlet_buffer_out": float(inlet_buffer_channel.volumes[-1].T),
            "T_outlet_1": float(t_outlet_list[0]),
            "T_outlet_2": float(t_outlet_list[1]),
            "T_outlet_3": float(t_outlet_list[2]),
            "T_outlet_buffer_out": float(outlet_buffer_channel.volumes[-1].T),
        }

        row["W_inlet_boundary"] = float(inlet_junction.W)
        row["W_outlet_boundary"] = float(outlet_junction.W)
        for idx, junc in enumerate(inlet_buffer_to_hot_leg, start=1):
            row[f"W_inlet_buffer_to_hotleg_{idx}"] = float(junc.W)
            row[f"W_inlet_buffer_to_hotleg_macro_{idx}"] = float(
                junc.get_mass_flow_for(inlet_buffer_channel.volumes[-1])
            )
        for idx, junc in enumerate(hot_leg_to_inlet_mix, start=1):
            row[f"W_hotleg_to_inlet_mix_{idx}"] = float(junc.W)
        for idx, junc in enumerate(inlet_mix_to_ring, start=1):
            row[f"W_inlet_mix_to_ring_{idx}"] = float(junc.W)
        for idx, junc in enumerate(ring_to_outlet_mix, start=1):
            row[f"W_ring_to_outlet_mix_{idx}"] = float(junc.W)
        for idx, junc in enumerate(outlet_mix_to_manifold, start=1):
            row[f"W_outlet_mix_to_manifold_{idx}"] = float(junc.W)
        for idx, junc in enumerate(manifold_to_outlet_buffer, start=1):
            row[f"W_manifold_to_outlet_buffer_{idx}"] = float(junc.W)
            row[f"W_manifold_to_outlet_buffer_macro_{idx}"] = float(
                junc.get_mass_flow_for(outlet_buffer_channel.volumes[0])
            )
        for idx, channel in enumerate(hot_legs, start=1):
            row[f"T_hotleg_{idx}_out"] = float(channel.volumes[-1].T)
        for idx, node in enumerate(inlet_mix_nodes, start=1):
            row[f"T_inlet_mix_{idx}"] = float(node.T)
            row[f"P_inlet_mix_{idx}"] = float(node.P)
        for idx, node in enumerate(outlet_mix_nodes, start=1):
            row[f"T_outlet_mix_{idx}"] = float(node.T)
            row[f"P_outlet_mix_{idx}"] = float(node.P)
        for idx, channel in enumerate(manifolds, start=1):
            row[f"T_manifold_{idx}_out"] = float(channel.volumes[-1].T)
        for idx in range(3):
            row[f"T_ring_inlet_node_{idx + 1}"] = float(
                ring_channel.volumes[INLET_NODE_INDICES[idx]].T
            )
            row[f"P_ring_inlet_node_{idx + 1}"] = float(
                ring_channel.volumes[INLET_NODE_INDICES[idx]].P
            )
            row[f"T_ring_outlet_node_{idx + 1}"] = float(
                ring_channel.volumes[OUTLET_NODE_INDICES[idx]].T
            )
            row[f"P_ring_outlet_node_{idx + 1}"] = float(
                ring_channel.volumes[OUTLET_NODE_INDICES[idx]].P
            )

        history.append(row)

        should_print = current_t >= t_end
        if next_print_time is not None and current_t >= next_print_time - 1.0e-12:
            should_print = True

        if should_print:
            print(
                f"t = {current_t:8.3f} s | "
                f"T_out_avg = {t_out_avg:.3f} K | "
                f"W_in_total = {w_in_total:.4f} kg/s | "
                f"W_ring_in_total = {w_ring_in_total:.4f} kg/s | "
                f"W_ring_closure = {ring_closure.W:.4f} kg/s"
            )
            while next_print_time is not None and current_t >= next_print_time - 1.0e-12:
                next_print_time += print_every_time

        if (
            restart_save_path is not None
            and restart_save_every > 0.0
            and next_restart_save_time is not None
            and current_t >= next_restart_save_time - 1.0e-12
        ):
            checkpoint_path = restart_checkpoint_path(
                restart_save_path,
                next_restart_save_time,
            )
            sys_mgr.save_global_state(checkpoint_path)
            print(f"Restart saved at t={current_t:.3f} s: {checkpoint_path}")
            while next_restart_save_time is not None and current_t >= next_restart_save_time - 1.0e-12:
                next_restart_save_time += restart_save_every

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}_history.csv")

    write_history_csv(csv_path, history)

    if restart_save_path is not None:
        sys_mgr.save_global_state(restart_save_path)
        print(f"Final restart saved: {restart_save_path}")

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print("=" * 70)

    return model, history


def print_model_summary(model):
    print("Model: geometry100hp potassium full-ring RingHP with inlet/outlet mixing nodes")
    print(f"Volumes: {len(model['all_vols'])}")
    print(f"Junctions: {len(model['all_juncs'])}")
    print(f"Nodes per 1/6    : {N_SECTOR}")
    print(f"HPs in whole ring : {sum(HP_MULTIPLIERS_RING)}")
    print(f"HP multipliers    : {HP_MULTIPLIERS_RING}")
    print(f"HP/fin emissivity : {HP_EMISSIVITY:.3f} / {FIN_EMISSIVITY:.3f}")
    print(f"Ring emissivity   : {RING_EMISSIVITY:.3f}, T_space = {T_SPACE:.1f} K")
    print("  Topology: inlet boundary -> inlet buffer -> MacroFlow x3 -> single hot legs")
    print("            -> inlet mix nodes -> ring")
    print("            ring -> outlet mix nodes -> single manifolds -> MacroFlow x3")
    print("            -> outlet buffer -> outlet boundary")
    print("  Inlet node indices :", INLET_NODE_INDICES)
    print("  Outlet node indices:", OUTLET_NODE_INDICES)
    print("  Ring closure       : last node -> node 0")


if __name__ == "__main__":
    run_case(
        case_name="collector_ring_full_ringhp_geometry100hp_potassium_mixed_50s",
        t_end=DEFAULT_T_END,
        print_every_time=DEFAULT_PRINT_EVERY_TIME,
        csv_path=os.path.join(current_dir, "collector_ring_full_ringhp_geometry100hp_potassium_mixed_50s_history.csv"),
        restart_save_path=os.path.join(current_dir, "collector_ring_full_ringhp_geometry100hp_potassium_mixed_50s_restart.npz"),
        restart_save_every=DEFAULT_RESTART_SAVE_EVERY,
        restart_from=None,
    )
