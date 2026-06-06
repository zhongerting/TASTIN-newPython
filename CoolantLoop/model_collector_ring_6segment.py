import csv
import io
import logging
import os
import sys
import time
from contextlib import redirect_stdout
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
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WickMaterial import WickMaterial
from Components.RingHP import RingHP
from profiler import TEASAProfiler


logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------
# 0. 几何与工况参数
# ---------------------------------------------------------
T_SPACE = 3.0
T_INLET = 843.0
P_OUTLET = 160000.0
T_INIT = 863.0

W_TOTAL = 2.2
W_BRANCH_TOTAL = W_TOTAL / 3.0
W_INLET_LEG_INIT = W_BRANCH_TOTAL

L_INLET_BUFFER = 0.20
N_INLET_BUFFER = 5

L_SECTOR = 0.793
N_SECTOR = 4
R_IN_RING = 0.020835
R_OUT_RING = 0.022835
DH_RING = 0.04167
AREA_RING = 0.0016065
PERIM_HEADER = 2.0 * np.pi * (DH_RING / 2.0)
HP_COUNT_PER_SECTOR = 26
HP_MULTIPLIERS_SECTOR = [6, 7, 6, 7]

L_HOT_LEG = 2.19632
R_IN_HOT_LEG = 0.0138
DH_HOT_LEG = 2.0 * R_IN_HOT_LEG
AREA_HOT_LEG = np.pi * R_IN_HOT_LEG**2
N_HOT_LEG = 28

AREA_INLET_BUFFER = 3.0 * AREA_HOT_LEG
DH_INLET_BUFFER = 2.0 * np.sqrt(AREA_INLET_BUFFER / np.pi)

L_OUTLET_BRANCH = 0.40911
R_IN_OUTLET_BRANCH = 0.009
DH_OUTLET_BRANCH = 2.0 * R_IN_OUTLET_BRANCH
AREA_OUTLET_BRANCH = np.pi * R_IN_OUTLET_BRANCH**2
N_OUTLET_BRANCH = 5

L_OUTLET_BUFFER = 0.20
AREA_OUTLET_BUFFER = 3.0 * AREA_OUTLET_BRANCH
DH_OUTLET_BUFFER = 2.0 * np.sqrt(AREA_OUTLET_BUFFER / np.pi)
N_OUTLET_BUFFER = 5

R_OUT_HP = 0.0085
R_IN_HP = 0.0081
R_VAPOR_HP = 0.0075
L_EVA = 0.0605
L_ABA = 0.0415
L_CON = 0.47
POROSITY = 0.966

THIN_FIN = 0.0003
FIN_HEIGHT = 22.65e-3
N_FIN_HEIGHT = 15

DEFAULT_T_END = 50.0
DEFAULT_PRINT_EVERY_TIME = 1.0
DEFAULT_RESTART_SAVE_EVERY = 10.0

PROFILER_KEY_FUNCTIONS = [
    "SystemManager.initialize_system",
    "SystemManager.compute_adaptive_dt",
    "SystemManager.step",
    "HydraulicNetwork.get_max_stable_dt",
    "HydraulicNetwork.step",
    "HydraulicNetwork.step_hydraulic",
    "HydraulicNetwork._step_energy_implicit",
    "HeatConduction2D.step",
    "FluidSolidCouple.execute",
]

if len(HP_MULTIPLIERS_SECTOR) != N_SECTOR:
    raise ValueError("HP_MULTIPLIERS_SECTOR length must equal N_SECTOR.")
if sum(HP_MULTIPLIERS_SECTOR) != HP_COUNT_PER_SECTOR:
    raise ValueError("HP_MULTIPLIERS_SECTOR sum must equal HP_COUNT_PER_SECTOR.")


# ---------------------------------------------------------
# 0.1 材料和组件
# ---------------------------------------------------------
nak = SodiumPotassium78()
mat_wall = SS316(name="SS316_wall")
mat_hp_fluid = SodiumHP(name="HP_Fluid_Na")
mat_wick = WickMaterial(
    name="WickMaterial",
    solid_mat=mat_wall,
    fluid_mat=mat_hp_fluid,
    porosity=POROSITY,
    r_vapor=R_VAPOR_HP,
    r_in_wall=R_IN_HP,
)


def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    pe = np.maximum(Re * Pr, 1.0)
    return 7.0 + 0.025 * (pe ** 0.8)


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
        hp_n_eva=1,
        hp_n_con=12,
        hp_n_aba=1,
        hp_n_wick=1,
        hp_n_wall=2,
        porosity_hp=POROSITY,
        HP_initial_temp=800.0,
        hp_wall_mat=mat_wall,
        hp_fluid_mat=mat_hp_fluid,
        hp_wick_mat=mat_wick,
        fin_thickness=THIN_FIN,
        fin_height=FIN_HEIGHT,
        n_fin_height=N_FIN_HEIGHT,
        fin_wrap_ratio=(2.0 * THIN_FIN) / (2.0 * np.pi * R_OUT_HP),
        emissivity=0.85,
        up_view_factor=0.0,
        down_view_factor=0.3,
        T_space=T_SPACE,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=None,
    )


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
    bare_area_array = solid.boundaries["right"].area / 2.0
    solid.boundaries["right"].add_dynamic_radiation_condition(
        emissivity=0.6,
        bare_area_array=bare_area_array,
        T_env=T_SPACE,
    )
    return solid


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
    # 2. Adiabatic inlet/outlet buffers + three inlet legs/outlet branches
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

    outlet_branches = [
        IncompressibleFluidChannel(
            name=f"OutletBranch_{i}",
            n_nodes=N_OUTLET_BRANCH,
            total_length=L_OUTLET_BRANCH,
            flow_area=AREA_OUTLET_BRANCH,
            hydraulic_diam=DH_OUTLET_BRANCH,
            initial_P=P_OUTLET,
            initial_T=T_INIT,
            material=nak,
        )
        for i in range(1, 4)
    ]

    # -----------------------------------------------------
    # 3. Build 6 stitched 1/6 RingHP sectors
    # -----------------------------------------------------
    sector_specs = [
        ("S1_I1_to_O1", "I1", "O1"),
        ("S2_O1_to_I2", "O1", "I2"),
        ("S3_I2_to_O2", "I2", "O2"),
        ("S4_O2_to_I3", "O2", "I3"),
        ("S5_I3_to_O3", "I3", "O3"),
        ("S6_O3_to_I1", "O3", "I1"),
    ]

    sectors = []
    solids = []
    ring_hps = []

    for sector_name, start_key, end_key in sector_specs:
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
        ring_hp = build_ring_hp(
            name=f"{sector_name}_RingHP",
            fluid_channel=channel,
            solid_header=solid,
            hp_multipliers=HP_MULTIPLIERS_SECTOR,
        )

        sectors.append(channel)
        solids.append(solid)
        ring_hps.append(ring_hp)

    # -----------------------------------------------------
    # 4. Use sector endpoint control volumes as shared ring nodes:
    #    I1=S1[0], O1=S2[0], I2=S3[0], O2=S4[0], I3=S5[0], O3=S6[0]
    # -----------------------------------------------------
    ring_nodes = {
        "I1": sectors[0].volumes[0],
        "O1": sectors[1].volumes[0],
        "I2": sectors[2].volumes[0],
        "O2": sectors[3].volumes[0],
        "I3": sectors[4].volumes[0],
        "O3": sectors[5].volumes[0],
    }

    sector_link_junctions = [
        FlowJunction(
            name="J_S1_O1_to_S2_O1",
            from_vol=sectors[0].volumes[-1],
            to_vol=ring_nodes["O1"],
            flow_area=AREA_RING,
            k_loss=ring_hps[0].outlet_k_loss,
            dynamic_loss_params=ring_hps[0].outlet_dynamic_loss_params,
        ),
        FlowJunction(
            name="J_S2_I2_to_S3_I2",
            from_vol=sectors[1].volumes[-1],
            to_vol=ring_nodes["I2"],
            flow_area=AREA_RING,
            k_loss=ring_hps[1].outlet_k_loss,
            dynamic_loss_params=ring_hps[1].outlet_dynamic_loss_params,
        ),
        FlowJunction(
            name="J_S3_O2_to_S4_O2",
            from_vol=sectors[2].volumes[-1],
            to_vol=ring_nodes["O2"],
            flow_area=AREA_RING,
            k_loss=ring_hps[2].outlet_k_loss,
            dynamic_loss_params=ring_hps[2].outlet_dynamic_loss_params,
        ),
        FlowJunction(
            name="J_S4_I3_to_S5_I3",
            from_vol=sectors[3].volumes[-1],
            to_vol=ring_nodes["I3"],
            flow_area=AREA_RING,
            k_loss=ring_hps[3].outlet_k_loss,
            dynamic_loss_params=ring_hps[3].outlet_dynamic_loss_params,
        ),
        FlowJunction(
            name="J_S5_O3_to_S6_O3",
            from_vol=sectors[4].volumes[-1],
            to_vol=ring_nodes["O3"],
            flow_area=AREA_RING,
            k_loss=ring_hps[4].outlet_k_loss,
            dynamic_loss_params=ring_hps[4].outlet_dynamic_loss_params,
        ),
        FlowJunction(
            name="J_S6_I1_to_S1_I1",
            from_vol=sectors[5].volumes[-1],
            to_vol=ring_nodes["I1"],
            flow_area=AREA_RING,
            k_loss=ring_hps[5].outlet_k_loss,
            dynamic_loss_params=ring_hps[5].outlet_dynamic_loss_params,
        ),
    ]

    # -----------------------------------------------------
    # 5. Connect inlet buffer -> three legs -> ring inlets,
    #    and ring outlets -> three outlet branches -> outlet buffer
    # -----------------------------------------------------
    inlet_junction = InletJunction(
        name="J_InletBoundary_InletBuffer",
        from_vol=inlet_boundary,
        to_vol=inlet_buffer_channel.volumes[0],
        W_initial=W_TOTAL,
    )

    inlet_buffer_to_hot_leg = []
    hot_leg_to_ring = []
    ring_to_outlet_branch = []
    outlet_branch_to_outlet_buffer = []
    outlet_junction = FlowJunction(
        name="J_OutletBuffer_OutletBoundary",
        from_vol=outlet_buffer_channel.volumes[-1],
        to_vol=outlet_boundary,
        flow_area=AREA_OUTLET_BUFFER,
        k_loss=0.0,
    )

    inlet_map = [("I1", 0), ("I2", 1), ("I3", 2)]
    outlet_map = [("O1", 0), ("O2", 1), ("O3", 2)]

    for interface_key, idx in inlet_map:
        inlet_buffer_to_hot_leg.append(
            FlowJunction(
                name=f"J_InletBuffer_HotLeg_{idx + 1}",
                from_vol=inlet_buffer_channel.volumes[-1],
                to_vol=hot_legs[idx].volumes[0],
                flow_area=AREA_HOT_LEG,
                k_loss=0.0,
            )
        )
        hot_leg_to_ring.append(
            MacroFlowJunction(
                name=f"J_HotLeg_{idx + 1}_{interface_key}",
                from_vol=hot_legs[idx].volumes[-1],
                to_vol=ring_nodes[interface_key],
                macro_vol=hot_legs[idx].volumes[-1],
                multiplier=2,
                flow_area=AREA_HOT_LEG,
                k_loss=0.0,
            )
        )

    for interface_key, idx in outlet_map:
        ring_to_outlet_branch.append(
            MacroFlowJunction(
                name=f"J_{interface_key}_OutletBranch_{idx + 1}",
                from_vol=ring_nodes[interface_key],
                to_vol=outlet_branches[idx].volumes[0],
                macro_vol=outlet_branches[idx].volumes[0],
                multiplier=2,
                flow_area=AREA_OUTLET_BRANCH,
                k_loss=0.0,
            )
        )
        outlet_branch_to_outlet_buffer.append(
            FlowJunction(
                name=f"J_OutletBranch_{idx + 1}_OutletBuffer",
                from_vol=outlet_branches[idx].volumes[-1],
                to_vol=outlet_buffer_channel.volumes[0],
                flow_area=AREA_OUTLET_BRANCH,
                k_loss=0.0,
            )
        )

    # -----------------------------------------------------
    # 6. Assemble network and system manager
    # -----------------------------------------------------
    all_vols = []
    all_vols.append(inlet_boundary)
    all_vols.append(outlet_boundary)
    all_vols.extend(inlet_buffer_channel.volumes)
    for channel in hot_legs:
        all_vols.extend(channel.volumes)
    for channel in outlet_branches:
        all_vols.extend(channel.volumes)
    all_vols.extend(outlet_buffer_channel.volumes)
    for channel in sectors:
        all_vols.extend(channel.volumes)

    all_juncs = []
    all_juncs.append(inlet_junction)
    all_juncs.extend(inlet_buffer_to_hot_leg)
    all_juncs.extend(hot_leg_to_ring)
    all_juncs.extend(ring_to_outlet_branch)
    all_juncs.extend(outlet_branch_to_outlet_buffer)
    all_juncs.append(outlet_junction)
    all_juncs.extend(sector_link_junctions)
    all_juncs.extend(inlet_buffer_channel.internal_junctions)
    for channel in hot_legs:
        all_juncs.extend(channel.internal_junctions)
    for channel in outlet_branches:
        all_juncs.extend(channel.internal_junctions)
    all_juncs.extend(outlet_buffer_channel.internal_junctions)
    for channel in sectors:
        all_juncs.extend(channel.internal_junctions)

    network = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    sys_mgr = SystemManager(fluid_network=network)
    for ring_hp in ring_hps:
        sys_mgr.add_component(ring_hp)

    return {
        "inlet_boundary": inlet_boundary,
        "outlet_boundary": outlet_boundary,
        "inlet_buffer_channel": inlet_buffer_channel,
        "hot_legs": hot_legs,
        "outlet_branches": outlet_branches,
        "outlet_buffer_channel": outlet_buffer_channel,
        "ring_nodes": ring_nodes,
        "sectors": sectors,
        "solids": solids,
        "ring_hps": ring_hps,
        "inlet_junction": inlet_junction,
        "inlet_buffer_to_hot_leg": inlet_buffer_to_hot_leg,
        "hot_leg_to_ring": hot_leg_to_ring,
        "ring_to_outlet_branch": ring_to_outlet_branch,
        "outlet_branch_to_outlet_buffer": outlet_branch_to_outlet_buffer,
        "outlet_junction": outlet_junction,
        "sector_link_junctions": sector_link_junctions,
        "all_vols": all_vols,
        "all_juncs": all_juncs,
        "network": network,
        "sys_mgr": sys_mgr,
        "sector_specs": sector_specs,
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


def install_profiler_hooks():
    def wrap_method(cls, method_name):
        method = getattr(cls, method_name)
        if getattr(method, "_teasa_profile_wrapped", False):
            return
        if hasattr(method, "__wrapped__"):
            return
        wrapped = TEASAProfiler.profile(method)
        wrapped._teasa_profile_wrapped = True
        setattr(cls, method_name, wrapped)

    wrap_method(SystemManager, "initialize_system")
    wrap_method(SystemManager, "compute_adaptive_dt")
    wrap_method(HydraulicNetwork, "get_max_stable_dt")


def reset_profiler_stats():
    TEASAProfiler.stats = {}


def get_profiler_snapshot_rows(sim_time, cpu_elapsed_s, wall_elapsed_s):
    rows = []
    for func_name in PROFILER_KEY_FUNCTIONS:
        data = TEASAProfiler.stats.get(func_name, {"count": 0, "time": 0.0})
        count = int(data["count"])
        total_time_s = float(data["time"])
        avg_time_ms = 1000.0 * total_time_s / count if count > 0 else 0.0
        rows.append(
            {
                "sim_time_s": float(sim_time),
                "cpu_elapsed_s": float(cpu_elapsed_s),
                "wall_elapsed_s": float(wall_elapsed_s),
                "function_name": func_name,
                "call_count": count,
                "total_time_s": total_time_s,
                "avg_time_ms": avg_time_ms,
            }
        )
    return rows


def get_profiler_summary_rows():
    rows = []
    for func_name in PROFILER_KEY_FUNCTIONS:
        data = TEASAProfiler.stats.get(func_name, {"count": 0, "time": 0.0})
        count = int(data["count"])
        total_time_s = float(data["time"])
        avg_time_ms = 1000.0 * total_time_s / count if count > 0 else 0.0
        rows.append(
            {
                "function_name": func_name,
                "call_count": count,
                "total_time_s": total_time_s,
                "avg_time_ms": avg_time_ms,
            }
        )
    return rows


def write_profiler_report(report_path):
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        TEASAProfiler.report()
    report_text = buffer.getvalue()
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_text)
    print(report_text, end="")


def get_model_statistics(model):
    solid_nodes = int(sum(np.size(solid.T) for solid in model["solids"]))
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
    case_name="collector_ring_6segment_buffered_50s",
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
    profiler_summary_path=None,
    profiler_snapshot_path=None,
    profiler_report_path=None,
):
    install_profiler_hooks()
    reset_profiler_stats()
    wall_start = time.perf_counter()
    cpu_start = time.process_time()

    model = build_model()
    sys_mgr = model["sys_mgr"]
    network = model["network"]
    inlet_boundary = model["inlet_boundary"]
    outlet_boundary = model["outlet_boundary"]
    inlet_buffer_channel = model["inlet_buffer_channel"]
    hot_legs = model["hot_legs"]
    outlet_branches = model["outlet_branches"]
    outlet_buffer_channel = model["outlet_buffer_channel"]
    ring_nodes = model["ring_nodes"]
    sectors = model["sectors"]
    inlet_junction = model["inlet_junction"]
    inlet_buffer_to_hot_leg = model["inlet_buffer_to_hot_leg"]
    hot_leg_to_ring = model["hot_leg_to_ring"]
    ring_to_outlet_branch = model["ring_to_outlet_branch"]
    outlet_branch_to_outlet_buffer = model["outlet_branch_to_outlet_buffer"]
    outlet_junction = model["outlet_junction"]

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
    profiler_snapshots = []
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
        cpu_elapsed_s = time.process_time() - cpu_start
        wall_elapsed_s = time.perf_counter() - wall_start

        w_in_total = float(inlet_junction.W)
        w_out_total = float(outlet_junction.W)
        w_ring_in_total = float(sum(j.W for j in hot_leg_to_ring))
        w_ring_out_total = float(sum(j.W for j in ring_to_outlet_branch))
        t_outlet_list = [channel.volumes[-1].T for channel in outlet_branches]
        t_out_avg = float(np.mean(t_outlet_list))

        row = {
            "time": current_t,
            "dt": dt,
            "cpu_elapsed_s": float(cpu_elapsed_s),
            "wall_elapsed_s": float(wall_elapsed_s),
            "W_in_total": w_in_total,
            "W_out_total": w_out_total,
            "W_ring_in_total": w_ring_in_total,
            "W_ring_out_total": w_ring_out_total,
            "T_out_avg": t_out_avg,
            "T_inlet_buffer_out": float(inlet_buffer_channel.volumes[-1].T),
            "T_outlet_1": float(t_outlet_list[0]),
            "T_outlet_2": float(t_outlet_list[1]),
            "T_outlet_3": float(t_outlet_list[2]),
            "T_outlet_buffer_out": float(outlet_buffer_channel.volumes[-1].T),
            "P_node_I1": float(ring_nodes["I1"].P),
            "P_node_O1": float(ring_nodes["O1"].P),
            "P_node_I2": float(ring_nodes["I2"].P),
            "P_node_O2": float(ring_nodes["O2"].P),
            "P_node_I3": float(ring_nodes["I3"].P),
            "P_node_O3": float(ring_nodes["O3"].P),
        }

        row["W_inlet_boundary"] = float(inlet_junction.W)
        row["W_outlet_boundary"] = float(outlet_junction.W)
        for idx, junc in enumerate(inlet_buffer_to_hot_leg, start=1):
            row[f"W_inlet_buffer_to_hotleg_{idx}"] = float(junc.W)
        for idx, junc in enumerate(hot_leg_to_ring, start=1):
            row[f"W_hotleg_to_ring_macro_{idx}"] = float(
                junc.get_mass_flow_for(hot_legs[idx - 1].volumes[-1])
            )
            row[f"W_hotleg_to_ring_{idx}"] = float(junc.W)
        for idx, junc in enumerate(ring_to_outlet_branch, start=1):
            row[f"W_ring_to_outlet_branch_{idx}"] = float(junc.W)
            row[f"W_ring_to_outlet_branch_macro_{idx}"] = float(
                junc.get_mass_flow_for(outlet_branches[idx - 1].volumes[0])
            )
        for idx, junc in enumerate(outlet_branch_to_outlet_buffer, start=1):
            row[f"W_outlet_branch_to_outlet_buffer_{idx}"] = float(junc.W)
        for idx, channel in enumerate(sectors, start=1):
            row[f"W_sector_{idx}_entry"] = float(channel.internal_junctions[0].W)
            row[f"T_sector_{idx}_out"] = float(channel.volumes[-1].T)
        for idx, channel in enumerate(hot_legs, start=1):
            row[f"T_hotleg_{idx}_out"] = float(channel.volumes[-1].T)
        for idx, channel in enumerate(outlet_branches, start=1):
            row[f"T_outlet_branch_{idx}_out"] = float(channel.volumes[-1].T)

        history.append(row)

        should_print = current_t >= t_end
        if next_print_time is not None and current_t >= next_print_time - 1.0e-12:
            should_print = True

        if should_print:
            profiler_snapshots.extend(
                get_profiler_snapshot_rows(
                    sim_time=current_t,
                    cpu_elapsed_s=cpu_elapsed_s,
                    wall_elapsed_s=wall_elapsed_s,
                )
            )
            print(
                f"t = {current_t:8.3f} s | "
                f"T_out_avg = {t_out_avg:.3f} K | "
                f"W_in_total = {w_in_total:.4f} kg/s | "
                f"W_ring_in_total = {w_ring_in_total:.4f} kg/s | "
                f"W_out_total = {w_out_total:.4f} kg/s | "
                f"CPU = {cpu_elapsed_s:.2f} s"
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
    if profiler_summary_path is None:
        profiler_summary_path = os.path.join(
            current_dir, f"{case_name}_profiler_summary.csv"
        )
    if profiler_snapshot_path is None:
        profiler_snapshot_path = os.path.join(
            current_dir, f"{case_name}_profiler_snapshots.csv"
        )
    if profiler_report_path is None:
        profiler_report_path = os.path.join(
            current_dir, f"{case_name}_profiler_report.txt"
        )

    write_history_csv(csv_path, history)
    write_history_csv(profiler_summary_path, get_profiler_summary_rows())
    write_history_csv(profiler_snapshot_path, profiler_snapshots)
    write_profiler_report(profiler_report_path)

    if restart_save_path is not None:
        sys_mgr.save_global_state(restart_save_path)
        print(f"Final restart saved: {restart_save_path}")

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print(f"Profiler summary : {profiler_summary_path}")
    print(f"Profiler snapshots: {profiler_snapshot_path}")
    print(f"Profiler report  : {profiler_report_path}")
    print("=" * 70)

    return model, history


def print_model_summary(model):
    print("Model: 6 stitched 1/6 RingHP sectors")
    print(f"Volumes: {len(model['all_vols'])}")
    print(f"Junctions: {len(model['all_juncs'])}")
    print(f"RingHP components: {len(model['ring_hps'])}")
    print(f"Nodes per 1/6    : {N_SECTOR}")
    print(f"HPs per 1/6 sector: {HP_COUNT_PER_SECTOR}")
    print(f"HPs in whole ring : {6 * HP_COUNT_PER_SECTOR}")
    print("  Topology: inlet boundary -> inlet buffer -> 3 hot legs -> ring inlets")
    print("            (branch-side flow is halved when entering the ring)")
    print("            ring outlets -> 3 outlet branches -> outlet buffer -> outlet boundary")
    for idx, (sector_name, start_key, end_key) in enumerate(model["sector_specs"], start=1):
        print(f"  Sector {idx}: {start_key} -> {end_key} ({sector_name})")
    print("  Inlet nodes: I1, I2, I3")
    print("  Outlet nodes: O1, O2, O3")


if __name__ == "__main__":
    run_case(
        case_name="collector_ring_6segment_buffered_50s",
        t_end=DEFAULT_T_END,
        print_every_time=DEFAULT_PRINT_EVERY_TIME,
        csv_path=os.path.join(current_dir, "collector_ring_6segment_buffered_50s_history.csv"),
        restart_save_path=os.path.join(current_dir, "collector_ring_6segment_buffered_50s_restart.npz"),
        restart_save_every=DEFAULT_RESTART_SAVE_EVERY,
        restart_from=None,
    )
