import csv
import logging
import os
import sys
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import (
    IncompressibleFluidChannel, IncompressibleFluidVolume,
    FlowJunction, MacroFlowJunction,
)
from Solvers.Hydrodynamics.BoundaryVolume import (
    IncompressibleBoundaryVolume, InletJunction,
)
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WickMaterial import WickMaterial
from Components.RingHP import RingHP

logging.basicConfig(level=logging.WARNING)

# ==========================================
# 辅助函数
# ==========================================
def get_boundary_radiation_by_node(boundary_region):
    q_node = np.zeros(boundary_region.shape, dtype=float)
    for condition in boundary_region.conditions:
        q_flux = getattr(condition, "last_q_flux", getattr(condition, "q_flux", None))
        if q_flux is not None:
            q_node += -np.array(q_flux, dtype=float)
    return q_node.reshape(-1)

def get_ring_wall_radiation_by_node(solid_ring):
    return get_boundary_radiation_by_node(solid_ring.boundaries["right"])

def get_hp_rejection_by_node(ring_hp, hp_multipliers):
    n_nodes = len(hp_multipliers)
    q_nodes = np.zeros(n_nodes, dtype=float)
    hp_idx = 0
    for i, multiplier in enumerate(hp_multipliers):
        if multiplier <= 0:
            continue
        hp_unit = ring_hp.hp_units[hp_idx]
        hp_idx += 1
        q_aba_dist, _ = hp_unit.get_heat_rejection_distribution()
        breakdown = hp_unit.get_heat_exchange_breakdown()
        q_single = (
            float(np.sum(q_aba_dist))
            + float(np.sum(breakdown["bare_radiation"]))
            + float(np.sum(breakdown["fin_radiation"]))
        )
        q_nodes[i] = q_single * float(multiplier)
    return q_nodes

def write_history_csv(csv_path, history):
    if not history:
        print("[WARN] No history rows to write.")
        return
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(history[0].keys()))
        writer.writeheader()
        writer.writerows(history)
    print(f"CSV saved: {csv_path}")

# ==========================================
# 0. 几何与工况参数
# ==========================================
T_space = 3.0
T_INLET = 843.0
W_TOTAL = 2.2
W_HOT_LEG = W_TOTAL / 3.0  # 单个进口流量
P_OUTLET = 160000.0
T_INIT = T_INLET

L_SECTOR = 0.793
N_SECTOR = 10
R_IN_RING = 0.020835
R_OUT_RING = 0.022835
DIAM_RING = 0.04167
AREA_RING = 0.0016065
PERIM_HEADER = 2.0 * np.pi * (DIAM_RING / 2.0)
HP_MULTIPLIERS = [3, 2, 3, 2, 3, 2, 3, 2, 3, 3]

L_HOT_LEG = 2.19632
R_IN_HOT_LEG = 0.0138
DH_HOT_LEG = 2.0 * R_IN_HOT_LEG
AREA_HOT_LEG = np.pi * R_IN_HOT_LEG**2
N_HOT_LEG = 28

L_MANIFOLD = 0.40911
R_IN_MANIFOLD = 0.009
DH_MANIFOLD = 2.0 * R_IN_MANIFOLD
AREA_MANIFOLD = np.pi * R_IN_MANIFOLD**2
N_MANIFOLD = 5

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

nak = SodiumPotassium78()
mat_fluid = nak
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

# ==========================================
# 1. 进出口边界
# ==========================================
inlet_boundary = IncompressibleBoundaryVolume(
    name="Full_Inlet_Boundary",
    material=nak,
    P=P_OUTLET + 5000.0,
    T=T_INLET,
)

outlet_boundary = IncompressibleBoundaryVolume(
    name="Full_Outlet_Boundary",
    material=nak,
    P=P_OUTLET,
    T=T_INIT,
)
outlet_boundary.is_pressure_boundary = True

# ==========================================
# 2. 构造辅助函数
# ==========================================
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
        header_dh=DIAM_RING,
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
        T_space=T_space,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=None,
    )

# ==========================================
# 3. 创建通道、实体和热管
# ==========================================
hot_legs = []
manifolds = []

# 3个热支路 和 3个汇流段
for i in range(1, 4):
    hl = IncompressibleFluidChannel(
        name=f"Hot_Leg_{i}",
        n_nodes=N_HOT_LEG,
        total_length=L_HOT_LEG,
        flow_area=AREA_HOT_LEG,
        hydraulic_diam=DH_HOT_LEG,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=mat_fluid,
    )
    hot_legs.append(hl)

    mf = IncompressibleFluidChannel(
        name=f"Manifold_{i}",
        n_nodes=N_MANIFOLD,
        total_length=L_MANIFOLD,
        flow_area=AREA_MANIFOLD,
        hydraulic_diam=DH_MANIFOLD,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=mat_fluid,
    )
    manifolds.append(mf)

sectors_fluid = []
sectors_solid = []
sectors_hp = []

# 6个 60度扇区
for i in range(1, 7):
    chan = IncompressibleFluidChannel(
        name=f"Sector_{i}_Channel",
        n_nodes=N_SECTOR,
        total_length=L_SECTOR,
        flow_area=AREA_RING,
        hydraulic_diam=DIAM_RING,
        initial_P=P_OUTLET,
        initial_T=T_INIT,
        material=mat_fluid,
    )
    sectors_fluid.append(chan)

    mesh = Mesh2D(
        x_dim=R_OUT_RING-R_IN_RING,
        n_x=1,
        y_dim=L_SECTOR,
        n_y=N_SECTOR,
        geometry_type="cylindrical",
        inner_radius=R_IN_RING,
    )
    solid = HeatConduction2D(
        mesh=mesh,
        material=SS316(),
        name=f"Solid_Sector_{i}",
        initial_temp=T_INIT,
    )
    bare_area_array = solid.boundaries['right'].area / 2
    solid.boundaries['right'].add_dynamic_radiation_condition(emissivity=0.6, bare_area_array=bare_area_array, T_env=T_space)
    sectors_solid.append(solid)

    ring_hp = build_ring_hp(
        name=f"Sector_{i}_RingHP",
        fluid_channel=chan,
        solid_header=solid,
        hp_multipliers=HP_MULTIPLIERS,
    )
    sectors_hp.append(ring_hp)

# ==========================================
# 4. 构建流网连接
# ==========================================
all_vols = [inlet_boundary, outlet_boundary]
for hl in hot_legs: all_vols.extend(hl.volumes)
for mf in manifolds: all_vols.extend(mf.volumes)
for sf in sectors_fluid: all_vols.extend(sf.volumes)

all_juncs = []

# 4.1 Inlet -> Hot Legs
for i, hl in enumerate(hot_legs):
    junc = InletJunction(
        name=f"J_Inlet_HotLeg_{i+1}",
        from_vol=inlet_boundary,
        to_vol=hl.volumes[0],
        W_initial=W_HOT_LEG,
    )
    all_juncs.append(junc)
    all_juncs.extend(hl.internal_junctions)

# 4.2 Manifolds -> Outlet
for i, mf in enumerate(manifolds):
    junc = FlowJunction(
        name=f"J_Manifold_{i+1}_Outlet",
        from_vol=mf.volumes[-1],
        to_vol=outlet_boundary,
        flow_area=AREA_MANIFOLD,
        k_loss=0.0,
    )
    all_juncs.append(junc)
    all_juncs.extend(mf.internal_junctions)

# 4.3 扇区内部连接
for sf in sectors_fluid:
    all_juncs.extend(sf.internal_junctions)

# 4.4 拓扑连接: 热支路 -> 扇区 (分配)
# Inlet 1 -> Sector 1, Sector 6
all_juncs.append(FlowJunction(name="J_HL1_S1", from_vol=hot_legs[0].volumes[-1], to_vol=sectors_fluid[0].volumes[0], flow_area=AREA_RING, k_loss=0.0))
all_juncs.append(FlowJunction(name="J_HL1_S6", from_vol=hot_legs[0].volumes[-1], to_vol=sectors_fluid[5].volumes[0], flow_area=AREA_RING, k_loss=0.0))

# Inlet 2 -> Sector 2, Sector 3
all_juncs.append(FlowJunction(name="J_HL2_S2", from_vol=hot_legs[1].volumes[-1], to_vol=sectors_fluid[1].volumes[0], flow_area=AREA_RING, k_loss=0.0))
all_juncs.append(FlowJunction(name="J_HL2_S3", from_vol=hot_legs[1].volumes[-1], to_vol=sectors_fluid[2].volumes[0], flow_area=AREA_RING, k_loss=0.0))

# Inlet 3 -> Sector 4, Sector 5
all_juncs.append(FlowJunction(name="J_HL3_S4", from_vol=hot_legs[2].volumes[-1], to_vol=sectors_fluid[3].volumes[0], flow_area=AREA_RING, k_loss=0.0))
all_juncs.append(FlowJunction(name="J_HL3_S5", from_vol=hot_legs[2].volumes[-1], to_vol=sectors_fluid[4].volumes[0], flow_area=AREA_RING, k_loss=0.0))

# 4.5 拓扑连接: 扇区 -> 汇流段 (汇集)
# Sector 1, 2 -> Outlet 1
all_juncs.append(FlowJunction(name="J_S1_MF1", from_vol=sectors_fluid[0].volumes[-1], to_vol=manifolds[0].volumes[0], flow_area=AREA_RING, k_loss=0.0))
all_juncs.append(FlowJunction(name="J_S2_MF1", from_vol=sectors_fluid[1].volumes[-1], to_vol=manifolds[0].volumes[0], flow_area=AREA_RING, k_loss=0.0))

# Sector 3, 4 -> Outlet 2
all_juncs.append(FlowJunction(name="J_S3_MF2", from_vol=sectors_fluid[2].volumes[-1], to_vol=manifolds[1].volumes[0], flow_area=AREA_RING, k_loss=0.0))
all_juncs.append(FlowJunction(name="J_S4_MF2", from_vol=sectors_fluid[3].volumes[-1], to_vol=manifolds[1].volumes[0], flow_area=AREA_RING, k_loss=0.0))

# Sector 5, 6 -> Outlet 3
all_juncs.append(FlowJunction(name="J_S5_MF3", from_vol=sectors_fluid[4].volumes[-1], to_vol=manifolds[2].volumes[0], flow_area=AREA_RING, k_loss=0.0))
all_juncs.append(FlowJunction(name="J_S6_MF3", from_vol=sectors_fluid[5].volumes[-1], to_vol=manifolds[2].volumes[0], flow_area=AREA_RING, k_loss=0.0))

# ==========================================
# 5. 系统初始化
# ==========================================
network = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
sys_mgr = SystemManager(fluid_network=network)

for shp in sectors_hp:
    sys_mgr.add_component(shp)

# ==========================================
# 6. 运行控制
# ==========================================
def run_case(
    case_name="Full_Ring_Model",
    t_end=50.0,
    min_dt=1e-3,
    max_dt=0.5,
    safety_factor=1.0,
    inner_iter=2,
    print_every=10,
    csv_path=None
):
    print("=" * 70)
    print(f"Running case: {case_name}")
    print("=" * 70)

    sys_mgr.initialize_system()
    print(f"System initialized. Total volumes: {len(all_vols)}, Total junctions: {len(all_juncs)}")

    history = []
    step_count = 0

    while sys_mgr.global_time < t_end:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=min_dt,
            max_dt=max_dt,
            safety_factor=safety_factor,
        )
        dt = min(dt, t_end - sys_mgr.global_time)

        sys_mgr.step(dt=dt, inner_iter=inner_iter)

        step_count += 1
        current_t = sys_mgr.global_time

        # 计算总散热
        q_wall_total = sum(float(np.sum(get_ring_wall_radiation_by_node(solid))) for solid in sectors_solid)
        q_hp_total = sum(float(np.sum(get_hp_rejection_by_node(hp, HP_MULTIPLIERS))) for hp in sectors_hp)
        q_total = q_wall_total + q_hp_total

        # 计算总流量
        w_in_total = sum(j.W for j in all_juncs if j.name.startswith("J_Inlet_HotLeg"))
        w_out_total = sum(j.W for j in all_juncs if j.name.endswith("Outlet"))

        # 平均出口温度
        T_out_avg = np.mean([mf.volumes[-1].T for mf in manifolds])

        row = {
            "time": current_t,
            "dt": dt,
            "W_in_total": w_in_total,
            "W_out_total": w_out_total,
            "T_out_avg": float(T_out_avg),
            "Q_wall_total": q_wall_total,
            "Q_hp_total": q_hp_total,
            "Q_system_total": q_total,
        }
        
        # 记录关键分支流量
        row["W_HL1_S1"] = next(j.W for j in all_juncs if j.name == "J_HL1_S1")
        row["W_HL1_S6"] = next(j.W for j in all_juncs if j.name == "J_HL1_S6")
        row["W_HL2_S2"] = next(j.W for j in all_juncs if j.name == "J_HL2_S2")
        row["W_HL2_S3"] = next(j.W for j in all_juncs if j.name == "J_HL2_S3")
        row["W_HL3_S4"] = next(j.W for j in all_juncs if j.name == "J_HL3_S4")
        row["W_HL3_S5"] = next(j.W for j in all_juncs if j.name == "J_HL3_S5")
        
        history.append(row)

        if step_count % print_every == 0 or current_t >= t_end:
            print(
                f"t = {current_t:8.3f} s | "
                f"T_out_avg = {T_out_avg:.3f} K | "
                f"W_in_total = {w_in_total:.4f} kg/s | "
                f"Q_total = {q_total:.3f} W"
            )

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}_history.csv")

    write_history_csv(csv_path, history)

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print(f"CSV: {csv_path}")
    print("=" * 70)

    return history

if __name__ == "__main__":
    history = run_case(
        case_name="Full_Ring_Model",
        t_end=50.0,
        csv_path=os.path.join(current_dir, "Full_Ring_Model_history.csv"),
    )
