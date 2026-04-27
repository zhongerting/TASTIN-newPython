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


# --- TASTIN 框架导入 ---
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

# 诊断辅助函数
def get_boundary_radiation_by_node(boundary_region):
    """
    返回某个固体边界上所有 Dynamic/Flux 类条件造成的向外散热量分布 [W/node]。

    约定：
    - Boundary 条件里的 q_flux 是流入固体为正；
    - 辐射散热时 q_flux 通常为负；
    - 因此这里返回 -q_flux，使“向外散热”为正。
    """
    q_node = np.zeros(boundary_region.shape, dtype=float)

    for condition in boundary_region.conditions:
        q_flux = getattr(condition, "last_q_flux", getattr(condition, "q_flux", None))
        if q_flux is not None:
            q_node += -np.array(q_flux, dtype=float)

    return q_node.reshape(-1)

def get_ring_wall_radiation_by_node(solid_ring):
    """
    集流环外表面直接向外辐射散热量，按节点返回 [W/node]。
    当前你的外表面辐射挂在 solid_ring.boundaries['right'] 上。
    """
    return get_boundary_radiation_by_node(solid_ring.boundaries["right"])

def get_hp_rejection_by_node(ring_hp, hp_multipliers):
    """
    返回每个流体节点对应的真实热管总散热量 [W/node]。
    这里显式按 `hp_multipliers` 折算，避免把代表热管单元误当成真实节点总量。
    """
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

def get_hp_radiation_breakdown_by_node(ring_hp, hp_multipliers):
    n_nodes = len(hp_multipliers)

    single_aba = np.zeros(n_nodes, dtype=float)
    single_con_bare = np.zeros(n_nodes, dtype=float)
    single_fin = np.zeros(n_nodes, dtype=float)

    hp_idx = 0
    for i, multiplier in enumerate(hp_multipliers):
        if multiplier <= 0:
            continue

        hp_unit = ring_hp.hp_units[hp_idx]
        hp_idx += 1

        q_aba_dist, _ = hp_unit.get_heat_rejection_distribution()
        breakdown = hp_unit.get_heat_exchange_breakdown()

        single_aba[i] = float(np.sum(q_aba_dist))
        single_con_bare[i] = float(np.sum(breakdown["bare_radiation"]))
        single_fin[i] = float(np.sum(breakdown["fin_radiation"]))

    multiplier_arr = np.asarray(hp_multipliers, dtype=float)
    single_total = single_aba + single_con_bare + single_fin

    return {
        "single_aba": single_aba,
        "single_con_bare": single_con_bare,
        "single_fin": single_fin,
        "single_total": single_total,
        "segment_aba": single_aba * multiplier_arr,
        "segment_con_bare": single_con_bare * multiplier_arr,
        "segment_fin": single_fin * multiplier_arr,
        "segment_total": single_total * multiplier_arr,
    }

def get_ring_total_rejection_by_node(solid_ring, ring_hp, hp_multipliers):
    """
    每个节点总散热量 = 集流环外壁直接辐射 + 热管/翅片散热。
    """
    q_wall = get_ring_wall_radiation_by_node(solid_ring)
    q_hp = get_hp_rejection_by_node(ring_hp, hp_multipliers)

    if len(q_wall) != len(q_hp):
        raise ValueError(
            f"Node mismatch: wall radiation has {len(q_wall)} nodes, "
            f"RingHP rejection has {len(q_hp)} nodes."
        )

    return q_wall + q_hp

# 输出诊断函数
def write_history_csv(csv_path, history):
    """
    history 是 list[dict] 时使用。
    """
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

# 运行辅助函数
def run_case(
    case_name="V5_up_ring_single_segment",
    t_end=50.0,
    min_dt=1e-3,
    max_dt=0.5,
    safety_factor=1.0,
    inner_iter=2,
    print_every=10,
    csv_path=None,
    restart_from=None,
    restart_save_path=None,
    restart_save_every=0.0,
):
    print("=" * 70)
    print(f"Running case: {case_name}")
    print("=" * 70)

    # 初始化或重启动
    if restart_from is not None:
        # 注意：必须先构建完全相同拓扑的 sys_mgr/network/up_ring，再 load。
        sys_mgr.load_global_state(restart_from)
        # 修复重启动后进口温度无法变化的bug
        inlet_boundary.set_boundary_state(T=T_INLET)
        inlet_idx = network.vol_to_idx[inlet_boundary]
        network.T_vec[inlet_idx] = inlet_boundary.T
        network.h_vec[inlet_idx] = inlet_boundary.h
        network.P_vec[inlet_idx] = inlet_boundary.P
        print(f"Restart loaded from: {restart_from}")
        print(f"Restart time: {sys_mgr.global_time:.6f} s")
    else:
        sys_mgr.initialize_system()
        print("System initialized from initial condition.")

    history = []
    step_count = 0
    next_restart_save_time = sys_mgr.global_time + restart_save_every if restart_save_every > 0 else None

    while sys_mgr.global_time < t_end:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=min_dt,
            max_dt=max_dt,
            safety_factor=safety_factor,
        )
        if next_restart_save_time is not None:
            if sys_mgr.global_time < next_restart_save_time < sys_mgr.global_time + dt:
                dt = next_restart_save_time - sys_mgr.global_time
        dt = min(dt, t_end - sys_mgr.global_time)

        sys_mgr.step(dt=dt, inner_iter=inner_iter)

        step_count += 1
        current_t = sys_mgr.global_time

        # ---------------------------------------------------------
        # 1. 出口冷却剂温度
        # ---------------------------------------------------------
        ring_exit_fluid_vol = manifold_channel.volumes[-1]
        outlet_fluid_vol = manifold_channel.volumes[-1]
        inlet_fluid_vol = inlet_boundary

        T_out = outlet_fluid_vol.T
        T_ring_exit = ring_exit_fluid_vol.T
        T_in = inlet_fluid_vol.T

        # 如果你希望入口焓严格使用入口边界温度，也可以改成 inlet_boundary.h。
        h_in = J_Inlet_HotLeg.from_vol.h
        h_in = inlet_boundary.h
        h_ring_exit = ring_exit_fluid_vol.h
        h_out = J_Manifold_Outlet.from_vol.h
        delta_h = h_in - h_out

        # 也可以额外记录按边界温度算的入口焓降
        h_in_boundary = inlet_boundary.h
        delta_h_boundary_to_out = h_in_boundary - h_out
        Q_coolant_enthalpy_drop = J_Manifold_Outlet.W * delta_h_boundary_to_out

        # ---------------------------------------------------------
        # 2. 集流环外表面散热量
        # ---------------------------------------------------------
        q_wall_nodes = get_ring_wall_radiation_by_node(solid_up_ring) + get_ring_wall_radiation_by_node(solid_down_ring)
        q_wall_total = float(np.sum(q_wall_nodes))

        # ---------------------------------------------------------
        # 3. 热管/翅片散热量与每节点总散热量
        # ---------------------------------------------------------
        hp_breakdown = get_hp_radiation_breakdown_by_node(up_ring, HP_MULTIPLIERS_UP)
        q_hp_nodes = hp_breakdown["segment_total"]
        q_hp_total = float(np.sum(q_hp_nodes))

        q_total_nodes = q_wall_nodes + q_hp_nodes
        q_total = float(np.sum(q_total_nodes))

        # ---------------------------------------------------------
        # 4. 写入一行 history
        # ---------------------------------------------------------
        row = {
            "time": current_t,
            "dt": dt,

            "W_in": J_Inlet_HotLeg.W,
            "W_ring_to_buffer": J_Manifold_Outlet.W,
            "W_out": J_Manifold_Outlet.W,

            "T_inlet_fluid": T_in,
            "T_ring_exit_fluid": T_ring_exit,
            "T_outlet_fluid": T_out,

            "h_inlet_fluid": h_in,
            "h_ring_exit_fluid": h_ring_exit,
            "h_outlet_fluid": h_out,
            "delta_h_inlet_to_outlet": delta_h,

            "h_inlet_boundary": h_in_boundary,
            "delta_h_boundary_to_outlet": delta_h_boundary_to_out,
            "Q_coolant_enthalpy_drop": Q_coolant_enthalpy_drop,

            "Q_ring_wall_rad_total": q_wall_total,
            "Q_hp_rad_total": q_hp_total,
            "Q_system_rad_total": q_total,

            "T_wall_avg": float(np.mean(solid_up_ring.T)),
            "T_fluid_avg": float(np.mean(chan_up_ring.temperature_vector)),
        }

        # 每个节点的散热量展开成 CSV 列
        for i, q in enumerate(q_wall_nodes):
            row[f"Q_wall_rad_node_{i+1:02d}"] = float(q)
            row[f"Q_ring_wall_rad_segment_{i+1:02d}"] = float(q)

        for i, q in enumerate(q_hp_nodes):
            row[f"Q_hp_rad_node_{i+1:02d}"] = float(q)
            row[f"Q_hp_rad_segment_{i+1:02d}"] = float(q)
            row[f"Q_hp_adiabatic_rad_segment_{i+1:02d}"] = float(hp_breakdown["segment_aba"][i])
            row[f"Q_hp_condenser_bare_rad_segment_{i+1:02d}"] = float(hp_breakdown["segment_con_bare"][i])
            row[f"Q_hp_fin_rad_segment_{i+1:02d}"] = float(hp_breakdown["segment_fin"][i])

        for i, q in enumerate(q_total_nodes):
            row[f"Q_total_rad_node_{i+1:02d}"] = float(q)

        for i, multiplier in enumerate(HP_MULTIPLIERS_UP):
            row[f"HP_{i+1:02d}_multiplier"] = float(multiplier)
            row[f"Q_HP_{i+1:02d}_single_adiabatic_rad"] = float(hp_breakdown["single_aba"][i])
            row[f"Q_HP_{i+1:02d}_single_condenser_bare_rad"] = float(hp_breakdown["single_con_bare"][i])
            row[f"Q_HP_{i+1:02d}_single_fin_rad"] = float(hp_breakdown["single_fin"][i])
            row[f"Q_HP_{i+1:02d}_single_total_rad"] = float(hp_breakdown["single_total"][i])

        # 可选：每个流体节点温度也存下来，方便排查沿程温降
        for i, vol in enumerate(chan_up_ring.volumes):
            row[f"T_fluid_node_{i+1:02d}"] = float(vol.T)
            row[f"P_fluid_node_{i+1:02d}"] = float(vol.P)
            row[f"h_fluid_node_{i+1:02d}"] = float(vol.h)

        for i, vol in enumerate(manifold_channel.volumes):
            row[f"T_buffer_node_{i+1:02d}"] = float(vol.T)
            row[f"P_buffer_node_{i+1:02d}"] = float(vol.P)
            row[f"h_buffer_node_{i+1:02d}"] = float(vol.h)

        history.append(row)

        if step_count % print_every == 0 or current_t >= t_end:
            print(
                f"t = {current_t:8.3f} s | "
                f"T_out = {T_out:.3f} K | "
                f"dh = {delta_h:.3f} J/kg | "
                f"Q_wall = {q_wall_total:.3f} W | "
                f"Q_hp = {q_hp_total:.3f} W | "
                f"Q_total = {q_total:.3f} W"
            )

        # ---------------------------------------------------------
        # 5. 可选：周期性保存重启动文件
        # ---------------------------------------------------------
        if (
            restart_save_path is not None
            and restart_save_every > 0.0
            and current_t >= next_restart_save_time
        ):
            checkpoint_path = restart_checkpoint_path(restart_save_path, next_restart_save_time)
            sys_mgr.save_global_state(checkpoint_path)
            print(f"Restart saved at t={current_t:.3f} s: {checkpoint_path}")
            next_restart_save_time += restart_save_every

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}_history.csv")

    write_history_csv(csv_path, history)

    if restart_save_path is not None:
        sys_mgr.save_global_state(restart_save_path)
        print(f"Final restart saved: {restart_save_path}")

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print(f"CSV: {csv_path}")
    print("=" * 70)

    return history



# ---------------------------------------------------------
# 0. 几何与工况参数
# ---------------------------------------------------------
T_space = 3.0       # K
# 0.1 ===================流道几何参数======================
# 入口与初值                                                                                                                                                             
T_INLET = 843.0                                                                                                                                                          
W_TOTAL = 2.2                                                                                                                                                            
W_BRANCH = W_TOTAL / 3.0                                                                                                                                                 
P_OUTLET = 160000.0                                                                                                                                                      
T_INIT = 863.0                                                                                                                                                           
                                                                                                                                                                        
# 1号集流环：沿用现有模型                                                                                                                                                
L_UP_RING = 0.793                                                                                                                                                        
N_UP_RING = 10                                                                                                                                                           
R_IN_RING = 0.020835                                                                                                                                                     
R_OUT_RING = 0.022835                                                                                                                                                    
DIAM_RING = 0.04167                                                                                                                                                      
AREA_RING = 0.0016065                                                                                                                                                    
PERIM_HEADER = 2.0 * np.pi * (DIAM_RING / 2.0)                                                                                                                           
HP_MULTIPLIERS_UP = [3, 2, 3, 2, 3, 2, 3, 2, 3, 3]                                                                                                                       
                                                                                                                                                                        
# 热支路：按半径重算面积                                                                                                                                                 
L_HOT_LEG = 2.19632                                                                                                                                                      
R_IN_HOT_LEG = 0.0138                                                                                                                                                    
DH_HOT_LEG = 2.0 * R_IN_HOT_LEG                                                                                                                                          
AREA_HOT_LEG = np.pi * R_IN_HOT_LEG**2   # 5.9828e-4 m^2                                                                                                                 
N_HOT_LEG = 28                                                                                                                                                           
                                                                                                                                                                        
# 2号集流环：长度不同，截面与水力直径沿用 1号集流环                                                                                                                      
L_DOWN_RING = 1.03386                                                                                                                                                    
N_DOWN_RING = 13                                                                                                                                                         
HP_MULTIPLIERS_DOWN = [2, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 2]  # sum = 31                                                                                                
                                                                                                                                                                        
# 汇流段：按半径重算面积                                                                                                                                                 
L_MANIFOLD = 0.40911                                                                                                                                                     
R_IN_MANIFOLD = 0.009                                                                                                                                                    
DH_MANIFOLD = 2.0 * R_IN_MANIFOLD
AREA_MANIFOLD = np.pi * R_IN_MANIFOLD**2   # 2.5447e-4 m^2                                                                                                               
N_MANIFOLD = 5

# 初始流量
W_UP_INIT = W_BRANCH / 2.0 / 3.0
W_DOWN_INIT = W_BRANCH / 2.0 / 3.0
W_HOT_LEG_INIT = 2.0 * (W_UP_INIT + W_DOWN_INIT)

# 0.3 热管尺寸
R_OUT_HP = 0.0085
R_IN_HP = 0.0081
R_VAPOR_HP = 0.0075
L_EVA = 0.0605
L_ABA = 0.0415
L_CON = 0.47
POROSITY = 0.966

# 0.4 翅片尺寸
THIN_FIN = 0.0003 # m
FIN_HEIGHT = 22.65e-3 # m
N_FIN_HEIGHT = 15 # 翅片高度数量

# 0.5 物性参数
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

# ---------------------------------------------------------
# 0. 进出口边界
# ---------------------------------------------------------
# 0.1 进口边界
inlet_boundary = IncompressibleBoundaryVolume(
    name="V5_Inlet_FlowBoundary",
    material=nak,
    P=P_OUTLET + 5000.0,
    T=T_INLET,
)

# 0.2 出口边界
outlet_boundary = IncompressibleBoundaryVolume(
    name="V5_Outlet_PressureBoundary",
    material=nak,
    P=P_OUTLET,
    T=T_INIT,
)
outlet_boundary.is_pressure_boundary = True

# ---------------------------------------------------------
# 1. 固体结构
# ---------------------------------------------------------
# 1.1 1号集流环固体壁面
mesh_ring_up = Mesh2D(
    x_dim=R_OUT_RING-R_IN_RING,
    n_x=1,
    y_dim=L_UP_RING,
    n_y=N_UP_RING,
    geometry_type="cylindrical",
    inner_radius=R_IN_RING,
)
solid_up_ring = HeatConduction2D(
    mesh=mesh_ring_up,
    material=SS316(),
    name="Solid_Up_Ring",
    initial_temp=T_INIT,
)
bare_area_array = solid_up_ring.boundaries['right'].area / 2
solid_up_ring.boundaries['right'].add_dynamic_radiation_condition(emissivity=0.6, bare_area_array=bare_area_array, T_env=T_space)

# 1.2 2号集流环固体壁面
mesh_ring_down = Mesh2D(
    x_dim=R_OUT_RING-R_IN_RING,
    n_x=1,
    y_dim=L_DOWN_RING,
    n_y=N_DOWN_RING,
    geometry_type="cylindrical",
    inner_radius=R_IN_RING,
)
solid_down_ring = HeatConduction2D(
    mesh=mesh_ring_down,
    material=SS316(),
    name="Solid_Down_Ring",
    initial_temp=T_INIT,
)
bare_area_array = solid_down_ring.boundaries['right'].area / 2
solid_down_ring.boundaries['right'].add_dynamic_radiation_condition(emissivity=0.6, bare_area_array=bare_area_array, T_env=T_space)

# ---------------------------------------------------------
# 2. 流体结构
# ---------------------------------------------------------
# 2.1 1号集流环流体
chan_up_ring = IncompressibleFluidChannel(
    name="Ring_Up_Channel",
    n_nodes=N_UP_RING,
    total_length=L_UP_RING,
    flow_area=AREA_RING,
    hydraulic_diam=DIAM_RING,
    initial_P=P_OUTLET,
    initial_T=T_INIT,
    material=mat_fluid,
)

# 2.2 2号集流环流体
chan_down_ring = IncompressibleFluidChannel(
    name="Ring_Down_Channel",
    n_nodes=N_DOWN_RING,
    total_length=L_DOWN_RING,
    flow_area=AREA_RING,
    hydraulic_diam=DIAM_RING,
    initial_P=P_OUTLET,
    initial_T=T_INIT,
    material=mat_fluid,
)

# 2.3 热支路
hot_leg_channel = IncompressibleFluidChannel(
    name="Hot_Leg_Channel",
    n_nodes=N_HOT_LEG,
    total_length=L_HOT_LEG,
    flow_area=AREA_HOT_LEG,
    hydraulic_diam=DH_HOT_LEG,
    initial_P=P_OUTLET,
    initial_T=T_INIT,
    material=mat_fluid,
)

# 2.4 集流环出口汇流管
manifold_channel = IncompressibleFluidChannel(
    name="Manifold_Channel",
    n_nodes=N_MANIFOLD,
    total_length=L_MANIFOLD,
    flow_area=AREA_MANIFOLD,
    hydraulic_diam=DH_MANIFOLD,
    initial_P=P_OUTLET,
    initial_T=T_INIT,
    material=mat_fluid,
)

# ---------------------------------------------------------
# 3. 接口结构
# ---------------------------------------------------------
# 3.1 入口->热支路
J_Inlet_HotLeg = InletJunction(
    name="J_Inlet_HotLeg",
    from_vol=inlet_boundary,
    to_vol=hot_leg_channel.volumes[0],
    W_initial=W_HOT_LEG_INIT,
)

# 3.2 热支路 -> 1号集流环                                                                                                                                                
# macro_vol 放在 hot_leg 侧，multiplier=2 表示“只建半环”的等效放大                                                                                                       
J_HotLeg_UpRing = MacroFlowJunction(                                                                                                                                     
    name="J_HotLeg_UpRing",                                                                                                                                              
    from_vol=hot_leg_channel.volumes[-1],                                                                                                                                
    to_vol=chan_up_ring.volumes[0],                                                                                                                                   
    macro_vol=hot_leg_channel.volumes[-1],                                                                                                                               
    multiplier=2,                                                                                                                                                        
    flow_area=AREA_HOT_LEG,                                                                                                                                              
    k_loss=0.0,                                                                                                                                                          
)

# 3.3 热支路 -> 2号集流环                                                                                                                                                
J_HotLeg_DownRing = MacroFlowJunction(                                                                                                                                   
    name="J_HotLeg_DownRing",                                                                                                                                            
    from_vol=hot_leg_channel.volumes[-1],                                                                                                                                
    to_vol=chan_down_ring.volumes[0],                                                                                                                                 
    macro_vol=hot_leg_channel.volumes[-1],                                                                                                                               
    multiplier=2,                                                                                                                                                        
    flow_area=AREA_HOT_LEG,                                                                                                                                              
    k_loss=0.0,                                                                                                                                                          
)

# 3.4 1号集流环 -> 汇流段                                                                                                                                                
# macro_vol 放在 manifold 侧，multiplier=2 表示“只建半环”的等效放大                                                                                                      
J_UpRing_Manifold = MacroFlowJunction(                                                                                                                                   
    name="J_UpRing_Manifold",                                                                                                                                            
    from_vol=chan_up_ring.volumes[-1],                                                                                                                                
    to_vol=manifold_channel.volumes[0],                                                                                                                                  
    macro_vol=manifold_channel.volumes[0],                                                                                                                               
    multiplier=2,                                                                                                                                                        
    flow_area=AREA_MANIFOLD,                                                                                                                                             
    k_loss=0.0,                                                                                                                                                          
)

# 3.5 2号集流环 -> 汇流段                                                                                                                                                
J_DownRing_Manifold = MacroFlowJunction(                                                                                                                                 
    name="J_DownRing_Manifold",                                                                                                                                          
    from_vol=chan_down_ring.volumes[-1],                                                                                                                              
    to_vol=manifold_channel.volumes[0],                                                                                                                                  
    macro_vol=manifold_channel.volumes[0],                                                                                                                               
    multiplier=2,                                                                                                                                                        
    flow_area=AREA_MANIFOLD,                                                                                                                                             
    k_loss=0.0,                                                                                                                                                          
)

# 3.6 汇流段 -> 出口边界                                                                                                                                                 
J_Manifold_Outlet = FlowJunction(                                                                                                                                        
    name="J_Manifold_Outlet",                                                                                                                                            
    from_vol=manifold_channel.volumes[-1],                                                                                                                               
    to_vol=outlet_boundary,                                                                                                                                              
    flow_area=AREA_MANIFOLD,                                                                                                                                             
    k_loss=0.0,                                                                                                                                                          
)

# ---------------------------------------------------------
# 4. 流网构建（这一步先不创建HydraulicNetwork）
# ---------------------------------------------------------
# 4.1 汇总所有控制体                                                                                                                                                     
all_vols = (                                                                                                                                                             
    [inlet_boundary, outlet_boundary]                                                                                                                                    
    + hot_leg_channel.volumes                                                                                                                                            
    + chan_up_ring.volumes                                                                                                                                            
    + chan_down_ring.volumes                                                                                                                                          
    + manifold_channel.volumes                                                                                                                                           
) 

# 4.2 汇总所有连接                                                                                                                                                       
# 先加入 section 3 中定义的“外部接口 junction”                                                                                                                           
all_juncs = [                                                                                                                                                            
    J_Inlet_HotLeg,                                                                                                                                                      
    J_HotLeg_UpRing,                                                                                                                                                     
    J_HotLeg_DownRing,                                                                                                                                                   
    J_UpRing_Manifold,                                                                                                                                                   
    J_DownRing_Manifold,                                                                                                                                                 
    J_Manifold_Outlet,                                                                                                                                                   
]  

# 再加入各流道内部的 junction                                                                                                                                            
all_juncs += hot_leg_channel.internal_junctions                                                                                                                          
all_juncs += chan_up_ring.internal_junctions                                                                                                                          
all_juncs += chan_down_ring.internal_junctions                                                                                                                        
all_juncs += manifold_channel.internal_junctions 

# ---------------------------------------------------------
# 5. 部件构建
# ---------------------------------------------------------
# 5.1 集流环内壁换热关系式
def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
        pe = np.maximum(Re * Pr, 1.0)
        return 7.0 + 0.025 * (pe ** 0.8)

# 5.2 集流环 RingHP 构建辅助函数                                                                                                                                         
def build_ring_hp(name, fluid_channel, solid_header, hp_multipliers):                                                                                                    
    return RingHP(                                                                                                                                                       
        name=name,                                                                                                                                                       
        fluid_channel=fluid_channel,                                                                                                                                     
        solid_header=solid_header,                                                                                                                                       
        hp_multipliers=hp_multipliers,                                                                                                                                   
                                                                                                                                                                        
        # 集流环流道几何：1号/2号集流环截面相同                                                                                                                          
        header_flow_area=AREA_RING,                                                                                                                                      
        header_dh=DIAM_RING,                                                                                                                                             
        header_heated_perimeter=PERIM_HEADER,                                                                                                                            
                                                                                                                                                                        
        # 热管几何                                                                                                                                                       
        hp_r_out=R_OUT_HP,                                                                                                                                               
        hp_r_in=R_IN_HP,                                                                                                                                                 
        hp_r_vapor=R_VAPOR_HP,                                                                                                                                           
        hp_L_eva=L_EVA,                                                                                                                                                  
        hp_L_con=L_CON,                                                                                                                                                  
        hp_L_aba=L_ABA,                                                                                                                                                  
                                                                                                                                                                        
        # 热管离散                                                                                                                                                       
        hp_n_eva=1,                                                                                                                                                      
        hp_n_con=12,                                                                                                                                                     
        hp_n_aba=1,                                                                                                                                                      
        hp_n_wick=1,                                                                                                                                                     
        hp_n_wall=2,                                                                                                                                                     
                                                                                                                                                                        
        # 材料与初值                                                                                                                                                     
        porosity_hp=POROSITY,                                                                                                                                            
        HP_initial_temp=800.0,                                                                                                                                           
        hp_wall_mat=mat_wall,                                                                                                                                            
        hp_fluid_mat=mat_hp_fluid,                                                                                                                                       
        hp_wick_mat=mat_wick,                                                                                                                                            
                                                                                                                                                                        
        # 翅片参数                                                                                                                                                       
        fin_thickness=THIN_FIN,                                                                                                                                          
        fin_height=FIN_HEIGHT,                                                                                                                                           
        n_fin_height=N_FIN_HEIGHT,                                                                                                                                       
        fin_wrap_ratio=(2.0 * THIN_FIN) / (2.0 * np.pi * R_OUT_HP),                                                                                                      
                                                                                                                                                                        
        # 辐射参数                                                                                                                                                       
        emissivity=0.85,                                                                                                                                                 
        up_view_factor=0.0,                                                                                                                                              
        down_view_factor=0.3,                                                                                                                                            
        T_space=T_space,                                                                                                                                                 
                                                                                                                                                                        
        # 关联式                                                                                                                                                         
        header_correlation_func=lyon_martinelli,                                                                                                                         
        hp_crossflow_base_func=lambda *args: 10.0,                                                                                                                       
        C_D=1.0,                                                                                                                                                         
                                                                                                                                                                        
        # 当前不考虑外热流                                                                                                                                               
        external_heat_config=None,                                                                                                                                       
    )

# 5.3 1号集流环                                                                                                                                                          
up_ring = build_ring_hp(                                                                                                                                                 
    name="up_ring",                                                                                                                                                      
    fluid_channel=chan_up_ring,                                                                                                                                       
    solid_header=solid_up_ring,                                                                                                                                          
    hp_multipliers=HP_MULTIPLIERS_UP,                                                                                                                                    
)

# 5.4 2号集流环                                                                                                                                                          
down_ring = build_ring_hp(                                                                                                                                               
    name="down_ring",                                                                                                                                                    
    fluid_channel=chan_down_ring,                                                                                                                                     
    solid_header=solid_down_ring,                                                                                                                                        
    hp_multipliers=HP_MULTIPLIERS_DOWN,                                                                                                                                  
)

# ---------------------------------------------------------
# 6. 系统构建
# ---------------------------------------------------------

# 6.1 创建流网                                                                                                                                                           
network = HydraulicNetwork(                                                                                                                                              
    all_vols,                                                                                                                                                            
    all_juncs,                                                                                                                                                           
    gravity_vector=0.0                                                                                                                                                   
)                                                                                                                                                                        
                                                                                                                                                                        
# 6.2 创建系统                                                                                                                                                           
sys_mgr = SystemManager(fluid_network=network)                                                                                                                           
                                                                                                                                                                        
# 6.3 添加组件                                                                                                                                                           
sys_mgr.add_component(up_ring)                                                                                                                                           
sys_mgr.add_component(down_ring)
sys_mgr.initialize_system()

# 6.5 基本检查输出                                                                                                                                                       
print(f"Volumes: {len(all_vols)}")
print(f"Junctions: {len(all_juncs)}")                                                                                                                                    
print(f"Registered solids: {len(sys_mgr.solid_components)}")                                                                                                             
print(f"Registered couplers: {len(sys_mgr.couplers)}")                                                                                                                   
                                                                                                                                                                        
print(f"Up ring outlet k_loss: {up_ring.outlet_k_loss}")                                                                                                                 
print(f"Down ring outlet k_loss: {down_ring.outlet_k_loss}")                                                                                                             
                                                                                                                                                                        
print(f"Inlet W: {J_Inlet_HotLeg.W}")                                                                                                                                    
print(f"HotLeg -> UpRing W: {J_HotLeg_UpRing.W}")                                                                                                                        
print(f"HotLeg -> DownRing W: {J_HotLeg_DownRing.W}")                                                                                                                    
print(f"UpRing -> Manifold W: {J_UpRing_Manifold.W}")                                                                                                                    
print(f"DownRing -> Manifold W: {J_DownRing_Manifold.W}")                                                                                                                
print(f"Manifold -> Outlet W: {J_Manifold_Outlet.W}")                                                                                                                    
                                                                                                                                                                        
print(f"Hot leg inlet T: {hot_leg_channel.volumes[0].T}")                                                                                                                
print(f"Up ring inlet T: {chan_up_ring.volumes[0].T}")                                                                                                                
print(f"Down ring inlet T: {chan_down_ring.volumes[0].T}")                                                                                                            
print(f"Manifold outlet T: {manifold_channel.volumes[-1].T}")

# ---------------------------------------------------------
# 7. 计算函数
# ---------------------------------------------------------
def get_hp_rejection_by_node(ring_hp, hp_multipliers):
    """
    返回每个流体节点对应的热管/翅片总散热量 [W/node]。
    返回长度严格等于 len(hp_multipliers)，可处理某些节点 multiplier=0 的情况。
    """
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


def get_ring_total_rejection_by_node(solid_ring, ring_hp, hp_multipliers):
    """
    每个节点总散热量 = 集流环外壁辐射 + 热管/翅片辐射 [W/node]
    """
    q_wall = get_ring_wall_radiation_by_node(solid_ring)
    q_hp = get_hp_rejection_by_node(ring_hp, hp_multipliers)

    if len(q_wall) != len(q_hp):
        raise ValueError(
            f"Node mismatch: wall radiation has {len(q_wall)} nodes, "
            f"HP rejection has {len(q_hp)} nodes."
        )

    return q_wall + q_hp


def run_case(
    case_name="V5_dual_ring_branch",
    t_end=50.0,
    min_dt=1e-3,
    max_dt=0.5,
    safety_factor=2.0,
    inner_iter=2,
    print_every=10,
    csv_path=None,
    restart_from=None,
    restart_save_path=None,
    restart_save_every=0.0,
):
    print("=" * 70)
    print(f"Running case: {case_name}")
    print("=" * 70)

    if restart_from is not None:
        sys_mgr.load_global_state(restart_from)
        inlet_boundary.set_boundary_state(T=T_INLET)

        inlet_idx = network.vol_to_idx[inlet_boundary]
        network.T_vec[inlet_idx] = inlet_boundary.T
        network.h_vec[inlet_idx] = inlet_boundary.h
        network.P_vec[inlet_idx] = inlet_boundary.P

        print(f"Restart loaded from: {restart_from}")
        print(f"Restart time: {sys_mgr.global_time:.6f} s")
    else:
        sys_mgr.initialize_system()
        print("System initialized from initial condition.")

    history = []
    step_count = 0
    next_restart_save_time = (
        sys_mgr.global_time + restart_save_every if restart_save_every > 0.0 else None
    )

    while sys_mgr.global_time < t_end:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=min_dt,
            max_dt=max_dt,
            safety_factor=safety_factor,
        )

        if next_restart_save_time is not None:
            if sys_mgr.global_time < next_restart_save_time < sys_mgr.global_time + dt:
                dt = next_restart_save_time - sys_mgr.global_time

        dt = min(dt, t_end - sys_mgr.global_time)
        sys_mgr.step(dt=dt, inner_iter=inner_iter)

        step_count += 1
        current_t = sys_mgr.global_time

        q_up_wall_nodes = get_ring_wall_radiation_by_node(solid_up_ring)
        q_down_wall_nodes = get_ring_wall_radiation_by_node(solid_down_ring)
        q_up_hp_nodes = get_hp_rejection_by_node(up_ring, HP_MULTIPLIERS_UP)
        q_down_hp_nodes = get_hp_rejection_by_node(down_ring, HP_MULTIPLIERS_DOWN)

        q_up_wall_total = float(np.sum(q_up_wall_nodes))
        q_down_wall_total = float(np.sum(q_down_wall_nodes))
        q_up_hp_total = float(np.sum(q_up_hp_nodes))
        q_down_hp_total = float(np.sum(q_down_hp_nodes))
        q_total = q_up_wall_total + q_down_wall_total + q_up_hp_total + q_down_hp_total

        w_hotleg_to_up_macro = 2.0 * J_HotLeg_UpRing.W
        w_hotleg_to_down_macro = 2.0 * J_HotLeg_DownRing.W
        w_up_to_manifold_macro = 2.0 * J_UpRing_Manifold.W
        w_down_to_manifold_macro = 2.0 * J_DownRing_Manifold.W

        row = {
            "time": current_t,
            "dt": dt,
            "W_inlet": J_Inlet_HotLeg.W,
            "W_hotleg_to_up_macro": w_hotleg_to_up_macro,
            "W_hotleg_to_down_macro": w_hotleg_to_down_macro,
            "W_up_to_manifold_macro": w_up_to_manifold_macro,
            "W_down_to_manifold_macro": w_down_to_manifold_macro,
            "W_outlet": J_Manifold_Outlet.W,
            "T_inlet": inlet_boundary.T,
            "T_hotleg_out": hot_leg_channel.volumes[-1].T,
            "T_up_ring_out": chan_up_ring.volumes[-1].T,
            "T_down_ring_out": chan_down_ring.volumes[-1].T,
            "T_manifold_out": manifold_channel.volumes[-1].T,
            "P_hotleg_out": hot_leg_channel.volumes[-1].P,
            "P_up_ring_out": chan_up_ring.volumes[-1].P,
            "P_down_ring_out": chan_down_ring.volumes[-1].P,
            "P_manifold_out": manifold_channel.volumes[-1].P,
            "Q_up_wall_rad_total": q_up_wall_total,
            "Q_up_hp_rad_total": q_up_hp_total,
            "Q_down_wall_rad_total": q_down_wall_total,
            "Q_down_hp_rad_total": q_down_hp_total,
            "Q_system_rad_total": q_total,
            "T_hotleg_avg": float(np.mean(hot_leg_channel.temperature_vector)),
            "T_up_ring_avg": float(np.mean(chan_up_ring.temperature_vector)),
            "T_down_ring_avg": float(np.mean(chan_down_ring.temperature_vector)),
            "T_manifold_avg": float(np.mean(manifold_channel.temperature_vector)),
            "T_up_wall_avg": float(np.mean(solid_up_ring.T)),
            "T_down_wall_avg": float(np.mean(solid_down_ring.T)),
        }
        history.append(row)

        if step_count % print_every == 0 or current_t >= t_end:
            print(
                f"t = {current_t:8.3f} s | "
                f"T_out = {manifold_channel.volumes[-1].T:.3f} K | "
                f"W_up = {w_hotleg_to_up_macro:.4f} kg/s | "
                f"W_down = {w_hotleg_to_down_macro:.4f} kg/s | "
                f"Q_total = {q_total:.3f} W"
            )

        if (
            restart_save_path is not None
            and restart_save_every > 0.0
            and current_t >= next_restart_save_time
        ):
            checkpoint_path = restart_checkpoint_path(
                restart_save_path,
                next_restart_save_time,
            )
            sys_mgr.save_global_state(checkpoint_path)
            print(f"Restart saved at t={current_t:.3f} s: {checkpoint_path}")
            next_restart_save_time += restart_save_every

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}_history.csv")

    write_history_csv(csv_path, history)

    if restart_save_path is not None:
        sys_mgr.save_global_state(restart_save_path)
        print(f"Final restart saved: {restart_save_path}")

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print(f"CSV: {csv_path}")
    print("=" * 70)

    return history


if __name__ == "__main__":
    history = run_case(
        case_name="V5_dual_ring_branch",
        t_end=1000.0,
        csv_path=os.path.join(current_dir, "V5_dual_ring_branch_history.csv"),
        restart_save_path=os.path.join(current_dir, "V5_dual_ring_branch_restart.npz"),
        restart_save_every=200.0,
        restart_from=None,
    )
