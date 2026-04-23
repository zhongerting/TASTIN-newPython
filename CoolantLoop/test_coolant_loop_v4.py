# ============================================================================
#  TOPAZ-II 冷却剂回路开环仿真 (V4)
# ----------------------------------------------------------------------------
#  日期：2026-04
#  基于：TASTIN-python 框架
# ============================================================================
#
#  【建模目标】
#      [入口边界 W,T] → 上腔室 → 热支路 → 集流环×2 → 汇流 → 辐射器
#                                                       → 主泵段 → 回流支路
#                                                       → [出口边界 P]
#
#  【边界条件设计】
#      • 入口边界: 流量 → 用 InletJunction 强制 W, 入口虚拟容器给定 T
#      • 出口边界: 压力 → 用 IncompressibleBoundaryVolume 锚定 P
#      堆芯被简化为一对边界条件: "送出 968K, 0.5445 kg/s 的流体"
#                                "回收时压力维持 161 kPa"
#
#  【主泵处理】
#  用管道代替泵，不实现泵特性曲线。流量由入口边界强制给定。
#
#  【拓扑结构】
#      上腔室 (T入口边界)
#          │
#          ▼
#      热支路 ×3 (用 multiplier=3 的 MacroFlowJunction 简化为单股)
#          │
#          ├──→ 集流环1 (Ring1, 6节点 + 78根代表辐射管)
#          │       │
#          ├──→ 集流环2 (Ring2, 6节点 + 78根代表辐射管)
#                  │
#                  ▼
#              汇流总管 → 辐射器内段 → 辐射器外段
#                  │
#                  ▼
#              泵段 A → 泵段 B → 泵分流段
#                  │
#                  ▼
#              回流支路 ×3 (multiplier=3 简化)
#                  │
#                  ▼
#              堆芯入口段 (P出口边界)
#
#  【热管/辐射管说明】
#  TOPAZ-II 真实辐射管是 NaK-78 单相强迫对流; 本模型借用 RingHP 组件
#  (含 Na 工质相变热管 + 翅片 + 太空辐射), 几何与散热能力等效。
#  每集流环上 78 根辐射管 (符合 TOPAZ-II 设计)。
#
#  【已知问题】
#  Q_rej (HP 算的辐射散热) 比理论值偏高 2 倍，已定位为 HPwithFin 源码 bug
#  (冷凝段裸壁辐射与翅片辐射双重计算)，详见 bug_report.md。
#  本仿真的流体温度场和回路能量平衡 (Q_loop) 不受此 bug 影响。
# ============================================================================

import os
import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

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
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WickMaterial import WickMaterial
from Components.RingHP import RingHP

logging.basicConfig(level=logging.WARNING)


# ============================================================================
#  HOT PATCH 区域: 修复 HeatPipe2D / HPwithFin 源码 bug (运行时修补，不动源码)
#  详见 bug_report.md 的 Bug 1, 2, 3
# ============================================================================
from Components.HPwithFin import HeatPipe2D as _HeatPipe2D_cls

_orig_hp2d_init = _HeatPipe2D_cls.__init__
_hp2d_counter = [0]   # 用列表包装实现闭包内可写

def _patched_hp2d_init(self, *args, **kwargs):
    """修复: HeatPipe2D 的 self.name 未赋值 bug + HPwithFin 不传 name 导致重名 bug"""
    name = kwargs.get('name', 'Unnamed_Solid')
    if name == 'Unnamed_Solid' or not name:
        _hp2d_counter[0] += 1
        name = f"HeatPipe2D_auto_{_hp2d_counter[0]}"
        kwargs['name'] = name
    self.name = name
    _orig_hp2d_init(self, *args, **kwargs)

def _patched_compute_fluxes(self, t):
    """修复: _compute_fluxes 用旧 R_*_inner 但父类已重构为 G_*_inner"""
    Q_net_2d = np.zeros(self.shape_nodes)
    T_2d = self.T.reshape(self.shape_nodes)

    flux_x = (T_2d[:-1, :] - T_2d[1:, :]) * self.G_x_inner
    Q_net_2d[:-1, :] -= flux_x
    Q_net_2d[1:, :]  += flux_x

    flux_y = (T_2d[:, :-1] - T_2d[:, 1:]) * self.G_y_inner
    Q_net_2d[:, :-1] -= flux_y
    Q_net_2d[:, 1:]  += flux_y

    if 'left' in self.boundaries:
        Q_net_2d[0, :] += self.boundaries['left'].compute_net_flux_for_solver()
    if 'bottom' in self.boundaries:
        Q_net_2d[:, 0] += self.boundaries['bottom'].compute_net_flux_for_solver()
    if 'top' in self.boundaries:
        Q_net_2d[:, -1] += self.boundaries['top'].compute_net_flux_for_solver()

    idx_eva = self.n_eva
    idx_aba = self.n_eva + self.n_aba
    idx_con = self.mesh.n_y
    if 'outer_eva' in self.boundaries:
        Q_net_2d[-1, 0:idx_eva] += self.boundaries['outer_eva'].compute_net_flux_for_solver()
    if 'outer_aba' in self.boundaries:
        Q_net_2d[-1, idx_eva:idx_aba] += self.boundaries['outer_aba'].compute_net_flux_for_solver()
    if 'outer_con' in self.boundaries:
        Q_net_2d[-1, idx_aba:idx_con] += self.boundaries['outer_con'].compute_net_flux_for_solver()

    return Q_net_2d.flatten()

_HeatPipe2D_cls.__init__ = _patched_hp2d_init
_HeatPipe2D_cls._compute_fluxes = _patched_compute_fluxes


# ============================================================================
#  Part 1. 边界条件参数
# ============================================================================
print("=" * 72)
print("  TOPAZ-II Coolant Loop - V4 (Open-Loop Final)")
print("=" * 72)

# --- 入口边界 (上腔室 = 堆芯出来的流体) ---
T_INLET = 968.0           # 入口温度 [K] - TOPAZ-II 堆芯出口设计温度
W_INLET_SINGLE = 0.1815   # 单股入口流量 [kg/s] (3 股并联)
W_TOTAL_DESIGN = 3 * W_INLET_SINGLE   # 总流量 0.5445 kg/s

# --- 出口边界 (堆芯入口段 = 流体即将进堆芯) ---
P_OUTLET = 1.61e5         # 出口压力 [Pa] - TOPAZ-II 一回路工作压力

# --- 系统初始温度 ---
T_INIT = 863.0            # 整个回路初始温度 [K]


# ============================================================================
#  Part 2. 几何参数 (TOPAZ-II 实测值或设计值)
# ============================================================================
# --- 热支路 (堆芯顶部到上腔室外的连接管) ---
L_HOT_LEG = 0.40          # 长度 [m]
N_HOT_LEG = 6             # 离散节点数
D_PIPE = 0.02             # 管内径 [m]
A_PIPE = np.pi * (D_PIPE / 2) ** 2   # 流通截面 [m²]

# --- 回流支路 (主泵到堆芯入口的连接管) ---
L_RET_LEG = 0.35
N_RET_LEG = 3

# --- 集流环 (上/下集流环结构相同) ---
L_RING = 0.793            # 集流环长度 [m]
N_RING = 6                # 集流环离散节点数
A_RING = 0.0016065        # 环形流道截面 [m²]
DH_RING = 0.04167         # 水力直径 [m]
WALL_THICKNESS_RING = 0.002    # 管壁厚 [m]
R_IN_RING = DH_RING / 2.0
R_OUT_RING = R_IN_RING + WALL_THICKNESS_RING
PERIM_RING = 2.0 * np.pi * R_IN_RING

# --- 辐射管阵列 (TOPAZ-II 设计: 每环 78 根) ---
HP_MULTIPLIERS = [13, 13, 13, 13, 13, 13]   # 每节点 13 根, 6 节点共 78 根
HP_R_OUT, HP_R_IN, HP_R_VAPOR = 0.0085, 0.0081, 0.0075   # 辐射管半径 [m]
HP_L_EVA, HP_L_CON = 0.06, 0.482            # 蒸发段/冷凝段长度 [m]
HP_N_EVA, HP_N_CON = 1, 12                  # 离散节点数
HP_N_WICK, HP_N_WALL = 1, 1
HP_POROSITY = 0.5
HP_INIT_TEMP = 800.0      # 辐射管初始温度 [K]

# --- 翅片 + 辐射散热参数 ---
FIN_THICKNESS = 0.0003                      # 翅片厚 [m]
FIN_HEIGHT = 22.65e-3                       # 翅片高 [m]
N_FIN_HEIGHT = 15
FIN_WRAP_RATIO = (2 * FIN_THICKNESS) / (2.0 * np.pi * HP_R_OUT)
EMISSIVITY = 0.93                           # 翅片表面发射率
UP_VF, DOWN_VF = 1.0, 0.675                 # 上下表面对深空角系数
T_SPACE = 3.0                               # 深空温度 [K]

# --- 汇流段 / 辐射器总管 / 泵段 ---
L_MANIFOLD = 0.30
L_RAD_INT = 0.50          # 辐射器内段长度 [m]
L_RAD_EXT = 0.60          # 辐射器外段长度 [m]
D_MANIFOLD = 0.025
A_MANIFOLD = np.pi * (D_MANIFOLD / 2) ** 2

# --- 局部阻力系数 ---
K_T_JUNCTION = 0.5        # 三通分流/合流
K_ELBOW = 0.3             # 弯头
K_CONTRACTION = 0.2       # 收缩
K_EXPANSION = 0.2         # 扩张


# ============================================================================
#  Part 3. 物性
# ============================================================================
# 主冷却剂: 用 Sodium 代用 NaK-78
nak = Sodium()

# 辐射管工质 (Na 相变模型)
hp_fluid = SodiumHP(name="HP_Fluid_Na")

# 管壁/集流环外壁: SS316 不锈钢
mat_wall = SS316(name="SS316_Wall")

# 辐射管吸液芯 (复合材料: SS316 骨架 + Na 工质)
mat_wick = WickMaterial(
    name="HP_Wick_Composite",
    solid_mat=SS316(), fluid_mat=hp_fluid,
    porosity=HP_POROSITY, r_vapor=HP_R_VAPOR, r_in_wall=HP_R_IN,
)


# ============================================================================
#  Part 4. 关联式
# ============================================================================
def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    """液态金属强迫对流 Nu 关联式 (Lyon-Martinelli)"""
    Pe = np.maximum(Re * Pr, 1.0)
    return 7.0 + 0.025 * (Pe ** 0.8)


# ============================================================================
#  Part 5. 流网组件构建
# ============================================================================
all_vols = []
all_juncs = []

# ----------------------------------------------------------------------------
#  5.1 ⭐ 入口边界: 上腔室 (堆芯出口流体)
#      物理: "堆芯送出 0.5445 kg/s, 968K 的流体"
#      实现: 用 IncompressibleBoundaryVolume 作为虚拟温度源
#            真正的流量约束由后续的 InletJunction 完成
# ----------------------------------------------------------------------------
upper_plenum = IncompressibleBoundaryVolume(
    name="UpperPlenum_Inlet",
    material=nak,
    P=P_OUTLET + 5000,   # 给个略高初始压力, 避免初始反流
    T=T_INLET,           # ⭐ 入口温度边界
)
# 注意: 不开 is_pressure_boundary, 压力锚在出口端
all_vols.append(upper_plenum)


# ----------------------------------------------------------------------------
#  5.2 热支路 (代表性单股, 用 multiplier=3 等效 3 股并联)
# ----------------------------------------------------------------------------
hot_leg = IncompressibleFluidChannel(
    name="HotLeg_Representative",
    n_nodes=N_HOT_LEG, total_length=L_HOT_LEG,
    flow_area=A_PIPE, hydraulic_diam=D_PIPE,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.extend(hot_leg.volumes)
all_juncs.extend(hot_leg.internal_junctions)


# ----------------------------------------------------------------------------
#  5.3 双集流环 (Ring1, Ring2)
#      每环组成: 流道 + 外壁固体 + RingHP复合组件 (热管+翅片+辐射)
# ----------------------------------------------------------------------------
# 5.3.1 集流环流道
ring1_channel = IncompressibleFluidChannel(
    name="Ring1_Channel", n_nodes=N_RING, total_length=L_RING,
    flow_area=A_RING, hydraulic_diam=DH_RING,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.extend(ring1_channel.volumes)
all_juncs.extend(ring1_channel.internal_junctions)

ring2_channel = IncompressibleFluidChannel(
    name="Ring2_Channel", n_nodes=N_RING, total_length=L_RING,
    flow_area=A_RING, hydraulic_diam=DH_RING,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.extend(ring2_channel.volumes)
all_juncs.extend(ring2_channel.internal_junctions)

# 5.3.2 集流环外壁固体 (SS316, 2D 圆柱网格)
def build_ring_wall_solid(name):
    """构建集流环的外壁固体导热模型 (外表面绝热, 散热全靠插在流道里的辐射管)"""
    mesh = Mesh2D(
        x_dim=WALL_THICKNESS_RING, n_x=1, y_dim=L_RING, n_y=N_RING,
        geometry_type='cylindrical', inner_radius=R_IN_RING,
    )
    solid = HeatConduction2D(mesh=mesh, material=mat_wall, initial_temp=T_INIT)
    solid.name = name
    solid.boundaries['right'].add_flux_condition(q_flux=0.0)
    solid.boundaries['top'].add_flux_condition(q_flux=0.0)
    solid.boundaries['bottom'].add_flux_condition(q_flux=0.0)
    return solid

ring1_wall_solid = build_ring_wall_solid("Ring1_Wall")
ring2_wall_solid = build_ring_wall_solid("Ring2_Wall")

# 5.3.3 RingHP 复合组件 (核心散热部件)
def build_ring_hp(name, channel, wall_solid):
    """组装 RingHP: 集流环 + 78根辐射管 + 翅片 + 太空辐射边界"""
    return RingHP(
        name=name, fluid_channel=channel, solid_header=wall_solid,
        hp_multipliers=HP_MULTIPLIERS,
        header_flow_area=A_RING, header_dh=DH_RING, header_heated_perimeter=PERIM_RING,
        hp_r_out=HP_R_OUT, hp_r_in=HP_R_IN, hp_r_vapor=HP_R_VAPOR,
        hp_L_eva=HP_L_EVA, hp_L_con=HP_L_CON,
        hp_n_eva=HP_N_EVA, hp_n_con=HP_N_CON,
        hp_n_wick=HP_N_WICK, hp_n_wall=HP_N_WALL,
        porosity_hp=HP_POROSITY, HP_initial_temp=HP_INIT_TEMP,
        fin_thickness=FIN_THICKNESS, fin_height=FIN_HEIGHT, n_fin_height=N_FIN_HEIGHT,
        fin_wrap_ratio=FIN_WRAP_RATIO, emissivity=EMISSIVITY,
        up_view_factor=UP_VF, down_view_factor=DOWN_VF, T_space=T_SPACE,
        hp_wall_mat=mat_wall, hp_fluid_mat=hp_fluid, hp_wick_mat=mat_wick,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
    )

ring1_hp = build_ring_hp("Ring1_HP", ring1_channel, ring1_wall_solid)
ring2_hp = build_ring_hp("Ring2_HP", ring2_channel, ring2_wall_solid)


# ----------------------------------------------------------------------------
#  5.4 汇流段 / 辐射器总管 / 泵段
# ----------------------------------------------------------------------------
manifold = IncompressibleFluidVolume(
    name="Manifold", volume=A_MANIFOLD * L_MANIFOLD, length=L_MANIFOLD,
    flow_area=A_MANIFOLD, hydraulic_diam=D_MANIFOLD,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.append(manifold)

radiator_internal = IncompressibleFluidVolume(
    name="Rad_Int", volume=A_MANIFOLD * L_RAD_INT, length=L_RAD_INT,
    flow_area=A_MANIFOLD, hydraulic_diam=D_MANIFOLD,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.append(radiator_internal)

radiator_external = IncompressibleFluidVolume(
    name="Rad_Ext", volume=A_MANIFOLD * L_RAD_EXT, length=L_RAD_EXT,
    flow_area=A_MANIFOLD, hydraulic_diam=D_MANIFOLD,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.append(radiator_external)

# 主泵段（"用管道代替泵"的指示, 用普通管道而非真实泵特性)
pump_A = IncompressibleFluidVolume(
    name="Pump_A", volume=A_MANIFOLD * 0.05, length=0.05,
    flow_area=A_MANIFOLD, hydraulic_diam=D_MANIFOLD,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.append(pump_A)

pump_B = IncompressibleFluidVolume(
    name="Pump_B", volume=A_MANIFOLD * 0.05, length=0.05,
    flow_area=A_MANIFOLD, hydraulic_diam=D_MANIFOLD,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.append(pump_B)

pump_dist = IncompressibleFluidVolume(
    name="Pump_Dist", volume=A_MANIFOLD * 0.10, length=0.10,
    flow_area=A_MANIFOLD, hydraulic_diam=D_MANIFOLD,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.append(pump_dist)


# ----------------------------------------------------------------------------
#  5.5 回流支路
# ----------------------------------------------------------------------------
return_leg = IncompressibleFluidChannel(
    name="ReturnLeg_Representative",
    n_nodes=N_RET_LEG, total_length=L_RET_LEG,
    flow_area=A_PIPE, hydraulic_diam=D_PIPE,
    initial_P=P_OUTLET, initial_T=T_INIT, material=nak,
)
all_vols.extend(return_leg.volumes)
all_juncs.extend(return_leg.internal_junctions)


# ----------------------------------------------------------------------------
#  5.6 ⭐ 出口边界: 堆芯入口段 (回路出口压力锚)
#      物理: "流体即将进入堆芯, 此处压力维持 161 kPa"
#      实现: IncompressibleBoundaryVolume + is_pressure_boundary=True
# ----------------------------------------------------------------------------
core_inlet_boundary = IncompressibleBoundaryVolume(
    name="CoreInlet_Outlet",
    material=nak, P=P_OUTLET, T=T_INIT,
)
core_inlet_boundary.is_pressure_boundary = True   # ⭐ 必须开启 (TASTIN 默认 False)
all_vols.append(core_inlet_boundary)


# ============================================================================
#  Part 6. 拓扑连接 (按拓扑图依次连接各段)
# ============================================================================

# 6.1 ⭐ 入口: 上腔室 → 热支路 (InletJunction 强制流量)
j_inlet = InletJunction(
    name="J_UpperPlenum2HotLeg",
    from_vol=upper_plenum, to_vol=hot_leg.volumes[0],
    W_initial=W_INLET_SINGLE,
)
all_juncs.append(j_inlet)

# 6.2 热支路 → 集流环1 / 集流环2 (并联分流, multiplier=3 把单股放大为 3 股)
all_juncs.append(MacroFlowJunction(
    name="J_HotLeg2Ring1",
    from_vol=hot_leg.volumes[-1], to_vol=ring1_channel.volumes[0],
    macro_vol=ring1_channel.volumes[0],
    multiplier=3, flow_area=A_PIPE, k_loss=K_T_JUNCTION,
))
all_juncs.append(MacroFlowJunction(
    name="J_HotLeg2Ring2",
    from_vol=hot_leg.volumes[-1], to_vol=ring2_channel.volumes[0],
    macro_vol=ring2_channel.volumes[0],
    multiplier=3, flow_area=A_PIPE, k_loss=K_T_JUNCTION,
))

# 6.3 集流环 → 汇流总管 (含辐射管阻塞产生的压损)
all_juncs.append(MacroFlowJunction(
    name="J_Ring1_2_Manifold",
    from_vol=ring1_channel.volumes[N_RING // 2], to_vol=manifold,
    macro_vol=ring1_channel.volumes[N_RING // 2],
    multiplier=3, flow_area=A_MANIFOLD,
    k_loss=K_T_JUNCTION + ring1_hp.outlet_k_loss,   # 叠加辐射管阵列阻力
))
all_juncs.append(MacroFlowJunction(
    name="J_Ring2_2_Manifold",
    from_vol=ring2_channel.volumes[N_RING // 2], to_vol=manifold,
    macro_vol=ring2_channel.volumes[N_RING // 2],
    multiplier=3, flow_area=A_MANIFOLD,
    k_loss=K_T_JUNCTION + ring2_hp.outlet_k_loss,
))

# 6.4 汇流 → 辐射器内段 → 辐射器外段
all_juncs.append(MacroFlowJunction(
    name="J_Manifold2RadInt",
    from_vol=manifold, to_vol=radiator_internal, macro_vol=radiator_internal,
    multiplier=3, flow_area=A_MANIFOLD, k_loss=K_T_JUNCTION,
))
all_juncs.append(FlowJunction(
    name="J_RadInt2Ext", from_vol=radiator_internal, to_vol=radiator_external,
    flow_area=A_MANIFOLD, k_loss=K_ELBOW,
))

# 6.5 辐射器外段 → 主泵段 → 泵分流段
all_juncs.append(FlowJunction(
    name="J_RadExt2PumpA", from_vol=radiator_external, to_vol=pump_A,
    flow_area=A_MANIFOLD, k_loss=K_CONTRACTION,
))
all_juncs.append(FlowJunction(
    name="J_PumpA2PumpB", from_vol=pump_A, to_vol=pump_B,
    flow_area=A_MANIFOLD, k_loss=0.1,
))
all_juncs.append(FlowJunction(
    name="J_PumpB2Dist", from_vol=pump_B, to_vol=pump_dist,
    flow_area=A_MANIFOLD, k_loss=K_EXPANSION,
))

# 6.6 泵分流段 → 回流支路
all_juncs.append(FlowJunction(
    name="J_PumpDist2RetLeg",
    from_vol=pump_dist, to_vol=return_leg.volumes[0],
    flow_area=A_PIPE, k_loss=K_CONTRACTION,
))

# 6.7 ⭐ 出口: 回流支路 → 堆芯入口段 (压力出口边界)
j_outlet = MacroFlowJunction(
    name="J_RetLeg2CoreInlet",
    from_vol=return_leg.volumes[-1], to_vol=core_inlet_boundary,
    macro_vol=return_leg.volumes[-1],
    multiplier=3, flow_area=A_PIPE, k_loss=K_T_JUNCTION,
)
all_juncs.append(j_outlet)


# ============================================================================
#  Part 7. 求解器组装 + 初始化
# ============================================================================
net = HydraulicNetwork(volumes=all_vols, junctions=all_juncs, gravity_vector=0.0)
sys_mgr = SystemManager(fluid_network=net)
sys_mgr.add_component(ring1_hp)   # 自动注册 RingHP 内部所有固体和耦合器
sys_mgr.add_component(ring2_hp)

print(f"  Volumes: {len(all_vols)}, Junctions: {len(all_juncs)}")
print(f"  Solids: {len(sys_mgr.solid_components)}, Couplers: {len(sys_mgr.couplers)}")
print(f"  Inlet boundary:  W = {W_TOTAL_DESIGN} kg/s, T = {T_INLET} K")
print(f"  Outlet boundary: P = {P_OUTLET/1000:.1f} kPa")

sys_mgr.initialize_system()
print("  System initialized\n")


# ============================================================================
#  Part 8. 稳态结果输出 (t=0)
# ============================================================================
print("[Part 8] Steady-State Baseline (t=0)")
print("-" * 72)
W_in = j_inlet.W * 3
W_out = j_outlet.W * 3
print(f"Mass Flow:")
print(f"  Inlet  (×3): {W_in:>8.4f} kg/s")
print(f"  Outlet (×3): {W_out:>8.4f} kg/s")
print(f"  Mass balance error: "
      f"{abs(W_in - W_out) / W_TOTAL_DESIGN * 100:>5.2f} %")

print(f"\nPressure key nodes:")
print(f"  Upper Plenum (in):  {upper_plenum.P:>10.1f} Pa")
print(f"  Hot Leg in:         {hot_leg.volumes[0].P:>10.1f} Pa")
print(f"  Ring1 in:           {ring1_channel.volumes[0].P:>10.1f} Pa")
print(f"  Ring1 out:          {ring1_channel.volumes[-1].P:>10.1f} Pa")
print(f"  Manifold:           {manifold.P:>10.1f} Pa")
print(f"  Pump_Dist:          {pump_dist.P:>10.1f} Pa")
print(f"  Core Inlet (out):   {core_inlet_boundary.P:>10.1f} Pa  (⭐ P边界锚)")

print(f"\nTemperature distribution:")
print(f"  Upper Plenum:  {upper_plenum.T:>7.2f} K  (⭐ T边界设定)")
print(f"  Hot Leg avg:   {np.mean(hot_leg.temperature_vector):>7.2f} K")
print(f"  Ring1 avg:     {np.mean(ring1_channel.temperature_vector):>7.2f} K")
print(f"  Ring2 avg:     {np.mean(ring2_channel.temperature_vector):>7.2f} K")
print(f"  Rad Ext:       {radiator_external.T:>7.2f} K")
print(f"  Core Inlet:    {core_inlet_boundary.T:>7.2f} K  (回路出口温度)")


# ============================================================================
#  Part 9. 瞬态仿真: 60 秒
# ============================================================================
T_END = 60.0
print(f"\n[Part 9] Transient (60 s)")
print("-" * 72)

history = {
    "t": [], "T_inlet": [], "T_outlet": [],
    "T_ring1_avg": [], "T_ring2_avg": [], "T_rad_ext": [],
    "Q_rej_total": [], "Q_loop_balance": [], "W_total": [],
}

t = 0.0
step_count = 0
print(f"  Running...")
while t < T_END:
    dt = sys_mgr.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=1.0)
    sys_mgr.step(dt=dt, inner_iter=2)
    t = sys_mgr.global_time
    step_count += 1

    # 计算两环总辐射散热: Σ (单根 HP 散热 × 该节点热管数)
    q_total_rej = 0.0
    for ring_hp in [ring1_hp, ring2_hp]:
        for i, hp in enumerate(ring_hp.hp_units):
            _, q_con = hp.get_heat_rejection_distribution()
            q_total_rej += np.sum(q_con) * HP_MULTIPLIERS[i]

    # 回路实际带走的热量 Q_loop = W × Cp × (T_in - T_out)
    Cp = nak.heat_capacity_liquid(T_INIT)
    Q_loop_balance = W_TOTAL_DESIGN * Cp * (T_INLET - core_inlet_boundary.T)

    history["t"].append(t)
    history["T_inlet"].append(T_INLET)
    history["T_outlet"].append(core_inlet_boundary.T)
    history["T_ring1_avg"].append(np.mean(ring1_channel.temperature_vector))
    history["T_ring2_avg"].append(np.mean(ring2_channel.temperature_vector))
    history["T_rad_ext"].append(radiator_external.T)
    history["Q_rej_total"].append(q_total_rej)
    history["Q_loop_balance"].append(Q_loop_balance)
    history["W_total"].append(j_outlet.W * 3)

    if step_count % 50 == 0:
        print(f"  t={t:>6.2f}s | T_outlet={core_inlet_boundary.T:7.2f}K "
              f"| ΔT_loop={T_INLET - core_inlet_boundary.T:6.2f}K "
              f"| Q_rej={q_total_rej/1000:.2f}kW "
              f"| Q_loop={Q_loop_balance/1000:.2f}kW")


# ============================================================================
#  Part 10. 末态总结 + 能量平衡分析
# ============================================================================
print(f"\n[Transient End]")
print(f"  Inlet  T:           {T_INLET:.2f} K  (恒定边界)")
print(f"  Outlet T:           {core_inlet_boundary.T:.2f} K")
print(f"  ΔT_loop (T_in - T_out): {T_INLET - core_inlet_boundary.T:.2f} K")
print(f"  Ring1 avg:          {np.mean(ring1_channel.temperature_vector):.2f} K")
print(f"  Ring2 avg:          {np.mean(ring2_channel.temperature_vector):.2f} K")
print()
print(f"  Q_loop  (W·Cp·ΔT):  {Q_loop_balance/1000:.2f} kW  (回路实际散热量)")
print(f"  Q_rej   (HP辐射):   {q_total_rej/1000:.2f} kW  (热管散到太空的热量)")
print(f"  能量平衡: Q_rej / Q_loop = {q_total_rej/Q_loop_balance:.2f}x")
print(f"  Note: 比值 ≈ 1.0 表示稳态; 比值偏高是 RingHP 源码 bug, 见 bug_report.md")
# 在 Part 10 末尾加这一行:
print(f"  Return Leg 末端 T (真实流体): {return_leg.volumes[-1].T:.2f} K")
print(f"  Hot Leg 入口 T (真实流体):   {hot_leg.volumes[0].T:.2f} K")


# ============================================================================
#  Part 11. 结果绘图
# ============================================================================
fig, axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

# 温度
axes[0].plot(history["t"], history["T_inlet"], 'r-', lw=2, label='Inlet (T edge)')
axes[0].plot(history["t"], history["T_outlet"], 'm-', lw=2, label='Outlet (T loop end)')
axes[0].plot(history["t"], history["T_ring1_avg"], 'b-', lw=1.5, label='Ring1 avg')
axes[0].plot(history["t"], history["T_ring2_avg"], 'c--', lw=1.5, label='Ring2 avg')
axes[0].plot(history["t"], history["T_rad_ext"], 'g-', lw=1.5, label='Rad Ext')
axes[0].axhline(y=T_INIT, color='gray', linestyle=':', alpha=0.5, label=f'T_init={T_INIT}K')
axes[0].set_ylabel('Temperature [K]')
axes[0].set_title(
    f'V4 Open-Loop Coolant System (T_in={T_INLET}K, W={W_TOTAL_DESIGN}kg/s)',
    fontweight='bold')
axes[0].legend(loc='best', fontsize=9)
axes[0].grid(True, linestyle='--', alpha=0.5)

# 能量平衡
axes[1].plot(history["t"], np.array(history["Q_rej_total"]) / 1000,
             'orange', lw=2, label='Q_rej (HP radiation)')
axes[1].plot(history["t"], np.array(history["Q_loop_balance"]) / 1000,
             'red', lw=2, linestyle='--', label='Q_loop (W·Cp·ΔT)')
axes[1].set_ylabel('Heat [kW]')
axes[1].set_title('Energy Balance: Q_rej should equal Q_loop at steady state',
                  fontweight='bold')
axes[1].legend(loc='best')
axes[1].grid(True, linestyle='--', alpha=0.5)

# 流量
axes[2].plot(history["t"], history["W_total"], 'purple', lw=2, label='Total Flow')
axes[2].axhline(y=W_TOTAL_DESIGN, color='gray', linestyle=':',
                label=f'W_target={W_TOTAL_DESIGN:.4f}')
axes[2].set_xlabel('Time [s]')
axes[2].set_ylabel('Mass Flow [kg/s]')
axes[2].legend(loc='best')
axes[2].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'transient_result_v4.png')
plt.savefig(out_path, dpi=120)
print(f"\n  Plot saved: {out_path}")
plt.close()

print("\n" + "=" * 72)
print("  V4 Open-Loop simulation completed.")
print("=" * 72)
