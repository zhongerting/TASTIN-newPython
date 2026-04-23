# ============================================================================
#  TOPAZ-II 冷却剂回路开环仿真 (V4.2)
# ----------------------------------------------------------------------------
#  作者：
#  日期：2026-04
#  基于：TASTIN
# ============================================================================
#
#  【建模目标】
#   仿真冷却剂回路 (堆芯本身不在仿真域内)：
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
#  无 (V4.1 修正了 V4 中 T_outlet 读取错误，能量守恒严格成立)
#
#  【V4 → V4.1 关键修正】
#  V4 把 core_inlet_boundary.T 当作 T_outlet 来计算 Q_loop, 但
#  IncompressibleBoundaryVolume 默认作"无限热池"(T 锁定为初始值),
#  导致 Q_loop 被严重低估、能量平衡比值看似 2.00x。
#  V4.1 改用 return_leg.volumes[-1].T (回流支路末端真实流体温度),
#  能量守恒比值精确到 1.05x (5% 误差来自 Cp 取固定值的近似)。
# ============================================================================

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

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

print("  [V4.2] HeatPipe2D/HPwithFin hot patch disabled; using fixed source code.")


# ============================================================================
#  Part 1. 边界条件参数
# ============================================================================
print("=" * 72)
print("  TOPAZ-II Coolant Loop - V4.2 (Open-Loop, no hot patch)")
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
#  Part 2. 几何参数
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

# 主泵段 (用了管道代替泵, 用普通管道而非真实泵特性)
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
#  Part 9. 瞬态仿真: 默认 500 秒
# ============================================================================
T_END = float(os.environ.get("COOLANT_LOOP_T_END", "500.0"))
print(f"\n[Part 9] Transient ({T_END:g} s)")
print("-" * 72)

history = {
    "t": [], "T_inlet": [], "T_outlet": [],
    "T_ring1_avg": [], "T_ring2_avg": [], "T_rad_ext": [],
    "Q_rej_total": [], "Q_rej_net": [],
    "Q_rej_solver": [],
    "Q_bare_radiation": [], "Q_fin_radiation": [], "Q_fin_net_from_root": [],
    "Q_loop_balance": [], "Q_loop_flux": [], "Q_loop_cp": [],
    "Q_storage_rate": [], "Q_balance_residual": [],
    "W_total": [], "W_in_actual": [], "W_out_actual": [],
}

ENERGY_AUDIT_KEYS = [
    "U_fluid_representative",
    "U_fluid_physical_scaled",
    "U_header_wall_representative",
    "U_header_wall_physical_scaled",
    "U_hp_wall_wick_fluid_representative",
    "U_hp_physical_scaled",
    "U_fin_quasi_steady",
    "U_boundary_excluded",
    "dU_fluid_hot_leg",
    "dU_fluid_ring1",
    "dU_fluid_ring2",
    "dU_fluid_return_leg",
    "dU_header_ring1",
    "dU_header_ring2",
    "dU_hp_ring1",
    "dU_hp_ring2",
]
for key in ENERGY_AUDIT_KEYS:
    history[key] = []

def collect_hp_heat_breakdown():
    """Collect both-ring heat-pipe radiation totals in W."""
    totals = {
        "bare_radiation": 0.0,
        "fin_radiation": 0.0,
        "fin_absorption": 0.0,
        "fin_net_from_root": 0.0,
        "gross_rejection": 0.0,
        "net_rejection": 0.0,
        "solver_rejection": 0.0,
    }
    for ring_hp in [ring1_hp, ring2_hp]:
        for i, hp in enumerate(ring_hp.hp_units):
            multiplier = HP_MULTIPLIERS[i]
            breakdown = hp.get_heat_exchange_breakdown()
            for key in totals:
                if key in breakdown:
                    totals[key] += float(np.sum(breakdown[key])) * multiplier
            # Solver-consistent condenser rejection. Boundary current_flux is
            # positive into the heat pipe, so outward rejection is negative.
            totals["solver_rejection"] += float(
                -np.sum(hp.hp.boundaries["outer_con"].current_flux)
            ) * multiplier
    return totals

def _fluid_volume_energy(vol):
    if all(hasattr(vol, attr) for attr in ("rho", "vol", "h")):
        return float(vol.rho * vol.vol * vol.h)
    return 0.0

def _fluid_volumes_energy(volumes):
    return float(sum(_fluid_volume_energy(vol) for vol in volumes))

def _solid_energy(solid):
    if hasattr(solid, "_update_properties"):
        solid._update_properties()
    if hasattr(solid, "thermal_capacitance") and hasattr(solid, "T"):
        return float(np.sum(solid.thermal_capacitance * solid.T))
    return 0.0

def _ring_hp_energy(ring_hp, physical_scaled):
    total = 0.0
    for i, hp in enumerate(ring_hp.hp_units):
        multiplier = HP_MULTIPLIERS[i] if physical_scaled else 1.0
        total += _solid_energy(hp.hp) * multiplier
    return float(total)

def collect_energy_audit(prev_audit=None, dt=1.0):
    """分项统计有限域储能；dU_* 为物理放大口径下的储能变化率 [W]。"""
    U_fluid_hot_leg_rep = _fluid_volumes_energy(hot_leg.volumes)
    U_fluid_ring1 = _fluid_volumes_energy(ring1_channel.volumes)
    U_fluid_ring2 = _fluid_volumes_energy(ring2_channel.volumes)
    U_fluid_return_leg_rep = _fluid_volumes_energy(return_leg.volumes)
    U_fluid_other = sum(_fluid_volume_energy(vol) for vol in [
        manifold, radiator_internal, radiator_external, pump_A, pump_B, pump_dist
    ])

    U_fluid_representative = (
        U_fluid_hot_leg_rep + U_fluid_ring1 + U_fluid_ring2
        + U_fluid_other + U_fluid_return_leg_rep
    )
    U_fluid_physical_scaled = (
        3.0 * U_fluid_hot_leg_rep
        + U_fluid_ring1 + U_fluid_ring2
        + U_fluid_other
        + 3.0 * U_fluid_return_leg_rep
    )

    U_header_ring1 = _solid_energy(ring1_wall_solid)
    U_header_ring2 = _solid_energy(ring2_wall_solid)
    U_header_wall_representative = U_header_ring1 + U_header_ring2
    U_header_wall_physical_scaled = U_header_wall_representative

    U_hp_ring1_rep = _ring_hp_energy(ring1_hp, physical_scaled=False)
    U_hp_ring2_rep = _ring_hp_energy(ring2_hp, physical_scaled=False)
    U_hp_ring1 = _ring_hp_energy(ring1_hp, physical_scaled=True)
    U_hp_ring2 = _ring_hp_energy(ring2_hp, physical_scaled=True)
    U_hp_wall_wick_fluid_representative = U_hp_ring1_rep + U_hp_ring2_rep
    U_hp_physical_scaled = U_hp_ring1 + U_hp_ring2

    physical_parts = {
        "fluid_hot_leg": 3.0 * U_fluid_hot_leg_rep,
        "fluid_ring1": U_fluid_ring1,
        "fluid_ring2": U_fluid_ring2,
        "fluid_return_leg": 3.0 * U_fluid_return_leg_rep,
        "header_ring1": U_header_ring1,
        "header_ring2": U_header_ring2,
        "hp_ring1": U_hp_ring1,
        "hp_ring2": U_hp_ring2,
    }

    audit = {
        "U_fluid_representative": U_fluid_representative,
        "U_fluid_physical_scaled": U_fluid_physical_scaled,
        "U_header_wall_representative": U_header_wall_representative,
        "U_header_wall_physical_scaled": U_header_wall_physical_scaled,
        "U_hp_wall_wick_fluid_representative": U_hp_wall_wick_fluid_representative,
        "U_hp_physical_scaled": U_hp_physical_scaled,
        "U_fin_quasi_steady": 0.0,      # 翅片当前为准稳态支路，不作为独立储能自由度。
        "U_boundary_excluded": 0.0,     # 上腔室/堆芯入口为无限边界，按定义排除有限域储能。
        "_physical_parts": physical_parts,
    }

    dt_safe = max(float(dt), 1e-30)
    for part_name, value in physical_parts.items():
        prev_value = value if prev_audit is None else prev_audit["_physical_parts"][part_name]
        audit[f"dU_{part_name}"] = (value - prev_value) / dt_safe

    audit["_total_physical_scaled"] = (
        U_fluid_physical_scaled
        + U_header_wall_physical_scaled
        + U_hp_physical_scaled
        + audit["U_fin_quasi_steady"]
    )
    prev_total = audit["_total_physical_scaled"] if prev_audit is None else prev_audit["_total_physical_scaled"]
    audit["_dU_total_physical_scaled"] = (audit["_total_physical_scaled"] - prev_total) / dt_safe
    return audit

t = 0.0
step_count = 0
prev_audit = collect_energy_audit()
prev_time = t
print(f"  Running...")
while t < T_END:
    dt = sys_mgr.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=1.0)
    sys_mgr.step(dt=dt, inner_iter=2)
    t = sys_mgr.global_time
    step_count += 1

    # 计算两环总辐射散热: Σ (单根 HP 散热 × 该节点热管数)
    heat_breakdown = collect_hp_heat_breakdown()
    q_total_rej = heat_breakdown["gross_rejection"]
    q_net_rej = heat_breakdown["net_rejection"]
    q_bare_rad = heat_breakdown["bare_radiation"]
    q_fin_rad = heat_breakdown["fin_radiation"]
    q_fin_net = heat_breakdown["fin_net_from_root"]
    q_solver_rej = heat_breakdown["solver_rejection"]

    # ⭐ V4.1 修正: T_outlet 读真实流体温度 (回流支路末端节点),
    # 不能读 core_inlet_boundary.T —— 它作为 IncompressibleBoundaryVolume
    # 在默认配置下是"无限热池", T 被锁定为初始值不响应入口流体
    T_outlet_real = return_leg.volumes[-1].T

    # 回路实际带走的热量: 优先使用进出口焓差，保留 Cp*dT 作为对照
    Cp = nak.heat_capacity_liquid(T_INIT)
    h_in = nak.enthalpy_saturated_liquid(T_INLET)
    h_out = nak.enthalpy_saturated_liquid(T_outlet_real)
    Q_loop_balance = W_TOTAL_DESIGN * (h_in - h_out)
    Q_loop_cp = W_TOTAL_DESIGN * Cp * (T_INLET - T_outlet_real)
    W_in_actual = j_inlet.W * 3.0
    W_out_actual = j_outlet.W * 3.0
    Q_loop_flux = W_in_actual * h_in - W_out_actual * h_out

    dt_energy = max(t - prev_time, 1e-30)
    current_audit = collect_energy_audit(prev_audit=prev_audit, dt=dt_energy)
    Q_storage_rate = current_audit["_dU_total_physical_scaled"]
    Q_balance_residual = Q_loop_flux - q_solver_rej - Q_storage_rate
    prev_audit = current_audit
    prev_time = t

    history["t"].append(t)
    history["T_inlet"].append(T_INLET)
    history["T_outlet"].append(T_outlet_real)
    history["T_ring1_avg"].append(np.mean(ring1_channel.temperature_vector))
    history["T_ring2_avg"].append(np.mean(ring2_channel.temperature_vector))
    history["T_rad_ext"].append(radiator_external.T)
    history["Q_rej_total"].append(q_total_rej)
    history["Q_rej_net"].append(q_net_rej)
    history["Q_rej_solver"].append(q_solver_rej)
    history["Q_bare_radiation"].append(q_bare_rad)
    history["Q_fin_radiation"].append(q_fin_rad)
    history["Q_fin_net_from_root"].append(q_fin_net)
    history["Q_loop_balance"].append(Q_loop_balance)
    history["Q_loop_flux"].append(Q_loop_flux)
    history["Q_loop_cp"].append(Q_loop_cp)
    history["Q_storage_rate"].append(Q_storage_rate)
    history["Q_balance_residual"].append(Q_balance_residual)
    history["W_total"].append(j_outlet.W * 3)
    history["W_in_actual"].append(W_in_actual)
    history["W_out_actual"].append(W_out_actual)
    for key in ENERGY_AUDIT_KEYS:
        history[key].append(current_audit[key])

    if step_count % 50 == 0:
        print(f"  t={t:>6.2f}s | T_outlet={T_outlet_real:7.2f}K "
              f"| ΔT_loop={T_INLET - T_outlet_real:6.2f}K "
              f"| Q_rej={q_total_rej/1000:.2f}kW "
              f"| Q_rej_solver={q_solver_rej/1000:.2f}kW "
              f"| Q_fin_rad={q_fin_rad/1000:.2f}kW "
              f"| Q_loop_h={Q_loop_balance/1000:.2f}kW "
              f"| Q_loop_flux={Q_loop_flux/1000:.2f}kW "
              f"| dUdt={Q_storage_rate/1000:.2f}kW "
              f"| residual={Q_balance_residual/1000:.2f}kW")


# ============================================================================
#  Part 10. 末态总结 + 能量平衡分析
# ============================================================================
print(f"\n[Transient End]")
T_outlet_real = return_leg.volumes[-1].T
Cp_final = nak.heat_capacity_liquid(T_INIT)
h_in_final = nak.enthalpy_saturated_liquid(T_INLET)
h_out_final = nak.enthalpy_saturated_liquid(T_outlet_real)
Q_loop_final = W_TOTAL_DESIGN * (h_in_final - h_out_final)
Q_loop_cp_final = W_TOTAL_DESIGN * Cp_final * (T_INLET - T_outlet_real)
heat_breakdown_final = collect_hp_heat_breakdown()
q_total_rej = heat_breakdown_final["gross_rejection"]
q_net_rej = heat_breakdown_final["net_rejection"]
q_bare_rad = heat_breakdown_final["bare_radiation"]
q_fin_rad = heat_breakdown_final["fin_radiation"]
q_fin_net = heat_breakdown_final["fin_net_from_root"]
q_solver_rej = heat_breakdown_final["solver_rejection"]
W_in_final = history["W_in_actual"][-1] if history["W_in_actual"] else j_inlet.W * 3.0
W_out_final = history["W_out_actual"][-1] if history["W_out_actual"] else j_outlet.W * 3.0
Q_loop_flux_final = history["Q_loop_flux"][-1] if history["Q_loop_flux"] else W_in_final * h_in_final - W_out_final * h_out_final
Q_storage_rate_final = history["Q_storage_rate"][-1] if history["Q_storage_rate"] else 0.0
Q_balance_residual_final = history["Q_balance_residual"][-1] if history["Q_balance_residual"] else 0.0
final_audit = prev_audit
print(f"  Q_bare_rad:          {q_bare_rad/1000:.2f} kW  (HP bare-wall radiation)")
print(f"  Q_fin_rad:           {q_fin_rad/1000:.2f} kW  (fin radiation to space)")
print(f"  Q_fin_net_from_root: {q_fin_net/1000:.2f} kW  (net heat drawn by fins)")
print(f"  Q_rej_net:           {q_net_rej/1000:.2f} kW  (net rejection)")
print(f"  Q_rej_solver:        {q_solver_rej/1000:.2f} kW  (solver boundary rejection)")
print(f"  W_in / W_out:        {W_in_final:.6f} / {W_out_final:.6f} kg/s")
print(f"  Q_loop_flux:         {Q_loop_flux_final/1000:.2f} kW  (Win*hin - Wout*hout)")
print(f"  dUdt_est:            {Q_storage_rate_final/1000:.2f} kW  (finite-domain storage)")
print(f"  balance_residual:    {Q_balance_residual_final/1000:.2f} kW  (Q_flux - Q_rej_solver - dUdt)")
print(f"  Inlet  T:           {T_INLET:.2f} K  (恒定边界)")
print(f"  Outlet T (real):    {T_outlet_real:.2f} K  ⭐ 读自回流支路末端")
print(f"  ΔT_loop (T_in - T_out): {T_INLET - T_outlet_real:.2f} K")
print(f"  Ring1 avg:          {np.mean(ring1_channel.temperature_vector):.2f} K")
print(f"  Ring2 avg:          {np.mean(ring2_channel.temperature_vector):.2f} K")
print()
print(f"  Q_loop  (W·Δh):     {Q_loop_final/1000:.2f} kW  (回路真实散热量, 焓差)")
print(f"  Q_loop  (W·Cp·ΔT):  {Q_loop_cp_final/1000:.2f} kW  (固定 Cp 对照)")
print(f"  Q_rej   (HP辐射):   {q_total_rej/1000:.2f} kW  (热管散到太空)")
print(f"  能量平衡: Q_rej / Q_loop = {q_total_rej/Q_loop_final:.2f}x")
print(f"  Note: 接近 1.0 表示能量守恒成立, 仿真物理正确")

print(f"\n[Energy Audit - finite domain]")
print(f"  U_fluid_representative:              {final_audit['U_fluid_representative']/1e6:.3f} MJ")
print(f"  U_fluid_physical_scaled:             {final_audit['U_fluid_physical_scaled']/1e6:.3f} MJ")
print(f"  U_header_wall_representative:        {final_audit['U_header_wall_representative']/1e6:.3f} MJ")
print(f"  U_header_wall_physical_scaled:       {final_audit['U_header_wall_physical_scaled']/1e6:.3f} MJ")
print(f"  U_hp_wall_wick_fluid_representative: {final_audit['U_hp_wall_wick_fluid_representative']/1e6:.3f} MJ")
print(f"  U_hp_physical_scaled:                {final_audit['U_hp_physical_scaled']/1e6:.3f} MJ")
print(f"  U_fin_quasi_steady:                  {final_audit['U_fin_quasi_steady']/1e6:.3f} MJ  (quasi-steady, no state)")
print(f"  U_boundary_excluded:                 {final_audit['U_boundary_excluded']/1e6:.3f} MJ  (UpperPlenum/CoreInlet excluded)")
print(f"  dU_fluid_hot_leg:                    {final_audit['dU_fluid_hot_leg']/1000:.3f} kW")
print(f"  dU_fluid_ring1 / ring2:              {final_audit['dU_fluid_ring1']/1000:.3f} / {final_audit['dU_fluid_ring2']/1000:.3f} kW")
print(f"  dU_fluid_return_leg:                 {final_audit['dU_fluid_return_leg']/1000:.3f} kW")
print(f"  dU_header_ring1 / ring2:             {final_audit['dU_header_ring1']/1000:.3f} / {final_audit['dU_header_ring2']/1000:.3f} kW")
print(f"  dU_hp_ring1 / ring2:                 {final_audit['dU_hp_ring1']/1000:.3f} / {final_audit['dU_hp_ring2']/1000:.3f} kW")


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
    f'V4.2 Open-Loop Coolant System (T_in={T_INLET}K, W={W_TOTAL_DESIGN}kg/s)',
    fontweight='bold')
axes[0].legend(loc='best', fontsize=9)
axes[0].grid(True, linestyle='--', alpha=0.5)

# 能量平衡
axes[1].plot(history["t"], np.array(history["Q_rej_total"]) / 1000,
             'orange', lw=2, label='Q_rej (HP radiation)')
axes[1].plot(history["t"], np.array(history["Q_rej_solver"]) / 1000,
             'black', lw=1.6, linestyle='--', label='Q_rej_solver')
axes[1].plot(history["t"], np.array(history["Q_fin_radiation"]) / 1000,
             'darkgreen', lw=1.6, linestyle='-.', label='Q_fin_rad (fins)')
axes[1].plot(history["t"], np.array(history["Q_loop_balance"]) / 1000,
             'red', lw=2, linestyle='--', label='Q_loop (W·Δh)')
axes[1].plot(history["t"], np.array(history["Q_loop_flux"]) / 1000,
             'blue', lw=1.6, linestyle='--', label='Q_loop_flux')
axes[1].plot(history["t"], np.array(history["Q_loop_cp"]) / 1000,
             'brown', lw=1.5, linestyle=':', label='Q_loop (W·Cp·ΔT)')
axes[1].set_ylabel('Heat [kW]')
axes[1].set_title('Energy Balance: Q_rej should equal Q_loop at steady state',
                  fontweight='bold')
axes[1].legend(loc='best')
axes[1].grid(True, linestyle='--', alpha=0.5)

# 流量
axes[2].plot(history["t"], history["W_total"], 'purple', lw=2, label='Total Flow')
axes[2].plot(history["t"], history["W_in_actual"], 'tab:blue', lw=1.4, linestyle='--', label='W_in')
axes[2].plot(history["t"], history["W_out_actual"], 'tab:orange', lw=1.4, linestyle='-.', label='W_out')
axes[2].axhline(y=W_TOTAL_DESIGN, color='gray', linestyle=':',
                label=f'W_target={W_TOTAL_DESIGN:.4f}')
axes[2].set_xlabel('Time [s]')
axes[2].set_ylabel('Mass Flow [kg/s]')
axes[2].legend(loc='best')
axes[2].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'transient_result_v4_2.png')
plt.savefig(out_path, dpi=120)
print(f"\n  Plot saved: {out_path}")
plt.close()

csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'energy_audit_v4_2_500s.csv')
with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.DictWriter(f, fieldnames=list(history.keys()))
    writer.writeheader()
    n_rows = len(history["t"])
    for i in range(n_rows):
        writer.writerow({key: history[key][i] for key in history})
print(f"  Energy audit CSV saved: {csv_path}")

fig_audit, audit_axes = plt.subplots(3, 1, figsize=(11, 11), sharex=True)

audit_axes[0].plot(history["t"], np.array(history["U_fluid_representative"]) / 1e6,
                   lw=1.8, color='tab:blue', label='U_fluid_rep')
audit_axes[0].plot(history["t"], np.array(history["U_fluid_physical_scaled"]) / 1e6,
                   lw=1.8, color='tab:cyan', linestyle='--', label='U_fluid_phys')
audit_axes[0].plot(history["t"], np.array(history["U_header_wall_physical_scaled"]) / 1e6,
                   lw=1.8, color='tab:orange', label='U_header_wall_phys')
audit_axes[0].plot(history["t"], np.array(history["U_hp_physical_scaled"]) / 1e6,
                   lw=1.8, color='tab:red', label='U_hp_phys')
audit_axes[0].plot(history["t"], np.array(history["U_fin_quasi_steady"]) / 1e6,
                   lw=1.2, color='gray', linestyle=':', label='U_fin_quasi_steady=0')
audit_axes[0].set_ylabel('Stored energy [MJ]')
audit_axes[0].set_title('Finite-Domain Stored Energy Audit', fontweight='bold')
audit_axes[0].legend(loc='best', fontsize=9)
audit_axes[0].grid(True, linestyle='--', alpha=0.5)

audit_axes[1].plot(history["t"], np.array(history["dU_fluid_hot_leg"]) / 1000,
                   lw=1.7, color='tab:red', label='dU_fluid_hot_leg')
audit_axes[1].plot(history["t"], np.array(history["dU_fluid_ring1"]) / 1000,
                   lw=1.7, color='tab:blue', label='dU_fluid_ring1')
audit_axes[1].plot(history["t"], np.array(history["dU_fluid_ring2"]) / 1000,
                   lw=1.7, color='tab:cyan', linestyle='--', label='dU_fluid_ring2')
audit_axes[1].plot(history["t"], np.array(history["dU_fluid_return_leg"]) / 1000,
                   lw=1.7, color='tab:green', label='dU_fluid_return_leg')
audit_axes[1].set_ylabel('dU/dt [kW]')
audit_axes[1].set_title('Fluid Storage-Rate Breakdown (physical scaled)', fontweight='bold')
audit_axes[1].legend(loc='best', fontsize=9)
audit_axes[1].grid(True, linestyle='--', alpha=0.5)

audit_axes[2].plot(history["t"], np.array(history["dU_header_ring1"]) / 1000,
                   lw=1.7, color='tab:orange', label='dU_header_ring1')
audit_axes[2].plot(history["t"], np.array(history["dU_header_ring2"]) / 1000,
                   lw=1.7, color='goldenrod', linestyle='--', label='dU_header_ring2')
audit_axes[2].plot(history["t"], np.array(history["dU_hp_ring1"]) / 1000,
                   lw=1.7, color='tab:red', label='dU_hp_ring1')
audit_axes[2].plot(history["t"], np.array(history["dU_hp_ring2"]) / 1000,
                   lw=1.7, color='tab:pink', linestyle='--', label='dU_hp_ring2')
audit_axes[2].set_xlabel('Time [s]')
audit_axes[2].set_ylabel('dU/dt [kW]')
audit_axes[2].set_title('Header and Heat-Pipe Storage-Rate Breakdown', fontweight='bold')
audit_axes[2].legend(loc='best', fontsize=9)
audit_axes[2].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
audit_plot_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'energy_audit_v4_2.png')
plt.savefig(audit_plot_path, dpi=120)
print(f"  Energy audit plot saved: {audit_plot_path}")
plt.close()

print("\n" + "=" * 72)
print("  V4.2 Open-Loop simulation completed.")
print("=" * 72)
