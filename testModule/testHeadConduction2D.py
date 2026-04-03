import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import logging

from MathSolvers.solver_module import NuclearODESolver
from Materials.Solids.MoNb import MoNb

from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Boundary import BoundaryRegion
from Solvers.Couplers import SolidSolidCouple2D

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("TestBasic")


# ==========================================
# 辅助类：耦合系统管理器
# ==========================================
class CoupledSystem:
    """
    用于管理两个导热对象及其耦合关系的封装类，
    提供给 ODE Solver 一个统一的接口。
    """

    def __init__(self, obj1, obj2, coupler):
        self.obj1 = obj1
        self.obj2 = obj2
        self.coupler = coupler

        # 记录切片位置，用于拆分状态向量
        self.n1 = self.obj1.N
        self.n2 = self.obj2.N
        self.total_N = self.n1 + self.n2

    def get_combined_derivatives(self, t, y_combined):
        """
        计算组合系统的导数。
        关键点：必须在计算 Flux 之前手动同步边界条件 (Sync)。
        """
        # 1. 拆分状态向量
        y1 = y_combined[:self.n1]
        y2 = y_combined[self.n1:]

        # 2. 更新两个对象的内部状态 (T, 物性, 内部热阻, 内部边界状态)
        #    注意：我们要手动展开 HeatConduction.get_derivatives 的前半部分

        # --- Object 1 ---
        self.obj1.T[:] = y1
        self.obj1._update_properties()
        self.obj1._compute_internal_resistance()
        self.obj1._update_boundaries_state()  # 更新 BoundaryRegion 的 T_adj_node 和 R_internal

        # --- Object 2 ---
        self.obj2.T[:] = y2
        self.obj2._update_properties()
        self.obj2._compute_internal_resistance()
        self.obj2._update_boundaries_state()

        # 3. 执行耦合同步 (Sync)
        #    利用刚刚更新的内部状态，交换边界信息 (T_ext, R_ext)
        self.coupler.sync()

        # 4. 计算热流和源项 (Flux & Source)
        #    注意：_compute_fluxes 会使用 Sync 步骤设定好的 T_ext
        Q_net1 = self.obj1._compute_fluxes(t)
        self.obj1._update_sources(t)
        dTdt1 = (Q_net1 + self.obj1.Q_source) / self.obj1.thermal_capacitance

        Q_net2 = self.obj2._compute_fluxes(t)
        self.obj2._update_sources(t)
        dTdt2 = (Q_net2 + self.obj2.Q_source) / self.obj2.thermal_capacitance

        # 5. 合并导数
        return np.concatenate((dTdt1, dTdt2))


# ==========================================
# 通用运行与绘图函数
# ==========================================

def run_comparison(case_name, setup_unified_func, setup_coupled_func, t_end=1000.0):
    """
    运行对比测试：Unified vs Coupled
    """
    logger.info(f"========== Running Case: {case_name} ==========")
    mat = MoNb()
    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)
    t_eval = np.linspace(0, t_end, 51)  # 50 steps

    # --- 1. 运行整体模型 (Unified) ---
    logger.info("  -> Running Unified Model...")
    phy_u = setup_unified_func(mat)
    y0_u = phy_u.T.copy()

    res_u = solver.solve(phy_u.get_derivatives, (0, t_end), y0_u, t_eval=t_eval)

    if not res_u['success']:
        logger.error("Unified simulation failed!")
        return

    # --- 2. 运行耦合模型 (Coupled) ---
    logger.info("  -> Running Coupled Model...")
    phy_a, phy_b, coupler = setup_coupled_func(mat)
    coupled_sys = CoupledSystem(phy_a, phy_b, coupler)
    y0_c = np.concatenate((phy_a.T, phy_b.T))

    res_c = solver.solve(coupled_sys.get_combined_derivatives, (0, t_end), y0_c, t_eval=t_eval)

    if not res_c['success']:
        logger.error("Coupled simulation failed!")
        return

    # --- 3. 数据后处理与对比 ---
    # 为了对比，我们需要提取特定位置的温度。
    # 简单起见，我们对比最终时刻的 全场平均温度 和 最高/最低温度。

    T_final_u = res_u['y'][:, -1]
    T_final_c = res_c['y'][:, -1]

    # 将耦合结果拼接以便对比 (注意顺序)
    # 对于 3.1 (左右) 和 3.2 (上下)，物理空间拼接顺序不同，但在 1D 数组中
    # Unified 的 mesh 生成顺序通常是行主序或列主序。
    # 我们这里通过统计量对比，并绘制沿特征线的分布。

    avg_u = np.mean(T_final_u)
    avg_c = np.mean(T_final_c)
    diff_avg = abs(avg_u - avg_c)

    logger.info(f"  [Result] Unified Avg T: {avg_u:.4f} K")
    logger.info(f"  [Result] Coupled Avg T: {avg_c:.4f} K")
    logger.info(f"  [Result] Abs Difference: {diff_avg:.6f} K")

    plot_comparison(case_name, res_u, res_c, phy_u, phy_a, phy_b)


def plot_comparison(case_name, res_u, res_c, mesh_u, mesh_a, mesh_b):
    """绘制对比图"""
    # 提取最终时刻
    T_u = res_u['y'][:, -1]
    T_c = res_c['y'][:, -1]
    T_a = T_c[:mesh_a.N]
    T_b = T_c[mesh_a.N:]

    plt.figure(figsize=(10, 5))

    # 根据案例类型选择绘图方式
    if "3.1" in case_name:
        # 笛卡尔左右拼接：取中间高度 (y=H/2) 沿 X 轴的分布
        # Mesh2D flatten order: i*ny + j (x index * ny + y index)
        # Unified
        ny_u = mesh_u.mesh.n_y
        mid_j = ny_u // 2
        # x indices: 0 to nx-1
        indices_u = [i * ny_u + mid_j for i in range(mesh_u.mesh.n_x)]
        x_u = mesh_u.mesh.x_centers

        # Coupled A (Left)
        ny_a = mesh_a.mesh.n_y
        mid_j_a = ny_a // 2
        indices_a = [i * ny_a + mid_j_a for i in range(mesh_a.mesh.n_x)]
        x_a = mesh_a.mesh.x_centers

        # Coupled B (Right) -> x 坐标要加上 A 的宽度
        ny_b = mesh_b.mesh.n_y
        mid_j_b = ny_b // 2
        indices_b = [i * ny_b + mid_j_b for i in range(mesh_b.mesh.n_x)]
        x_b = mesh_b.mesh.x_centers + 0.05  # A width

        plt.plot(x_u, T_u[indices_u], 'k-', lw=3, label='Unified Model')
        plt.plot(x_a, T_a[indices_a], 'r--', lw=2, marker='o', label='Coupled Part A')
        plt.plot(x_b, T_b[indices_b], 'b--', lw=2, marker='x', label='Coupled Part B')
        plt.xlabel("X (m)")

    elif "3.2" in case_name:
        # 柱坐标上下拼接：取中间半径 (mid R) 沿 Z 轴 (Y dim) 的分布
        # Unified
        nx_u = mesh_u.mesh.n_x
        mid_i = nx_u // 2
        # y indices: 0 to ny-1. Index = mid_i * ny + j
        ny_u = mesh_u.mesh.n_y
        indices_u = [mid_i * ny_u + j for j in range(ny_u)]
        z_u = mesh_u.mesh.y_centers

        # Coupled A (Bottom)
        nx_a = mesh_a.mesh.n_x
        mid_i_a = nx_a // 2
        ny_a = mesh_a.mesh.n_y
        indices_a = [mid_i_a * ny_a + j for j in range(ny_a)]
        z_a = mesh_a.mesh.y_centers

        # Coupled B (Top) -> Z 坐标加 A 的高度
        nx_b = mesh_b.mesh.n_x
        mid_i_b = nx_b // 2
        ny_b = mesh_b.mesh.n_y
        indices_b = [mid_i_b * ny_b + j for j in range(ny_b)]
        z_b = mesh_b.mesh.y_centers + 0.1  # A height

        plt.plot(z_u, T_u[indices_u], 'k-', lw=3, label='Unified Model')
        plt.plot(z_a, T_a[indices_a], 'r--', lw=2, marker='o', label='Coupled Part A (Bottom)')
        plt.plot(z_b, T_b[indices_b], 'b--', lw=2, marker='x', label='Coupled Part B (Top)')
        plt.xlabel("Z (Axial) (m)")

    elif "3.3" in case_name:
        # 柱坐标内外拼接：取中间高度 (mid Z) 沿 R 轴 (X dim) 的分布
        # Unified
        ny_u = mesh_u.mesh.n_y
        mid_j = ny_u // 2
        indices_u = [i * ny_u + mid_j for i in range(mesh_u.mesh.n_x)]
        r_u = mesh_u.mesh.x_centers

        # Coupled A (Inner)
        ny_a = mesh_a.mesh.n_y
        mid_j_a = ny_a // 2
        indices_a = [i * ny_a + mid_j_a for i in range(mesh_a.mesh.n_x)]
        r_a = mesh_a.mesh.x_centers

        # Coupled B (Outer)
        ny_b = mesh_b.mesh.n_y
        mid_j_b = ny_b // 2
        indices_b = [i * ny_b + mid_j_b for i in range(mesh_b.mesh.n_x)]
        r_b = mesh_b.mesh.x_centers  # Mesh B is created with inner_radius=0.35, so coords are absolute

        plt.plot(r_u, T_u[indices_u], 'k-', lw=3, label='Unified Model')
        plt.plot(r_a, T_a[indices_a], 'r--', lw=2, marker='o', label='Coupled Part A (Inner)')
        plt.plot(r_b, T_b[indices_b], 'b--', lw=2, marker='x', label='Coupled Part B (Outer)')
        plt.xlabel("R (Radial) (m)")

    plt.ylabel("Temperature (K)")
    plt.title(f"Temperature Profile Comparison: {case_name}")
    plt.legend()
    plt.grid(True)
    filename = f"Compare_{case_name}.png"
    plt.savefig(filename)
    logger.info(f"  [Plot] Saved to {filename}")
    plt.close()


# ==========================================
def run_case_3_1():
    # 物理参数: 左=500K, 右=Conv(300K, 20), 上下=Adia

    def set_common_bcs(phy, is_unified=False, is_left_part=False, is_right_part=False):
        # Top/Bottom Adiabatic
        phy.boundaries['top'].add_resistance_condition(T_ext=300, R_ext=1e15)
        phy.boundaries['bottom'].add_resistance_condition(T_ext=300, R_ext=1e15)

        # Left Boundary
        if is_unified or is_left_part:
            phy.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)

        # Right Boundary
        if is_unified or is_right_part:
            phy.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=20.0)

    def setup_unified(mat):
        # 宽0.1, 高0.1
        mesh = Mesh2D(x_dim=0.1, n_x=20, y_dim=0.1, n_y=10, geometry_type='cartesian')
        phy = HeatConduction2D(mesh, mat, initial_temp=500.0)
        phy.boundaries['left'].clear_conditions()
        phy.boundaries['right'].clear_conditions()
        set_common_bcs(phy, is_unified=True)
        return phy

    def setup_coupled(mat):
        # A: 宽0.05 (Left)
        mesh_a = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cartesian')
        phy_a = HeatConduction2D(mesh_a, mat, initial_temp=500.0)
        phy_a.boundaries['left'].clear_conditions();
        phy_a.boundaries['right'].clear_conditions()
        set_common_bcs(phy_a, is_left_part=True)  # Right is open for coupling

        # B: 宽0.05 (Right)
        mesh_b = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cartesian')
        phy_b = HeatConduction2D(mesh_b, mat, initial_temp=500.0)
        phy_b.boundaries['left'].clear_conditions();
        phy_b.boundaries['right'].clear_conditions()
        set_common_bcs(phy_b, is_right_part=True)  # Left is open for coupling

        # Coupling: A(Right) <-> B(Left)
        # Couplers.py: direction 'right' means obj1.right connects to obj2.left
        coupler = SolidSolidCouple2D(phy_a, phy_b, direction='right')
        return phy_a, phy_b, coupler

    run_comparison("3.1_Cartesian_LeftRight", setup_unified, setup_coupled)


# ==========================================
# 3.2 柱坐标管道换热 (上下/轴向拼接)
# ==========================================
def run_case_3_2():
    # 参数: R_in=0.3, Thick=0.05 (R_out=0.35).
    # Left(Inner)=500K, Right(Outer)=Conv(300, 20).
    # Top/Bottom=Adia (Unified ends)

    def set_common_bcs(phy, is_unified=False, is_bottom_part=False, is_top_part=False):
        # Inner (Left in Mesh2D) = 500K
        phy.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)

        # Outer (Right in Mesh2D) = Conv
        phy.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=20.0)

        # Bottom Boundary (z=min)
        if is_unified or is_bottom_part:
            phy.boundaries['bottom'].add_resistance_condition(T_ext=300, R_ext=1e15)  # Adiabatic

        # Top Boundary (z=max)
        if is_unified or is_top_part:
            phy.boundaries['top'].add_resistance_condition(T_ext=300, R_ext=1e15)  # Adiabatic

    def setup_unified(mat):
        # Length 0.2
        mesh = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.2, n_y=20, geometry_type='cylindrical', inner_radius=0.3)
        phy = HeatConduction2D(mesh, mat, initial_temp=500.0)
        for k in phy.boundaries: phy.boundaries[k].clear_conditions()
        set_common_bcs(phy, is_unified=True)
        return phy

    def setup_coupled(mat):
        # A: Length 0.1 (Bottom)
        mesh_a = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cylindrical', inner_radius=0.3)
        phy_a = HeatConduction2D(mesh_a, mat, initial_temp=500.0)
        for k in phy_a.boundaries: phy_a.boundaries[k].clear_conditions()
        set_common_bcs(phy_a, is_bottom_part=True)  # Top is coupling

        # B: Length 0.1 (Top)
        mesh_b = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cylindrical', inner_radius=0.3)
        phy_b = HeatConduction2D(mesh_b, mat, initial_temp=500.0)
        for k in phy_b.boundaries: phy_b.boundaries[k].clear_conditions()
        set_common_bcs(phy_b, is_top_part=True)  # Bottom is coupling

        # Coupling: A(Top) <-> B(Bottom)
        # Couplers.py: direction 'top' means obj1.top connects to obj2.bottom
        coupler = SolidSolidCouple2D(phy_a, phy_b, direction='top')
        return phy_a, phy_b, coupler

    run_comparison("3.2_Cylindrical_Axial", setup_unified, setup_coupled)


# ==========================================
# 3.3 柱坐标管道换热 (内外/径向拼接)
# ==========================================
def run_case_3_3():
    # 参数: Length=0.1.
    # Unified: Inner=0.3, Thick=0.1 (Out=0.4). BC: Inner=Adia, Outer=600K.
    # Coupled: A(In=0.3, Th=0.05, Out=0.35), B(In=0.35, Th=0.05, Out=0.4).
    # Coupling: A.Out <-> B.In.

    def set_common_bcs(phy, is_unified=False, is_inner_part=False, is_outer_part=False):
        # Top/Bottom Adiabatic
        phy.boundaries['top'].add_resistance_condition(T_ext=300, R_ext=1e15)
        phy.boundaries['bottom'].add_resistance_condition(T_ext=300, R_ext=1e15)

        # Inner Boundary (Left)
        if is_unified or is_inner_part:
            # Adiabatic
            phy.boundaries['left'].add_resistance_condition(T_ext=300, R_ext=1e15)

        # Outer Boundary (Right)
        if is_unified or is_outer_part:
            # Fixed 600K
            phy.boundaries['right'].add_resistance_condition(T_ext=600.0, R_ext=0.0)

    def setup_unified(mat):
        # Thick 0.1 (0.3 -> 0.4)
        mesh = Mesh2D(x_dim=0.1, n_x=20, y_dim=0.1, n_y=10, geometry_type='cylindrical', inner_radius=0.3)
        phy = HeatConduction2D(mesh, mat, initial_temp=500.0)
        for k in phy.boundaries: phy.boundaries[k].clear_conditions()
        set_common_bcs(phy, is_unified=True)
        return phy

    def setup_coupled(mat):
        # A: Inner Part (0.3 -> 0.35)
        mesh_a = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cylindrical', inner_radius=0.3)
        phy_a = HeatConduction2D(mesh_a, mat, initial_temp=500.0)
        for k in phy_a.boundaries: phy_a.boundaries[k].clear_conditions()
        set_common_bcs(phy_a, is_inner_part=True)  # Right is coupling

        # B: Outer Part (0.35 -> 0.40)
        mesh_b = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cylindrical', inner_radius=0.35)
        phy_b = HeatConduction2D(mesh_b, mat, initial_temp=500.0)
        for k in phy_b.boundaries: phy_b.boundaries[k].clear_conditions()
        set_common_bcs(phy_b, is_outer_part=True)  # Left is coupling

        # Coupling: A(Right/Outer) <-> B(Left/Inner)
        coupler = SolidSolidCouple2D(phy_a, phy_b, direction='right')
        return phy_a, phy_b, coupler

    run_comparison("3.3_Cylindrical_Radial", setup_unified, setup_coupled)


def run_case(case_name, mesh_params, bc_setup_func, t_end=1000.0):
    """
    运行单个测试用例的通用流程

    :param case_name: 算例名称
    :param mesh_params: Mesh2D 初始化字典
    :param bc_setup_func: 回调函数，用于设置物理对象的边界条件 func(physics)
    :param t_end: 模拟时长
    """
    logger.info(f"=== Starting Case: {case_name} ===")

    # 1. 准备材料与网格
    mat = MoNb()
    try:
        mesh = Mesh2D(**mesh_params)
    except Exception as e:
        logger.error(f"Mesh generation failed: {e}")
        return None

    # 2. 初始化物理求解器
    # 初始温度 500K
    physics = HeatConduction2D(mesh=mesh, material=mat, initial_temp=500.0)

    # 3. 设置边界条件
    # 先清除默认设置，应用自定义设置
    for key in physics.boundaries:
        physics.boundaries[key].clear_conditions()

    bc_setup_func(physics)

    # 4. 配置 ODE 求解器
    # 状态向量是 flattened 1D array
    y0 = physics.T.copy()

    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)

    # 5. 定义 RHS 包装器
    def rhs(t, y):
        return physics.get_derivatives(t, y)

    # 6. 求解
    t_eval = np.linspace(0, t_end, 101)
    res = solver.solve(rhs, (0, t_end), y0, t_eval=t_eval)

    if res['success']:
        logger.info(f"Case {case_name} Successful. Max Temp: {np.max(res['y'][:, -1]):.2f} K")
        return res, mesh
    else:
        logger.error(f"Case {case_name} Failed: {res['message']}")
        return None, None

def run_case(case_name, mesh_params, bc_setup_func, t_end=1000.0):
    """
    运行单个测试用例的通用流程

    :param case_name: 算例名称
    :param mesh_params: Mesh2D 初始化字典
    :param bc_setup_func: 回调函数，用于设置物理对象的边界条件 func(physics)
    :param t_end: 模拟时长
    """
    logger.info(f"=== Starting Case: {case_name} ===")

    # 1. 准备材料与网格
    mat = MoNb()
    try:
        mesh = Mesh2D(**mesh_params)
    except Exception as e:
        logger.error(f"Mesh generation failed: {e}")
        return None

    # 2. 初始化物理求解器
    # 初始温度 500K
    physics = HeatConduction2D(mesh=mesh, material=mat, initial_temp=500.0)

    # 3. 设置边界条件
    # 先清除默认设置，应用自定义设置
    for key in physics.boundaries:
        physics.boundaries[key].clear_conditions()

    bc_setup_func(physics)

    # 4. 配置 ODE 求解器
    # 状态向量是 flattened 1D array
    y0 = physics.T.copy()

    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)

    # 5. 定义 RHS 包装器
    def rhs(t, y):
        return physics.get_derivatives(t, y)

    # 6. 求解
    t_eval = np.linspace(0, t_end, 101)
    res = solver.solve(rhs, (0, t_end), y0, t_eval=t_eval)

    if res['success']:
        logger.info(f"Case {case_name} Successful. Max Temp: {np.max(res['y'][:, -1]):.2f} K")
        return res, mesh
    else:
        logger.error(f"Case {case_name} Failed: {res['message']}")
        return None, None


def plot_temperature_history(res, mesh, title, filename):
    """绘制关键点的温度随时间变化"""
    if res is None: return

    t = res['t']
    y_history = res['y']  # shape (N_nodes, N_time)

    # 选取几个特征点：内(左)、中、外(右)
    # 网格重塑逻辑：Mesh2D 中 flatten 顺序为 i*ny + j (X主序? 需检查 Mesh.py)
    # Mesh.py: flatten_index return i * self.n_y + j -> X 变化慢，Y 变化快 (Row-Major if X is Row?)
    # 通常: Index = i * ny + j.
    # Left (i=0, j=ny/2), Right (i=nx-1, j=ny/2), Center (i=nx/2, j=ny/2)

    nx, ny = mesh.n_x, mesh.n_y
    mid_y = ny // 2

    idx_left = mesh.flatten_index(0, mid_y)
    idx_center = mesh.flatten_index(nx // 2, mid_y)
    idx_right = mesh.flatten_index(nx - 1, mid_y)

    plt.figure(figsize=(10, 6))
    plt.plot(t, y_history[idx_left, :], label='Left/Inner (Mid-Y)', linewidth=2)
    plt.plot(t, y_history[idx_center, :], label='Center (Mid-Y)', linestyle='--')
    plt.plot(t, y_history[idx_right, :], label='Right/Outer (Mid-Y)', linewidth=2)

    plt.title(title)
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(filename)
    plt.close()
    logger.info(f"Plot saved to {filename}")


def plot_2d_contour(res, mesh, title, filename):
    """绘制最终时刻的 2D 温度云图"""
    if res is None:
        return

    T_final_flat = res['y'][:, -1]
    # 重塑为 2D 数组 (nx, ny)
    # 注意 Mesh.py 中的 stored logic: x_centers 是 (nx,), y_centers 是 (ny,)
    # flatten order 是 i*ny + j. 所以 reshape((nx, ny)) 是正确的
    T_2d = T_final_flat.reshape((mesh.n_x, mesh.n_y))

    # 获取网格中心用于绘图
    X, Y = np.meshgrid(mesh.x_centers, mesh.y_centers, indexing='ij')

    plt.figure(figsize=(8, 6))
    cp = plt.contourf(X, Y, T_2d, levels=20, cmap='inferno')
    plt.colorbar(cp, label='Temperature (K)')
    plt.title(title + " (Final State)")
    if mesh.geometry_type == 'cartesian':
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
    else:
        plt.xlabel('R (m)')
        plt.ylabel('Z (m)')
    plt.savefig(filename)
    plt.close()


# ==========================================
# 主程序
# ==========================================
if __name__ == "__main__":
    logger.info("Starting Coupling Verification Tests...")
    # run_case_3_1()
    run_case_3_2()
    run_case_3_3()
    logger.info("All tests completed.")

    # ---------------------------------------------------------
    # 案例 1: 左边界热流 30000, 右边界定温 550K, 上下绝热
    # ---------------------------------------------------------
    def setup_bc_case1(physics):
        # Left/Inner: Heat Flux 30000
        # 注意: 您的 HeatConduction 中 FluxBC 直接返回 q_flux。
        # 物理上通常定义流入为正。
        physics.boundaries['left'].add_flux_condition(q_flux=30000.0)

        # Right/Outer: Fixed T 550K (Resistance=0)
        physics.boundaries['right'].add_resistance_condition(T_ext=550.0, R_ext=0.0)

        # Top/Bottom: Adiabatic (默认绝热，或者显式加个大热阻)
        physics.boundaries['top'].add_resistance_condition(T_ext=300.0, R_ext=1e15)
        physics.boundaries['bottom'].add_resistance_condition(T_ext=300.0, R_ext=1e15)


    # --- 1.1 笛卡尔 ---
    mesh_params_1_cart = {
        'x_dim': 0.1, 'n_x': 20,
        'y_dim': 0.1, 'n_y': 10,
        'geometry_type': 'cartesian'
    }
    res1c, mesh1c = run_case("1.1_Cartesian", mesh_params_1_cart, setup_bc_case1)
    plot_temperature_history(res1c, mesh1c, "Case 1.1 (Cartesian): Flux -> Fixed T", "Result_1_1_History.png")
    plot_2d_contour(res1c, mesh1c, "Case 1.1 Temperature Field", "Result_1_1_Field.png")

    # --- 1.2 柱坐标 ---
    # 内径 0.3, 厚度 0.05 => 外径 0.35
    mesh_params_1_cyl = {
        'x_dim': 0.05, 'n_x': 20,  # x_dim 在圆柱中是厚度
        'y_dim': 0.1, 'n_y': 10,
        'geometry_type': 'cylindrical',
        'inner_radius': 0.3
    }
    res1r, mesh1r = run_case("1.2_Cylindrical", mesh_params_1_cyl, setup_bc_case1)
    plot_temperature_history(res1r, mesh1r, "Case 1.2 (Cylindrical): Flux -> Fixed T", "Result_1_2_History.png")
    plot_2d_contour(res1r, mesh1r, "Case 1.2 Temperature Field", "Result_1_2_Field.png")

    # ---------------------------------------------------------
    # 案例 2: 左边界绝热, 右边界对流 (300K, h=50)
    # ---------------------------------------------------------
    def setup_bc_case2(physics):
        # Left/Inner: Adiabatic
        physics.boundaries['left'].add_resistance_condition(T_ext=300.0, R_ext=1e15)

        # Right/Outer: Convection
        physics.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=50.0)

        # Top/Bottom: Adiabatic
        physics.boundaries['top'].add_resistance_condition(T_ext=300.0, R_ext=1e15)
        physics.boundaries['bottom'].add_resistance_condition(T_ext=300.0, R_ext=1e15)


    # --- 2.1 笛卡尔 ---
    res2c, mesh2c = run_case("2.1_Cartesian", mesh_params_1_cart, setup_bc_case2)
    plot_temperature_history(res2c, mesh2c, "Case 2.1 (Cartesian): Adiabatic -> Convection", "Result_2_1_History.png")

    # --- 2.2 柱坐标 ---
    res2r, mesh2r = run_case("2.2_Cylindrical", mesh_params_1_cyl, setup_bc_case2)
    plot_temperature_history(res2r, mesh2r, "Case 2.2 (Cylindrical): Adiabatic -> Convection", "Result_2_2_History.png")

    print("\nAll basic cases completed. Please check the generated PNG files.")
