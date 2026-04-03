import sys
import os
import numpy as np
import logging
import matplotlib.pyplot as plt

from MathSolvers.solver_module import NuclearODESolver
from Materials.Solids.MoNb import MoNb
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.Couplers import SolidSolidCouple2D

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger("TestCase3.1")


# --- 2. 耦合系统封装类 (核心逻辑) ---
class CoupledModelWrapper:
    """
    专门用于将两个 HeatConduction2D 对象和一个 Couple 对象
    封装成 ODE Solver 能够识别的单一右端项函数 (RHS)。
    """

    def __init__(self, obj_A, obj_B, coupler):
        self.obj_A = obj_A
        self.obj_B = obj_B
        self.coupler = coupler

        # 记录切分点，用于将长向量 y 拆分给 A 和 B
        self.len_A = self.obj_A.N
        self.len_B = self.obj_B.N

    def dydt(self, t, y_combined):
        """
        计算组合系统的导数 dy/dt。
        步骤：拆分 -> 更新状态 -> 同步(Sync) -> 计算热流 -> 合并导数
        """
        # 1. 拆分状态向量
        y_A = y_combined[:self.len_A]
        y_B = y_combined[self.len_A:]

        # 2. 分别更新 A 和 B 的内部状态 (温度、物性、热阻)
        # -------------------------------------------------
        # 对象 A
        self.obj_A.T[:] = y_A
        self.obj_A._update_properties()
        self.obj_A._compute_internal_resistance()
        self.obj_A._update_boundaries_state()  # 关键：算出 A 的边界 T 和 R

        self.obj_A._compute_fluxes(t)

        # 对象 B
        self.obj_B.T[:] = y_B
        self.obj_B._update_properties()
        self.obj_B._compute_internal_resistance()
        self.obj_B._update_boundaries_state()  # 关键：算出 B 的边界 T 和 R

        self.obj_B._compute_fluxes(t)

        # 3. 【核心】执行耦合器同步
        # -------------------------------------------------
        # 这一步会交换 A 和 B 接触面的 T_ext 和 R_ext
        self.coupler.sync()

        # 4. 计算各自的热流和源项
        # -------------------------------------------------
        # 对象 A
        Q_net_A = self.obj_A._compute_fluxes(t)  # 此时用的是同步后的边界条件
        self.obj_A._update_sources(t)
        dTdt_A = (Q_net_A + self.obj_A.Q_source) / self.obj_A.thermal_capacitance

        # 对象 B
        Q_net_B = self.obj_B._compute_fluxes(t)
        self.obj_B._update_sources(t)
        dTdt_B = (Q_net_B + self.obj_B.Q_source) / self.obj_B.thermal_capacitance

        # 5. 合并导数并返回
        return np.concatenate((dTdt_A, dTdt_B))


# --- 3. 测试主程序 ---
if __name__ == "__main__":
    logger.info("=== 开始执行测试案例 3.1: 笛卡尔平板拼接验证 ===")

    # 公共参数
    mat = MoNb()
    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)
    t_end = 5000.0
    t_eval = np.linspace(0, t_end, 5001)

    # ==========================================
    # 场景 A: 整体模型 (Unified)
    # 尺寸: 0.1m x 0.1m
    # ==========================================
    logger.info("\n--- [1/2] 正在计算整体模型 (Unified Model) ---")

    # 1. 网格 (Cartesian)
    mesh_u = Mesh2D(x_dim=0.1, n_x=20, y_dim=0.1, n_y=10, geometry_type='cartesian')

    # 2. 物理对象
    phy_u = HeatConduction2D(mesh_u, mat, initial_temp=500.0)

    # 3. 边界条件
    # 清除默认
    for k in phy_u.boundaries:
        phy_u.boundaries[k].clear_conditions()

    # 左: 定温 500K
    phy_u.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)
    # 右: 对流 300K, h=20
    phy_u.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=800.0)
    # 上下: 绝热 (给个极大热阻)
    phy_u.boundaries['top'].add_flux_condition(q_flux=0.0)
    phy_u.boundaries['bottom'].add_flux_condition(q_flux=0.0)

    # 给定一个恒定热源，大小为900000W/m3
    q_volumetric = 900000
    source_array = np.full(phy_u.N, q_volumetric) * mesh_u.geom_data.volumes
    phy_u.link_source_buffer(source_array)

    # 4. 求解
    res_u = solver.solve(phy_u.get_derivatives, (0, t_end), phy_u.T.copy(), t_eval=t_eval)

    if not res_u['success']:
        logger.error("整体模型计算失败!")
        sys.exit(1)
    logger.info(f"整体模型计算完成。最终平均温度: {np.mean(res_u['y'][:, -1]):.4f} K")

    # ==========================================
    # 场景 B: 拼接模型 (Coupled)
    # 左块 A (0.05x0.1) + 右块 B (0.05x0.1)
    # ==========================================
    logger.info("\n--- [2/2] 正在计算耦合模型 (Coupled Model) ---")

    # 1. 网格 (注意 x_dim 是 0.05)
    mesh_a = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cartesian')
    mesh_b = Mesh2D(x_dim=0.05, n_x=10, y_dim=0.1, n_y=10, geometry_type='cartesian')

    # 2. 物理对象
    phy_a = HeatConduction2D(mesh_a, mat, initial_temp=500.0)
    phy_b = HeatConduction2D(mesh_b, mat, initial_temp=500.0)

    # 3. 边界条件
    # --- 块 A (左边) ---
    for k in phy_a.boundaries:
        phy_a.boundaries[k].clear_conditions()
    # A左: 定温 500K
    phy_a.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)
    # A上下: 绝热
    phy_a.boundaries['top'].add_flux_condition(q_flux=0.0)
    phy_a.boundaries['bottom'].add_flux_condition(q_flux=0.0)
    # A右: 待耦合 (暂空或默认)

    # 给定一个恒定热源，大小为900000W/m3
    source_array_a = np.full(phy_a.N, q_volumetric) * mesh_a.geom_data.volumes
    phy_a.link_source_buffer(source_array_a)

    # --- 块 B (右边) ---
    for k in phy_b.boundaries:
        phy_b.boundaries[k].clear_conditions()
    # B左: 待耦合
    # B右: 对流 300K, h=20
    phy_b.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=800.0)
    # B上下: 绝热
    phy_b.boundaries['top'].add_flux_condition(q_flux=0.0)
    phy_b.boundaries['bottom'].add_flux_condition(q_flux=0.0)

    # 给定一个恒定热源，大小为900000W/m3
    source_array_b = np.full(phy_b.N, q_volumetric) * mesh_b.geom_data.volumes
    phy_b.link_source_buffer(source_array_b)

    # 4. 建立耦合
    # A 的右边 (Right/Outer) 连接 B 的左边 (Left/Inner)
    # 笛卡尔坐标下 direction='right' 表示 obj1 在左，obj2 在右，连接 obj1.right 和 obj2.left
    coupler = SolidSolidCouple2D(phy_a, phy_b, direction='right')

    # 5. 封装求解
    wrapper = CoupledModelWrapper(phy_a, phy_b, coupler)
    y0_coupled = np.concatenate((phy_a.T.copy(), phy_b.T.copy()))

    res_c = solver.solve(wrapper.dydt, (0, t_end), y0_coupled, t_eval=t_eval)

    if not res_c['success']:
        logger.error("耦合模型计算失败!")
        sys.exit(1)
    logger.info(f"耦合模型计算完成。最终平均温度: {np.mean(res_c['y'][:, -1]):.4f} K")

    # ==========================================
    # 结果对比与绘图
    # ==========================================
    logger.info("\n--- 结果验证 ---")

    # 提取最终时刻温度场
    T_u = res_u['y'][:, -1]  # 整体
    T_c = res_c['y'][:, -1]  # 耦合总和

    # 计算全场平均误差
    # 注意：如果网格划分完全一致（整体20格，分开各10格），节点位置应该是一一对应的
    # 直接对比数组值
    if len(T_u) == len(T_c):
        max_diff = np.max(np.abs(T_u - T_c))
        avg_diff = np.mean(np.abs(T_u - T_c))
        logger.info(f"最大节点温差: {max_diff:.6e} K")
        logger.info(f"平均节点温差: {avg_diff:.6e} K")

        if avg_diff < 1e-3:
            logger.info("✅ 验证通过：结果一致！")
        else:
            logger.warning("⚠️ 验证警告：存在较大误差，请检查网格或边界设置。")
    else:
        logger.warning("网格节点数不一致，无法进行逐点对比。")

    # 绘图：取中轴线 (y=0.05m) 处的温度分布
    plt.figure(figsize=(10, 6))

    # 1. 整体模型曲线
    # Reshape为 (nx, ny)
    T_u_2d = T_u.reshape((mesh_u.n_x, mesh_u.n_y))
    mid_j = mesh_u.n_y // 2
    plt.plot(mesh_u.x_centers, T_u_2d[:, mid_j], 'k-', linewidth=3, label='Unified (0.1m)')

    # 2. 耦合模型曲线
    # 拆分数据
    T_a = T_c[:phy_a.N].reshape((mesh_a.n_x, mesh_a.n_y))
    T_b = T_c[phy_a.N:].reshape((mesh_b.n_x, mesh_b.n_y))

    # A 的坐标 (0 ~ 0.05)
    plt.plot(mesh_a.x_centers, T_a[:, mid_j], 'r--', marker='o', markersize=4, label='Block A (0~0.05m)')

    # B 的坐标 (0.05 ~ 0.1) -> 需要加上偏移量 0.05
    plt.plot(mesh_b.x_centers + 0.05, T_b[:, mid_j], 'b--', marker='x', markersize=4, label='Block B (0.05~0.1m)')

    plt.title("Case 3.1: Cartesian Plate Coupling Verification")
    plt.xlabel("X Position (m)")
    plt.ylabel("Temperature (K)")
    plt.legend()
    plt.grid(True)

    img_name = "Case_3_1_Result.png"
    plt.savefig(img_name)
    logger.info(f"结果对比图已保存为: {img_name}")
    plt.show()
