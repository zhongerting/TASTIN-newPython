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

# ==============================================================================
# 测试配置
# ==============================================================================
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("TestCase3.3")


def run_test_3_3():
    """
    测试用例 3.3: 柱坐标系下管道换热 (径向耦合验证)

    对比:
    Case A (Whole): 整体厚壁圆管 (内径0.3m, 厚度0.1m)
    Case B (Split): 双层套筒结构 (A层厚0.05m + B层厚0.05m) 内外连接

    参数:
    - 长度: 0.1m
    - 初始温度: 500K
    - 边界:
        - 内边界 (r=0.3): 绝热
        - 外边界 (r=0.4): 定温 600K
        - 上下边界: 绝热
    """
    logger.info(">>> 启动测试 3.3: 柱坐标径向(内外)耦合验证 <<<")

    # --- 1. 物理参数 ---
    material = MoNb()
    length = 0.1
    r_inner = 0.3

    # 厚度参数
    thick_A = 0.05
    thick_B = 0.05
    thick_total = thick_A + thick_B  # 0.1m

    # 网格参数 (保持几何重合)
    # Z方向: 统一10个网格
    # R方向: 整体20个 -> 分层各10个
    n_z = 10
    n_r_split = 10
    n_r_total = 2 * n_r_split  # 20

    t_end = 1000.0
    t_eval = np.linspace(0, t_end, 101)

    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)

    # ==========================================================================
    # Case A: 整体模型 (Whole Model)
    # ==========================================================================
    logger.info("Building Case A (Whole Model)...")

    # Mesh2D: x_dim=厚度, inner_radius=内径
    mesh_whole = Mesh2D(
        x_dim=thick_total, n_x=n_r_total,
        y_dim=length, n_y=n_z,
        geometry_type='cylindrical',
        inner_radius=r_inner
    )

    model_whole = HeatConduction2D(mesh_whole, material, initial_temp=500.0)

    # --- 边界条件 ---
    # 1. 内边界 (Left, r=0.3): 固定热流为-6000W/m2
    model_whole.boundaries['left'].add_flux_condition(q_flux=-6000.0)
    # HeatConduction2D 默认绝热，无需操作 (或者显式 Flux=0)

    # 2. 外边界 (Right, r=0.4): 定温 600K
    model_whole.boundaries['right'].add_resistance_condition(T_ext=600.0, R_ext=0.0)

    # 3. 上下边界: 绝热 (默认)

    # --- 求解 Case A ---
    logger.info("Solving Case A...")
    res_whole = solver.solve(fun=model_whole.get_derivatives, t_span=(0, t_end), y0=model_whole.T.copy(), t_eval=t_eval)

    if not res_whole['success']:
        logger.error(f"Case A failed: {res_whole['message']}")
        return

    # ==========================================================================
    # Case B: 分层耦合模型 (Coupled Model)
    # ==========================================================================
    logger.info("Building Case B (Radial Split Model)...")

    # --- Inner Cylinder (Block A) ---
    # r: 0.3 -> 0.35
    mesh_in = Mesh2D(
        x_dim=thick_A, n_x=n_r_split,
        y_dim=length, n_y=n_z,
        geometry_type='cylindrical',
        inner_radius=r_inner
    )
    model_in = HeatConduction2D(mesh_in, material, initial_temp=500.0)

    # 边界: Left= 固定热流为-6000W/m2, Right=Coupled(待定)
    model_in.boundaries['left'].add_flux_condition(q_flux=-6000.0)

    # --- Outer Cylinder (Block B) ---
    # r: 0.35 -> 0.40
    # 注意: inner_radius 必须是 Inner Cylinder 的外径
    r_mid = r_inner + thick_A
    mesh_out = Mesh2D(
        x_dim=thick_B, n_x=n_r_split,
        y_dim=length, n_y=n_z,
        geometry_type='cylindrical',
        inner_radius=r_mid
    )
    model_out = HeatConduction2D(mesh_out, material, initial_temp=500.0)

    # 边界: Left=Coupled(待定), Right=定温 600K
    model_out.boundaries['right'].add_resistance_condition(T_ext=600.0, R_ext=0.0)

    # --- 建立耦合 ---
    # Inner 的 Right 连接 Outer 的 Left
    # direction='right' 表示 obj2 在 obj1 的右侧 (外侧)
    coupler = SolidSolidCouple2D(
        obj1=model_in,
        obj2=model_out,
        direction='right',
        contact_resistance=0.0
    )

    # --- 定义耦合系统 RHS ---
    def coupled_rhs(t, y):
        # 1. 拆分状态
        n_in = model_in.N
        y_in = y[:n_in]
        y_out = y[n_in:]

        # 2. 更新内部状态 (赋值 + 算物性 + 算内部热阻 + 推送边界快照)
        model_in.T[:] = y_in
        model_in._update_properties()
        model_in._compute_internal_resistance()
        model_in._update_boundaries_state()

        model_out.T[:] = y_out
        model_out._update_properties()
        model_out._compute_internal_resistance()
        model_out._update_boundaries_state()

        # 3. 同步边界 (交换 T_surf 和 R_int)
        coupler.sync()

        # 4. 计算导数
        dy_in = model_in.get_derivatives(t, y_in)
        dy_out = model_out.get_derivatives(t, y_out)

        return np.concatenate([dy_in, dy_out])

    # --- 求解 Case B ---
    logger.info("Solving Case B...")
    y0_coupled = np.concatenate([model_in.T, model_out.T])

    res_coupled = solver.solve(fun=coupled_rhs, t_span=(0, t_end), y0=y0_coupled, t_eval=t_eval)

    if not res_coupled['success']:
        logger.error(f"Case B failed: {res_coupled['message']}")
        return

    # ==========================================================================
    # 结果验证
    # ==========================================================================
    logger.info("Comparing results...")

    y_w = res_whole['y']
    y_c = res_coupled['y']

    # 获取稳态(最后时刻)结果
    T_final_w = y_w[:, -1]
    T_final_c = y_c[:, -1]

    # 计算误差
    # 由于两个 Mesh 生成顺序一致 (都是 Row-major 或 Column-major，且径向分层顺序自然衔接)
    # 直接拼接对比即可
    diff = np.abs(T_final_w - T_final_c)
    max_error = np.max(diff)
    l2_error = np.linalg.norm(diff) / np.linalg.norm(T_final_w)

    logger.info(f"Final Max Diff: {max_error:.6e} K")
    logger.info(f"Relative L2 Error: {l2_error:.6e}")

    if l2_error < 1e-4:
        logger.info(">>> TEST PASSED: Radial Coupling matches Whole Model! <<<")
    else:
        logger.warning(">>> TEST FAILED: Significant discrepancy detected! <<<")

    # ==========================================================================
    # 可视化
    # ==========================================================================
    plot_radial_profile(y_w, y_c, n_r_total, n_z, mesh_whole)


def plot_radial_profile(y_whole, y_coupled, n_r, n_z, mesh):
    """
    绘制径向温度分布对比 (取中间高度)
    """
    plt.figure(figsize=(10, 6))

    # 取 Z 方向中间层索引
    z_idx = n_z // 2

    # 提取径向节点温度
    # Mesh2D 的 flatten 顺序: i * n_y + j (i是x/r索引, j是y/z索引)
    # 我们需要固定 j=z_idx, 遍历所有 i=0..n_r-1
    indices = [i * n_z + z_idx for i in range(n_r)]

    r_coords = mesh.geom_data.node_centers_x.reshape((n_r, n_z))[:, z_idx]

    T_w_profile = y_whole[indices, -1]  # 最后时刻
    T_c_profile = y_coupled[indices, -1]

    plt.plot(r_coords, T_w_profile, 'o-', label='Whole Model', markersize=6)
    plt.plot(r_coords, T_c_profile, 'x--', label='Coupled Model', markersize=8)

    plt.title(f"Radial Temperature Profile at t=1000s (Mid-Height)")
    plt.xlabel("Radius (m)")
    plt.ylabel("Temperature (K)")
    plt.legend()
    plt.grid(True)
    plt.savefig("Test3.3_Radial_Result.png")
    logger.info("Plot saved to Test3.3_Radial_Result.png")


if __name__ == "__main__":
    run_test_3_3()
