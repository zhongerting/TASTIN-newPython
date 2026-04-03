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
logger = logging.getLogger("TestCase3.2")


def run_test_3_2():
    """
    测试用例 3.2: 柱坐标系下管道换热 (轴向耦合验证)

    对比:
    Case A (Whole): 整体长管 (L=0.2m)
    Case B (Split): 两段短管 (L=0.1m + L=0.1m) 上下连接

    参数:
    - 内径: 0.3m
    - 厚度: 0.05m (外径 0.35m)
    - 初始温度: 500K
    - 左边界(内): 定温 500K
    - 右边界(外): 对流 T_inf=300K, h=20 W/m2K
    - 上下边界: 绝热
    """
    logger.info(">>> 启动测试 3.2: 柱坐标轴向耦合验证 <<<")

    # --- 1. 物理参数定义 ---
    material = MoNb()
    r_in = 0.3
    thickness = 0.05  # x_dim
    length_total = 0.2

    # 网格参数
    # 为了精确对比，必须保证几何节点完全重合
    # R方向: 5个网格
    # Z方向: 整体20个网格 -> 分段各10个网格
    n_r = 5
    n_z_total = 20
    n_z_split = 10

    t_end = 2000.0
    t_eval = np.linspace(0, t_end, 2001)

    # 求解器
    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)

    # ==========================================================================
    # Case A: 整体模型 (Whole Model)
    # ==========================================================================
    logger.info("Building Case A (Whole Model)...")

    # Mesh2D(x_dim, n_x, y_dim, n_y, geometry_type, inner_radius)
    # 注意: x_dim 在圆柱坐标下代表径向厚度 (R_out - R_in)
    mesh_whole = Mesh2D(
        x_dim=thickness, n_x=n_r,
        y_dim=length_total, n_y=n_z_total,
        geometry_type='cylindrical',
        inner_radius=r_in
    )

    model_whole = HeatConduction2D(mesh_whole, material, initial_temp=500.0)

    # --- 边界条件设置 ---
    # HeatConduction2D 边界键: 'left'(内), 'right'(外), 'bottom', 'top'

    # 1. 左边界(内表面): 定温 500K
    # 使用 add_resistance_condition(T_ext=500, R_ext=0) 模拟定温
    model_whole.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)

    # 2. 右边界(外表面): 对流 T=300K, h=20
    model_whole.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=200.0)

    # 3. 上下边界: 绝热 (HeatConduction2D 默认无 Flux，即绝热，无需显式添加)
    model_whole.boundaries['top'].add_flux_condition(q_flux=0.0)
    model_whole.boundaries['bottom'].add_flux_condition(q_flux=0.0)

    # 给定一个恒定热源，大小为1000000W/m3
    q_volumetric = 1000000
    source_array_whole = np.full(model_whole.N, q_volumetric) * mesh_whole.geom_data.volumes
    model_whole.link_source_buffer(source_array_whole)

    # --- 求解 Case A ---
    logger.info("Solving Case A...")
    res_whole = solver.solve(fun=model_whole.get_derivatives, t_span=(0, t_end), y0=model_whole.T.copy(), t_eval=t_eval)

    if not res_whole['success']:
        logger.error(f"Case A failed: {res_whole['message']}")
        return

    # ==========================================================================
    # Case B: 分段耦合模型 (Coupled Model)
    # ==========================================================================
    logger.info("Building Case B (Split Coupled Model)...")

    # --- Bottom Pipe (Z: 0.0 -> 0.1) ---
    mesh_bot = Mesh2D(
        x_dim=thickness, n_x=n_r,
        y_dim=0.1, n_y=n_z_split,
        geometry_type='cylindrical',
        inner_radius=r_in
    )
    model_bot = HeatConduction2D(mesh_bot, material, initial_temp=500.0)

    # 边界
    model_bot.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)
    model_bot.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=200.0)
    # Top 边界留给耦合器处理
    # 下边界绝热
    model_bot.boundaries['bottom'].add_flux_condition(q_flux=0.0)

    # 给定一个恒定热源，大小为1000000W/m3
    q_volumetric = 1000000
    source_array = np.full(model_bot.N, q_volumetric) * mesh_bot.geom_data.volumes
    model_bot.link_source_buffer(source_array)

    # --- Top Pipe (Z: 0.1 -> 0.2) ---
    mesh_top = Mesh2D(
        x_dim=thickness, n_x=n_r,
        y_dim=0.1, n_y=n_z_split,
        geometry_type='cylindrical',
        inner_radius=r_in
    )
    model_top = HeatConduction2D(mesh_top, material, initial_temp=500.0)

    # 边界
    model_top.boundaries['left'].add_resistance_condition(T_ext=500.0, R_ext=0.0)
    model_top.boundaries['right'].add_convection_condition(T_fluid=300.0, h_coeff=200.0)
    # Bottom 边界留给耦合器处理
    # 上边界绝热
    model_top.boundaries['top'].add_flux_condition(q_flux=0.0)

    # 给定一个恒定热源，大小为1000000W/m3
    q_volumetric = 1000000
    source_array_top = np.full(model_top.N, q_volumetric) * mesh_top.geom_data.volumes
    model_top.link_source_buffer(source_array_top)

    # --- 建立耦合 ---
    # Bottom 的 Top 连接 Top 的 Bottom
    # SolidSolidCouple2D(obj1, obj2, direction='top')
    # 这里的 direction='top' 意味着 obj2 在 obj1 的上方 (obj1.top <-> obj2.bottom)
    coupler = SolidSolidCouple2D(
        obj1=model_bot,
        obj2=model_top,
        direction='top',
        contact_resistance=0.0
    )

    # --- 定义耦合系统的 RHS ---
    # 需要手动拼接状态向量，并在计算导数前同步边界条件

    def coupled_rhs(t, y):
        # 1. 拆分状态向量
        n_bot = model_bot.N
        y_bot = y[:n_bot]
        y_top = y[n_bot:]

        # 2. 更新模型内部状态 (为了计算物性和热阻)
        # 注意: 这里手动赋值是为了确保 sync() 能读取到最新的 T_surface
        # get_derivatives 内部也会赋值，但 sync 需要在 get_derivatives 之前拿到最新 T
        model_bot.T[:] = y_bot
        model_top.T[:] = y_top

        # 3. 必须先更新物性，才能计算出正确的内部热阻 R_int，供 sync 使用
        model_bot._update_properties()
        model_bot._compute_internal_resistance()
        model_bot._update_boundaries_state()  # 更新 T_surf, R_int 到 BoundaryRegion

        model_bot._compute_fluxes(t)

        model_top._update_properties()
        model_top._compute_internal_resistance()
        model_top._update_boundaries_state()

        model_top._compute_fluxes(t)

        # 4. 执行耦合同步 (交换边界条件)
        coupler.sync()

        # 5. 计算导数 (此时边界条件已是最新)
        # get_derivatives 会再次调用 update_properties 等，但开销可接受
        dy_bot = model_bot.get_derivatives(t, y_bot)
        dy_top = model_top.get_derivatives(t, y_top)

        return np.concatenate([dy_bot, dy_top])

    # --- 求解 Case B ---
    logger.info("Solving Case B...")
    y0_coupled = np.concatenate([model_bot.T, model_top.T])

    res_coupled = solver.solve(fun=coupled_rhs, t_span=(0, t_end), y0=y0_coupled, t_eval=t_eval)

    if not res_coupled['success']:
        logger.error(f"Case B failed: {res_coupled['message']}")
        return

    # ==========================================================================
    # 结果对比与验证
    # ==========================================================================
    logger.info("Comparing results...")

    y_w = res_whole['y']  # Shape: (N_total, Time)
    y_c = res_coupled['y']  # Shape: (N_bot + N_top, Time)

    # 获取最后时刻的温度场
    T_end_whole = y_w[:, -1]
    T_end_coupled = y_c[:, -1]

    # 计算误差
    diff = np.abs(T_end_whole - T_end_coupled)
    max_error = np.max(diff)
    l2_error = np.linalg.norm(diff) / np.linalg.norm(T_end_whole)

    logger.info(f"Final Max Diff: {max_error:.6e} K")
    logger.info(f"Relative L2 Error: {l2_error:.6e}")

    if l2_error < 1e-4:
        logger.info(">>> TEST PASSED: Coupling matches whole model! <<<")
    else:
        logger.warning(">>> TEST FAILED: Discrepancy detected! <<<")

    # ==========================================================================
    # 可视化
    # ==========================================================================
    plot_results(t_eval, y_w, y_c, n_r, n_z_total)


def plot_results(t, y_whole, y_coupled, n_r, n_z):
    """
    绘制特定节点的温度随时间变化
    """
    plt.figure(figsize=(10, 6))

    # 选取几个特征点
    # 1. 底部-内侧 (Bottom-Inner) -> Index 0
    # 2. 中部-外侧 (Middle-Outer) -> Interface 附近的外表面
    #    Middle z index approx n_z // 2. Outer x index is n_r - 1.
    #    Flatten index = i * n_y + j  (注意 Mesh2D flatten 顺序)
    #    Mesh2D 逻辑: n_x (rows), n_y (cols).
    #    Flatten: row-major? Numpy default is 'C' style (last axis fastest).
    #    So index = i * n_y + j.
    #    Middle Z (j = n_z // 2), Outer R (i = n_r - 1)

    idx_bot_in = 0  # (0, 0)

    j_mid = n_z // 2
    i_out = n_r - 1
    idx_mid_out = i_out * n_z + j_mid

    j_top = n_z - 1
    i_in = 0
    idx_top_in = i_in * n_z + j_top

    points = [
        (idx_bot_in, 'Bottom-Inner (Inlet)'),
        (idx_mid_out, 'Mid-Outer (Coupling Interface)'),
        (idx_top_in, 'Top-Inner (Outlet)')
    ]

    for idx, label in points:
        line, = plt.plot(t, y_whole[idx, :], '-', lw=2, label=f'Whole: {label}')
        plt.plot(t, y_coupled[idx, :], '--', color=line.get_color(), lw=2, label=f'Coupled: {label}')

    plt.title("Test 3.2: Cylindrical Axial Coupling Verification")
    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("Test3.2_Result.png")
    logger.info("Plot saved to Test3.2_Result.png")


if __name__ == "__main__":
    run_test_3_2()