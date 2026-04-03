import numpy as np
import matplotlib.pyplot as plt
import logging
from typing import Tuple, Literal

from MathSolvers.solver_module import NuclearODESolver
from Solvers.HeatConduction.Mesh import Mesh1D
from Solvers.HeatConduction.Boundary import BoundaryRegion
from Solvers.HeatConduction.HeatConduction import HeatConduction1D

from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.MoNb import MoNb

from Solvers.Couplers import SolidSolidCouple1D

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# ==============================================================================
# 绘图辅助：合并两个网格的数据
# ==============================================================================
def merge_results_for_plot(mesh1, mesh2, y1_res, y2_res, geo_type):
    """
    将两个物体的网格和结果合并，以便绘制统一的云图
    """
    # 1. 坐标合并
    x1 = mesh1.node_centers

    # 笛卡尔坐标下，Mesh2 的坐标是从 0 开始的，需要加上 Mesh1 的长度平移
    # 柱坐标下，Mesh2 的坐标是绝对半径，直接使用即可
    if geo_type == 'cartesian':
        offset = mesh1.L_or_R
        x2 = mesh2.node_centers + offset
    else:
        x2 = mesh2.node_centers

    x_combined = np.concatenate([x1, x2])

    # 2. 温度合并
    # y_res shape: (N, TimeSteps)
    # 按节点维度拼接
    y_combined = np.concatenate([y1_res, y2_res], axis=0)

    return x_combined, y_combined


def plot_coupled_spacetime(ax, x_coords, times, temps, title_prefix):
    """绘制合并后的时空云图"""
    T_grid, X_grid = np.meshgrid(times, x_coords)

    # 使用 pcolormesh 绘图
    c = ax.pcolormesh(T_grid, X_grid, temps, cmap='plasma', shading='auto')

    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label('Temperature [K]')

    ax.set_title(f"{title_prefix}\nSpacetime Map (Coupled)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [m]")

    # 标出界面位置
    interface_pos = x_coords[len(x_coords) // 2]  # 粗略取中间，仅作参考
    # 实际上界面在 mesh1 的最后一个面
    # 这里我们只画图，不画线了，云图能看出来

# ==============================================================================
# 辅助求解函数
# ==============================================================================
def solve_system_every_second(hc_object, t_end):
    """
    通用求解驱动函数：强制每隔 1s 输出一次数据
    """

    # 1. 定义 RHS 包装器
    def rhs_wrapper(t, y):
        return hc_object.get_derivatives(t, y)

    # 2. 配置求解器 (使用 BDF 处理刚性)
    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)

    # 3. 设置时间步: 0, 1, 2, ..., t_end (共 t_end+1 个点)
    t_eval = np.arange(0, t_end + 1.0, 1.0)

    # 4. 执行求解
    res = solver.solve(fun=rhs_wrapper, t_span=(0, t_end), y0=hc_object.T.copy(), t_eval=t_eval)

    if not res['success']:
        raise RuntimeError(f"Solver failed: {res['message']}")

    return res['t'], res['y']


# ==============================================================================
# 绘图函数：二维时空云图
# ==============================================================================
def plot_spacetime_contour(ax, mesh, times, temps, title_prefix):
    """
    绘制 (X=时间, Y=位置, Color=温度) 的平滑云图
    """
    # 准备网格数据
    # times: (T,) -> X轴
    # mesh.node_centers: (N,) -> Y轴
    # temps: (N, T) -> Z轴 (颜色)

    # 创建网格
    T_grid, X_grid = np.meshgrid(times, mesh.node_centers)

    # 绘制平滑云图
    # 修改点：使用 contourf 代替 pcolormesh
    # levels=100 表示将颜色分为100个层级，视觉上会形成连续的渐变
    c = ax.contourf(T_grid, X_grid, temps, levels=100, cmap='turbo')

    # 备选修改点：如果你必须用 pcolormesh，可以使用 gouraud 插值
    # c = ax.pcolormesh(T_grid, X_grid, temps, cmap='turbo', shading='gouraud')

    # 添加颜色条
    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label('Temperature [K]')

    # (可选) 消除 contourf 生成的极细微白色间隔线（用于导出PDF/SVG时）
    # for collection in c.collections:
    #     collection.set_edgecolor("face")

    # 标签与标题
    ax.set_title(f"{title_prefix}\nSpacetime Temperature Map")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [m]")

    # 设置显示范围
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(mesh.node_centers[0], mesh.node_centers[-1])


# ==============================================================================
# 工况 1: 固定温度 + 绝热
# ==============================================================================
def run_case_1(geo_type: Literal['cartesian', 'cylindrical']):
    logger.info(f"--- Running Case 1: {geo_type} ---")

    # --- 1. 参数定义 ---
    # 题目要求：计算域尺寸 0.1m
    # 柱坐标：内径 0.05m -> 外径 = 0.05 + 0.1 = 0.15m
    length = 0.1
    r_in = 0.05 if geo_type == 'cylindrical' else 0.0

    # Mesh 定义: total_dim 在圆柱下为外径，在笛卡尔下为总长
    total_dim = r_in + length if geo_type == 'cylindrical' else length

    N_cells = 50  # 网格数 (增加一点以提高云图细腻度)

    # --- 2. 构建对象 ---
    mesh = Mesh1D(total_dim=total_dim, n_volumes=N_cells,
                  geometry_type=geo_type, inner_radius=r_in)
    mat = AusteniticStainlessSteel()
    hc = HeatConduction1D(mesh, mat, initial_temp=300.0)

    # --- 3. 设置边界条件 ---
    # 左边界 (Inner): 固定温度 500K (使用极小热阻模拟)
    b_left = BoundaryRegion(shape=(1,), area_array=np.array([hc.mesh.face_areas[0]]))
    b_left.add_resistance_condition(T_ext=500.0, R_ext=1e-7)
    hc.attach_boundary('inner', b_left)

    # 右边界 (Outer): 绝热
    b_right = BoundaryRegion(shape=(1,), area_array=np.array([hc.mesh.face_areas[-1]]))
    b_right.add_flux_condition(q_flux=0.0)
    hc.attach_boundary('outer', b_right)

    # --- 4. 求解 ---
    t_final = 500.0
    times, temps = solve_system_every_second(hc, t_final)

    return mesh, times, temps


# ==============================================================================
# 工况 2: 热流 (Left) + 对流 (Right)
# ==============================================================================
def run_case_2(geo_type: Literal['cartesian', 'cylindrical']):
    logger.info(f"--- Running Case 2: {geo_type} ---")

    # --- 1. 参数定义 ---
    # 计算域尺寸 0.1m
    length = 0.02

    # 几何设置
    if geo_type == 'cylindrical':
        r_in = 0.05
        # 柱坐标 total_dim 是外径
        total_dim = r_in + length
    else:
        r_in = 0.0
        # 笛卡尔 total_dim 是总长
        total_dim = length

    N_cells = 150  # 网格数

    # --- 2. 构建对象 ---
    mesh = Mesh1D(total_dim=total_dim, n_volumes=N_cells,
                  geometry_type=geo_type, inner_radius=r_in)

    mat = AusteniticStainlessSteel()

    # 初始温度 300K
    hc = HeatConduction1D(mesh, mat, initial_temp=300.0)

    # --- 3. 设置边界条件 ---

    # 左边界 (Inner): 热流密度 10000 W/m^2 (加热)
    # FluxBC 正值代表流入节点 (Input/Heating)
    b_left = BoundaryRegion(shape=(1,), area_array=np.array([hc.mesh.face_areas[0]]))
    b_left.add_flux_condition(q_flux=10000.0)
    hc.attach_boundary('inner', b_left)

    # 右边界 (Outer): 对流 (T=300K, h=500)
    # Boundary 类会自动将 h 和 A 转换为热阻 R_conv
    b_right = BoundaryRegion(shape=(1,), area_array=np.array([hc.mesh.face_areas[-1]]))
    b_right.add_convection_condition(T_fluid=300.0, h_coeff=500.0)
    hc.attach_boundary('outer', b_right)

    # --- 4. 求解 ---
    t_final = 10000.0
    times, temps = solve_system_every_second(hc, t_final)

    return mesh, times, temps


# ==============================================================================
# 工况 3: 接触耦合 (Coupling)
# ==============================================================================
def run_case_3(geo_type: Literal['cartesian', 'cylindrical']):
    logger.info(f"--- Running Case 3: {geo_type} (Coupled) ---")

    # --- 1. 几何定义 ---
    len1 = 0.08
    len2 = 0.02

    if geo_type == 'cylindrical':
        # 柱坐标: 定义绝对半径
        r0 = 0.05  # 内径
        r1 = 0.05 + len1  # 界面
        r2 = 0.05 + len1 + len2  # 外径

        # Mesh1D(total_dim=外径, inner_radius=内径)
        mesh1 = Mesh1D(total_dim=r1, n_volumes=40, geometry_type='cylindrical', inner_radius=r0)
        mesh2 = Mesh1D(total_dim=r2, n_volumes=20, geometry_type='cylindrical', inner_radius=r1)
    else:
        # 笛卡尔: 定义各段长度 (Mesh2 也是从 0 开始，后续绘图时平移)
        mesh1 = Mesh1D(total_dim=len1, n_volumes=40, geometry_type='cartesian')
        mesh2 = Mesh1D(total_dim=len2, n_volumes=20, geometry_type='cartesian')

    # --- 2. 物体初始化 ---
    # 物体 1: SS, 初始 300K
    obj1 = HeatConduction1D(mesh1, AusteniticStainlessSteel(), initial_temp=300.0)

    # 设置体积热源 50000 W/m3
    # HeatConduction 需要的是总功率 [W]，所以需要乘以体积
    q_volumetric = 50000.0
    source_array = np.full(obj1.N, q_volumetric) * mesh1.volumes
    obj1.link_source_buffer(source_array)  # 静态绑定

    # 物体 2: MoNb, 初始 300K
    obj2 = HeatConduction1D(mesh2, MoNb(), initial_temp=300.0)

    # --- 3. 边界条件 ---

    # Obj1 左 (Inner): 绝热
    b1_in = BoundaryRegion(shape=(1,),area_array=np.array([obj1.mesh.face_areas[0]]))
    b1_in.add_flux_condition(0.0)
    obj1.attach_boundary('inner', b1_in)

    # Obj1 右 (Outer): 接触面 (初始设为绝热占位，RHS中更新)
    # 使用 ResistanceBC 模拟接触，R_ext 和 T_ext 将动态来自 Obj2
    b1_out = BoundaryRegion(shape=(1,),area_array=np.array([obj1.mesh.face_areas[-1]]))
    # 初始: T_ext=300, R_ext=0 (无阻), R_add=0 (无接触热阻)
    contact_bc_1 = b1_out.add_resistance_condition(T_ext=300.0, R_ext=0.0, R_add=0.0)
    obj1.attach_boundary('outer', b1_out)

    # Obj2 左 (Inner): 接触面
    b2_in = BoundaryRegion(shape=(1,),area_array=np.array([obj1.mesh.face_areas[0]]))
    contact_bc_2 = b2_in.add_resistance_condition(T_ext=300.0, R_ext=0.0, R_add=0.0)
    obj2.attach_boundary('inner', b2_in)

    # Obj2 右 (Outer): 复合边界 (对流 + 固定散热)
    b2_out = BoundaryRegion(shape=(1,),area_array=np.array([obj1.mesh.face_areas[-1]]))
    # 1. 对流: 300K, h=300
    b2_out.add_convection_condition(T_fluid=300.0, h_coeff=300.0)
    # 2. 固定散热: 2000 W/m2 (流出)
    # FluxBC 正值=流入，负值=流出
    b2_out.add_flux_condition(q_flux=-3000.0)
    obj2.attach_boundary('outer', b2_out)

    # --- 4. 定义耦合微分方程 ---

    N1 = obj1.N
    N2 = obj2.N

    def coupled_rhs(t, y_combined):
        # A. 拆分状态向量
        y1 = y_combined[:N1]
        y2 = y_combined[N1:]

        # B. 赋值给对象 (触发内部 T 更新)
        # 注意: get_derivatives 会再次赋值，但我们需要先赋值以更新物性 R_int
        # 为了避免重复计算，这里直接调用 get_derivatives 的前置步骤是最高效的，
        # 但为了代码清晰，我们显式手动调用更新方法。

        obj1.T[:] = y1
        obj2.T[:] = y2

        # C. 更新物性和内部热阻 (计算 R_internal)
        obj1._update_properties()
        obj1._compute_internal_resistance()
        obj1._update_boundaries_state()  # 更新边界上的 T_node, R_int

        obj2._update_properties()
        obj2._compute_internal_resistance()
        obj2._update_boundaries_state()

        # D. === 交换边界信息 (Coupling) ===
        # Obj1 Outer 看到的是 -> Obj2 Inner
        # 获取 Obj2 Inner 处的信息
        # T_ext 是 Obj2 第一个节点的温度
        # R_ext 是 Obj2 从表面到第一个节点的内部热阻
        T2_node0 = obj2.boundaries['inner'].T_adj_node[0]  # 标量
        R2_int0 = obj2.boundaries['inner'].R_internal[0]

        # 更新 Obj1 的边界条件
        contact_bc_1.update_params(T_ext=T2_node0, R_ext=R2_int0)

        # Obj2 Inner 看到的是 -> Obj1 Outer
        T1_nodeN = obj1.boundaries['outer'].T_adj_node[0]
        R1_intN = obj1.boundaries['outer'].R_internal[0]

        # 更新 Obj2 的边界条件
        contact_bc_2.update_params(T_ext=T1_nodeN, R_ext=R1_intN)

        # E. 计算导数
        # 由于我们已经更新了 BC 参数，现在调用 get_derivatives 会使用最新的边界条件计算 Flux
        dy1 = obj1.get_derivatives(t, y1)
        dy2 = obj2.get_derivatives(t, y2)

        return np.concatenate([dy1, dy2])

    # --- 5. 求解 ---
    t_end = 3000.0
    y0 = np.concatenate([obj1.T, obj2.T])

    # 设置求解器
    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)
    t_eval = np.arange(0, t_end + 1.0, 1.0)

    res = solver.solve(coupled_rhs, (0, t_end), y0, t_eval=t_eval)

    if not res['success']:
        raise RuntimeError(f"Solver failed: {res['message']}")

    # --- 6. 结果拆分 ---
    y_total = res['y']  # (N1+N2, Time)
    y1_res = y_total[:N1, :]
    y2_res = y_total[N1:, :]
    times = res['t']

    return mesh1, mesh2, times, y1_res, y2_res

# ==============================================================================
# 工况 3: 接触耦合 (Coupling)
# 采用Couplers.py 中的耦合方法
# ==============================================================================
def run_case_3_new(geo_type: Literal['cartesian', 'cylindrical']):
    logger.info(f"--- Running Case 3 (Coupled) : {geo_type} ---")

    # --- 1. 几何与网格 ---
    len1 = 0.08
    len2 = 0.12

    if geo_type == 'cylindrical':
        r0 = 0.05
        r1 = 0.05 + len1
        r2 = 0.05 + len1 + len2
        mesh1 = Mesh1D(total_dim=r1, n_volumes=40, geometry_type='cylindrical', inner_radius=r0)
        mesh2 = Mesh1D(total_dim=r2, n_volumes=20, geometry_type='cylindrical', inner_radius=r1)
    else:
        mesh1 = Mesh1D(total_dim=len1, n_volumes=40, geometry_type='cartesian')
        mesh2 = Mesh1D(total_dim=len2, n_volumes=20, geometry_type='cartesian')

    # --- 2. 物体初始化 ---
    # Obj1: 不锈钢, 初始300K, 内热源 5000 W/m3
    obj1 = HeatConduction1D(mesh1, AusteniticStainlessSteel(), initial_temp=300.0)
    # [关键] 设置体积热源 (转换为总功率数组)
    q_vol_source = np.full(obj1.N, 50000.0)
    obj1.link_source_buffer(q_vol_source * mesh1.volumes)

    # Obj2: MoNb, 初始300K
    obj2 = HeatConduction1D(mesh2, MoNb(), initial_temp=300.0)

    # --- 3. 边界条件 (非耦合面) ---
    # Obj1 左侧: 绝热
    b1_in = BoundaryRegion(shape=(1,), area_array=np.array([obj1.mesh.face_areas[0]]))
    b1_in.add_flux_condition(0.0)
    obj1.attach_boundary('inner', b1_in)

    # Obj2 右侧: 对流 + 散热
    b2_out = BoundaryRegion(shape=(1,), area_array=np.array([obj2.mesh.face_areas[-1]]))
    b2_out.add_convection_condition(T_fluid=300.0, h_coeff=300.0)
    b2_out.add_flux_condition(q_flux=-3000.0)  # 负值代表流出
    obj2.attach_boundary('outer', b2_out)

    # --- 4. 【关键】建立耦合 ---
    # 使用您在 Couplers.py 中定义的类
    # Obj2 在 Obj1 的 "Outer" (右/外) 方向
    couple = SolidSolidCouple1D(obj1, obj2, direction='outer', contact_resistance=0.0)

    # --- 5. 定义 RHS 函数 ---
    N1 = obj1.N

    def coupled_rhs(t, y_combined):
        # A. 拆分并赋值状态
        y1, y2 = y_combined[:N1], y_combined[N1:]
        obj1.T[:] = y1
        obj2.T[:] = y2

        # B. 更新各自内部状态 (使用 HeatConduction.py 中新加的 pre_step_update)
        obj1.pre_step_update()
        obj2.pre_step_update()

        # C. 同步边界条件 (使用 Couplers.py 中的 sync)
        couple.sync()

        # D. 计算导数
        dy1 = obj1.get_derivatives(t, y1)
        dy2 = obj2.get_derivatives(t, y2)

        return np.concatenate([dy1, dy2])

    # --- 6. 求解 ---
    t_end = 3000.0
    y0 = np.concatenate([obj1.T, obj2.T])
    solver = NuclearODESolver(method='BDF', rtol=1e-6, atol=1e-8)

    res = solver.solve(coupled_rhs, (0, t_end), y0, t_eval=np.arange(0, t_end + 1, 1.0))

    if not res['success']:
        raise RuntimeError(f"Solver failed: {res['message']}")

    return mesh1, mesh2, res['t'], res['y'][:N1, :], res['y'][N1:, :]


# ==============================================================================
if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # --- 运行笛卡尔 ---
    m1_c, m2_c, t_c, y1_c, y2_c = run_case_3_new('cartesian')
    x_comb_c, y_comb_c = merge_results_for_plot(m1_c, m2_c, y1_c, y2_c, 'cartesian')
    plot_coupled_spacetime(axes[0], x_comb_c, t_c, y_comb_c, "Cartesian (Coupled)")

    # --- 运行柱坐标 ---
    m1_r, m2_r, t_r, y1_r, y2_r = run_case_3_new('cylindrical')
    x_comb_r, y_comb_r = merge_results_for_plot(m1_r, m2_r, y1_r, y2_r, 'cylindrical')
    plot_coupled_spacetime(axes[1], x_comb_r, t_r, y_comb_r, "Cylindrical (Coupled)")

    plt.tight_layout()
    plt.show()

    # --- 验证打印 ---
    print("\n--- 结果验证 (t=300s) ---")
    print(f"[Cartesian] Max Temp (Obj1 Left): {y1_c[0, -1]:.2f} K")
    print(f"[Cartesian] Interface Temp: {y1_c[-1, -1]:.2f} K (Obj1) / {y2_c[0, -1]:.2f} K (Obj2)")
    print(f"[Cartesian] Min Temp (Obj2 Right): {y2_c[-1, -1]:.2f} K")

    print(f"\n[Cylindrical] Max Temp (Obj1 Left): {y1_r[0, -1]:.2f} K")
    print(f"[Cylindrical] Interface Temp: {y1_r[-1, -1]:.2f} K / {y2_r[0, -1]:.2f} K")

    # 物理合理性检查:
    # 界面处由于无接触热阻，温度应连续 (Numerical上可能会有极小差异，取决于网格密度)
    diff = abs(y1_c[-1, -1] - y2_c[0, -1])
    print(f"\nInterface Delta T (Cartesian): {diff:.4e} K (Should be close to 0)")

    # # 创建画布：1行2列
    # fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    #
    # # --- 1. 运行并绘制笛卡尔坐标算例 ---
    # mesh_cart, t_cart, y_cart = run_case_2('cartesian')
    # plot_spacetime_contour(axes[0], mesh_cart, t_cart, y_cart, "Cartesian (Plate)")
    #
    # # --- 2. 运行并绘制柱坐标算例 ---
    # mesh_cyl, t_cyl, y_cyl = run_case_2('cylindrical')
    # plot_spacetime_contour(axes[1], mesh_cyl, t_cyl, y_cyl, "Cylindrical (Pipe/Rod)")
    #
    # plt.tight_layout()
    # plt.show()
    #
    # # 简单的数值验证打印
    # print("\n--- 结果验证 ---")
    # print(f"[Cartesian] Time steps: {len(t_cart)} (Expected 101)")
    # print(f"[Cartesian] T_end distribution (Left -> Right):")
    # print(f"  Left Node:  {y_cart[0, -1]:.2f} K")
    # print(f"  Right Node: {y_cart[-1, -1]:.2f} K")
    # print(f"[Cylindrical] T_end distribution (Left -> Right):")
    # print(f"  Left Node:  {y_cyl[0, -1]:.2f} K")
    # print(f"  Right Node: {y_cyl[-1, -1]:.2f} K")
