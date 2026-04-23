import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

# from profiler import TEASAProfiler

# =============================================================================
# 0. 路径与系统导入配置 (与您的 test_TEC 保持一致)
# =============================================================================
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

# 导入 TASTIN 系统组件
from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import IncompressibleFluidChannel, FlowJunction
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume

# 导入您原生材料库
from Materials.Fluids.Sodium import Sodium

# 导入我们刚刚编写的带翅片热管组件
from Components.HPwithFin import HPwithFin

# =============================================================================
# [核心修改] 导入您真实的物性材料库
# 请根据您真实的包路径(如 Materials.WallMaterial) 进行调整
# =============================================================================
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WickMaterial import WickMaterial
from Materials.Fluids.Sodium import Sodium  # 仅用于 Dummy 液力网络的常规钠物性

# =============================================================================
# 1. 哑流体网络构造器 (专为 SystemManager 启动服务)
# =============================================================================
def create_dummy_fluid_network() -> HydraulicNetwork:
    """构建极简哑流体网络，只为了满足 SystemManager 的底层依赖"""
    mat_fluid = Sodium()

    L_channel = 0.1
    D_inner = 0.010
    N_fluid = 3
    T_init = 600.0
    P_out = 1.58e5
    dP_drive = 100.0
    P_in = P_out + dP_drive
    area_flow = np.pi * (D_inner / 2) ** 2

    inlet_plenum = IncompressibleBoundaryVolume("Dummy_Inlet", mat_fluid, P=P_in, T=T_init)
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum = IncompressibleBoundaryVolume("Dummy_Outlet", mat_fluid, P=P_out, T=T_init)
    outlet_plenum.is_pressure_boundary = True

    dummy_chan = IncompressibleFluidChannel(
        name="Dummy_Chan", n_nodes=N_fluid, total_length=L_channel, flow_area=area_flow,
        hydraulic_diam=D_inner, initial_P=P_out, initial_T=T_init, material=mat_fluid
    )

    j_in = FlowJunction("J_In", inlet_plenum, dummy_chan.volumes[0], flow_area=area_flow)
    j_out = FlowJunction("J_Out", dummy_chan.volumes[-1], outlet_plenum, flow_area=area_flow)

    all_vols = [inlet_plenum, outlet_plenum] + dummy_chan.volumes
    all_juncs = [j_in, j_out] + dummy_chan.internal_junctions

    return HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)

# =============================================================================
# 2. 主测试程序
# =============================================================================
def main():
    print("=== TASTIN HPwithFin 真实物性级测试 ===")

    # ---------------------------------------------------------
    # 2.1 物理与几何参数设定
    # ---------------------------------------------------------
    T_init = 800.0  # 热管初始均匀温度
    T_eva_ext = 800.0  # 蒸发段强制加热边界
    T_env = 3.0  # 太空环境温度
    up_vf = 0.0  # 内侧向上角系数
    down_vf = 0.3  # 内侧向下角系数
    emissivity = 0.93  # 发射率

    # 轴向尺寸划分
    L_eva = 0.06
    n_eva = 1
    L_aba = 0.0
    n_aba = 0
    L_con = 0.482
    n_con = 12

    # 径向尺寸划分
    r_out_wall = 0.0085  # 8.5 mm
    r_in_wall = r_out_wall - 0.0004  # 壁厚 0.4 mm
    r_vapor = r_in_wall - 0.0006  # 吸液芯厚 0.6 mm -> 7.5mm
    porosity = 0.5  # 孔隙率
    n_wick = 1  # 吸液芯 1 格
    n_wall = 1  # 管壁 1 格

    fin_height = 22.65e-3  # 22.65 mm
    fin_thickness = 0.0003  # 0.3 mm
    n_fin = 2  # 对称双翅片

    # 几何包覆率折算
    fin_wrap_ratio = (n_fin * fin_thickness) / (2.0 * np.pi * r_out_wall)

    # ---------------------------------------------------------
    # 2.2 [核心] 实例化真实的固体材料和工质
    # ---------------------------------------------------------
    # 1. 热管管壁使用 316 不锈钢
    mat_wall = SS316(name="HP_Wall_SS316")

    # 2. 热管内部相变工质使用钠 (NaHP)
    mat_fluid = SodiumHP(name="HP_Fluid_Na")

    # 3. 复合吸液芯：316 不锈钢骨架 + 钠工质
    mat_wick = WickMaterial(
        name="HP_Wick_Composite",
        solid_mat=SS316(),  # 丝网骨架材质
        fluid_mat=mat_fluid,  # 灌注的工质
        porosity=porosity,  # 0.5
        r_vapor=r_vapor,  # 传入内部蒸汽腔半径计算赝导热
        r_in_wall=r_in_wall  # 传入管壁内半径计算赝导热
    )

    # ---------------------------------------------------------
    # 2.3 实例化 HPwithFin 并注入极小热阻加热边界
    # ---------------------------------------------------------
    print("\n[*] 正在构建装载真实物性的 HPwithFin 辐射器组件...")
    hp_radiator = HPwithFin(
        name="Main_Radiator",
        r_out_wall=r_out_wall, r_in_wall=r_in_wall, r_vapor=r_vapor,
        L_eva=L_eva, L_aba=L_aba, L_con=L_con,
        n_eva=n_eva, n_aba=n_aba, n_con=n_con,
        n_wick=n_wick, n_wall=n_wall,
        wall_mat=mat_wall, fluid_mat=mat_fluid, wick_struct_mat=mat_wick, porosity=porosity,
        fin_thickness=fin_thickness, fin_height=fin_height, n_fin_height=15,
        fin_wrap_ratio=fin_wrap_ratio,
        emissivity=emissivity,
        up_view_factor=up_vf, down_view_factor=down_vf,
        T_env=T_env,
        initial_temp=T_init
    )

    # 用接触热阻边界强制定温
    hp_radiator.hp.boundaries['outer_eva'].add_resistance_condition(
        T_ext=T_eva_ext,
        R_ext=1e-8
    )

    # ---------------------------------------------------------
    # 2.4 注册到 SystemManager 并执行积分
    # ---------------------------------------------------------
    dummy_net = create_dummy_fluid_network()
    sys_manager = SystemManager(fluid_network=dummy_net, start_time=0.0)
    sys_manager.add_component(hp_radiator)
    sys_manager.initialize_system(dt_init=0.05)

    t_end = 10.0
    dt = 0.05

    time_history = []
    q_con_history = []
    temp_contour_history = []  # <--- [新增] 初始化空列表用于存储壁温演化
    temp_line_history = []

    print("\n" + "=" * 50)
    print("🚀 开始基于真实物性的瞬态隐式积分...")
    print("=" * 50)
    start_cpu = time.time()

    # 为了避免画太多线，我们可以设置一个记录间隔，例如每 10 步记录一次温度场
    record_interval = 200
    step_count = 0

    while sys_manager.global_time <= t_end + 1e-6:
        sys_manager.step(dt=dt, inner_iter=1)

        current_t = sys_manager.global_time
        _, q_con_total_dist = hp_radiator.get_heat_rejection_distribution()
        total_q = np.sum(q_con_total_dist)

        # 每隔指定步数记录一次温度分布，避免生成过多的图
        step_count += 1
        if step_count % record_interval == 0:
            T_2d = hp_radiator.get_temperature_distribution()
            T_outer_wall = T_2d[-1, :]  # 提取最外层管壁温度
            temp_line_history.append(T_outer_wall)
            time_history.append(current_t)

        time_history.append(current_t)
        q_con_history.append(total_q)

        if len(time_history) % 10 == 0:
            print(f"  时间: {current_t:.3f} s | 冷凝段总散热: {total_q:.2f} W")

    # ---------------------------------------------------------
    # 2.5 可视化后处理
    # ---------------------------------------------------------
    _, q_con_final_dist = hp_radiator.get_heat_rejection_distribution()
    z_nodes = np.linspace(1, n_con, n_con)

    plt.figure(figsize=(10, 5))
    plt.bar(z_nodes, q_con_final_dist, color='coral', edgecolor='black')
    plt.xlabel('Condenser Node Index')
    plt.ylabel('Heat Rejection (W)')
    plt.title(f'Heat Rejection Distribution with Real SS316 & NaHP at t={t_end}s')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

    # 2. 绘制轴向温度演化线图
    print("\n[*] 正在生成外壁温轴向分布线图...")
    z_centers = hp_radiator.hp_mesh.y_centers  # 获取轴向坐标
    plot_axial_temp_lines(time_history, z_centers, temp_line_history)

def plot_axial_temp_lines(time_history, z_centers, temp_history_2d):
    """
    绘制热管轴向外壁面温度随时间演化的线图。
    每 6 个时刻画一张图，超过则新建图。

    :param time_history: 记录的时间点列表
    :param z_centers: 轴向网格中心坐标列表
    :param temp_history_2d: 对应的外壁温度二维列表 (N_t, N_z)
    """
    n_times = len(time_history)
    lines_per_plot = 6

    # 计算需要画几张图
    n_plots = int(np.ceil(n_times / lines_per_plot))

    for i in range(n_plots):
        start_idx = i * lines_per_plot
        end_idx = min((i + 1) * lines_per_plot, n_times)

        plt.figure(figsize=(10, 6))

        for j in range(start_idx, end_idx):
            t = time_history[j]
            T_axial = temp_history_2d[j]
            plt.plot(z_centers, T_axial, marker='.', label=f't = {t:.2f} s')

        plt.xlabel('Axial Position Z (m)', fontsize=12)
        plt.ylabel('Outer Wall Temperature (K)', fontsize=12)
        plt.title(f'Axial Temperature Distribution (Plot {i + 1}/{n_plots})', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.show()


def plot_axial_temp_contour(time_history, z_centers, temp_history_2d):
    """
    绘制热管轴向外壁面温度随时间演化的云图
    :param time_history: 一维时间数组 (N_t,)
    :param z_centers: 一维轴向坐标数组 (N_z,)
    :param temp_history_2d: 二维温度数组 (N_t, N_z)
    """
    # 将列表转换为 NumPy 数组
    time_array = np.array(time_history)
    z_array = np.array(z_centers)
    temp_array = np.array(temp_history_2d)

    # 创建网格 (Z, T)
    Z, T_time = np.meshgrid(z_array, time_array)

    plt.figure(figsize=(10, 6))
    # 使用 contourf 绘制填充等高线云图，levels=100 保证色彩过渡极其平滑
    contour = plt.contourf(Z, T_time, temp_array, levels=100, cmap='inferno')

    # 添加颜色条
    cbar = plt.colorbar(contour)
    cbar.set_label('Outer Wall Temperature (K)', fontsize=12)

    plt.xlabel('Axial Position Z (m)', fontsize=12)
    plt.ylabel('Time (s)', fontsize=12)
    plt.title('Transient Spatiotemporal Contour of HP Outer Wall Temperature', fontsize=14)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
    # TEASAProfiler.report()

