import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys

from profiler import TEASAProfiler

# 1. 获取当前测试脚本所在的绝对路径 (.../TASTIN_Project/Tests)
current_dir = os.path.dirname(os.path.abspath(__file__))

# 2. 推导出 ThermoCalc 模块所在的绝对路径
thermo_calc_dir = os.path.join(current_dir, '..', 'ThermoCalc')
thermo_calc_dir = os.path.abspath(thermo_calc_dir)

# 3. 将该路径加入系统搜索路径
if thermo_calc_dir not in sys.path:
    sys.path.insert(0, thermo_calc_dir)

# =============================================================================
# 导入您的 TASTIN 系统组件 (请确保路径与您的工程一致)
# =============================================================================
from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import FluidVolume, FlowJunction, IncompressibleFluidChannel
from Solvers.Hydrodynamics.BoundaryVolume import BoundaryVolume

from Components.basicComponents.TECPair import TECPair
from Components.TECCircuitManager import TECCircuitManager

# 导入随便一个流体和固体物性 (用于测试)
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.MoNb import MoNb

from Materials.Fluids.Sodium import Sodium
from Solvers.Hydrodynamics.Components import IncompressibleFluidVolume, FlowJunction
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume

# =============================================================================
# 1. 测试数据准备 (50 个轴向节点的固定边界温度)
# =============================================================================
TE_data = np.array([
    1491.29, 1500.11, 1516.00, 1537.34, 1562.78, 1591.18, 1621.63, 1653.35, 1685.71, 1718.16,
    1750.24, 1781.57, 1811.82, 1840.71, 1868.00, 1893.49, 1917.00, 1938.40, 1957.56, 1974.40,
    1988.83, 2000.80, 2010.27, 2017.21, 2021.62, 2023.48, 2022.83, 2019.69, 2014.09, 2006.08,
    1995.73, 1983.11, 1968.29, 1951.39, 1932.50, 1911.77, 1889.32, 1865.34, 1840.01, 1813.57,
    1786.27, 1758.44, 1730.45, 1702.77, 1675.95, 1650.68, 1627.80, 1608.39, 1593.75, 1585.55
])

TC_data = np.array([
    760.774, 762.065, 763.823, 765.810, 767.951, 770.213, 772.576, 775.027, 777.551, 780.139,
    782.779, 785.462, 788.178, 790.920, 793.677, 796.445, 799.211, 801.974, 804.720, 807.451,
    810.151, 812.824, 815.452, 818.044, 820.578, 823.068, 825.488, 827.858, 830.144, 832.379,
    834.517, 836.603, 838.580, 840.509, 842.317, 844.084, 845.715, 847.323, 848.779, 850.236,
    851.526, 852.851, 853.995, 855.227, 856.262, 857.463, 858.453, 859.715, 860.732, 862.062
])


def create_dummy_fluid_network() -> HydraulicNetwork:
    """
    构建哑流体网络，完全隔离，只为了满足 SystemManager 的依赖
    [修改后]：采用不可压缩体积模型，并明确指定压力边界。
    """

    mat_fluid = Sodium()

    # --- 1. 几何与工况参数 (严格对齐 test_coupled_heating.py) ---
    L_channel = 0.375
    D_inner = 0.010
    N_fluid = 15
    T_init = 600.0
    P_out = 1.58e5
    dP_drive = 150.0
    P_in = P_out + dP_drive
    area_flow = np.pi * (D_inner / 2) ** 2

    # --- 2. 构建压力边界 ---
    inlet_plenum = IncompressibleBoundaryVolume("Dummy_Inlet", mat_fluid, P=P_in, T=T_init)
    inlet_plenum.is_pressure_boundary = True  # 声明为定压点

    outlet_plenum = IncompressibleBoundaryVolume("Dummy_Outlet", mat_fluid, P=P_out, T=T_init)
    outlet_plenum.is_pressure_boundary = True  # 声明为定压点

    # --- 3. 构建流体通道 (包含 30 个控制体) ---
    dummy_chan = IncompressibleFluidChannel(
        name="Dummy_Chan",
        n_nodes=N_fluid,
        total_length=L_channel,
        flow_area=area_flow,
        hydraulic_diam=D_inner,
        initial_P=P_out,
        initial_T=T_init,
        material=mat_fluid
    )

    # --- 4. 构建进出口连接 ---
    j_in = FlowJunction("J_In", inlet_plenum, dummy_chan.volumes[0], flow_area=area_flow)
    j_out = FlowJunction("J_Out", dummy_chan.volumes[-1], outlet_plenum, flow_area=area_flow)

    # --- 5. 组装网络 ---
    all_vols = [inlet_plenum, outlet_plenum] + dummy_chan.volumes
    # 注意：Channel 内部会自动生成 internal_junctions
    all_juncs = [j_in, j_out] + dummy_chan.internal_junctions

    # 重力设置为 0.0 (与测试用例保持一致)
    return HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)


def main():
    print("=== TASTIN Thermo-Electric Coupling System Test ===")

    # ---------------------------------------------------------
    # 1. 组装 6 个热离子电极对并强加温度边界
    # ---------------------------------------------------------
    n_elem = 6
    n_node = 50
    total_L = 0.377  # 真实物理长度 [m]
    L_node = total_L / n_node

    tfe_list = []
    mat_e = MoNb()
    mat_c = Molybdenum()

    for i in range(n_elem):
        # 实例化单根 TFE
        tfe = TECPair(
            name=f"TFE_{i}",
            L_node=L_node, n_node=n_node,
            R_e_in=0.008, delta_e=0.0015, delta_gap=0.0003, delta_c=0.0015,  # 示例几何
            n_rad_e=3, n_rad_c=3,
            mat_emitter=mat_e, mat_collector=mat_c,
            T_init_e=1800.0, T_init_c=800.0
        )

        # [核心测试逻辑]：利用极其微小的接触热阻 (R=1e-10) 强制固定内外边界温度
        tfe.inner_boundary.add_resistance_condition(T_ext=TE_data, R_ext=1e-10)
        tfe.outer_boundary.add_resistance_condition(T_ext=TC_data, R_ext=1e-10)

        tfe_list.append(tfe)
        print(f"[*] TECPair '{tfe.name}' initialized with forced boundary conditions.")

    # ---------------------------------------------------------
    # 2. 挂载到电热宏观电路管理器
    # ---------------------------------------------------------
    circuit_manager = TECCircuitManager(
        name="Core_Circuit",
        tfe_list=tfe_list,
        Tcs_init=610.0,  # 铯池温度设定
        V_init=0.2
    )

    # ---------------------------------------------------------
    # 3. 初始化 SystemManager (搭载 Dummy 流体)
    # ---------------------------------------------------------
    dummy_net = create_dummy_fluid_network()
    sys_manager = SystemManager(fluid_network=dummy_net, start_time=0.0)

    # 将包含 12 个固体、6 个间隙的电路代理器一次性注册进去
    sys_manager.add_component(circuit_manager)

    # 初始化液力 (Dummy 管道会瞬间收敛)
    sys_manager.initialize_system(dt_init=0.1)

    # ---------------------------------------------------------
    # 4. 瞬态主循环 (Transient Integration)
    # ---------------------------------------------------------
    t_end = 2.0  # 仿真总时长 [s]
    dt = 0.05  # 显式步进推荐使用小步长

    time_history = []
    I_history = []
    U_history = []

    print("\nStarting transient integration...")
    start_cpu = time.time()

    while sys_manager.global_time <= t_end + 1e-6:
        # 系统总调度：推进流体、固体、同步电热、点堆(如果有)
        sys_manager.step(dt=dt, inner_iter=1)  # 建议开启 2 次 Picard 迭代收敛导热

        current_t = sys_manager.global_time

        # 提取系统级宏观电学参数记录
        macro_res = circuit_manager.circuit.get_global_results()
        if macro_res:
            time_history.append(current_t)
            I_history.append(macro_res['Iout'])
            U_history.append(macro_res['Uout'])

        if len(time_history) % 10 == 0:
            print(f"  Time: {current_t:.3f} s | I_out: {macro_res['Iout']:.2f} A | U_out: {macro_res['Uout']:.3f} V")

    print(f"\nSimulation finished in {time.time() - start_cpu:.2f} seconds.")

    # ---------------------------------------------------------
    # 5. 可视化后处理 (Post-Processing)
    # ---------------------------------------------------------
    # A. 宏观电学曲线
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(time_history, I_history, 'b-', label='Total Current (I)')
    plt.xlabel('Time [s]')
    plt.ylabel('Current [A]')
    plt.grid(True);
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(time_history, U_history, 'r-', label='Output Voltage (U)')
    plt.xlabel('Time [s]')
    plt.ylabel('Voltage [V]')
    plt.grid(True);
    plt.legend()
    plt.tight_layout()
    plt.show()

    # B. 提取最后一根 TFE (TFE_5) 的轴向微观参数
    res_tfe5 = circuit_manager.circuit.get_tec_results(5)
    if res_tfe5:
        z_axis = np.linspace(0, total_L, n_node)

        plt.figure(figsize=(10, 5))
        plt.plot(z_axis, res_tfe5['J'], 'g.-', label='Current Density J [A/cm2]')
        plt.xlabel('Axial Position [m]')
        plt.ylabel('Current Density')
        plt.title('TFE_5 Axial Current Density at t_end')
        plt.grid(True);
        plt.legend()
        plt.show()


if __name__ == '__main__':
    main()
    TEASAProfiler.report()
