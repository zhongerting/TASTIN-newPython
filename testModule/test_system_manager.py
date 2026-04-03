import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

# --- 模拟项目路径引入 ---
# 假设脚本运行在项目根目录下
try:
    from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
    from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, IncompressibleFluidChannel
    from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
    from Solvers.SystemManager import SystemManager
    from Materials.Fluids.Sodium import Sodium
except ImportError as e:
    print(f"Import Error: {e}")
    sys.exit(1)

# 配置日志输出，方便观察求解器收敛情况
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- 基础参数 ---
L = 0.5  # 通道长度 [m]
D_in = 0.02  # 水力直径 [m]
P_out = 1.6e5  # 出口压力 [Pa]
P_drop = 80.0  # 目标驱动压差 [Pa]
P_in = P_out + P_drop
T_in = 600.0  # 初始温度 [K]
n_axi = 25  # 节点数


# --- 辅助验证函数: 理论摩擦系数计算 ---
def calc_theoretical_friction(Re):
    """复刻 HydraulicNetwork 中的摩擦系数逻辑进行验证"""
    if Re <= 1000.0:
        return 64.0 / max(Re, 1e-5)
    elif Re < 100000.0:
        return 0.3164 / (max(Re, 1e-5) ** 0.25)
    else:
        # 简化处理，高雷诺数仅作近似对比
        return 0.3164 / (max(Re, 1e-5) ** 0.25)


def AdiabaticFlow_test():
    print("==========================================================")
    print("      TEST 1: Adiabatic Flow (Momentum Verification)      ")
    print("==========================================================")

    mat_fluid = Sodium()
    flow_area = np.pi * (D_in / 2) ** 2

    # 1. 构建流体通道
    channel1 = IncompressibleFluidChannel(
        name="channel1",
        n_nodes=n_axi,
        total_length=L,
        flow_area=flow_area,
        hydraulic_diam=D_in,
        initial_P=P_out,  # 初始压力给出口压力，等待求解器建立梯度
        initial_T=T_in,
        material=mat_fluid
    )

    # 2. 构建边界 Plenum
    inlet = IncompressibleBoundaryVolume("inlet", mat_fluid, P=P_in, T=T_in)
    inlet.is_pressure_boundary = True  # 标记为定压边界

    outlet = IncompressibleBoundaryVolume("outlet", mat_fluid, P=P_out, T=T_in)
    outlet.is_pressure_boundary = True  # 标记为定压边界

    # 3. 构建连接 (Junctions)
    # 注意：流道入口/出口连接如果不指定 custom_length，会自动计算惯性长度
    j_in_a = FlowJunction("J_In_A", inlet, channel1.volumes[0], flow_area=flow_area)
    j_out_a = FlowJunction("J_Out_A", channel1.volumes[-1], outlet, flow_area=flow_area)

    # 4. 组装网络
    all_vols = [inlet, outlet] + channel1.volumes
    all_juncs = [j_in_a, j_out_a] + channel1.internal_junctions

    # 绝热测试：重力设为 0，单纯验证压差与摩擦的平衡
    net = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)

    # 5. 系统管理器
    sys_mgr = SystemManager(fluid_network=net)

    # --- 关键步骤：初始化 ---
    # 这会先运行稳态水力迭代，消除初始时刻 W=0 带来的非物理震荡
    print(">> Initializing System (Hydraulic Steady State)...")
    sys_mgr.initialize_system(dt_init=0.1, tol=1e-8, max_iter=500)

    # --- 瞬态循环 ---
    t_end = 50.0  # 绝热流动很快稳定，不需要跑太久
    dt = 0.01  # 时间步
    current_time = 0.0

    # 数据记录
    history = {
        'time': [],
        'W_in': [],  # 入口流量
        'W_out': [],  # 出口流量
        'P_in_node': [],  # 通道第一个节点压力
        'P_out_node': [],  # 通道最后一个节点压力
        'T_out': []  # 出口温度
    }

    print(f">> Starting Transient Loop (Target dP = {P_drop} Pa)...")

    while current_time < t_end:
        current_time += dt

        # 执行一步计算
        # 绝热算例不需要耦合器和固体，只运行流体即可
        # 但 SystemManager.step 会自动处理
        sys_mgr.step(dt)

        # 记录
        history['time'].append(current_time)
        history['W_in'].append(j_in_a.W)
        history['W_out'].append(j_out_a.W)
        history['P_in_node'].append(channel1.volumes[0].P)
        history['P_out_node'].append(channel1.volumes[-1].P)
        history['T_out'].append(channel1.volumes[-1].T)

    print(">> Simulation Completed.")

    # ==========================================================================
    # 专家验证：物理自洽性检查
    # ==========================================================================
    print("\n[VERIFICATION REPORT]")

    # 1. 获取最终稳定状态
    W_final = history['W_in'][-1]
    T_final = history['T_out'][-1]

    # 2. 计算流体物性 (基于 T_final)
    rho = mat_fluid.density(T_final, P_out)
    mu = mat_fluid.viscosity(T_final, P_out)

    # 3. 计算流速与无量纲数
    vel = W_final / (rho * flow_area)
    Re = (rho * vel * D_in) / mu

    # 4. 计算理论压降
    # Darcy-Weisbach: dP = f * (L/D) * (rho * v^2 / 2)
    f_theory = calc_theoretical_friction(Re)

    # 动态压头 (Dynamic Head)
    dyn_head = 0.5 * rho * vel ** 2

    # 理论摩擦压降
    dP_friction_theory = f_theory * (L / D_in) * dyn_head

    # 5. 误差分析
    # 注意：Code 中的压降是施加在 Plenum 上的，中间经过了 Inlet Junction, Channel, Outlet Junction
    # 总压降 = P_in - P_out = 600 Pa
    # 理论计算应该尽量接近这个值 (假设局部阻力为0)
    dP_boundary = P_in - P_out
    error_pct = abs(dP_friction_theory - dP_boundary) / dP_boundary * 100.0

    print(f"1. Flow Stability Check:")
    print(f"   Final Mass Flow : {W_final:.6f} kg/s")
    print(f"   Inlet vs Outlet : {abs(history['W_in'][-1] - history['W_out'][-1]):.2e} kg/s (Should be ~0)")

    print(f"2. Temperature Stability Check:")
    print(f"   Initial Temp    : {T_in:.2f} K")
    print(f"   Final Temp      : {T_final:.2f} K")
    print(f"   Drift           : {abs(T_final - T_in):.4f} K (Should be 0 for adiabatic)")

    print(f"3. Momentum Balance Check (The 'Expert' Test):")
    print(f"   Reynolds Number : {Re:.2f}")
    print(f"   Friction Factor : {f_theory:.6f} (Blasius/Theoretical)")
    print(f"   Velocity        : {vel:.4f} m/s")
    print(f"   ------------------------------------------------")
    print(f"   Target Boundary dP : {dP_boundary:.2f} Pa")
    print(f"   Calc. Friction dP  : {dP_friction_theory:.2f} Pa (Based on simulated W)")
    print(f"   ------------------------------------------------")
    print(f"   Error              : {error_pct:.2f}%")

    if error_pct < 1.0:
        print(">> RESULT: PASS (Momentum Equation Verified)")
    else:
        print(">> RESULT: WARNING (Check friction correlation or junction losses)")

    # --- 绘图 ---
    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    # Plot 1: Flow Rate
    axes[0].plot(history['time'], history['W_in'], 'b-', label='Inlet Flow')
    axes[0].plot(history['time'], history['W_out'], 'r--', label='Outlet Flow')
    axes[0].set_title('Mass Flow Rate Stability')
    axes[0].set_ylabel('Mass Flow [kg/s]')
    axes[0].legend()
    axes[0].grid(True)

    # Plot 2: Pressure Distribution (Snapshot at end)
    # 重构压力分布数组
    P_dist = [inlet.P] + [v.P for v in channel1.volumes] + [outlet.P]
    x_dist = np.linspace(0, L, len(P_dist))  # 简化坐标

    axes[1].plot(x_dist, P_dist, 'k-o')
    axes[1].set_title('Pressure Distribution (End of Simulation)')
    axes[1].set_xlabel('Position (approx) [m]')
    axes[1].set_ylabel('Pressure [Pa]')
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    AdiabaticFlow_test()