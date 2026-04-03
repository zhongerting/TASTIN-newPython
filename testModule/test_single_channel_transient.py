import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Materials.Fluids.Sodium import Sodium
from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, IncompressibleFluidChannel
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume


def run_test_case():
    print("==========================================================")
    print("   TEST: Single Heated Channel Transient (Sodium)         ")
    print("==========================================================")

    # 1. 物理参数定义
    # ----------------------------------------------------------
    mat = Sodium()

    P_out = 1.6e5  # 出口压力 150 kPa
    dP_drive = 480.0  # 驱动压差 120 Pa
    P_in = P_out + dP_drive

    T_init = 600.0  # 初始温度 600 K

    # 几何参数 (水平管)
    L_channel = 0.5  # m
    D_channel = 0.01  # m (10 mm)
    N_nodes = 50
    A_flow = np.pi * (D_channel / 2) ** 2

    print(f"Fluid: Sodium")
    print(f"Geometry: L={L_channel}m, D={D_channel}m, Nodes={N_nodes}")
    print(f"Boundary: Pin={P_in:.1f} Pa, Pout={P_out:.1f} Pa, dT={dP_drive} Pa")

    # 2. 组件创建
    # ----------------------------------------------------------

    # A. 边界 (Boundaries)
    # 使用 IncompressibleBoundaryVolume 作为 Plenum
    # 关键: 设置 is_pressure_boundary = True，告诉 Solver 这是 Dirichlet 边界
    inlet_plenum = IncompressibleBoundaryVolume("Inlet", mat, P=P_in, T=T_init)
    inlet_plenum.is_pressure_boundary = True

    outlet_plenum = IncompressibleBoundaryVolume("Outlet", mat, P=P_out, T=T_init)
    outlet_plenum.is_pressure_boundary = True

    # B. 通道 (Channel)
    channel = IncompressibleFluidChannel(
        name="TestChan",
        n_nodes=N_nodes,
        total_length=L_channel,
        flow_area=A_flow,
        hydraulic_diam=D_channel,
        initial_P=P_out,  # 初始压力给出口压力，稍后修正
        initial_T=T_init,
        material=mat
    )

    # 3. 功能演示: 根据给定流量推导压力分布 (Initialization Helper)
    # ----------------------------------------------------------
    print("\n[Step 1] Demonstration: Calculate Pressure from Guess Flow")

    # 假设一个猜测流量 (例如 0.05 kg/s)
    W_guess = 0.05
    print(f"   -> Guessing flow rate: {W_guess} kg/s")

    # 手动将通道内部连接的流量设为猜测值
    for junc in channel.internal_junctions:
        junc.W = W_guess

    # 调用 FluidChannel 的方法从上游向下游递推压力
    # 这会更新 channel.volumes 中每个节点的 P
    channel.update_pressure_distribution_downstream(P_inlet=P_in)

    # 检查末端压力是否与 P_out 一致 (由于流量是瞎猜的，肯定不一致)
    P_end_calc = channel.volumes[-1].P
    print(f"   -> Calculated End P: {P_end_calc:.1f} Pa")
    print(f"   -> Target End P:     {P_out:.1f} Pa")
    print(f"   -> Mismatch indicates that W={W_guess} is not the hydraulic solution.")

    # 4. 网络组装 (Topology Assembly)
    # ----------------------------------------------------------
    # 手动创建连接 Channel 与 Boundary 的 Junction
    # 注意: Inlet -> Channel[0]
    j_inlet = FlowJunction("J_Inlet", inlet_plenum, channel.volumes[0], flow_area=A_flow)

    # 注意: Channel[-1] -> Outlet
    j_outlet = FlowJunction("J_Outlet", channel.volumes[-1], outlet_plenum, flow_area=A_flow)

    # 收集所有的 Volume 和 Junction
    all_volumes = [inlet_plenum] + channel.volumes + [outlet_plenum]
    all_junctions = [j_inlet] + channel.internal_junctions + [j_outlet]

    # 实例化求解器
    # 假设为水平管，重力加速度设为 0 以验证纯摩擦压降
    net = HydraulicNetwork(all_volumes, all_junctions, gravity_vector=0.0)

    # 5. 稳态初始化 (Hydraulic Initialization)
    # ----------------------------------------------------------
    print("\n[Step 2] Solver Initialization (Finding True Flow)")

    # 调用我们在 Phase 3 完成的初始化方法
    # 它会自动调整流量和压力，直到满足动量方程 (压差 = 阻力)
    success = net.initialize_hydraulics(dt=1, tol=1e-6, max_iter=100)

    if not success:
        print("Error: Initialization failed!")
        return

    W_steady = j_inlet.W
    print(f"   -> Converged Flow Rate: {W_steady:.6f} kg/s")

    # 简单验证: dP = f * (L/D) * (v^2/2) * rho
    # 取中间节点物性
    vol_mid = channel.volumes[N_nodes // 2]
    rho = vol_mid.rho
    mu = vol_mid.mu
    vel = W_steady / (rho * A_flow)
    Re = rho * vel * D_channel / mu

    # 获取摩擦系数 (调用内部静态方法)
    f = HydraulicNetwork._calc_friction_factor_static(Re)

    # 估算压降 (Channel 长度)
    # 注意：总压降还包括 inlet/outlet junction 的惯性长度部分，
    # 但由于 plenum 很大，惯性长度主要由 channel 贡献，这里做近似验证
    dP_est = f * (L_channel / D_channel) * 0.5 * rho * vel ** 2

    print(f"   -> Check Physics:")
    print(f"      Re = {Re:.1f}, f = {f:.5f}")
    print(f"      Est. Friction dP = {dP_est:.2f} Pa (Approx)")
    print(f"      Actual Drive dP  = {dP_drive:.2f} Pa")

    # 6. 瞬态计算 (Transient Loop)
    # ----------------------------------------------------------
    print("\n[Step 3] Transient Run (0-5s, Heating at t=2s)")

    dt = 0.005
    t_end = 10.0
    current_time = 0.0

    # 记录数据
    history = {
        'time': [],
        'W_in': [],
        'W_out': [],
        'T_out': [],
        'P_mid': []
    }

    # 加热功率: 100 kW 总功率 -> 分配到 10 个节点 -> 10 kW/node
    # Q_wall 正值代表加热
    heat_power_per_node = 50.0  # W
    heat_dist = np.ones(N_nodes) * heat_power_per_node

    step_count = 0

    while current_time < t_end:
        step_count += 1
        current_time += dt

        # A. 施加边界条件/源项
        # 必须在每步开始前清除上一布的源项 (重要!)
        channel.clear_sources()

        if current_time > 2.0:
            # t > 2s 施加阶跃加热
            channel.add_heat_source_distribution(heat_dist)

        # B. 执行一步计算
        # 使用 Picard 迭代保证非线性收敛
        net.step_Picard(dt, max_iter=100, tol=1e-5)

        # C. 记录数据
        history['time'].append(current_time)
        history['W_in'].append(j_inlet.W)
        history['W_out'].append(j_outlet.W)
        # 记录最后一个流体节点的温度
        history['T_out'].append(channel.volumes[-1].T)
        history['P_mid'].append(channel.volumes[N_nodes // 2].P)

        if step_count % 100 == 0:
            print(f"   t = {current_time:.2f} s | T_out = {channel.volumes[-1].T:.2f} K | W = {j_inlet.W:.5f} kg/s")

    print("Transient simulation completed.")

    # 7. 结果绘图
    # ----------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:red'
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Temperature Out (K)', color=color)
    ax1.plot(history['time'], history['T_out'], color=color, linewidth=2, label='T_out')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:blue'
    ax2.set_ylabel('Mass Flow Rate (kg/s)', color=color)
    # 稍微偏移一下线型以便区分入口和出口流量
    ax2.plot(history['time'], history['W_in'], color=color, linestyle='--', label='W_in')
    ax2.plot(history['time'], history['W_out'], color='green', linestyle=':', label='W_out', alpha=0.7)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.legend(loc='upper right')

    plt.title('Single Channel Transient: Step Heating 100kW at t=2.0s')
    fig.tight_layout()
    plt.show()

    # 额外：绘制最终的压力分布
    # 手动收集 z 坐标用于绘图 (Inlet/Outlet 设为虚拟坐标)
    z_coords = [v.z_coordinate for v in channel.volumes]
    z_plot = [-0.1] + z_coords + [L_channel + 0.1]
    p_plot = [inlet_plenum.P] + [v.P for v in channel.volumes] + [outlet_plenum.P]

    plt.figure()
    plt.plot(z_plot, p_plot, 'o-', markersize=4)
    plt.xlabel('Position (m)')
    plt.ylabel('Pressure (Pa)')
    plt.title('Final Pressure Distribution along Channel')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    run_test_case()
