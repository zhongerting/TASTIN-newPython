import numpy as np
import matplotlib.pyplot as plt
import logging
import sys
import os

# 设置日志显示
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# --- 引入项目模块 ---
# try:
#     from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
#     from Materials.Fluids.Sodium import Sodium
#     from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, IncompressibleFluidChannel
#     from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
# except ImportError as e:
#     print(f"Import Error: {e}")
#     print("请确保脚本运行在包含 Materials 和 Solvers 文件夹的根目录下。")
#     sys.exit(1)

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
from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WickMaterial import WickMaterial
from Components.RingHP import RingHP


def run_parallel_test():
    print("==========================================================")
    print("   TEST: Parallel Channels Flow Distribution (Sodium)     ")
    print("==========================================================")

    # 1. 物理参数与几何定义
    # ----------------------------------------------------------
    mat = SodiumPotassium78()

    # 边界条件: 恒定压差驱动
    P_out = 1.5e5  # 出口压力 150 kPa
    dP_drive = 775  # 驱动压差 775 Pa
    P_in = P_out + dP_drive

    T_init = 600.0  # 初始温度 600 K

    # 通道 A (低阻通道 - 主流道)
    L_A = 1.0  # m
    D_A = 0.01  # 10 mm
    N_A = 10
    A_A = np.pi * (D_A / 2) ** 2

    # 通道 B (高阻通道 - 细管/旁路)
    L_B = 1.0  # m
    D_B = 0.02  # 8 mm (更细 -> 阻力更大)
    N_B = 10
    A_B = np.pi * (D_B / 2) ** 2

    print(f"Boundary Conditions:")
    print(f"  Pin  = {P_in:.1f} Pa")
    print(f"  Pout = {P_out:.1f} Pa")
    print(f"  dP   = {dP_drive:.1f} Pa")
    print(f"Channel Configuration:")
    print(f"  Channel A (Low R) : L={L_A}m, D={D_A * 1000}mm")
    print(f"  Channel B (High R): L={L_B}m, D={D_B * 1000}mm")

    # 2. 组件构建
    # ----------------------------------------------------------

    # A. 边界 (Boundaries) - 公用的 Inlet 和 Outlet Plenum
    inlet_plenum = IncompressibleBoundaryVolume("Inlet", mat, P=P_in, T=T_init)
    inlet_plenum.is_pressure_boundary = True  # 标记为定压边界

    outlet_plenum = IncompressibleBoundaryVolume("Outlet", mat, P=P_out, T=T_init)
    outlet_plenum.is_pressure_boundary = True  # 标记为定压边界

    # B. 通道 (Channels)
    channel_a = IncompressibleFluidChannel(
        name="ChanA", n_nodes=N_A, total_length=L_A,
        flow_area=A_A, hydraulic_diam=D_A,
        initial_P=P_out, initial_T=T_init, material=mat
    )

    channel_b = IncompressibleFluidChannel(
        name="ChanB", n_nodes=N_B, total_length=L_B,
        flow_area=A_B, hydraulic_diam=D_B,
        initial_P=P_out, initial_T=T_init, material=mat
    )

    # C. 连接 (Junctions) - 手动构建拓扑
    # 需要创建 4 个特定的 Junction 来连接 Plenum 和 Channel

    # Channel A 的连接
    j_in_a = FlowJunction("J_In_A", inlet_plenum, channel_a.volumes[0], flow_area=A_A)
    j_out_a = FlowJunction("J_Out_A", channel_a.volumes[-1], outlet_plenum, flow_area=A_A)

    # Channel B 的连接
    j_in_b = FlowJunction("J_In_B", inlet_plenum, channel_b.volumes[0], flow_area=A_B)
    j_out_b = FlowJunction("J_Out_B", channel_b.volumes[-1], outlet_plenum, flow_area=A_B)

    # 3. 网络组装
    # ----------------------------------------------------------
    # 收集所有对象
    all_volumes = [inlet_plenum, outlet_plenum] + channel_a.volumes + channel_b.volumes

    # 注意 Junction 的顺序不影响计算，但包含了内部和外部连接
    all_junctions = [j_in_a, j_out_a, j_in_b, j_out_b] + \
                    channel_a.internal_junctions + \
                    channel_b.internal_junctions

    # 实例化求解器 (水平布置 g=0 以专注摩擦分配)
    net = HydraulicNetwork(all_volumes, all_junctions, gravity_vector=0.0)

    # 4. 稳态初始化 (Stage 1)
    # ----------------------------------------------------------
    print("\n[Step 1] Steady-State Initialization...")

    # 使用较大的 dt 进行伪瞬态松弛，快速找到稳态
    success = net.initialize_hydraulics(dt=0.5, tol=1e-7, max_iter=200)

    if not success:
        print("Initialization Failed! Check input parameters.")
        return

    W_A_init = j_in_a.W
    W_B_init = j_in_b.W

    print(f"   -> Converged Flow A: {W_A_init:.6f} kg/s")
    print(f"   -> Converged Flow B: {W_B_init:.6f} kg/s")
    print(f"   -> Flow Ratio (A/B): {W_A_init / W_B_init:.2f}")

    # 理论检查: 粗略估算流量比
    # 假设湍流 dP ~ W^2 / D^4.75，则 W ~ D^2.375
    # (D_A/D_B)^(2.375) = (10/8)^2.375 ≈ 1.7
    # 实际层流/湍流混合可能略有不同，但 A 应该显著大于 B
    if W_A_init <= W_B_init:
        print("WARNING: Flow distribution logic might be wrong (A should be > B)")

    # 5. 非对称加热瞬态 (Stage 2)
    # ----------------------------------------------------------
    print("\n[Step 2] Transient: Asymmetric Heating (Heat A ONLY)")
    print("   -> t = 0~2s:  Adiabatic")
    print("   -> t = 2~10s: Channel A heated with 5 kW (Channel B adiabatic)")

    dt = 0.01
    t_end = 10.0
    time = 0.0

    # 加热配置
    power_A = 5000.0  # W
    heat_dist_A = np.ones(N_A) * (power_A / N_A)

    history = {
        'time': [],
        'W_A': [], 'W_B': [],
        'T_out_A': [], 'T_out_B': []
    }

    step_count = 0
    while time < t_end:
        step_count += 1
        time += dt

        # A. 清除上一部源项
        channel_a.clear_sources()
        channel_b.clear_sources()

        # B. 施加非对称加热
        if time > 2.0:
            channel_a.add_heat_source_distribution(heat_dist_A)
            # Channel B 保持绝热 (无操作)

        # C. 求解一步
        # Picard 迭代处理非线性阻力系数变化
        net.step_Picard(dt, max_iter=20, tol=1e-5)

        # D. 记录数据
        history['time'].append(time)
        history['W_A'].append(j_in_a.W)
        history['W_B'].append(j_in_b.W)
        history['T_out_A'].append(channel_a.volumes[-1].T)
        history['T_out_B'].append(channel_b.volumes[-1].T)

        if step_count % 100 == 0:
            print(f"   t={time:.2f}s | WA={j_in_a.W:.5f} | WB={j_in_b.W:.5f} "
                  f"| TA_out={channel_a.volumes[-1].T:.1f} K | Wall={j_in_a.W+j_in_b.W:.5f}")

    # 6. 绘图分析
    # ----------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # 流量图
    ax1.set_title('Parallel Channels Flow Redistribution Test')
    ax1.plot(history['time'], history['W_A'], label='Channel A (Heated, D=10mm)', color='tab:red', linewidth=2)
    ax1.plot(history['time'], history['W_B'], label='Channel B (Adiabatic, D=8mm)', color='tab:blue', linewidth=2)
    ax1.set_ylabel('Mass Flow Rate (kg/s)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 温度图
    ax2.plot(history['time'], history['T_out_A'], label='T_out Channel A', color='tab:red', linestyle='--')
    ax2.plot(history['time'], history['T_out_B'], label='T_out Channel B', color='tab:blue', linestyle='--')
    ax2.set_ylabel('Outlet Temperature (K)')
    ax2.set_xlabel('Time (s)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_parallel_test()
