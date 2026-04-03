import numpy as np
import matplotlib.pyplot as plt
import logging
import sys

# --- 模拟项目路径引入 ---
try:
    from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
    from Solvers.Hydrodynamics.Components import FluidChannel, FlowJunction, IncompressibleFluidChannel
    from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume

    from Solvers.HeatConduction.HeatConduction import HeatConduction2D
    from Solvers.HeatConduction.Mesh import Mesh2D
    from Solvers.HeatConduction.Boundary import FluxBC, ResistanceBC

    from Solvers.Couplers import FluidSolidCouple
    from Solvers.SystemManager import SystemManager

    from Materials.Fluids.Sodium import Sodium
    from Materials.Solids.StainlessSteel import AusteniticStainlessSteel

except ImportError as e:
    print(f"Import Error: {e}")
    sys.exit(1)

# 配置日志输出，方便观察 dt 的变化
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# ==============================================================================
# 0. 辅助类定义
# ==============================================================================

def lyon_martinelli_correlation(Re, Pr, P_D_ratio):
    """液态金属 Nusselt 数经验关联式 (Lyon-Martinelli)"""
    Pe = Re * Pr
    Pe = np.maximum(Pe, 1.0)
    return 7.0 + 0.025 * (Pe ** 0.8)


# ==============================================================================
# 1. 测试主程序
# ==============================================================================

def run_coupled_test():
    print("==========================================================")
    print("   TEST: Asymmetric Dual-Pipe Coupled Heating (Adaptive DT) ")
    print("==========================================================")

    # --- A. 几何与工况参数 ---
    L_channel = 0.375
    D_inner = 0.010
    D_outer = 0.014
    N_fluid = 15
    N_solid_r = 1
    N_solid_z = N_fluid

    T_init = 600.0
    P_out = 1.58e5
    dP_drive = 150.0
    P_in = P_out + dP_drive

    mat_fluid = Sodium()
    mat_solid = AusteniticStainlessSteel()

    # --- B. 构建流体网络 ---
    print(">> Building Hydraulic Network...")
    inlet_plenum = IncompressibleBoundaryVolume("Inlet", mat_fluid, P=P_in, T=T_init)
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum = IncompressibleBoundaryVolume("Outlet", mat_fluid, P=P_out, T=T_init)
    outlet_plenum.is_pressure_boundary = True

    area_flow = np.pi * (D_inner / 2) ** 2

    chan_a = IncompressibleFluidChannel(
        name="Chan_A", n_nodes=N_fluid, total_length=L_channel,
        flow_area=area_flow, hydraulic_diam=D_inner,
        initial_P=P_out, initial_T=T_init, material=mat_fluid
    )

    chan_b = IncompressibleFluidChannel(
        name="Chan_B", n_nodes=N_fluid, total_length=L_channel,
        flow_area=area_flow, hydraulic_diam=D_inner,
        initial_P=P_out, initial_T=T_init, material=mat_fluid
    )

    j_in_a = FlowJunction("J_In_A", inlet_plenum, chan_a.volumes[0], flow_area=area_flow)
    j_out_a = FlowJunction("J_Out_A", chan_a.volumes[-1], outlet_plenum, flow_area=area_flow)
    j_in_b = FlowJunction("J_In_B", inlet_plenum, chan_b.volumes[0], flow_area=area_flow)
    j_out_b = FlowJunction("J_Out_B", chan_b.volumes[-1], outlet_plenum, flow_area=area_flow)

    all_vols = [inlet_plenum, outlet_plenum] + chan_a.volumes + chan_b.volumes
    all_juncs = [j_in_a, j_out_a, j_in_b, j_out_b] + chan_a.internal_junctions + chan_b.internal_junctions

    net = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)

    print("   -> Initializing flow...")
    net.initialize_hydraulics(dt=0.5, tol=1e-6, max_iter=100)

    # --- C. 构建固体域 ---
    print(">> Building Solid Domains...")
    radial_thickness = (D_outer - D_inner) / 2.0
    inner_radius = D_inner / 2.0

    mesh_a = Mesh2D(x_dim=radial_thickness, n_x=N_solid_r, y_dim=L_channel, n_y=N_solid_z,
                    geometry_type='cylindrical', inner_radius=inner_radius)
    mesh_b = Mesh2D(x_dim=radial_thickness, n_x=N_solid_r, y_dim=L_channel, n_y=N_solid_z,
                    geometry_type='cylindrical', inner_radius=inner_radius)

    solid_a = HeatConduction2D(mesh=mesh_a, material=mat_solid, initial_temp=T_init)
    solid_a.name = "Pipe_A"
    solid_b = HeatConduction2D(mesh=mesh_b, material=mat_solid, initial_temp=T_init)
    solid_b.name = "Pipe_B"

    # [修改点 1] 强制初始化固体状态，并提取热容数据
    # -------------------------------------------------------------------------
    print("   -> Pre-calculating solid thermal capacitance...")
    # 这里的 initialize_state 已经在 __init__ 中被调用过一次，但为了保险起见，
    # 同时也为了显式展示逻辑，我们确保状态是最新的。
    solid_a.initialize_state()
    solid_b.initialize_state()

    # 提取与流体接触面 ('left' = 内表面) 的节点热容
    # 这将用于计算耦合稳定性限制: dt < C / (h*A)
    cap_a = solid_a.get_boundary_node_capacitance('left')
    cap_b = solid_b.get_boundary_node_capacitance('left')
    # -------------------------------------------------------------------------

    # 设置边界条件
    # Pipe A: 50 kW/m2, Pipe B: 2 kW/m2
    solid_a.boundaries['right'].add_flux_condition(q_flux=50.0)
    solid_b.boundaries['right'].add_flux_condition(q_flux=20.0)

    solid_a.boundaries['top'].add_flux_condition(0.0)
    solid_a.boundaries['bottom'].add_flux_condition(0.0)
    solid_b.boundaries['top'].add_flux_condition(0.0)
    solid_b.boundaries['bottom'].add_flux_condition(0.0)

    # --- D. 构建耦合 ---
    print(">> Building Couplers...")
    heated_perimeter = np.pi * D_inner

    # [修改点 2] 传入 solid_node_capacitance
    # -------------------------------------------------------------------------
    coupler_a = FluidSolidCouple(
        name="Coup_A", fluid=chan_a, solid_boundary_region=solid_a.boundaries['left'],
        heated_perimeter=heated_perimeter, correlation_func=lyon_martinelli_correlation,
        solid_node_capacitance=cap_a  # <--- 注入热容数据
    )
    coupler_b = FluidSolidCouple(
        name="Coup_B", fluid=chan_b, solid_boundary_region=solid_b.boundaries['left'],
        heated_perimeter=heated_perimeter, correlation_func=lyon_martinelli_correlation,
        solid_node_capacitance=cap_b  # <--- 注入热容数据
    )
    # -------------------------------------------------------------------------

    # --- E. 系统组装 ---
    print(">> Assembling System...")
    sys_mgr = SystemManager(fluid_network=net)
    sys_mgr.add_solid_component(solid_a)
    sys_mgr.add_solid_component(solid_b)
    # sys_mgr.add_coupler(coupler_a)
    # sys_mgr.add_coupler(coupler_b)

    sys_mgr.initialize_system()

    # --- F. 瞬态循环 (自适应步长) ---
    print("\n[Start Adaptive Transient Simulation]")

    t_end = 10.0
    current_time = 0.0

    # 历史记录
    history = {'time': [], 'dt': [],
               'T_f_out_A': [], 'T_f_out_B': [],
               'T_w_out_A': [], 'T_w_out_B': [],
               'W_A': [], 'W_B': []}

    step_count = 0

    while current_time < t_end:
        step_count += 1

        # [修改点 3] 计算自适应时间步长
        # -------------------------------------------------------------------------
        # min_dt: 防止步长过小卡死
        # max_dt: 防止步长过大导致物理过程失真 (例如 0.1s 内功率突变)
        # safety_factor: 0.8 是比较稳健的选择
        dt = sys_mgr.compute_adaptive_dt(min_dt=1e-5, max_dt=0.1, safety_factor=1)

        # dt = 0.01
        # -------------------------------------------------------------------------

        # [修改点 4] 执行带预估-校正的步进
        # -------------------------------------------------------------------------
        # inner_iter=2: 开启 Predictor-Corrector
        # 这对于大幅度动态步长至关重要，能保证二阶精度和稳定性
        sys_mgr.step(dt, inner_iter=1)
        # -------------------------------------------------------------------------

        # 更新当前时间 (直接读取 SystemManager 的时钟，保证同步)
        current_time = sys_mgr.global_time

        # 记录数据
        history['time'].append(current_time)
        history['dt'].append(dt)  # 记录 dt 变化曲线

        history['T_f_out_A'].append(chan_a.volumes[-1].T)
        history['T_f_out_B'].append(chan_b.volumes[-1].T)

        T_surf_A = np.mean(solid_a.boundaries['right'].T_surface)
        T_surf_B = np.mean(solid_b.boundaries['right'].T_surface)
        history['T_w_out_A'].append(T_surf_A)
        history['T_w_out_B'].append(T_surf_B)

        history['W_A'].append(j_in_a.W)
        history['W_B'].append(j_in_b.W)

        # 每隔一定步数打印状态
        if step_count % 50 == 0:
            print(f"   t={current_time:.4f}s | dt={dt:.1e}s | "
                  f"TA_out={chan_a.volumes[-1].T:.1f}K | "
                  f"TB_out={chan_b.volumes[-1].T:.1f}K | "
                  f"WA={j_in_a.W:.4f}")

    print("Simulation Completed.")

    # --- G. 结果绘图 ---
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # Plot 1: Time Step Evolution (新图表)
    ax = axes[0]
    ax.set_title('Adaptive Time Step Evolution')
    ax.plot(history['time'], history['dt'], 'g-', linewidth=1.5, label='dt [s]')
    ax.set_ylabel('dt [s]')
    ax.set_yscale('log')  # 对数坐标看数量级变化
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Fluid Temperatures
    ax = axes[1]
    ax.set_title('Asymmetric Heating: Fluid Outlet Temperatures')
    ax.plot(history['time'], history['T_f_out_A'], 'r-', linewidth=2, label='Channel A (50 kW/m2)')
    ax.plot(history['time'], history['T_f_out_B'], 'b--', linewidth=2, label='Channel B (2 kW/m2)')
    ax.set_ylabel('Fluid Temp [K]')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Solid Wall Temperatures
    ax = axes[2]
    ax.set_title('Solid Outer Wall Temperatures')
    ax.plot(history['time'], history['T_w_out_A'], 'r-', linewidth=2, label='Wall A')
    ax.plot(history['time'], history['T_w_out_B'], 'b--', linewidth=2, label='Wall B')
    ax.set_ylabel('Wall Temp [K]')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Mass Flow Rates
    ax = axes[3]
    ax.set_title('Mass Flow Rate')
    ax.plot(history['time'], history['W_A'], 'r-', label='Flow A')
    ax.plot(history['time'], history['W_B'], 'b--', label='Flow B')
    ax.set_ylabel('Flow [kg/s]')
    ax.set_xlabel('Time [s]')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_coupled_test()
