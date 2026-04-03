import numpy as np
import matplotlib.pyplot as plt

# 导入底层组件
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import IncompressibleFluidChannel, FlowJunction
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Solvers.SystemManager import SystemManager
from Components.AnnularPipe import AnnularPipe

# 导入时间统计组件
from profiler import TEASAProfiler

# 导入换热关联式
from Correlations.Correlations import nu_ringpipe

@TEASAProfiler.profile
def run_annular_pipes_test():
    print("==========================================================")
    print("   TEST: 5 Parallel Annular Pipes (Fixed dP, K_loss, FFTF)   ")
    print("==========================================================")

    # ---------------------------------------------------------
    # 1. 几何与工况参数
    # ---------------------------------------------------------
    L_pipe = 0.375  # 管道长度 [m]

    # 内套管 (Inner Solid)
    r_in_inner = 0.0119  # 内半径
    r_out_inner = 0.01225  # 外半径
    thick_inner = r_out_inner - r_in_inner

    # 外套管 (Outer Solid)
    r_in_outer = 0.01295  # 内半径
    r_out_outer = 0.01430  # 外半径-原0.01330
    thick_outer = r_out_outer - r_in_outer

    # 环形流道 (Fluid Channel)
    gap_width = r_in_outer - r_out_inner
    area_flow = np.pi * (r_in_outer ** 2 - r_out_inner ** 2)
    D_hydraulic = 2.0 * gap_width

    # 加热/润湿周长
    perim_inner = 2 * np.pi * r_out_inner  # 流体接触的内侧周长
    perim_outer = 2 * np.pi * r_in_outer  # 流体接触的外侧周长

    N_z = 30  # 轴向网格数
    N_r = 1  # 径向厚度网格数

    T_init = 600.0  # 初始温度 [K]
    P_out = 1.58e5  # 出口压力 [Pa]
    P_in = P_out + 3000.0  # 入口压力 (压差 200 Pa)

    # 并联通道非对称参数
    powers = [10, 20, 30, 40, 50]  # 内侧加热总功率 [W] (直接使用)
    # powers = [0.0, 0.0, 0.0, 0.0, 0.0]
    k_loss_list = [1.5, 2.0, 2.5, 3.0, 3.5]  # 出口局部阻力系数

    n_pipes = 5

    # ==========================================================
    # 【新增代码】定义一个接口适配器 (Adapter Function)
    # 它接收底层的 (Re, Pr, P_D)，内部调用您真实的 4 参数函数
    # ==========================================================
    def nu_ringpipe_adapter(Re, Pr, pd_dummy=1.1):
        # 将局部几何参数 R_out 和 R_in 闭包捕获，传给真实的换热关系式
        return nu_ringpipe(R_out=r_in_outer, R_in=r_out_inner, Re=Re, Pr=Pr)

    mat_fluid = Sodium()
    mat_solid = AusteniticStainlessSteel()

    # ---------------------------------------------------------
    # 2. 构建公共进出口联箱 (Plenums)
    # ---------------------------------------------------------
    inlet_plenum = IncompressibleBoundaryVolume("Inlet_Plenum", mat_fluid, P=P_in, T=T_init)
    outlet_plenum = IncompressibleBoundaryVolume("Outlet_Plenum", mat_fluid, P=P_out, T=T_init)

    # 硬性锚定为定压边界，避免奇异矩阵
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum.is_pressure_boundary = True

    all_vols = [inlet_plenum, outlet_plenum]
    all_juncs = []
    inlet_junctions = []
    annular_pipes = []

    # ---------------------------------------------------------
    # 3. 批量生成 5 组环形管道装配体
    # ---------------------------------------------------------
    for i in range(n_pipes):
        idx = i + 1

        # --- A. 流体支路 ---
        chan = IncompressibleFluidChannel(
            name=f"Chan_{idx}",
            n_nodes=N_z, total_length=L_pipe,
            flow_area=area_flow, hydraulic_diam=D_hydraulic,
            initial_P=P_out, initial_T=T_init, material=mat_fluid
        )
        all_vols.extend(chan.volumes)
        all_juncs.extend(chan.internal_junctions)

        # --- B. 进出口 Junction (出口添加局部阻力) ---
        j_in = FlowJunction(f"J_in_{idx}", inlet_plenum, chan.volumes[0], flow_area=area_flow, k_loss=0.0)
        j_out = FlowJunction(f"J_out_{idx}", chan.volumes[-1], outlet_plenum, flow_area=area_flow,
                             k_loss=k_loss_list[i])

        all_juncs.extend([j_in, j_out])
        inlet_junctions.append(j_in)

        # --- C. 内套管 (Inner Solid) ---
        mesh_in = Mesh2D(x_dim=thick_inner, n_x=N_r, y_dim=L_pipe, n_y=N_z,
                         geometry_type='cylindrical', inner_radius=r_in_inner)
        solid_in = HeatConduction2D(mesh=mesh_in, material=mat_solid, initial_temp=T_init)
        solid_in.name = f"Solid_Inner_{idx}"

        # 内侧 (left) 加热：直接给定总功率
        solid_in.boundaries['left'].add_flux_condition(q_flux=powers[i])
        solid_in.boundaries['top'].add_flux_condition(q_flux=0.0)
        solid_in.boundaries['bottom'].add_flux_condition(q_flux=0.0)

        # --- D. 外套管 (Outer Solid) ---
        mesh_out = Mesh2D(x_dim=thick_outer, n_x=N_r, y_dim=L_pipe, n_y=N_z,
                          geometry_type='cylindrical', inner_radius=r_in_outer)
        solid_out = HeatConduction2D(mesh=mesh_out, material=mat_solid, initial_temp=T_init)
        solid_out.name = f"Solid_Outer_{idx}"

        # 外侧 (right) 绝热
        solid_out.boundaries['right'].add_flux_condition(q_flux=0.0)
        solid_out.boundaries['top'].add_flux_condition(q_flux=0.0)
        solid_out.boundaries['bottom'].add_flux_condition(q_flux=0.0)

        # --- E. 使用 AnnularPipe 宏观组件一键打包 ---
        annular_comp = AnnularPipe(
            name=f"AnnularPipe_{idx}",
            solid_inner=solid_in,
            solid_outer=solid_out,
            fluid_channel=chan,
            heated_perimeter_inner=perim_inner,
            heated_perimeter_outer=perim_outer,
            correlation_func=nu_ringpipe_adapter,
            boundary_inner_solid='right',  # 内套管的外表面接触流体
            boundary_outer_solid='left'  # 外套管的内表面接触流体
        )
        annular_pipes.append(annular_comp)

    # ---------------------------------------------------------
    # 4. 系统组装与初始化
    # ---------------------------------------------------------
    net = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    net.initialize_hydraulics(dt=0.01, tol=1e-5, max_iter=200)

    sys_mgr = SystemManager(fluid_network=net)
    for p in annular_pipes:
        sys_mgr.add_component(p)

    # 系统初始化
    sys_mgr.initialize_system()

    # ---------------------------------------------------------
    # 5. 瞬态执行
    # ---------------------------------------------------------
    print("\n>> Starting Transient Simulation...")
    t_end = 1.0

    # 1. 循环外：初始化 history 字典
    history = {'time': [], 'flows': [], 'temps': []}

    while sys_mgr.global_time < t_end:
        if sys_mgr.global_time > 25.0:
            inlet_plenum.set_boundary_state(P=P_out+3500)
        elif sys_mgr.global_time > 20.0:
            inlet_plenum.set_boundary_state(P=P_out+3400)
        elif sys_mgr.global_time > 15.0:
            inlet_plenum.set_boundary_state(P=P_out+3300)
        elif sys_mgr.global_time > 10.0:
            inlet_plenum.set_boundary_state(P=P_out+3200)
        elif sys_mgr.global_time > 5.0:
            inlet_plenum.set_boundary_state(P=P_out+3100)
        dt = sys_mgr.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=2.0)
        sys_mgr.step(dt=dt, inner_iter=2)

        # 2. 循环内：每次推演后记录当前时刻的数据
        history['time'].append(sys_mgr.global_time)
        history['flows'].append([j.W for j in inlet_junctions])
        history['temps'].append([p.fluid_channel.volumes[-1].T for p in annular_pipes])

        if int(sys_mgr.global_time / dt) % 20 == 0:
            current_total_W = sum([j.W for j in inlet_junctions])
            print(f" t = {sys_mgr.global_time:>5.2f}s | Delta_P = {inlet_plenum.P-outlet_plenum.P} Pa | Total Flow = {current_total_W:.5f} kg/s")

    # ---------------------------------------------------------
    # 6. 打印最终稳态结果
    # ---------------------------------------------------------
    print("\n==========================================================")
    print("   Final Steady-State Results at t = {:.2f} s".format(sys_mgr.global_time))
    print("==========================================================")
    for i, p in enumerate(annular_pipes):
        T_f_out = p.fluid_channel.volumes[-1].T
        T_inner_wall_avg = np.mean(p.solid_inner.boundaries['right'].T_surface)
        T_outer_wall_avg = np.mean(p.solid_outer.boundaries['left'].T_surface)
        j_in = inlet_junctions[i]

        print(f"[{p.name}] - Input Power: {powers[i]:.1f} W | Outlet K_loss: {k_loss_list[i]:.1f}")
        print(f"  - Mass Flow Rate         : {j_in.W:.5f} kg/s")
        print(f"  - Fluid Outlet Temp      : {T_f_out:.2f} K")
        print(f"  - Inner Wall (wet) Temp  : {T_inner_wall_avg:.2f} K")
        print(f"  - Outer Wall (wet) Temp  : {T_outer_wall_avg:.2f} K\n")

    plot_transient_results(history)

def plot_transient_results(history):
    """
    绘制并联通道的瞬态流量与温度变化曲线

    :param history: 包含历史数据的字典，必须包含以下键：
        - 'time': 一维列表，记录仿真时间点 [s]
        - 'flows': 二维列表或数组，形状为 (时间步数, 通道数)，记录各通道质量流量 [kg/s]
        - 'temps': 二维列表或数组，形状为 (时间步数, 通道数)，记录各通道出口流体温度 [K]
    """
    times = np.array(history['time'])
    flows = np.array(history['flows'])
    temps = np.array(history['temps'])

    # 获取通道数量
    n_pipes = flows.shape[1] if len(flows.shape) > 1 else 0

    # 创建上下两个子图，共享横坐标
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # ==========================================
    # 子图 1: 流量变化曲线
    # ==========================================
    ax1 = axes[0]
    ax1.set_title('Transient Mass Flow Rate Distribution', fontsize=14, fontweight='bold')

    # 为不同通道使用不同的颜色和标记
    for i in range(n_pipes):
        ax1.plot(times, flows[:, i], linewidth=2, label=f'Pipe {i + 1}')

    ax1.set_ylabel('Mass Flow Rate [kg/s]', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))  # 将图例放在图表右侧

    # ==========================================
    # 子图 2: 出口温度变化曲线
    # ==========================================
    ax2 = axes[1]
    ax2.set_title('Transient Fluid Outlet Temperature', fontsize=14, fontweight='bold')

    for i in range(n_pipes):
        ax2.plot(times, temps[:, i], linewidth=2, linestyle='-', label=f'Pipe {i + 1}')

    ax2.set_xlabel('Time [s]', fontsize=12)
    ax2.set_ylabel('Fluid Outlet Temp [K]', fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    # 调整布局并保存/显示
    plt.tight_layout()
    plt.savefig('annular_pipes_transient.png', dpi=300, bbox_inches='tight')
    print(">> 瞬态结果曲线已成功保存为 'annular_pipes_transient.png'")

    TEASAProfiler.report()


if __name__ == "__main__":
    run_annular_pipes_test()

