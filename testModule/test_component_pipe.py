import numpy as np
import matplotlib.pyplot as plt

# 导入底层组件
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import IncompressibleFluidChannel, FlowJunction
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Fluids.Sodium import Sodium         # 导入 Sodium 作为 NaK
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Solvers.SystemManager import SystemManager
from Components.Pipe import Pipe

# 换热关联式 (Lyon-Martinelli for Liquid Metals)

from Correlations.Correlations import nu_fftf, nu_aoki

def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
    Pe = np.maximum(Re * Pr, 1.0)
    return 7.0 + 0.025 * (Pe ** 0.8)

def run_parallel_pipes_test():
    print("==========================================================")
    print("   TEST: 5 Parallel Pipes with Flow Distribution & CHT    ")
    print("==========================================================")

    # ---------------------------------------------------------
    # 1. 几何与工况参数
    # ---------------------------------------------------------
    L_pipe = 0.375                           # 管道长度 [m]
    D_inner = 0.01                           # 内径 [m]
    D_outer = 0.018                          # 外径 [m]
    thickness = (D_outer - D_inner) / 2.0    # 壁厚 [m] (0.003m)
    area_flow = np.pi * (D_inner / 2)**2     # 流通面积
    heated_perim = np.pi * D_inner           # 润湿/加热周长

    N_z = 30                                 # 轴向网格数
    N_r = 1                                  # 径向网格数

    T_init = 600.0                           # 初始温度 [K]
    P_out = 1.58e5                           # 出口压力 [Pa]
    W_target = 0.25                          # 目标总流量 [kg/s]

    # 并联通道非对称参数
    T_outer_list = [610.0, 620.0, 630.0, 640.0, 650.0]
    k_loss_list = [0.3, 0.25, 0.2, 0.15, 0.1]
    n_pipes = 5

    # ---------------------------------------------------------
    # 2. 实例化材料
    # ---------------------------------------------------------
    mat_fluid = Sodium()  # NaK 替代
    mat_solid = AusteniticStainlessSteel()  # 不锈钢管壁

    # ---------------------------------------------------------
    # 3. 构建公共进出口联箱 (Plenums)
    # ---------------------------------------------------------
    # 给定一个初始的猜测入口压力 (稍后会由系统自动调节以满足 0.25kg/s)
    P_in_guess = P_out + 500.0

    inlet_plenum = IncompressibleBoundaryVolume("Inlet_Plenum", mat_fluid, P=P_in_guess, T=T_init)
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum = IncompressibleBoundaryVolume("Outlet_Plenum", mat_fluid, P=P_out, T=T_init)
    inlet_plenum.is_pressure_boundary = True

    all_vols = [inlet_plenum, outlet_plenum]
    all_juncs = []
    inlet_junctions = []  # 记录入口 Junction 用于统计总流量
    pipes = []

    # ---------------------------------------------------------
    # 4. 批量生成 5 组并联通道与 Pipe 组件
    # ---------------------------------------------------------
    for i in range(n_pipes):
        idx = i + 1

        # --- A. 流体支路 ---
        chan = IncompressibleFluidChannel(
            name=f"Channel_{idx}",
            n_nodes=N_z, total_length=L_pipe,
            flow_area=area_flow, hydraulic_diam=D_inner,
            initial_P=P_out, initial_T=T_init, material=mat_fluid
        )
        all_vols.extend(chan.volumes)
        all_juncs.extend(chan.internal_junctions)

        # --- B. 进出口 Junction (添加局部阻力) ---
        j_in = FlowJunction(f"J_in_{idx}", inlet_plenum, chan.volumes[0], flow_area=area_flow, k_loss=0.0)
        j_out = FlowJunction(f"J_out_{idx}", chan.volumes[-1], outlet_plenum, flow_area=area_flow,
                             k_loss=k_loss_list[i])

        all_juncs.extend([j_in, j_out])
        inlet_junctions.append(j_in)

        # --- C. 固体管壁 ---
        mesh = Mesh2D(x_dim=thickness, n_x=N_r, y_dim=L_pipe, n_y=N_z,
                      geometry_type='cylindrical', inner_radius=D_inner / 2.0)
        solid = HeatConduction2D(mesh=mesh, material=mat_solid, initial_temp=T_init)

        # 外边界 (right) 定温：采用极小热阻模拟
        # solid.boundaries['right'].add_resistance_condition(T_ext=T_outer_list[i], R_ext=1e-5)
        solid.boundaries['right'].add_flux_condition(q_flux=50.0)
        # 上下绝热 (Top / Bottom)
        solid.boundaries['top'].add_flux_condition(q_flux=0.0)
        solid.boundaries['bottom'].add_flux_condition(q_flux=0.0)

        # --- D. 使用 Pipe 宏观组件一键打包 ---
        pipe_comp = Pipe(
            name=f"Pipe_{idx}",
            solid_wall=solid,
            fluid_channel=chan,
            heated_perimeter=heated_perim,
            correlation_func=lyon_martinelli,
            coupled_boundary_name='left'  # 流体接触固体内壁
        )

        pipes.append(pipe_comp)

    # ---------------------------------------------------------
    # 5. 系统组装与初始化
    # ---------------------------------------------------------
    net = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    net.initialize_hydraulics(dt=0.5, tol=1e-5, max_iter=200)

    sys_mgr = SystemManager(fluid_network=net)
    for p in pipes:
        sys_mgr.add_component(p)

    # ---------------------------------------------------------
    # 6. 瞬态执行 (引入虚拟 P-Controller 控制总流量)
    # ---------------------------------------------------------
    print("\n>> Starting Transient Simulation with Flow Control...")
    t_end = 15.0  # 跑 15 秒虚拟瞬态以达到稳态

    while sys_mgr.global_time < t_end:
        # [流控黑科技] P-Controller: 调整入口压力以锁定 0.25 kg/s （不启用，还是使用）
        current_total_W = sum([j.W for j in inlet_junctions])
        # error_W = W_target - current_total_W
        # 比例增益 K_p，视系统敏感度可调
        # inlet_plenum.P += 100000.0 * error_W

        # 计算自适应步长并执行推演
        dt = sys_mgr.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=1.0)
        sys_mgr.step(dt=dt)

        # 打印监控日志 (每模拟 1 秒打印一次)
        if int(sys_mgr.global_time / dt) % 20 == 0:
            print(f" t = {sys_mgr.global_time:>5.2f}s | P_in = {inlet_plenum.P:>10.1f} Pa | Total Flow = {current_total_W:.4f} kg/s")

    # ---------------------------------------------------------
    # 7. 打印最终稳态结果
    # ---------------------------------------------------------
    print("\n==========================================================")
    print("   Final Steady-State Results at t = {:.2f} s".format(sys_mgr.global_time))
    print("==========================================================")
    for i, p in enumerate(pipes):
        T_f_out = p.fluid_channel.volumes[-1].T
        # 提取整根管子的内壁面平均温度
        T_w_in_avg = np.mean(p.solid.boundaries['left'].T_surface)
        j_in = inlet_junctions[i]

        print(f"[{p.name}]")
        print(f"  - Resistance (K_loss)  : {k_loss_list[i]:.2f}")
        print(f"  - Outer Temp Boundary  : {T_outer_list[i]:.1f} K")
        print(f"  - Mass Flow Rate       : {j_in.W:.4f} kg/s")
        print(f"  - Fluid Outlet Temp    : {T_f_out:.1f} K")
        print(f"  - Wall Inner Temp(Avg) : {T_w_in_avg:.1f} K\n")


if __name__ == "__main__":
    run_parallel_pipes_test()
