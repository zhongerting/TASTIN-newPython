import numpy as np
import matplotlib.pyplot as plt

# --- 导入底层求解器与组件 ---
from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidChannel

# --- 导入材料库 ---
from Materials.Solids.UO2 import UO2
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.ZrH import ZirconiumHydride
from Materials.Solids.GasGaps import Xenon, Cesium, CarbonDioxide, Helium
from Materials.Fluids.Sodium import Sodium  # 假设液态金属为 NaK78

# --- 导入 TFE 组件与数据结构 ---
from Components.TFEUnit import TFEUnit, TFEGeometry, TFEMeshParams, GapConfig
from Components.ReactorCore import ReactorCore
from profiler import TEASAProfiler


def run_test():
    # =========================================================================
    # 1. 核心参数设置 (请填写您的实际几何数据)
    # =========================================================================
    H = 0.375  # TFE 有效高度 [m]
    N_axial = 30  # 轴向网格数

    # 1.1 几何参数 (需填写真实半径)
    geometry = TFEGeometry(
        r_pellet_inner=4e-3,
        r_pellet_outer=8.5e-3,
        r_fission_gas_outer=8.65e-3,
        r_emitter_outer=9.8e-3,
        r_collector_inner=10.3e-3,
        r_collector_outer=11.85e-3,
        r_inner_clad_inner=11.9e-3,
        r_inner_clad_outer=12.25e-3,
        r_coolant_inner=12.25e-3,
        r_coolant_outer=12.95e-3,
        r_outer_clad_outer=13.30e-3,
        r_moderator_inner=13.52e-3,
        r_moderator_outer=14.52e-3,
        height=H
    )

    # 1.2 网格参数
    mesh_params = TFEMeshParams(
        n_axial=N_axial,
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=2,
        n_r_outer_clad=2,
        n_r_moderator=3
    )

    # 1.3 固体物性字典
    materials_dict = {
        'UO2': UO2(),
        'MoNb': MoNb(),
        'Molybdenum': Molybdenum(),
        'StainlessSteel': AusteniticStainlessSteel(),
        'ZrH': ZirconiumHydride()
    }

    # =========================================================================
    # 2. 气隙配置 (仅测试热源与对流，CO2外侧绝热)
    # =========================================================================
    fg_config = GapConfig(mode='simplified', h_eq=5678.0, material=Xenon(), emissivity_inner=0.95,
                          emissivity_outer=0.95)
    tec_config = GapConfig(mode='simplified', h_eq=1.45, material=Cesium(), emissivity_inner=0.40,
                           emissivity_outer=0.40)
    he_config = GapConfig(mode='simplified', h_eq=5678.0, material=Helium())
    # [关键] 将 CO2 的换热系数设为极小值以模拟绝热外边界
    co2_config = GapConfig(mode='simplified', h_eq=1e-20, material=CarbonDioxide())

    # =========================================================================
    # 3. 构建严格符合 HydraulicNetwork 接口的流体支路
    # =========================================================================
    fluid_mat = Sodium()
    A_flow = np.pi * (geometry.r_coolant_outer ** 2 - geometry.r_coolant_inner ** 2)
    D_h = 2.0 * (geometry.r_coolant_outer - geometry.r_coolant_inner)

    # 初始化创建流体网络中节点和接管
    all_vols = []
    all_juncs = []

    # 3.1 实例化边界体积 (严格遵循 BoundaryVolume 签名，传入 material)
    inlet = IncompressibleBoundaryVolume(name="Inlet", material=fluid_mat, P=1.162e5, T=743.0)
    outlet = IncompressibleBoundaryVolume(name="Outlet", material=fluid_mat, P=1.16e5, T=743.0)

    inlet.is_pressure_boundary = True
    outlet.is_pressure_boundary = True

    all_vols.extend([inlet, outlet])

    # 3.2 实例化冷却剂 (采用IncompressibleFluidChannel创建)
    chan = IncompressibleFluidChannel(
        name=f"Chan_1",
        n_nodes=mesh_params.n_axial, total_length=geometry.height,
        flow_area=A_flow, hydraulic_diam=D_h,
        initial_P=1.16e5, initial_T=743.0, material=fluid_mat
    )

    all_vols.extend(chan.volumes)
    all_juncs.extend(chan.internal_junctions)

    # 3.3 实例化 Junctions (严格包含 InletJunction)
    # 进口流量约束 0.0365 kg/s
    j_in = InletJunction(name="J_Inlet", from_vol=inlet, to_vol=chan.volumes[0], W_initial=0.0365)
    all_juncs.append(j_in)

    j_out = FlowJunction(name="J_Outlet", from_vol=chan.volumes[-1], to_vol=outlet, flow_area=A_flow, k_loss=0)
    all_juncs.append(j_out)

    # =========================================================================
    # 4. 实例化 TFEUnit 与 ReactorCore
    # =========================================================================

    # 4.1 轴向功率分布计算 (Cos 分布)
    z_centers = np.linspace(-np.pi / 2, np.pi / 2, N_axial)
    cos_profile = np.cos(z_centers)
    axial_profile = cos_profile / np.sum(cos_profile)

    # 4.2 实例化单根 TFEUnit 对象
    tfe = TFEUnit(
        name="Test_TFE",
        geometry=geometry,
        mesh_params=mesh_params,
        materials=materials_dict,
        coolant_channel=chan,
        fission_gas_config=fg_config,
        tec_gap_config=tec_config,
        he_gap_config=he_config,
        co2_gap_config=co2_config,
        power_fraction=1.0,
        axial_power_profile=axial_profile
    )

    # =========================================================================
    # [新增] 4.3 实例化 ReactorCore 容器，将 TFE 包装进去
    # =========================================================================
    core = ReactorCore(name="MyCore", tfe_dict={"TFE_Center": tfe}, tfe_multipliers={"TFE_Center": 6}, enable_tec_coupled=False)

    # =========================================================================
    # 5. 组装流体网络
    # =========================================================================
    fluid_network = HydraulicNetwork(volumes=all_vols, junctions=all_juncs, gravity_vector=0.0)

    # =========================================================================
    # 6. 系统组装与初始化 (严格遵循 SystemManager 签名)
    # =========================================================================
    system = SystemManager(fluid_network=fluid_network, start_time=0.0)

    # 【核心修改】将 ReactorCore 注册进系统，SystemManager会自动将其内部网格榨取出来
    system.add_component(core)

    # =========================================================================
    # 7. 瞬态测试循环与重启动逻辑 (Restart Logic)
    # =========================================================================
    dt = 0.5
    t_end = 1500.0

    # --- 重启动开关 ---
    # 如果想从头算，设为 None；如果想续算，设为文件路径如 "restart_t100.npz"
    # TEASA提示：首次运行建议设为 None，跑出存档后再设为相应文件名
    restart_file = "test_core_assemble_v2_restart_t600.npz"

    if restart_file:
        system.load_global_state(restart_file)
        current_time = system.global_time
        print(f"\n成功从存档恢复！当前仿真时间跳跃至: {current_time:.2f} s")

        # 【关键修正】：强制覆盖文件中的旧标志，开启热电耦合计算！
        # core.enable_tec_coupled = True
        core.enable_tec_coupled = False

        # 【TEASA 黑科技】：手动触发一次 post_step，把刚刚读进来的 1600K 热态温度推进给 C++ 电路！
        # 否则 C++ 第一步会带着 600K 的初始冷态温度去算电流，导致数值崩溃。
        core.post_step(0.0, current_time)
        print("已强制激活 TEC 电计算，并同步初始高温场至 C++ 模型！")
    else:
        # [分支 B: 全新启动]
        current_time = 0.0
        system.initialize_system()
        print("\n系统冷态初始化完成。")

    # [新增] 设置下一个自动保存的时间节点 (每隔 100s 保存一次)
    next_save_time = current_time + 100.0

    # 数据记录列表
    history_time = []
    history_T_pellet_max = []
    history_T_coolant_out = []
    history_I = []
    history_U = []

    print("开始瞬态测试...")
    while current_time <= t_end:
        # A. 施加边界条件 (均匀分配 3000 W 总功率)
        # TEASA说明：因为 tfe 是由对象的内存地址引用的，即使被包进了 core，
        # 在这里直接调用 tfe.update_neutronic_power 依然是有效且极其高效的！
        tfe.update_neutronic_power(p_total=3000.0)

        # 获取自适应步长
        dt = system.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=50.0)

        # B. 推进求解
        system.step(dt)

        # 提取温度
        T_pellet_max = np.max(tfe.solids['pellet'].T)
        T_coolant_out = chan.volumes[-1].T

        # 提取电学结果
        tec_res = core.thermo_calc.get_global_results()
        if tec_res is not None and core.enable_tec_coupled:
            I_tot = tec_res['Iout']
            U_tot = tec_res['Uout']
        else:
            I_tot, U_tot = 0.0, 0.0

        history_time.append(current_time)
        history_T_pellet_max.append(T_pellet_max)
        history_T_coolant_out.append(T_coolant_out)
        history_I.append(I_tot)
        history_U.append(U_tot)

        if len(history_time) % 10 == 0:
            print(
                f"Time: {current_time:6.2f} s | Max Pellet T: {T_pellet_max:7.2f} K | Coolant Out T: {T_coolant_out:7.2f} K")

        current_time += dt

        if current_time >= next_save_time:
            # 格式化文件名，如 restart_t100.npz, restart_t200.npz
            save_path = f"test_core_assemble_v2_restart_t{int(next_save_time)}.npz"
            print(f"\n[Checkpoint] 达到保存节点，正在写入: {save_path}")
            system.save_global_state(save_path)

            # 靶标向后推 100 秒
            next_save_time += 100.0

    print("测试完成。")

    # =========================================================================
    # 7. 结果绘图
    # =========================================================================
    plt.figure(figsize=(15, 4))

    plt.subplot(1, 3, 1)
    plt.plot(history_time, history_T_pellet_max, label='Max Pellet Temp', color='red')
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.title('Max Fuel Temp Transient')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(history_time, history_T_coolant_out, label='Coolant Outlet Temp', color='blue')
    plt.xlabel('Time [s]')
    plt.title('Coolant Outlet Temp Transient')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(history_time, history_I, label='Current (A)', color='green')
    plt.plot(history_time, np.array(history_U) * 100, label='Voltage (V x100)', color='purple', linestyle='--')
    plt.xlabel('Time [s]')
    plt.title('TEC Electrical Output')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_test()
    TEASAProfiler.report()
