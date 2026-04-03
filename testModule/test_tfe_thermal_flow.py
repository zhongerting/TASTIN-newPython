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
# 请确保按您实际的文件路径导入
from Components.TFEUnit import TFEUnit, TFEGeometry, TFEMeshParams, GapConfig


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
    fg_config = GapConfig(mode='simplified', h_eq=5678.0, material=Xenon(), emissivity_inner=0.15, emissivity_outer=0.35)
    # fg_config = GapConfig(mode='meshed', material=Xenon(), n_radial_nodes=1, emissivity_inner=0.15, emissivity_outer=0.4)
    tec_config = GapConfig(mode='simplified', h_eq=168.0, material=Cesium(), emissivity_inner=0.35, emissivity_outer=0.6)
    he_config = GapConfig(mode='simplified', h_eq=168.0, material=Helium())
    # [关键] 将 CO2 的换热系数设为极小值以模拟绝热外边界
    co2_config = GapConfig(mode='simplified', h_eq=1e-20, material=CarbonDioxide())

    # =========================================================================
    # 3. 构建严格符合 HydraulicNetwork 接口的流体支路
    # =========================================================================
    fluid_mat = Sodium()
    A_flow = np.pi * (geometry.r_coolant_outer ** 2 - geometry.r_coolant_inner ** 2)
    D_h = 2.0 * (geometry.r_coolant_outer - geometry.r_coolant_inner)
    # V_node = A_flow * (H / N_axial)
    # L_node = H / N_axial

    # volumes = []
    # junctions = []
    # coolant_channel_vols = []  # 用于传递给 TFEUnit 的流体引用

    # 初始化创建流体网络中节点和接管
    all_vols = []
    all_juncs = []

    # 3.1 实例化边界体积 (严格遵循 BoundaryVolume 签名，传入 material)
    inlet = IncompressibleBoundaryVolume(name="Inlet", material=fluid_mat, P=1.162e5, T=743.0)
    outlet = IncompressibleBoundaryVolume(name="Outlet", material=fluid_mat, P=1.16e5, T=743.0)

    inlet.is_pressure_boundary = True
    outlet.is_pressure_boundary = True

    # volumes.append(inlet)
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

    # 3.2 实例化冷却剂内部网格 (IncompressibleFluidVolume)
    # for i in range(N_axial):
    #     vol = IncompressibleFluidVolume(
    #         name=f"Coolant_Vol_{i}", volume=V_node, length=L_node,
    #         flow_area=A_flow, hydraulic_diam=D_h,
    #         initial_P=1.0e5, initial_T=743.0, material=fluid_mat
    #     )
    #     volumes.append(vol)
    #     coolant_channel_vols.append(vol)
    #
    # volumes.append(outlet)

    # 3.3 实例化 Junctions (严格包含 InletJunction)
    # 进口流量约束 0.0365 kg/s

    j_in = InletJunction(name="J_Inlet", from_vol=inlet, to_vol=chan.volumes[0], W_initial=0.0365)
    all_juncs.append(j_in)

    # for i in range(N_axial - 1):
    #     j_int = FlowJunction(name=f"J_int_{i}", from_vol=coolant_channel_vols[i], to_vol=coolant_channel_vols[i + 1])
    #     junctions.append(j_int)

    j_out = FlowJunction(name="J_Outlet", from_vol=chan.volumes[-1], to_vol=outlet, flow_area=A_flow, k_loss=0)
    all_juncs.append(j_out)

    # =========================================================================
    # 4. 实例化 TFEUnit
    # =========================================================================

    # 4.1 轴向功率分布计算 (Cos 分布)
    z_centers = np.linspace(-np.pi / 2, np.pi / 2, N_axial)
    cos_profile = np.cos(z_centers)
    axial_profile = cos_profile / np.sum(cos_profile)

    # 4.2 实例化 TFEUnit 对象
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
    # 5. 组装流体网络 TFEUnit
    # =========================================================================
    fluid_network = HydraulicNetwork(volumes=all_vols, junctions=all_juncs, gravity_vector=0.0)

    # =========================================================================
    # 6. 系统组装与初始化 (严格遵循 SystemManager 签名)
    # =========================================================================
    system = SystemManager(fluid_network=fluid_network, start_time=0.0)

    system.add_component(tfe)

    # for solid in tfe.get_solids():
    #     system.add_solid(solid)  # 假设您的 SystemManager 具备标准的 add_solid 方法
    # for coupler in tfe.get_couplers():
    #     system.add_coupler(coupler)  # 假设您的 SystemManager 具备标准的 add_coupler 方法

    # =========================================================================
    # 7. 瞬态测试循环与重启动逻辑 (Restart Logic)
    # =========================================================================
    dt = 0.5
    t_end = 600.0

    # --- [新增] 重启动开关 ---
    # 如果想从头算，设为 None；如果想续算，设为文件路径如 "restart_t100.npz"
    restart_file = "restart_t500.npz"

    if restart_file:
        # [分支 A: 断点续算]
        system.load_global_state(restart_file)

        # 【核心修正 1】将系统的全局时间同步给局部循环变量
        current_time = system.global_time
        print(f"\n成功从存档恢复！当前仿真时间跳跃至: {current_time:.2f} s")
    else:
        # [分支 B: 全新启动]
        current_time = 0.0
        system.initialize_system()
        print("\n系统冷态初始化完成。")

    # [新增] 设置下一个自动保存的时间节点 (每隔 100s 保存一次)
    next_save_time = current_time + 100.0

    history_time = []
    history_T_pellet_max = []
    history_T_coolant_out = []

    # system.initialize_system()
    #
    # system.load_global_state("restart_t300.npz")

    print("开始瞬态测试...")
    while current_time <= t_end:
        # A. 施加边界条件 (均匀分配 3000 W 总功率)
        tfe.update_neutronic_power(p_total=3000.0)
        dt = system.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=10.0)
        # tfe.pre_step(dt, current_time)

        # B. 推进求解
        system.step(dt)

        # C. 后处理反馈
        # tfe.post_step(dt, current_time)

        # 提取关键物理量用于监控
        T_pellet_max = np.max(tfe.solids['pellet'].T)
        T_coolant_out = chan.volumes[-1].T  # 获取流道最后一个体积的温度

        history_time.append(current_time)
        history_T_pellet_max.append(T_pellet_max)
        history_T_coolant_out.append(T_coolant_out)

        if len(history_time) % 10 == 0:
            print(
                f"Time: {current_time:6.2f} s | Max Pellet T: {T_pellet_max:7.2f} K | Coolant Out T: {T_coolant_out:7.2f} K")

        current_time += dt

        # --- [新增] 自动保存断点续算文件 ---
        if current_time >= next_save_time:
            # 格式化文件名，如 restart_t100.npz, restart_t200.npz
            save_path = f"restart_t{int(next_save_time)}.npz"
            print(f"\n[Checkpoint] 达到保存节点，正在写入: {save_path}")
            system.save_global_state(save_path)

            # 靶标向后推 100 秒
            next_save_time += 100.0

    print("测试完成。")

    # =========================================================================
    # 8. 结果绘图
    # =========================================================================
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(history_time, history_T_pellet_max, label='Max Pellet Temp', color='red')
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.title('Max Fuel Temp Transient')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(history_time, history_T_coolant_out, label='Coolant Outlet Temp', color='blue')
    plt.xlabel('Time [s]')
    plt.title('Coolant Outlet Temp Transient')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_test()
