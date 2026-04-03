import numpy as np
import matplotlib.pyplot as plt

from Solvers.HeatConduction.Mesh import Mesh2D
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
        r_moderator_outer=13.90e-3,
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
    fg_config = GapConfig(mode='simplified', h_eq=5678.0, material=Xenon(), emissivity_inner=0.15,
                          emissivity_outer=0.30)
    tec_config = GapConfig(mode='simplified', h_eq=1.45, material=Cesium(), emissivity_inner=0.30,
                           emissivity_outer=0.30)
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

    # 3.2 轴向功率分布计算 (Cos 分布)
    # z_centers = np.linspace(-np.pi / 2, np.pi / 2, N_axial)
    # cos_profile = np.cos(z_centers)
    # axial_profile = cos_profile / np.sum(cos_profile)

    # 3.2 轴向功率分布计算
    raw_power_data = np.array([
        0.000540698, 0.000674943, 0.000790854, 0.000894659, 0.000986677,
        0.001067197, 0.001136471, 0.001194716, 0.001242119, 0.001278831,
        0.001304965, 0.001320617, 0.001325828, 0.001320613, 0.001304968,
        0.001278831, 0.001242123, 0.001194720, 0.001136475, 0.001067202,
        0.000986683, 0.000894664, 0.000790859, 0.000674947, 0.000546582
    ])

    # 2. 为原始数据生成归一化的相对坐标 (从 -1 到 1，代表堆芯底部到顶部)
    z_raw = np.linspace(-1, 1, len(raw_power_data))

    # 3. 为你当前的计算网格生成相对坐标
    # 不管你的 N_axial 是 10, 50 还是 100，这里都映射到 -1 到 1 的区间
    z_centers = np.linspace(-1, 1, N_axial)

    # 4. 使用 numpy 的线性插值，将25个点的数据映射到你需要的 N_axial 个点上
    interpolated_profile = np.interp(z_centers, z_raw, raw_power_data)

    # 5. 归一化 (保证所有节点的功率因子总和为1)
    axial_profile = interpolated_profile / np.sum(interpolated_profile)

    # 3.3 堆芯配置数据
    ring_names = ["TFE_Center", "TFE_Ring1", "TFE_Ring2", "TFE_Ring3"]
    multipliers = [1, 6, 12, 15]
    tfes = {}

    # 3.4 循环构建流道与 TFEUnit 元件
    for name, mult in zip(ring_names, multipliers):
        # 1. 为每一种 TFE 构建专属的单根冷却剂流道
        chan = IncompressibleFluidChannel(
            name=f"Chan_{name}", n_nodes=mesh_params.n_axial, total_length=geometry.height,
            flow_area=A_flow, hydraulic_diam=D_h, initial_P=1.16e5, initial_T=743.0, material=fluid_mat
        )
        all_vols.extend(chan.volumes)
        all_juncs.extend(chan.internal_junctions)

        # 进出口连接 (单根管的流量)
        j_in = InletJunction(name=f"J_In_{name}", from_vol=inlet, to_vol=chan.volumes[0], W_initial=0.0365)
        j_out = FlowJunction(name=f"J_Out_{name}", from_vol=chan.volumes[-1], to_vol=outlet, flow_area=A_flow, k_loss=0)
        all_juncs.extend([j_in, j_out])

        # 2. 实例化独立的 TFEUnit
        tfes[name] = TFEUnit(
            name=name, geometry=geometry, mesh_params=mesh_params, materials=materials_dict,
            coolant_channel=chan, fission_gas_config=fg_config, tec_gap_config=tec_config,
            he_gap_config=he_config, co2_gap_config=co2_config,
            power_fraction=0.95,  # [物理修正] 燃料内发热占比 95%
            axial_power_profile=axial_profile
        )

    # 3.5 定义全局慢化剂二维网格
    r_bounds = [0, 21e-3, 53.25e-3, 86.75e-3, 109e-3]  # 4 个圆环的半径边界 [m]
    mod_meshes = []
    for i in range(4):
        delta_r = r_bounds[i + 1] - r_bounds[i]
        mesh = Mesh2D(x_dim=delta_r, n_x=3, y_dim=H, n_y=N_axial,
                      geometry_type='cylindrical', inner_radius=r_bounds[i])
        mod_meshes.append(mesh)

    # 3.6 装配映射字典与 ReactorCore
    tfe_multipliers = {name: mult for name, mult in zip(ring_names, multipliers)}
    ring_mapping = {"TFE_Center": 0, "TFE_Ring1": 1, "TFE_Ring2": 2, "TFE_Ring3": 3}

    core = ReactorCore(
        name="MyCore",
        tfe_dict=tfes,
        tfe_multipliers=tfe_multipliers,
        enable_tec_coupled=False,
        mod_meshes=mod_meshes,
        mod_material=ZirconiumHydride(),
        ring_mapping=ring_mapping,
        T_space=250.0  # 太空背景温度 250K
    )

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
    t_end = 1540.0

    # --- 重启动开关 ---
    # 如果想从头算，设为 None；如果想续算，设为文件路径如 "restart_t100.npz"
    # TEASA提示：首次运行建议设为 None，跑出存档后再设为相应文件名
    restart_file = "test_core_assemble_v3_restart_t1500.npz"

    if restart_file:
        system.load_global_state(restart_file)
        current_time = system.global_time
        print(f"\n成功从存档恢复！当前仿真时间跳跃至: {current_time:.2f} s")

        # 【关键修正】：强制覆盖文件中的旧标志，开启热电耦合计算！
        # core.enable_tec_coupled = False
        core.enable_tec_coupled = True

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
    next_save_time = current_time + 200.0

    # --- 初始化数据记录字典 ---
    history_time = []
    history_I = []
    history_U = []

    # 为每一种 TFE 创建独立的温度记录列表
    history_T_pellet = {name: [] for name in tfes.keys()}
    history_T_coolant = {name: [] for name in tfes.keys()}

    # [新增] 为每一圈慢化剂创建温度记录列表
    history_T_mod = {f"ModRing_{i}": [] for i in range(4)}

    print("开始瞬态测试...")
    while current_time <= t_end:
        # A. 施加边界条件 (均匀分配 3000 W 总功率)
        # TEASA说明：因为 tfe 是由对象的内存地址引用的，即使被包进了 core，
        # 在这里直接调用 tfe.update_neutronic_power 依然是有效且极其高效的！

        for name, tfe in tfes.items():
            tfe.update_neutronic_power(p_total=3100.0)

        # 获取自适应步长
        dt = system.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=10.0)

        # B. 推进求解
        system.step(dt)

        # C. 数据提取
        history_time.append(current_time)

        # C.1 提取每种 TFE 的最高燃料温度和冷却剂出口温度
        global_max_pellet = 0.0
        for name, tfe in tfes.items():
            T_pellet_max = np.max(tfe.solids['pellet'].T)
            # tfe.coolant 保存了流道对象的引用，直接获取最后一个 volume 的温度
            T_coolant_out = tfe.coolant.volumes[-1].T

            history_T_pellet[name].append(T_pellet_max)
            history_T_coolant[name].append(T_coolant_out)

            if T_pellet_max > global_max_pellet:
                global_max_pellet = T_pellet_max

        # C.2 提取每一圈全局慢化剂的平均温度
        for i, ring in enumerate(core.mod_rings):
            T_mod_avg = np.mean(ring.T)
            history_T_mod[f"ModRing_{i}"].append(T_mod_avg)

        # C.3 提取全局电学结果
        tec_res = core.thermo_calc.get_global_results()
        if tec_res is not None and core.enable_tec_coupled:
            I_tot = tec_res['Iout']
            U_tot = tec_res['Uout']
        else:
            I_tot, U_tot = 0.0, 0.0

        history_I.append(I_tot)
        history_U.append(U_tot)

        # 控制台打印监控 (仅打印全局最高温度和总电流)
        if len(history_time) % 10 == 0:
            print(f"Time: {current_time:6.2f} s | Global Max Pellet T: {global_max_pellet:7.2f} K | "
                  f"I: {I_tot:6.1f} A | U: {U_tot:6.2f} V")

        current_time += dt

        if current_time >= next_save_time:
            # 格式化文件名，如 restart_t100.npz, restart_t200.npz
            save_path = f"test_core_assemble_v3_restart_t{int(next_save_time)}.npz"
            print(f"\n[Checkpoint] 达到保存节点，正在写入: {save_path}")
            system.save_global_state(save_path)

            # 靶标向后推 100 秒
            next_save_time += 100.0

    print("测试完成。")

    # =========================================================================
    # 8. 结果绘图 (2x2 全局视图)
    # =========================================================================
    # 使用更大的画布来容纳 4 张子图
    plt.figure(figsize=(14, 10))

    # 图 1: 各圈 TFE 芯块最高温度
    plt.subplot(2, 2, 1)
    for name in tfes.keys():
        plt.plot(history_time, history_T_pellet[name], label=name)
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.title('Max Fuel Temp Transient (By Ring)')
    plt.grid(True)
    plt.legend()

    # 图 2: 各圈 TFE 冷却剂出口温度
    plt.subplot(2, 2, 2)
    for name in tfes.keys():
        plt.plot(history_time, history_T_coolant[name], label=name)
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.title('Coolant Outlet Temp Transient (By Ring)')
    plt.grid(True)
    plt.legend()

    # 图 3: [新增] 全局慢化剂各圈平均温度分布
    plt.subplot(2, 2, 3)
    for i in range(4):
        plt.plot(history_time, history_T_mod[f"ModRing_{i}"], label=f"Moderator Ring {i}")
    plt.xlabel('Time [s]')
    plt.ylabel('Temperature [K]')
    plt.title('Moderator Average Temp Transient')
    plt.grid(True)
    plt.legend()

    # 图 4: 热电转换输出 (全堆芯总输出)
    plt.subplot(2, 2, 4)
    plt.plot(history_time, history_I, label='Total Current (A)', color='green')
    plt.plot(history_time, np.array(history_U) * 100, label='Voltage (V x100)', color='purple', linestyle='--')
    plt.xlabel('Time [s]')
    plt.ylabel('Electrical Output')
    plt.title('TEC Core Electrical Output')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_test()
    TEASAProfiler.report()
