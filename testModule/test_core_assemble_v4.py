import os
import numpy as np
import matplotlib.pyplot as plt

import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, '..'))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from profiler import TEASAProfiler

# --- 导入底层求解器与组件 ---
from Solvers.SystemManager import SystemManager
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidVolume
# [TEASA 新增] 导入非均匀流道与魔法接口
from Solvers.Hydrodynamics.Components import NonUniformIncompressibleFluidChannel, MacroFlowJunction

# --- 导入材料库 ---
from Materials.Solids.UO2 import UO2
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum
from Materials.Solids.StainlessSteel import AusteniticStainlessSteel
from Materials.Solids.ZrH import ZirconiumHydride
from Materials.Solids.BerylliumOxide import BerylliumOxide  # [新增] 反射层材料
from Materials.Solids.GasGaps import Xenon, Cesium, CarbonDioxide, Helium
from Materials.Fluids.Sodium import Sodium

# --- 导入 TFE 组件与数据结构 ---
from Components.TFEUnit import TFEUnit, TFEGeometry, TFEMeshParams, GapConfig
from Components.ReactorCore import (
    ReactorCore,
    GlobalAnnulusStructureConfig,
    GlobalGapStructureConfig,
)
from profiler import TEASAProfiler


def run_test_v4(t_end: float = 3000.0,
                save_interval: float = 200.0,
                restart_file: str = "test_core_assemble_v4_restart_t1200.npz",
                enable_plot: bool = False):
    """
    全堆芯非均匀反射层 + 虚拟代表通道流量归集 瞬态测试
    """
    print("=== TASTIN System Test V4: Non-Uniform Core & Flow Hubs ===")

    # =========================================================================
    # Part 1: 核心参数设置 (Core Parameter Setup)
    # =========================================================================
    print("1. Setting up core parameters (Geometry, Mesh arrays, Power profiles)...")
    # --- 1.1 轴向非均匀网格与几何分配 ---
    H_active = 0.375  # 活性区高度 [m]
    L_reflector = 0.065  # 上下反射层高度 [m]
    N_active = 30  # 活性区轴向网格数
    dz_active = H_active / N_active

    # 构造流道用的 node_lengths 数组 (长度为 32)
    node_lengths = [L_reflector] + [dz_active] * N_active + [L_reflector]
    node_lengths_arr = np.array(node_lengths)

    # 构造 TFEUnit 用的分段指令
    axial_length_alloc = [L_reflector, H_active, L_reflector]
    axial_node_alloc = [1, N_active, 1]
    N_total = sum(axial_node_alloc)  # 总网格数为 32

    # --- 1.2 径向几何参数设置 ---
    # 【专家修正】: height 必须是真实的物理总高度，即活性区 + 两个反射层
    geom_params = TFEGeometry(
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
        height=H_active + 2 * L_reflector  # 0.505 m
    )

    # --- 1.3 径向与轴向网格密度设置 ---
    mesh_params = TFEMeshParams(
        n_axial=N_total,  # 这里的总节点数传入 32
        n_r_pellet=5,
        n_r_emitter=1,
        n_r_collector=1,
        n_r_inner_clad=2,
        n_r_outer_clad=2,
        n_r_moderator=3
    )

    # --- 1.4 构造带反射层补零的轴向功率分布 (截断余弦) ---
    # 先计算活性区 (长度为 30) 的截断余弦分布
    z_centers = np.linspace(-1, 1, N_active)
    active_power_profile = np.cos(z_centers * np.pi / 2 * 0.8)  # 0.8为截断因子

    # 【陷阱规避】: 首尾各拼接一个 0.0，使数组长度严格等于 32
    padded_power_profile = np.concatenate(([0.0], active_power_profile, [0.0]))

    # --- 1.5 实例化物性字典 ---
    materials_dict = {
        'UO2': UO2(),
        'MoNb': MoNb(),
        'Molybdenum': Molybdenum(),
        'StainlessSteel': AusteniticStainlessSteel(),
        'ZrH': ZirconiumHydride(),
        'BerylliumOxide': BerylliumOxide(),  # [新增] 必须包含反射层材料
        'Sodium': Sodium()
    }

    # =========================================================================
    # Part 2: 气隙配置 (Gas Gap Configuration)
    # =========================================================================
    print("2. Configuring gas gaps...")
    # 提取自 v3 版本，保持 Simplified 模式下的等效导热系数和黑度
    cfg_fg = GapConfig(
        mode='simplified', h_eq=5678.0,
        emissivity_inner=0.15, emissivity_outer=0.3
    )

    # TEC 极间隙：由 TECCouple2D 专门处理电热耦合和辐射，气体本身对流/导热设为0
    cfg_tec = GapConfig(
        mode='simplified', h_eq=1.45,
        emissivity_inner=0.3, emissivity_outer=0.3
    )

    # 氦气隙
    cfg_he = GapConfig(
        mode='simplified', h_eq=5678.0,
        emissivity_inner=0.8, emissivity_outer=0.8
    )

    # CO2 气隙
    cfg_co2 = GapConfig(
        mode='simplified', h_eq=53.6,
        emissivity_inner=0.8, emissivity_outer=0.8
    )

    # =========================================================================
    # Part 3: 流体支路构建 (Fluid Branch Construction)
    # =========================================================================
    print("3. Building fluid network branches and macro junctions...")

    # --- 3.1 提取流体物性与系统边界条件 (严格继承 v3) ---
    sodium = materials_dict['Sodium']
    P_inlet_sys = 1.162e5  # 系统入口压力 [Pa]
    P_outlet_sys = 1.16e5  # 系统出口压力 [Pa]
    T_inlet_sys = 743.0  # 冷却剂入口初始温度 [K]

    # --- 3.2 严密的几何换算 (单根流道) ---
    # 严格照搬 v3 的几何计算逻辑
    r_c_in = geom_params.r_coolant_inner
    r_c_out = geom_params.r_coolant_outer
    A_flow = np.pi * (r_c_out ** 2 - r_c_in ** 2)  # 单管流通面积 [m^2]
    D_h = 2.0 * (r_c_out - r_c_in)  # 单管水力直径 [m]

    # --- 3.3 定义堆芯的径向布局与并联倍数 (严格继承 v3) ---
    ring_names = ['Center', 'Ring1', 'Ring2', 'Ring3']
    multipliers = [1, 6, 12, 15]
    W_single_design = 0.0365

    all_fluid_vols = []
    all_fluid_juncs = []
    fluid_channels = {}

    # --- 3.4 建立全局定压边界联箱 (System Pressure References) ---
    # 严格遵循 v3 签名，传入 material
    inlet_plenum = IncompressibleBoundaryVolume(
        name="Global_Inlet", material=sodium, P=P_inlet_sys, T=T_inlet_sys
    )
    outlet_plenum = IncompressibleBoundaryVolume(
        name="Global_Outlet", material=sodium, P=P_outlet_sys, T=T_inlet_sys
    )
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum.is_pressure_boundary = True
    all_fluid_vols.extend([inlet_plenum, outlet_plenum])

    # --- 3.5 逐环构建并联流道与宏观乘子拓扑 ---
    for name, mult in zip(ring_names, multipliers):
        # [A] 实例化非均匀单根物理流道 (使用我们之前定义的 NonUniform 类型)
        chan = NonUniformIncompressibleFluidChannel(
            name=f"Chan_{name}",
            node_lengths=node_lengths_arr,  # 长度为 32 的非均匀数组
            flow_area=A_flow,  # 严格单管面积
            hydraulic_diam=D_h,  # 严格单管水力直径
            initial_P=P_inlet_sys,
            initial_T=T_inlet_sys,
            material=sodium
        )
        fluid_channels[name] = chan
        all_fluid_vols.extend(chan.volumes)
        all_fluid_juncs.extend(chan.internal_junctions)

        # [B] 建立中间控制体 (Intermediate Hubs)
        # 面积放大 mult 倍，用于质量归集
        L_buffer = 0.01
        inter_in = IncompressibleFluidVolume(
            name=f"InterIn_{name}",
            volume=A_flow * mult * L_buffer, length=L_buffer,
            flow_area=A_flow * mult, hydraulic_diam=D_h, material=sodium,
            initial_P=P_inlet_sys, initial_T=T_inlet_sys
        )
        inter_out = IncompressibleFluidVolume(
            name=f"InterOut_{name}",
            volume=A_flow * mult * L_buffer, length=L_buffer,
            flow_area=A_flow * mult, hydraulic_diam=D_h, material=sodium,
            initial_P=P_outlet_sys, initial_T=T_inlet_sys
        )
        all_fluid_vols.extend([inter_in, inter_out])

        # [C] 全局联箱 -> 中间缓冲 (强迫循环驱动)
        j_in = InletJunction(
            name=f"J_In_{name}",
            from_vol=inlet_plenum,
            to_vol=inter_in,
            W_initial=W_single_design * mult
        )

        # [D] 进口端魔法接口：符合你的 MacroFlowJunction 定义
        # 从 inter_in (宏观端) 到 chan.volumes[0] (微观端)
        j_macro_in = MacroFlowJunction(
            name=f"J_MacroIn_{name}",
            from_vol=inter_in,
            to_vol=chan.volumes[0],
            macro_vol=inter_in,  # 指定 inter_in 是需要缩放的宏观端
            multiplier=mult,  # 缩放倍数
            flow_area=A_flow  # 动量计算强制使用单管面积
        )

        # [E] 出口端魔法接口：符合你的 MacroFlowJunction 定义
        # 从 chan.volumes[-1] (微观端) 到 inter_out (宏观端)
        j_macro_out = MacroFlowJunction(
            name=f"J_MacroOut_{name}",
            from_vol=chan.volumes[-1],
            to_vol=inter_out,
            macro_vol=inter_out,  # 指定 inter_out 是需要缩放的宏观端
            multiplier=mult,  # 缩放倍数
            flow_area=A_flow  # 动量计算强制使用单管面积
        )

        # [F] 中间缓冲 -> 全局出口联箱
        j_out = FlowJunction(
            name=f"J_Out_{name}",
            from_vol=inter_out,
            to_vol=outlet_plenum,
            flow_area=A_flow * mult
        )

        all_fluid_juncs.extend([j_in, j_macro_in, j_macro_out, j_out])

    # TODO: 建立全局 InletPlenum 和 OutletPlenum (Pressure Reference)
    # TODO: 循环各圈 (Ring 1~4)，建立 NonUniformIncompressibleFluidChannel
    # TODO: 建立微小缓冲控制体 (InterIn, InterOut)
    # TODO: 建立 MacroFlowJunction (魔法接口) 和普通 FlowJunction/InletJunction
    # TODO: 将生成的 volumes 和 junctions 收集到 all_fluid_vols 和 all_fluid_juncs 中

    # =========================================================================
    # Part 4: 堆芯组件构建 (Core Component Construction)
    # =========================================================================
    print("4. Building TFE Units and Reactor Core...")

    tfes = {}

    # --- 4.1 组装多物理场 TFE 单管 ---
    for name in ring_names:
        tfes[name] = TFEUnit(
            name=name,
            geometry=geom_params,
            mesh_params=mesh_params,
            materials=materials_dict,
            coolant_channel=fluid_channels[name],  # 挂载专属流道
            fission_gas_config=cfg_fg,
            tec_gap_config=cfg_tec,
            he_gap_config=cfg_he,
            co2_gap_config=cfg_co2,
            power_fraction=1.0,  # 归一化在 Core 统一处理
            axial_power_profile=padded_power_profile,  # 补零后的功率数组 (Len=32)
            axial_length_allocation=axial_length_alloc,
            axial_node_allocation=axial_node_alloc,
            axial_contact_resistance=0.05  # 注入 0.05 K/W 轴向接触热阻
        )

    # --- 4.2 构建全局慢化剂网格 (The Mesh Master 对齐) ---
    tfe_multipliers = {name: mult for name, mult in zip(ring_names, multipliers)}
    # 功率因子按“单根真实燃料元件的功率份额”定义。
    # 这里先采用最简单的等功率假设：每根真实燃料元件占全堆总功率的 1 / N_total。
    total_real_elements = float(sum(multipliers))
    tfe_power_factors = {name: 1.0 / total_real_elements for name in ring_names}
    ring_mapping = {name: i for i, name in enumerate(ring_names)}

    # 【核心对齐】: 累加求出所有网格界面的绝对坐标 (Len=33)
    common_y_faces = np.insert(np.cumsum(node_lengths_arr), 0, 0.0)

    mod_meshes = []
    r_in = geom_params.r_moderator_outer

    # 这里按 Fortran 版本“慢化剂 -> 间隙 -> 筒体 -> 间隙 -> 反射层”的思路，
    # 给出一组可直接运行的外层结构示例参数。
    # 后续若接 CoreInput 解析器，只需要把这些数值替换成输入文件读入值即可。
    r_mod_outer_sys = 57.0e-3
    delta_r = (r_mod_outer_sys - r_in) / 4

    for i in range(4):
        # 强制传入 y_faces=common_y_faces，确保轴向高度和非均匀节点分布绝对吻合
        mesh = Mesh2D(x_dim=delta_r, n_x=3,
                      y_dim=geom_params.height, n_y=N_total,
                      y_faces=common_y_faces,  # <--- 锁定轴向拓扑
                      geometry_type='cylindrical', inner_radius=r_in)
        mod_meshes.append(mesh)
        r_in += delta_r

    # --- 4.3 构建堆芯外层结构配置 ---
    # 第一层间隙：全局慢化剂外表面 -> 筒体内表面
    moderator_barrel_gap_cfg = GlobalGapStructureConfig(
        mode='simplified',
        width=5.0e-3,
        h_eq=0.0,
        emissivity_inner=0.8,
        emissivity_outer=0.8
    )

    # 筒体：当前先用不锈钢建模，径向控制体数参考旧版输入风格取 3
    barrel_cfg = GlobalAnnulusStructureConfig(
        material=materials_dict['StainlessSteel'],
        outer_radius=65.0e-3,
        n_radial=3,
        initial_temp=743.0,
        outer_surface_emissivity=0.05
    )

    # 第二层间隙：筒体外表面 -> 反射层内表面
    barrel_reflector_gap_cfg = GlobalGapStructureConfig(
        mode='simplified',
        width=2.0e-3,
        h_eq=0.0,
        emissivity_inner=0.8,
        emissivity_outer=0.8
    )

    # 外部反射层：按你的要求默认保留 8 个径向控制体
    reflector_cfg = GlobalAnnulusStructureConfig(
        material=materials_dict['BerylliumOxide'],
        outer_radius=102.0e-3,
        n_radial=8,
        initial_temp=743.0,
        outer_surface_emissivity=0.6
    )

    # --- 4.4 组装宏观 ReactorCore ---
    core = ReactorCore(
        name="TASTIN_Core",
        tfe_dict=tfes,
        tfe_multipliers=tfe_multipliers,
        tfe_power_factors=tfe_power_factors,
        mod_meshes=mod_meshes,
        mod_material=materials_dict['ZrH'],
        ring_mapping=ring_mapping,
        barrel_config=barrel_cfg,
        reflector_config=reflector_cfg,
        moderator_barrel_gap_config=moderator_barrel_gap_cfg,
        barrel_reflector_gap_config=barrel_reflector_gap_cfg,
        T_space=250.0,  # 继承自 v3 的太空辐射环境温度
        alpha_tec=0.5,
        enable_tec_coupled=False  # 稳态/热工水力测试验证期，暂关电耦合
    )

    # =========================================================================
    # Part 5: 流体网络组装 (Fluid Network Assembly)
    # =========================================================================
    print("5. Assembling Hydraulic Network...")

    # 【严谨性校验】: 收集所有控制体和连接器。
    # 得益于 MacroFlowJunction，HydraulicNetwork 将自动识别分并并联的质量/能量缩放。
    hydraulic_net = HydraulicNetwork(
        volumes=all_fluid_vols,
        junctions=all_fluid_juncs,
        gravity_vector=0.0  # 严格继承 v3，维持重力压降计算一致性
    )

    # =========================================================================
    # Part 6: 系统组装与初始化 (System Assembly & Initialization)
    # =========================================================================
    print("6. Initializing System Manager...")

    system = SystemManager(fluid_network=hydraulic_net, start_time=0.0)
    system.add_component(core)

    # 【核心防御】: 执行初始化。
    # 该步骤会完成 FluidSolidCouple 的网格匹配、物性对象绑定以及 ODE 索引分配。
    system.initialize_system()

    # =========================================================================
    # Part 7: 瞬态测试循环与重启动逻辑 (Transient Test Loop & Restart Logic)
    # =========================================================================
    # --- 7.1 仿真时间参数 ---
    t_end = float(t_end)  # 仿真结束时间 [s]
    save_interval = float(save_interval)  # 自动保存时间间隔 [s]

    # --- 7.2 重启动逻辑 (严格使用 v3 的 load_global_state) ---
    # 如果想从头算，设为 None；如果想续算，设为文件路径如 "restart_tXXX.npz"
    if restart_file and os.path.exists(restart_file):
        system.load_global_state(restart_file)
        current_time = system.global_time
        print(f"   => [Restart] 成功从存档恢复！当前仿真时间跳跃至: {current_time:.2f} s")

        core.enable_tec_coupled = True

        # 如果开启了热电耦合，可以在这里进行一次同步
        core.post_step(0.0, current_time)
    else:
        current_time = 0.0
        system.initialize_system()
        # 由 ReactorCore 统一持有点堆模块，并在初始化时把稳态总功率同步到 Fuel 内热源。
        core.initialize_point_reactor(total_power_initial=3100.0)
        print("   => [Init] 系统冷态初始化完成。")

    # 设置下一个自动保存的时间节点
    next_save_time = current_time + save_interval

    # --- 7.3 初始化数据记录字典 (严格遵循 v3 结构) ---
    history_time = []
    history_I = []
    history_U = []
    # 为每一种 TFE 创建独立的记录列表
    history_T_pellet = {name: [] for name in ring_names}
    history_T_coolant = {name: [] for name in ring_names}
    # 为每一圈慢化剂创建记录列表
    history_T_mod = {f"ModRing_{i}": [] for i in range(4)}
    # [新增] 流量验证
    history_W_flow = {name: [] for name in ring_names}

    # --- 7.4 瞬态主循环 ---
    print("开始执行求解...")
    while current_time <= t_end:

        # A. 设置外加反应性。当前测试先采用零外加反应性，只验证耦合链路本身。
        current_rho_control = 0.0

        # B. 获取自适应步长 (严格调用 v3 方法)
        dt = system.compute_adaptive_dt(min_dt=1e-3, max_dt=0.2, safety_factor=20.0)

        # C. 推进求解。中子状态更新与功率下发由 ReactorCore 在 SystemManager.step 内部完成。
        system.step(dt, reactivity_control=current_rho_control)

        current_p_total = 0.0 if core.point_reactor is None else core.point_reactor.total_power

        # D. 数据提取与记录
        history_time.append(current_time)

        global_max_pellet = 0.0
        for name, tfe in tfes.items():
            # 提取芯块最高温度
            T_pellet_max = np.max(tfe.solids['pellet'].T)
            # 从 tfe 挂载的流道对象中提取最后一个 volume 的温度
            T_coolant_out = tfe.coolant.volumes[-1].T

            history_T_pellet[name].append(T_pellet_max)
            history_T_coolant[name].append(T_coolant_out)

            # 记录流道流量 (验证 MacroFlowJunction)
            history_W_flow[name].append(tfe.coolant.internal_junctions[0].W)

            if T_pellet_max > global_max_pellet:
                global_max_pellet = T_pellet_max

        # 提取慢化剂各圈平均温度
        for i, ring in enumerate(core.mod_rings):
            history_T_mod[f"ModRing_{i}"].append(np.mean(ring.T))

        # 提取全局电学结果 (如果开启了 TEC)
        tec_res = core.thermo_calc.get_global_results()
        if tec_res is not None and core.enable_tec_coupled:
            history_I.append(tec_res.get('Iout', 0.0))
            history_U.append(tec_res.get('Uout', 0.0))
        else:
            history_I.append(0.0)
            history_U.append(0.0)

        # 控制台打印监控
        if len(history_time) % 20 == 0:
            print(f"Time: {current_time:6.2f} s | Power: {current_p_total / 1000:4.0f} kW | "
                  f"Max Fuel T: {global_max_pellet:7.1f} K | Flow[Center]: {history_W_flow['Center'][-1]:.4f}")

        # 时间推进
        current_time += dt

        # E. 断点自动保存 (严格调用 v3 方法)
        if current_time >= next_save_time:
            save_path = f"test_core_assemble_v4_restart_t{int(next_save_time)}.npz"
            print(f"\n[Checkpoint] 正在写入存档: {save_path}")
            system.save_global_state(save_path)
            next_save_time += save_interval

    print("\n瞬态测试循环完成。")

    # =========================================================================
    # Part 8: 结果绘图 (Result Plotting)
    # =========================================================================
    print("8. Plotting results...")
    if enable_plot:
        # TODO: 创建 plt.figure()
        # TODO: 绘制芯块温度、流道出口温度、全局流量等瞬态曲线
        # TODO: plt.tight_layout(), plt.savefig(), plt.show()
        pass

    # 返回关键信息，便于自动化验证脚本复用
    return {
        'final_time': current_time,
        'history_time': history_time,
        'history_I': history_I,
        'history_U': history_U,
        'history_T_pellet': history_T_pellet,
        'history_T_coolant': history_T_coolant,
        'history_T_mod': history_T_mod,
        'history_W_flow': history_W_flow,
        'core': core,
        'system': system,
        'tfes': tfes,
    }


if __name__ == "__main__":
    run_test_v4()
    TEASAProfiler.report()
