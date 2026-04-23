import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, '..'))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from profiler import TEASAProfiler

thermo_calc_dir = os.path.abspath(os.path.join(current_dir, '..', 'ThermoCalc'))
if thermo_calc_dir not in sys.path:
    sys.path.insert(0, thermo_calc_dir)

from Components.ExternalHeatSources import (
    AlbedoHeatSource,
    CompositeHeatSource,
    EarthIRHeatSource,
    ExternalHeatFluxBC,
    OrbitalHeatSource,
    OrbitalTableHeatSource,
)
from Components.HPwithFin import HPwithFin
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidChannel
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager


def create_dummy_fluid_network() -> HydraulicNetwork:
    """Build the smallest hydraulic network required by SystemManager."""
    mat_fluid = Sodium()

    L_channel = 0.1
    D_inner = 0.010
    N_fluid = 3
    T_init = 600.0
    P_out = 1.58e5
    dP_drive = 100.0
    P_in = P_out + dP_drive
    area_flow = np.pi * (D_inner / 2) ** 2

    inlet_plenum = IncompressibleBoundaryVolume("Dummy_Inlet", mat_fluid, P=P_in, T=T_init)
    inlet_plenum.is_pressure_boundary = True
    outlet_plenum = IncompressibleBoundaryVolume("Dummy_Outlet", mat_fluid, P=P_out, T=T_init)
    outlet_plenum.is_pressure_boundary = True

    dummy_chan = IncompressibleFluidChannel(
        name="Dummy_Chan",
        n_nodes=N_fluid,
        total_length=L_channel,
        flow_area=area_flow,
        hydraulic_diam=D_inner,
        initial_P=P_out,
        initial_T=T_init,
        material=mat_fluid
    )

    j_in = FlowJunction("J_In", inlet_plenum, dummy_chan.volumes[0], flow_area=area_flow)
    j_out = FlowJunction("J_Out", dummy_chan.volumes[-1], outlet_plenum, flow_area=area_flow)

    all_vols = [inlet_plenum, outlet_plenum] + dummy_chan.volumes
    all_juncs = [j_in, j_out] + dummy_chan.internal_junctions

    return HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)


def build_external_heat_source(shape: tuple, external_heat_config: dict) -> CompositeHeatSource:
    """
    根据配置构造轨道外热流源。

    当前支持两大类路线：
    1. 查表法：
       - use_embedded_table=True
       - 使用 OrbitalTableHeatSource
       - 轨道热流随时间按表插值
    2. 解析/半解析法：
       - 通过 add_solar / add_albedo / add_earth_ir 组合
       - 其中太阳热流目前由 OrbitalHeatSource 显式计算
       - 地球反照和地球红外当前仍是简化常值模型
    """
    composite_source = CompositeHeatSource(shape)

    if external_heat_config.get('use_embedded_table', False):
        composite_source.add_source(
            OrbitalTableHeatSource(
                shape=shape,
                table_ids=external_heat_config.get('table_ids', 1),
                scale_factor=external_heat_config.get('table_scale_factor', 1.0),
                offset=external_heat_config.get('table_offset', 0.0),
                periodic=external_heat_config.get('table_periodic', True)
            )
        )
        return composite_source

    if external_heat_config.get('add_solar', False):
        composite_source.add_source(
            OrbitalHeatSource(
                shape=shape,
                solar_constant=external_heat_config.get('solar_constant', 1361.0),
                orbit_height=external_heat_config.get('orbit_height', 800.0),
                orbit_period=external_heat_config.get('orbit_period', 7644.0),
                orbit_inclination=external_heat_config.get('orbit_inclination', 0.0),
                surface_normal_angles=external_heat_config.get('surface_normal_angles', (0.0, 0.0))
            )
        )

    if external_heat_config.get('add_albedo', False):
        composite_source.add_source(
            AlbedoHeatSource(
                shape=shape,
                albedo_factor=external_heat_config.get('albedo_factor', 0.3),
                solar_constant=external_heat_config.get('solar_constant', 1361.0)
            )
        )

    if external_heat_config.get('add_earth_ir', False):
        composite_source.add_source(
            EarthIRHeatSource(
                shape=shape,
                earth_ir_flux=external_heat_config.get('earth_ir_flux', 237.0)
            )
        )

    return composite_source


def create_hp_radiator_with_external_heat(case_name: str, external_heat_config: dict) -> HPwithFin:
    """
    创建一个带轨道外热流的 HPwithFin 散热器。

    这里额外负责处理“翅片外热流如何加载”的两种方式：
    1. lumped_root_area：
       把翅片受照折算成等效吸热面积，直接并入冷凝段外壁边界。
    2. distributed_fin_absorption：
       只有管壁受照部分挂在外壁边界上，翅片受照部分直接送入 HPwithFin 的降维翅片方程。
    """
    T_init = 800.0
    T_eva_ext = 800.0
    T_env = 3.0
    up_vf = 0.0
    down_vf = 0.675
    emissivity = 0.93

    L_eva = 0.06
    n_eva = 1
    L_aba = 0.0
    n_aba = 0
    L_con = 0.482
    n_con = 12

    r_out_wall = 0.0085
    r_in_wall = r_out_wall - 0.0004
    r_vapor = r_in_wall - 0.0006
    porosity = 0.5
    n_wick = 1
    n_wall = 1

    fin_height = 22.65e-3
    fin_thickness = 0.0003
    n_fin = 2
    fin_wrap_ratio = (n_fin * fin_thickness) / (2.0 * np.pi * r_out_wall)

    mat_wall = SS316(name=f"{case_name}_HP_Wall_SS316")
    mat_fluid = SodiumHP(name=f"{case_name}_HP_Fluid_Na")
    mat_wick = WickMaterial(
        name=f"{case_name}_HP_Wick_Composite",
        solid_mat=SS316(),
        fluid_mat=mat_fluid,
        porosity=porosity,
        r_vapor=r_vapor,
        r_in_wall=r_in_wall
    )

    hp_radiator = HPwithFin(
        name=f"{case_name}_Main_Radiator",
        r_out_wall=r_out_wall,
        r_in_wall=r_in_wall,
        r_vapor=r_vapor,
        L_eva=L_eva,
        L_aba=L_aba,
        L_con=L_con,
        n_eva=n_eva,
        n_aba=n_aba,
        n_con=n_con,
        n_wick=n_wick,
        n_wall=n_wall,
        wall_mat=mat_wall,
        fluid_mat=mat_fluid,
        wick_struct_mat=mat_wick,
        porosity=porosity,
        fin_thickness=fin_thickness,
        fin_height=fin_height,
        n_fin_height=15,
        fin_wrap_ratio=fin_wrap_ratio,
        emissivity=emissivity,
        up_view_factor=up_vf,
        down_view_factor=down_vf,
        T_env=T_env,
        initial_temp=T_init
    )

    hp_radiator.hp.boundaries['outer_eva'].add_resistance_condition(
        T_ext=T_eva_ext,
        R_ext=1e-8
    )

    outer_con_boundary = hp_radiator.hp.boundaries['outer_con']
    area_con = outer_con_boundary.area
    shape = (n_con,)

    composite_source = build_external_heat_source(shape, external_heat_config)

    # wall_illumination_factor:
    #   冷凝段圆管外壁真正能接收到轨道外热流的面积比例。
    # fin_illuminated_area_scale:
    #   翅片受照投影面积的缩放系数，便于后续引入遮挡/姿态修正。
    # fin_loading_mode:
    #   选择翅片外热流加载方法。
    wall_illumination_factor = external_heat_config.get('wall_illumination_factor', 0.5)
    fin_illuminated_area_scale = external_heat_config.get('fin_illuminated_area_scale', 1.0)
    fin_loading_mode = external_heat_config.get('fin_loading_mode', 'lumped_root_area')

    wall_absorption_area = area_con * wall_illumination_factor
    fin_absorption_area = hp_radiator.get_fin_illuminated_area_array(fin_illuminated_area_scale)

    if fin_loading_mode == 'lumped_root_area':
        # 简化法：
        # 将“外壁吸热 + 翅片吸热”统一折算成冷凝段外壁边界上的等效输入功率。
        effective_boundary_area = wall_absorption_area + fin_absorption_area
        external_bc = ExternalHeatFluxBC(
            heat_source=composite_source,
            area_array=effective_boundary_area
        )
        outer_con_boundary.conditions.append(external_bc)
    elif fin_loading_mode == 'distributed_fin_absorption':
        # 更严格的方法：
        # 外壁受照仍然走 BoundaryRegion 边界，
        # 翅片受照则直接进入 HPwithFin.pre_step() 内部的翅片方程。
        external_bc = ExternalHeatFluxBC(
            heat_source=composite_source,
            area_array=wall_absorption_area
        )
        outer_con_boundary.conditions.append(external_bc)
        hp_radiator.set_fin_external_heat_source(
            composite_source,
            illuminated_area_scale=fin_illuminated_area_scale
        )
    else:
        raise ValueError(
            "Unsupported fin_loading_mode. Use 'lumped_root_area' or 'distributed_fin_absorption'."
        )

    hp_radiator.configure_external_heat_accounting(
        composite_source,
        wall_area_array=wall_absorption_area,
        fin_area_array=fin_absorption_area
    )

    return hp_radiator


def plot_comparison_results(case_results: dict):
    """Plot transient and final-state comparisons for all cases."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    case_names = list(case_results.keys())
    colors = [plt.get_cmap('tab10')(i % 10) for i in range(len(case_names))]

    for idx, (case_name, result) in enumerate(case_results.items()):
        time_history = result['time_history']
        temp_contour = result['temp_contour']

        if len(time_history) == 0:
            continue

        ax1 = axes[0, 0]
        ax1.plot(result['z_centers'], temp_contour[-1], label=case_name, color=colors[idx], linewidth=2)
        ax1.set_xlabel('Axial Position Z (m)')
        ax1.set_ylabel('Outer Wall Temperature (K)')
        ax1.set_title('Final Outer Wall Temperature Distribution')
        ax1.grid(True, linestyle='--', alpha=0.7)

        ax2 = axes[0, 1]
        ax2.plot(time_history, result['q_total_history'], label=case_name, color=colors[idx])
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Gross Heat Rejection (W)')
        ax2.set_title('Gross Heat Rejection vs Time')
        ax2.grid(True, linestyle='--', alpha=0.7)

        ax3 = axes[1, 0]
        ax3.plot(time_history, result['q_net_history'], label=case_name, color=colors[idx])
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Net Heat Rejection (W)')
        ax3.set_title('Net Heat Rejection vs Time')
        ax3.grid(True, linestyle='--', alpha=0.7)

        ax4 = axes[1, 1]
        final_q = result['q_dist_history'][-1]
        ax4.plot(np.arange(len(final_q)), final_q, marker='o', label=case_name, color=colors[idx], alpha=0.85)
        ax4.set_xlabel('Condenser Node Index')
        ax4.set_ylabel('Gross Heat Rejection (W)')
        ax4.set_title('Final Condenser Heat Rejection Distribution')
        ax4.grid(True, linestyle='--', alpha=0.7)

    for ax in axes.flat:
        ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig('external_heat_source_comparison.png', dpi=150)
    plt.show()


def run_single_case(case_name: str,
                    external_heat_config: dict,
                    t_end: float = 10.0,
                    dt: float = 0.05) -> dict:
    """
    运行单个算例，并记录：
    1. gross heat rejection：总向外散热量
    2. absorbed orbital heat：吸收的轨道外热流
    3. net heat rejection：净排热 = gross - absorbed
    """
    print(f"\n{'=' * 70}")
    print(f"Running Case: {case_name}")
    print(f"External Heat Config: {external_heat_config}")
    print(f"{'=' * 70}")

    hp_radiator = create_hp_radiator_with_external_heat(case_name, external_heat_config)

    dummy_net = create_dummy_fluid_network()
    sys_manager = SystemManager(fluid_network=dummy_net, start_time=0.0)
    sys_manager.add_component(hp_radiator)
    sys_manager.initialize_system(dt_init=dt)

    time_history = []
    temp_contour_history = []
    q_total_history = []
    q_net_history = []
    q_abs_history = []
    q_dist_history = []

    record_interval = 50
    step_count = 0

    print(f"\n[*] Starting transient simulation for {case_name}...")
    start_cpu = time.time()

    while sys_manager.global_time <= t_end + 1e-6:
        sys_manager.step(dt=dt, inner_iter=1)

        current_t = sys_manager.global_time
        _, q_con_total_dist = hp_radiator.get_heat_rejection_distribution()
        _, _, q_abs_total_dist = hp_radiator.get_external_heat_absorption_distribution(current_t)

        total_q = np.sum(q_con_total_dist)
        absorbed_q = np.sum(q_abs_total_dist)
        net_q = total_q - absorbed_q

        step_count += 1
        if step_count % record_interval == 0:
            T_2d = hp_radiator.get_temperature_distribution()
            T_outer_wall = T_2d[-1, :]
            temp_contour_history.append(T_outer_wall.copy())
            time_history.append(current_t)
            q_total_history.append(total_q)
            q_net_history.append(net_q)
            q_abs_history.append(absorbed_q)
            q_dist_history.append(q_con_total_dist.copy())

        if step_count % 20 == 0:
            print(
                f"  Time: {current_t:.3f} s | Gross Heat Rejection: {total_q:.2f} W | "
                f"Absorbed Orbital Heat: {absorbed_q:.2f} W"
            )

    elapsed = time.time() - start_cpu
    print(f"\n[*] {case_name} completed in {elapsed:.2f} seconds")

    return {
        'time_history': time_history,
        'temp_contour': temp_contour_history,
        'q_total_history': q_total_history,
        'q_net_history': q_net_history,
        'q_abs_history': q_abs_history,
        'q_dist_history': q_dist_history,
        'z_centers': hp_radiator.hp_mesh.y_centers,
        'config': dict(external_heat_config)
    }


def _tail_mean(values, n_tail: int = 3) -> float:
    if not values:
        return 0.0
    tail = values[-min(len(values), n_tail):]
    return float(np.mean(tail))


def summarize_pairwise_results(case_results: dict, comparison_pairs: list):
    """输出“简化法 vs 直接法”的对比结果。"""
    print("\nPairwise comparison: Lumped root-area model vs distributed fin absorption")
    print("-" * 90)

    for label, lumped_name, direct_name in comparison_pairs:
        lumped = case_results[lumped_name]
        direct = case_results[direct_name]

        gross_lumped = _tail_mean(lumped['q_total_history'])
        gross_direct = _tail_mean(direct['q_total_history'])
        net_lumped = _tail_mean(lumped['q_net_history'])
        net_direct = _tail_mean(direct['q_net_history'])
        abs_lumped = _tail_mean(lumped['q_abs_history'])
        abs_direct = _tail_mean(direct['q_abs_history'])

        gross_rel_err = 0.0
        if abs(gross_direct) > 1e-12:
            gross_rel_err = abs(gross_lumped - gross_direct) / abs(gross_direct) * 100.0

        net_rel_err = 0.0
        if abs(net_direct) > 1e-12:
            net_rel_err = abs(net_lumped - net_direct) / abs(net_direct) * 100.0

        print(f"{label}:")
        print(f"  Lumped gross rejection      : {gross_lumped:.3f} W")
        print(f"  Direct gross rejection      : {gross_direct:.3f} W")
        print(f"  Gross rejection rel. diff   : {gross_rel_err:.3f} %")
        print(f"  Lumped net rejection        : {net_lumped:.3f} W")
        print(f"  Direct net rejection        : {net_direct:.3f} W")
        print(f"  Net rejection rel. diff     : {net_rel_err:.3f} %")
        print(f"  Lumped absorbed orbital heat: {abs_lumped:.3f} W")
        print(f"  Direct absorbed orbital heat: {abs_direct:.3f} W")
        print()


def main():
    print("=" * 70)
    print("TASTIN HPwithFin orbital heat loading comparison")
    print("=" * 70)

    base_table_config = {
        'use_embedded_table': True,
        'table_scale_factor': 1.0,
        'table_offset': 0.0,
        'table_periodic': True,
        'wall_illumination_factor': 0.5,
        'fin_illuminated_area_scale': 1.0
    }

    # 这里选 3 个不同的轨道热流表号。
    # 每个表号都各跑两遍：
    # 1. Lumped  : 翅片吸热折算回根部边界
    # 2. Direct  : 翅片吸热直接进入翅片方程
    table_ids = [4, 18, 32]
    cases = {}
    comparison_pairs = []

    for idx, table_id in enumerate(table_ids, start=1):
        lumped_name = f"Case{idx}_Table_{table_id:02d}_Lumped"
        direct_name = f"Case{idx}_Table_{table_id:02d}_Direct"

        cases[lumped_name] = {
            **base_table_config,
            'table_ids': table_id,
            'fin_loading_mode': 'lumped_root_area'
        }
        cases[direct_name] = {
            **base_table_config,
            'table_ids': table_id,
            'fin_loading_mode': 'distributed_fin_absorption'
        }
        comparison_pairs.append((f"Table {table_id:02d}", lumped_name, direct_name))

    all_results = {}
    for case_name, config in cases.items():
        all_results[case_name] = run_single_case(case_name, config, t_end=20.0, dt=0.05)

    print("\n[*] Generating comparison plots...")
    plot_comparison_results(all_results)

    print("\n" + "=" * 70)
    print("Simulation completed")
    print("=" * 70)

    print("\nFinal case summary")
    print("-" * 70)
    for case_name, result in all_results.items():
        final_gross = _tail_mean(result['q_total_history'], n_tail=1)
        final_net = _tail_mean(result['q_net_history'], n_tail=1)
        final_abs = _tail_mean(result['q_abs_history'], n_tail=1)
        final_T = result['temp_contour'][-1] if result['temp_contour'] else []
        avg_T = float(np.mean(final_T)) if len(final_T) > 0 else 0.0
        print(f"{case_name}:")
        print(f"  Gross heat rejection : {final_gross:.3f} W")
        print(f"  Net heat rejection   : {final_net:.3f} W")
        print(f"  Absorbed orbital heat: {final_abs:.3f} W")
        print(f"  Average wall temp    : {avg_T:.3f} K")
        print()

    summarize_pairwise_results(all_results, comparison_pairs)


if __name__ == '__main__':
    main()
    TEASAProfiler.report()
