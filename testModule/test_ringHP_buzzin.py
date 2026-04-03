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

from Components.RingHP import RingHP
from Materials.Fluids.Sodium import Sodium
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Solvers.Hydrodynamics.BoundaryVolume import IncompressibleBoundaryVolume, InletJunction
from Solvers.Hydrodynamics.Components import FlowJunction, IncompressibleFluidChannel
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.SystemManager import SystemManager


def finalize_figure(fig, filename: str):
    """
    统一处理测试绘图输出。

    在交互式后端下直接 `show()`；
    在 `Agg` 等非交互式后端下自动保存到当前测试目录，避免运行级验证时出现无意义告警。
    """
    backend = plt.get_backend().lower()
    if 'agg' in backend:
        out_path = os.path.join(current_dir, filename)
        fig.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"[*] 图像已保存: {out_path}")
    else:
        plt.show()


def build_ringhp_case(external_heat_config=None):
    """
    构造一个单节点 RingHP 算例。

    参数 external_heat_config 直接透传给 RingHP，用于控制：
    1. 是否启用轨道外热流
    2. 使用查表法还是解析/半解析法
    3. 翅片外热流采用折算法还是直接受照法
    """
    # ---------------------------------------------------------
    # 1. 几何与工况参数
    # ---------------------------------------------------------
    T_in = 863.0
    P_out = 161000.0
    W_in = 0.1815
    T_env = 3.0

    L_header = 0.793
    area_flow = 0.0016065
    D_h = 0.04167
    wall_thickness = 0.002

    r_in_header = D_h / 2.0
    perim_header = 2.0 * np.pi * r_in_header

    r_out_hp = 0.0085
    r_in_hp = 0.0081
    r_vapor = 0.0075
    L_eva = 0.06
    n_eva = 1
    L_con = 0.482
    n_con = 12
    n_wick = 1
    n_wall = 1
    porosity = 0.5

    fin_thickness = 0.0003
    fin_height = 22.65e-3
    n_fin_height = 15
    fin_wrap_ratio = (2 * fin_thickness) / (2.0 * np.pi * r_out_hp)
    emissivity = 0.93
    up_vf = 1.0
    down_vf = 0.675

    # ---------------------------------------------------------
    # 2. 物性
    # ---------------------------------------------------------
    mat_fluid = Sodium()
    mat_wall = SS316(name="SS316_Wall")
    mat_hp_fluid = SodiumHP(name="HP_Fluid_Na")
    mat_wick = WickMaterial(
        name="HP_Wick_Composite",
        solid_mat=SS316(),
        fluid_mat=mat_hp_fluid,
        porosity=porosity,
        r_vapor=r_vapor,
        r_in_wall=r_in_hp
    )

    # ---------------------------------------------------------
    # 3. 水力网络
    # ---------------------------------------------------------
    inlet_plenum = IncompressibleBoundaryVolume("Inlet_Plenum", mat_fluid, P=P_out + 1000, T=T_in)
    outlet_plenum = IncompressibleBoundaryVolume("Outlet_Plenum", mat_fluid, P=P_out, T=T_in)
    outlet_plenum.is_pressure_boundary = True

    chan = IncompressibleFluidChannel(
        name="Header_Channel",
        n_nodes=1,
        total_length=L_header,
        flow_area=area_flow,
        hydraulic_diam=D_h,
        initial_P=P_out + 500,
        initial_T=T_in,
        material=mat_fluid
    )

    j_in = InletJunction("J_in", inlet_plenum, chan.volumes[0], W_initial=W_in)
    j_out = FlowJunction("J_out", chan.volumes[-1], outlet_plenum, flow_area=area_flow, k_loss=0.0)

    # ---------------------------------------------------------
    # 4. 集流环固体壁面
    # ---------------------------------------------------------
    mesh_header = Mesh2D(
        x_dim=wall_thickness,
        n_x=1,
        y_dim=L_header,
        n_y=1,
        geometry_type='cylindrical',
        inner_radius=r_in_header
    )
    solid_header = HeatConduction2D(mesh=mesh_header, material=mat_wall, initial_temp=T_in)
    solid_header.name = "Header_Outer_Wall"

    solid_header.boundaries['right'].add_flux_condition(q_flux=0.0)
    solid_header.boundaries['top'].add_flux_condition(q_flux=0.0)
    solid_header.boundaries['bottom'].add_flux_condition(q_flux=0.0)

    # ---------------------------------------------------------
    # 5. 集流环内壁换热相关式
    # ---------------------------------------------------------
    def lyon_martinelli(Re, Pr, P_D_ratio=1.0):
        pe = max(Re * Pr, 1.0)
        return 7.0 + 0.025 * (pe ** 0.8)

    ring_hp = RingHP(
        name="Main_RingHP",
        fluid_channel=chan,
        solid_header=solid_header,
        hp_multipliers=[26.0],
        header_flow_area=area_flow,
        header_dh=D_h,
        header_heated_perimeter=perim_header,
        hp_r_out=r_out_hp,
        hp_r_in=r_in_hp,
        hp_r_vapor=r_vapor,
        hp_L_eva=L_eva,
        hp_L_con=L_con,
        hp_n_eva=n_eva,
        hp_n_con=n_con,
        hp_n_wick=n_wick,
        hp_n_wall=n_wall,
        porosity_hp=porosity,
        HP_initial_temp=800.0,
        fin_thickness=fin_thickness,
        fin_height=fin_height,
        n_fin_height=n_fin_height,
        fin_wrap_ratio=fin_wrap_ratio,
        emissivity=emissivity,
        up_view_factor=up_vf,
        down_view_factor=down_vf,
        T_space=T_env,
        hp_wall_mat=mat_wall,
        hp_fluid_mat=mat_hp_fluid,
        hp_wick_mat=mat_wick,
        header_correlation_func=lyon_martinelli,
        hp_crossflow_base_func=lambda *args: 10.0,
        C_D=1.0,
        external_heat_config=external_heat_config,
    )

    j_out.k_loss = ring_hp.outlet_k_loss

    all_vols = [inlet_plenum, outlet_plenum] + chan.volumes
    all_juncs = [j_in, j_out]
    net = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)

    sys_mgr = SystemManager(fluid_network=net)
    sys_mgr.add_component(ring_hp)

    return {
        'system_manager': sys_mgr,
        'ring_hp': ring_hp,
        'channel': chan,
        'solid_header': solid_header,
        'inlet_junction': j_in,
        'outlet_junction': j_out,
    }


def run_ringhp_case(case_name: str, external_heat_config=None, t_end: float = 50.0):
    """运行一个 RingHP 算例，并记录与轨道外热流相关的关键量。"""
    print("=" * 70)
    print(f"Running RingHP case: {case_name}")
    print(f"External heat config: {external_heat_config}")
    print("=" * 70)

    case = build_ringhp_case(external_heat_config=external_heat_config)
    sys_mgr = case['system_manager']
    ring_hp = case['ring_hp']
    chan = case['channel']
    solid_header = case['solid_header']
    j_in = case['inlet_junction']

    print(f"[*] 注入阻塞阻力系数: K_loss = {case['outlet_junction'].k_loss:.4f}")

    sys_mgr.initialize_system()

    history = {
        'time': [],
        'W_in': [],
        'T_fluid': [],
        'T_header_wall': [],
        'Q_reject_gross': [],
        'Q_external_absorb': [],
        'Q_reject_net': [],
    }

    print("\n" + "=" * 50)
    print(f"开始瞬态积分: {case_name}")
    print("=" * 50)

    start_cpu = time.time()
    while sys_mgr.global_time <= t_end:
        dt = sys_mgr.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, safety_factor=0.5)
        sys_mgr.step(dt=dt, inner_iter=2)

        current_t = sys_mgr.global_time
        current_W = j_in.W
        current_T_f = chan.volumes[0].T
        current_T_wall = np.mean(solid_header.boundaries['left'].T_surface)

        gross_q = ring_hp.get_total_heat_rejection()
        external_q = ring_hp.get_total_external_heat_absorption(current_t)
        net_q = gross_q - external_q

        history['time'].append(current_t)
        history['W_in'].append(current_W)
        history['T_fluid'].append(current_T_f)
        history['T_header_wall'].append(current_T_wall)
        history['Q_reject_gross'].append(gross_q)
        history['Q_external_absorb'].append(external_q)
        history['Q_reject_net'].append(net_q)

        if len(history['time']) % 10 == 0:
            print(
                f" t = {current_t:>5.2f}s | W = {current_W:.4f} kg/s | "
                f"T_fluid = {current_T_f:.2f} K | Q_gross = {gross_q:.2f} W | "
                f"Q_abs = {external_q:.2f} W | Q_net = {net_q:.2f} W"
            )

    elapsed = time.time() - start_cpu
    print(f"[*] {case_name} completed in {elapsed:.2f} s")
    print()

    return history


def plot_case_results(case_name: str, history: dict):
    """绘制单个 RingHP 算例结果。"""
    times = np.array(history['time'])

    fig, axes = plt.subplots(4, 1, figsize=(10, 14), sharex=True)

    axes[0].plot(times, history['W_in'], 'b-', linewidth=2)
    axes[0].set_title(f'{case_name}: Inlet Mass Flow Rate', fontweight='bold')
    axes[0].set_ylabel('Mass Flow [kg/s]')
    axes[0].grid(True, linestyle='--')

    axes[1].plot(times, history['T_fluid'], 'r-', linewidth=2, label='Fluid Temp')
    axes[1].plot(times, history['T_header_wall'], 'g--', linewidth=2, label='Header Wall Temp')
    axes[1].set_title(f'{case_name}: Temperatures Evolution', fontweight='bold')
    axes[1].set_ylabel('Temperature [K]')
    axes[1].legend()
    axes[1].grid(True, linestyle='--')

    axes[2].plot(times, history['Q_reject_gross'], color='orange', linewidth=2, label='Gross Rejection')
    axes[2].plot(times, history['Q_external_absorb'], color='purple', linewidth=2, label='Absorbed External Heat')
    axes[2].plot(times, history['Q_reject_net'], color='black', linewidth=2, linestyle='--', label='Net Rejection')
    axes[2].set_title(f'{case_name}: RingHP Heat Balance', fontweight='bold')
    axes[2].set_ylabel('Heat Power [W]')
    axes[2].legend()
    axes[2].grid(True, linestyle='--')

    axes[3].plot(times, history['T_fluid'], 'r-', linewidth=2.5, label='Fluid Outlet Temperature')
    final_t = history['T_fluid'][-1]
    axes[3].axhline(y=final_t, color='gray', linestyle='--', alpha=0.7)
    axes[3].text(times[-1] * 0.65, final_t + 0.5, f'Steady State: {final_t:.2f} K', fontsize=10)
    axes[3].set_title(f'{case_name}: Fluid Outlet Temperature', fontweight='bold')
    axes[3].set_xlabel('Time [s]')
    axes[3].set_ylabel('Temperature [K]')
    axes[3].legend()
    axes[3].grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    safe_case_name = case_name.replace(' ', '_')
    finalize_figure(fig, f'{safe_case_name}_ringhp_history.png')


def plot_comparison(case_histories: dict):
    """绘制基线算例与外热流算例的对比图。"""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    colors = [plt.get_cmap('tab10')(i % 10) for i in range(len(case_histories))]

    for idx, (case_name, history) in enumerate(case_histories.items()):
        times = np.array(history['time'])

        axes[0, 0].plot(times, history['T_fluid'], linewidth=2, color=colors[idx], label=case_name)
        axes[0, 1].plot(times, history['Q_reject_gross'], linewidth=2, color=colors[idx], label=case_name)
        axes[1, 0].plot(times, history['Q_external_absorb'], linewidth=2, color=colors[idx], label=case_name)
        axes[1, 1].plot(times, history['Q_reject_net'], linewidth=2, color=colors[idx], label=case_name)

    axes[0, 0].set_title('Fluid Temperature')
    axes[0, 0].set_ylabel('Temperature [K]')
    axes[0, 0].grid(True, linestyle='--', alpha=0.7)

    axes[0, 1].set_title('Gross Heat Rejection')
    axes[0, 1].set_ylabel('Heat Power [W]')
    axes[0, 1].grid(True, linestyle='--', alpha=0.7)

    axes[1, 0].set_title('Absorbed External Heat')
    axes[1, 0].set_xlabel('Time [s]')
    axes[1, 0].set_ylabel('Heat Power [W]')
    axes[1, 0].grid(True, linestyle='--', alpha=0.7)

    axes[1, 1].set_title('Net Heat Rejection')
    axes[1, 1].set_xlabel('Time [s]')
    axes[1, 1].set_ylabel('Heat Power [W]')
    axes[1, 1].grid(True, linestyle='--', alpha=0.7)

    for ax in axes.flat:
        ax.legend()

    plt.tight_layout()
    finalize_figure(fig, 'ringhp_case_comparison.png')


@TEASAProfiler.profile
def main():
    print("==========================================================")
    print("   TEST: RingHP Baseline & External Orbital Heat Cases    ")
    print("==========================================================")

    baseline_history = run_ringhp_case(
        case_name='Baseline_No_External_Heat',
        external_heat_config=None,
        t_end=50.0,
    )

    orbital_external_heat_config = {
        'use_embedded_table': True,
        'table_ids': 4,
        'table_scale_factor': 1.0,
        'table_offset': 0.0,
        'table_periodic': True,
        'wall_illumination_factor': 0.5,
        'fin_illuminated_area_scale': 1.0,
        'fin_loading_mode': 'distributed_fin_absorption',
    }
    external_history = run_ringhp_case(
        case_name='Case_With_Orbital_Heat_Table04',
        external_heat_config=orbital_external_heat_config,
        t_end=50.0,
    )

    print("Final summary")
    print("-" * 70)
    for case_name, history in {
        'Baseline_No_External_Heat': baseline_history,
        'Case_With_Orbital_Heat_Table04': external_history,
    }.items():
        print(f"{case_name}:")
        print(f"  Fluid outlet temperature : {history['T_fluid'][-1]:.3f} K")
        print(f"  Gross heat rejection     : {history['Q_reject_gross'][-1]:.3f} W")
        print(f"  Absorbed external heat   : {history['Q_external_absorb'][-1]:.3f} W")
        print(f"  Net heat rejection       : {history['Q_reject_net'][-1]:.3f} W")
        print()

    plot_case_results('Baseline_No_External_Heat', baseline_history)
    plot_case_results('Case_With_Orbital_Heat_Table04', external_history)
    plot_comparison({
        'Baseline_No_External_Heat': baseline_history,
        'Case_With_Orbital_Heat_Table04': external_history,
    })


if __name__ == '__main__':
    main()
    TEASAProfiler.report()
