import argparse
import csv
import os
import sys
import time

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

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
    """Same dummy network as testHPwithFin.py."""
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
        material=mat_fluid,
    )

    j_in = FlowJunction("J_In", inlet_plenum, dummy_chan.volumes[0], flow_area=area_flow)
    j_out = FlowJunction("J_Out", dummy_chan.volumes[-1], outlet_plenum, flow_area=area_flow)

    all_vols = [inlet_plenum, outlet_plenum] + dummy_chan.volumes
    all_juncs = [j_in, j_out] + dummy_chan.internal_junctions

    return HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)


def build_hp_radiator() -> HPwithFin:
    """Reproduce the HPwithFin setup from testHPwithFin.py."""
    T_init = 800.0
    T_eva_ext = 800.0
    T_env = 3.0
    up_vf = 0.0
    down_vf = 0.3
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

    mat_wall = SS316(name="HP_Wall_SS316")
    mat_fluid = SodiumHP(name="HP_Fluid_Na")
    mat_wick = WickMaterial(
        name="HP_Wick_Composite",
        solid_mat=SS316(),
        fluid_mat=mat_fluid,
        porosity=porosity,
        r_vapor=r_vapor,
        r_in_wall=r_in_wall,
    )

    hp_radiator = HPwithFin(
        name="Main_Radiator",
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
        initial_temp=T_init,
    )

    hp_radiator.hp.boundaries["outer_eva"].add_resistance_condition(
        T_ext=T_eva_ext,
        R_ext=1e-8,
    )

    return hp_radiator


def write_csv(path: str, rows):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def save_restart(path: str, hp_radiator: HPwithFin, rows):
    if not path or not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    final = rows[-1]
    np.savez_compressed(
        path,
        hp_T=np.array(hp_radiator.hp.T, copy=True),
        hp_current_time=np.array([float(hp_radiator.hp.current_time)]),
        audit_time=np.array([float(final["t"])]),
        Q_con_out=np.array([float(final["Q_con_out"])]),
        T_hp_min=np.array([float(final["T_hp_min"])]),
        T_hp_mean=np.array([float(final["T_hp_mean"])]),
        T_hp_max=np.array([float(final["T_hp_max"])]),
        last_fin_temperature=np.array(hp_radiator.last_fin_temperature, copy=True),
        last_fin_radiation_distribution=np.array(
            hp_radiator.last_fin_radiation_distribution, copy=True
        ),
        last_fin_absorption_distribution=np.array(
            hp_radiator.last_fin_absorption_distribution, copy=True
        ),
        last_fin_net_from_root_distribution=np.array(
            hp_radiator.last_fin_net_from_root_distribution, copy=True
        ),
    )
    print(f"Restart written to: {path}")


def load_restart(path: str, hp_radiator: HPwithFin) -> float:
    with np.load(path, allow_pickle=False) as data:
        saved_T = data["hp_T"]
        if saved_T.shape != hp_radiator.hp.T.shape:
            raise ValueError(
                f"Restart HP shape {saved_T.shape} does not match current "
                f"shape {hp_radiator.hp.T.shape}."
            )
        hp_radiator.hp.T[:] = saved_T
        if "hp_current_time" in data:
            hp_radiator.hp.current_time = float(data["hp_current_time"][0])
        if "last_fin_temperature" in data:
            hp_radiator.last_fin_temperature[:] = data["last_fin_temperature"]
        if "last_fin_radiation_distribution" in data:
            hp_radiator.last_fin_radiation_distribution[:] = data[
                "last_fin_radiation_distribution"
            ]
        if "last_fin_absorption_distribution" in data:
            hp_radiator.last_fin_absorption_distribution[:] = data[
                "last_fin_absorption_distribution"
            ]
        if "last_fin_net_from_root_distribution" in data:
            hp_radiator.last_fin_net_from_root_distribution[:] = data[
                "last_fin_net_from_root_distribution"
            ]
        hp_radiator.hp._update_properties()
        hp_radiator.hp._compute_internal_resistance()
        hp_radiator.hp._update_boundaries_state(current_time=0.0)
        return float(data["Q_con_out"][0])


def clear_boundary_conditions(boundary):
    boundary.conditions.clear()
    boundary.clear_boundary_conditions()


def attach_evaporator_power_boundary(hp_radiator: HPwithFin, total_power: float):
    """Attach total evaporator input power in W, distributed by boundary area."""
    boundary = hp_radiator.hp.boundaries["outer_eva"]
    clear_boundary_conditions(boundary)
    area = np.asarray(boundary.area, dtype=float)
    weights = area / np.sum(area)
    # BoundaryRegion.add_flux_condition expects power per boundary cell [W],
    # not heat flux density [W/m^2].
    q_power_by_cell = total_power * weights
    boundary.add_flux_condition(q_flux=q_power_by_cell)


def run_case(
    t_end: float,
    dt: float,
    csv_path: str = None,
    restart_path: str = None,
    print_every: int = 10,
):
    print("=== Test 3 Phase 0: reproduce testHPwithFin.py case ===")

    hp_radiator = build_hp_radiator()
    dummy_net = create_dummy_fluid_network()
    sys_manager = SystemManager(fluid_network=dummy_net, start_time=0.0)
    sys_manager.add_component(hp_radiator)
    sys_manager.initialize_system(dt_init=dt)

    rows = []
    step_count = 0
    start_cpu = time.time()

    print("\n" + "=" * 50)
    print("Start transient integration with the testHPwithFin.py setup")
    print("=" * 50)

    while sys_manager.global_time <= t_end + 1.0e-6:
        sys_manager.step(dt=dt, inner_iter=1)

        current_t = sys_manager.global_time
        _, q_con_total_dist = hp_radiator.get_heat_rejection_distribution()
        total_q = float(np.sum(q_con_total_dist))
        T_2d = hp_radiator.get_temperature_distribution()
        T_outer_wall = T_2d[-1, :]

        step_count += 1
        row = {
            "step": step_count,
            "t": float(current_t),
            "Q_con_out": total_q,
            "T_hp_min": float(np.min(hp_radiator.hp.T)),
            "T_hp_mean": float(np.mean(hp_radiator.hp.T)),
            "T_hp_max": float(np.max(hp_radiator.hp.T)),
            "T_outer_min": float(np.min(T_outer_wall)),
            "T_outer_mean": float(np.mean(T_outer_wall)),
            "T_outer_max": float(np.max(T_outer_wall)),
        }
        rows.append(row)

        if step_count == 1 or step_count % max(print_every, 1) == 0:
            print(
                f"  time: {current_t:.3f} s | "
                f"condenser heat rejection: {total_q:.2f} W | "
                f"Tmin: {row['T_hp_min']:.2f} K"
            )

    elapsed = time.time() - start_cpu
    print("\n" + "-" * 50)
    print(f"Completed {step_count} steps in {elapsed:.2f} s")
    if rows:
        final = rows[-1]
        print(f"Final Q_con_out: {final['Q_con_out']:.6f} W")
        print(f"Final Tmin:      {final['T_hp_min']:.6f} K")

    if csv_path:
        write_csv(csv_path, rows)
        print(f"CSV written to: {csv_path}")

    if restart_path:
        save_restart(restart_path, hp_radiator, rows)

    return rows


def run_power_restart_case(
    restart_path: str,
    input_power: float,
    t_end: float,
    dt: float,
    csv_path: str = None,
    restart_out_path: str = None,
    print_every: int = 10,
):
    print("=== Test 3 Phase 1: restart with evaporator fixed power input ===")
    hp_radiator = build_hp_radiator()
    restart_power = load_restart(restart_path, hp_radiator)
    attach_evaporator_power_boundary(hp_radiator, input_power)

    dummy_net = create_dummy_fluid_network()
    sys_manager = SystemManager(fluid_network=dummy_net, start_time=0.0)
    sys_manager.add_component(hp_radiator)
    sys_manager.initialize_system(dt_init=dt)

    rows = []
    step_count = 0
    start_cpu = time.time()

    print(f"Restart source Q_con_out: {restart_power:.6f} W")
    print(f"Applied evaporator power: {input_power:.6f} W")
    print("Note: evaporator power boundary is in W per boundary cell, not W/m^2.")
    print("\n" + "=" * 50)
    print("Start fixed-power restart integration")
    print("=" * 50)

    while sys_manager.global_time <= t_end + 1.0e-6:
        sys_manager.step(dt=dt, inner_iter=1)

        current_t = sys_manager.global_time
        _, q_con_total_dist = hp_radiator.get_heat_rejection_distribution()
        total_q = float(np.sum(q_con_total_dist))
        eva_power = float(np.sum(hp_radiator.hp.boundaries["outer_eva"].current_flux))
        T_2d = hp_radiator.get_temperature_distribution()
        T_outer_wall = T_2d[-1, :]

        step_count += 1
        row = {
            "step": step_count,
            "t": float(current_t),
            "Q_eva_input_boundary": eva_power,
            "Q_con_out": total_q,
            "Q_balance_simple": eva_power - total_q,
            "T_hp_min": float(np.min(hp_radiator.hp.T)),
            "T_hp_mean": float(np.mean(hp_radiator.hp.T)),
            "T_hp_max": float(np.max(hp_radiator.hp.T)),
            "T_outer_min": float(np.min(T_outer_wall)),
            "T_outer_mean": float(np.mean(T_outer_wall)),
            "T_outer_max": float(np.max(T_outer_wall)),
        }
        rows.append(row)

        if step_count == 1 or step_count % max(print_every, 1) == 0:
            print(
                f"  time: {current_t:.3f} s | "
                f"Qeva: {eva_power:.2f} W | "
                f"Qout: {total_q:.2f} W | "
                f"diff: {eva_power - total_q:.3e} W | "
                f"Tmin: {row['T_hp_min']:.2f} K"
            )

    elapsed = time.time() - start_cpu
    print("\n" + "-" * 50)
    print(f"Completed {step_count} steps in {elapsed:.2f} s")
    if rows:
        final = rows[-1]
        print(f"Final Q_eva:    {final['Q_eva_input_boundary']:.6f} W")
        print(f"Final Q_con_out:{final['Q_con_out']:.6f} W")
        print(f"Final diff:     {final['Q_balance_simple']:.6e} W")
        print(f"Final Tmin:     {final['T_hp_min']:.6f} K")

    if csv_path:
        write_csv(csv_path, rows)
        print(f"CSV written to: {csv_path}")

    if restart_out_path:
        save_restart(restart_out_path, hp_radiator, rows)

    return rows


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reproduce the testHPwithFin.py case and record condenser heat rejection."
    )
    parser.add_argument("--t-end", type=float, default=10.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--print-every", type=int, default=10)
    parser.add_argument(
        "--csv-path",
        default=os.path.join(current_dir, "single_hp_fin_testhpwithfin_repro.csv"),
    )
    parser.add_argument(
        "--restart-path",
        default=os.path.join(current_dir, "single_hp_fin_testhpwithfin_repro_restart.npz"),
    )
    parser.add_argument("--restart-out-path", default=None)
    parser.add_argument("--no-csv", action="store_true")
    parser.add_argument("--no-restart", action="store_true")
    parser.add_argument(
        "--mode",
        choices=("isothermal", "power-restart"),
        default="isothermal",
    )
    parser.add_argument(
        "--input-power",
        type=float,
        default=1449.687439,
        help="Total evaporator input power [W] for power-restart mode.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.mode == "power-restart":
        run_power_restart_case(
            restart_path=args.restart_path,
            input_power=args.input_power,
            t_end=args.t_end,
            dt=args.dt,
            csv_path=None if args.no_csv else args.csv_path,
            restart_out_path=None if args.no_restart else args.restart_out_path,
            print_every=args.print_every,
        )
    else:
        run_case(
            t_end=args.t_end,
            dt=args.dt,
            csv_path=None if args.no_csv else args.csv_path,
            restart_path=None if args.no_restart else args.restart_path,
            print_every=args.print_every,
        )
