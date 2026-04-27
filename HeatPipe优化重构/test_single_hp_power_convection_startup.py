import argparse
import csv
import os
import sys
import time

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from Components.basicComponents.HeatPipe2D import HeatPipe2D
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Solvers.HeatConduction.Mesh import Mesh2D


def clear_boundary_conditions(boundary):
    boundary.conditions.clear()
    boundary.clear_boundary_conditions()


def attach_evaporator_power_boundary(hp: HeatPipe2D, total_power: float):
    boundary = hp.boundaries["outer_eva"]
    clear_boundary_conditions(boundary)
    area = np.asarray(boundary.area, dtype=float)
    weights = area / np.sum(area)
    q_power_by_cell = total_power * weights
    boundary.add_flux_condition(q_flux=q_power_by_cell)


def attach_condenser_convection_boundary(hp: HeatPipe2D, coolant_temp: float, h_coeff: float):
    boundary = hp.boundaries["outer_con"]
    clear_boundary_conditions(boundary)
    boundary.add_convection_condition(T_fluid=coolant_temp, h_coeff=h_coeff)


def attach_adiabatic_boundary(boundary, reference_temp: float = 293.15):
    clear_boundary_conditions(boundary)
    boundary.add_resistance_condition(T_ext=reference_temp, R_ext=1.0e15)


def build_startup_hp(initial_temp: float,
                     input_power: float,
                     condenser_h: float,
                     condenser_temp: float,
                     n_eva: int = 6,
                     n_aba: int = 5,
                     n_con: int = 10,
                     n_wick: int = 1,
                     n_wall: int = 1) -> HeatPipe2D:
    """
    直接构造裸热管启动算例。
    轴向离散默认按约 0.05 m 一个单元：
    - eva: 0.265 m -> 6 cells
    - aba: 0.235 m -> 5 cells
    - con: 0.500 m -> 10 cells
    """
    L_eva = 0.265
    L_aba = 0.235
    L_con = 0.500

    r_out_wall = 0.025 / 2.0
    wall_thickness = 0.0015
    wick_thickness = 0.00064
    r_in_wall = r_out_wall - wall_thickness
    r_vapor = r_in_wall - wick_thickness
    porosity = 0.966

    faces_eva = np.linspace(0.0, L_eva, n_eva + 1)
    faces_aba = np.linspace(L_eva, L_eva + L_aba, n_aba + 1)[1:]
    faces_con = np.linspace(L_eva + L_aba, L_eva + L_aba + L_con, n_con + 1)[1:]
    y_faces_custom = np.concatenate([faces_eva, faces_aba, faces_con])

    faces_wick = np.linspace(r_vapor, r_in_wall, n_wick + 1)
    faces_wall = np.linspace(r_in_wall, r_out_wall, n_wall + 1)[1:]
    x_faces_custom = np.concatenate([faces_wick, faces_wall])

    mesh = Mesh2D(
        x_dim=0.0,
        n_x=n_wick + n_wall,
        y_dim=0.0,
        n_y=n_eva + n_aba + n_con,
        geometry_type='cylindrical',
        inner_radius=r_vapor,
        y_faces=y_faces_custom,
        x_faces=x_faces_custom
    )

    mat_wall = SS316(name="Startup_HP_Wall_SS316")
    mat_fluid = SodiumHP(name="Startup_HP_Fluid_Na")
    mat_wick_struct = SS316(name="Startup_HP_Wick_Structure_SS316")

    hp = HeatPipe2D(
        mesh=mesh,
        solid1=mat_wall,
        solid2=mat_fluid,
        solid3=mat_wick_struct,
        n_wick=n_wick,
        porosity=porosity,
        n_eva=n_eva,
        n_aba=n_aba,
        n_con=n_con,
        name="Single_HP_Power_Convection_Startup",
        initial_temp=initial_temp
    )

    attach_evaporator_power_boundary(hp, input_power)
    attach_condenser_convection_boundary(hp, condenser_temp, condenser_h)
    attach_adiabatic_boundary(hp.boundaries["outer_aba"], reference_temp=condenser_temp)
    attach_adiabatic_boundary(hp.boundaries["left"], reference_temp=condenser_temp)
    attach_adiabatic_boundary(hp.boundaries["top"], reference_temp=condenser_temp)
    attach_adiabatic_boundary(hp.boundaries["bottom"], reference_temp=condenser_temp)

    hp.T.fill(initial_temp)
    hp.current_time = 0.0
    hp.initialize_state()
    return hp


def resolve_axial_cell_counts(axial_refine_factor: float,
                              n_eva: int = None,
                              n_aba: int = None,
                              n_con: int = None):
    if axial_refine_factor <= 0.0:
        raise ValueError(f"axial_refine_factor must be positive, got {axial_refine_factor}")

    base_n_eva = 6
    base_n_aba = 5
    base_n_con = 10

    resolved_n_eva = int(np.ceil(base_n_eva * axial_refine_factor)) if n_eva is None else int(n_eva)
    resolved_n_aba = int(np.ceil(base_n_aba * axial_refine_factor)) if n_aba is None else int(n_aba)
    resolved_n_con = int(np.ceil(base_n_con * axial_refine_factor)) if n_con is None else int(n_con)

    if min(resolved_n_eva, resolved_n_aba, resolved_n_con) <= 0:
        raise ValueError(
            f"Invalid axial mesh counts: n_eva={resolved_n_eva}, "
            f"n_aba={resolved_n_aba}, n_con={resolved_n_con}"
        )

    return resolved_n_eva, resolved_n_aba, resolved_n_con


def write_csv(path: str, rows):
    if not path or not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_npz(path: str, payload: dict):
    if not path:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, **payload)


def run_case(t_end: float,
             dt: float,
             initial_temp: float = 293.15,
             input_power: float = 1800.0,
             condenser_h: float = 147.67,
             condenser_temp: float = 293.15,
             axial_refine_factor: float = 1.0,
             n_eva: int = None,
             n_aba: int = None,
             n_con: int = None,
             print_every: int = 100,
             record_every: int = 1,
             csv_path: str = None,
             npz_path: str = None):
    resolved_n_eva, resolved_n_aba, resolved_n_con = resolve_axial_cell_counts(
        axial_refine_factor=axial_refine_factor,
        n_eva=n_eva,
        n_aba=n_aba,
        n_con=n_con,
    )

    hp = build_startup_hp(
        initial_temp=initial_temp,
        input_power=input_power,
        condenser_h=condenser_h,
        condenser_temp=condenser_temp,
        n_eva=resolved_n_eva,
        n_aba=resolved_n_aba,
        n_con=resolved_n_con,
    )

    rows = []
    time_history = []
    temperature_history = []
    outer_wall_temperature_history = []
    wick_k_history = []
    q_eva_input_distribution_history = []
    q_aba_loss_distribution_history = []
    q_con_rejection_distribution_history = []
    q_con_boundary_inflow_distribution_history = []
    step_count = 0
    start_cpu = time.time()

    print("=" * 78)
    print("Single HP startup with evaporator power and condenser convection")
    print(f"initial_temp = {initial_temp:.3f} K")
    print(f"input_power = {input_power:.3f} W")
    print(f"condenser_h = {condenser_h:.5f} W/m^2/K")
    print(f"condenser_temp = {condenser_temp:.3f} K")
    print(
        f"axial mesh = n_eva {resolved_n_eva}, "
        f"n_aba {resolved_n_aba}, n_con {resolved_n_con}"
    )
    print(f"t_end = {t_end:.3f} s, dt = {dt:.3f} s")
    print("=" * 78)

    while hp.current_time + 0.5 * dt <= t_end:
        ok = hp.step(dt=dt)
        if not ok:
            raise RuntimeError(f"HeatPipe2D step failed at t={hp.current_time:.6f} s")

        current_t = hp.current_time
        T_2d = hp.T.reshape(hp.shape_nodes)
        k_2d = hp.k_node.reshape(hp.shape_nodes)
        wick_k = k_2d[:hp.n_wick, :]

        q_eva_in_dist = np.array(hp.boundaries["outer_eva"].current_flux, copy=True)
        q_aba_loss_dist = -np.array(hp.boundaries["outer_aba"].current_flux, copy=True)
        q_con_boundary_inflow_dist = np.array(hp.boundaries["outer_con"].current_flux, copy=True)
        q_con_rejection_dist = -q_con_boundary_inflow_dist

        row = {
            "step": step_count + 1,
            "t": float(current_t),
            "Q_eva_in": float(np.sum(q_eva_in_dist)),
            "Q_aba_loss": float(np.sum(q_aba_loss_dist)),
            "Q_con_rejection": float(np.sum(q_con_rejection_dist)),
            "T_hp_min": float(np.min(hp.T)),
            "T_hp_mean": float(np.mean(hp.T)),
            "T_hp_max": float(np.max(hp.T)),
            "T_eva_wall_mean": float(np.mean(T_2d[-1, hp._slice_eva])),
            "T_aba_wall_mean": float(np.mean(T_2d[-1, hp._slice_aba])),
            "T_con_wall_mean": float(np.mean(T_2d[-1, hp._slice_con])),
            "k_wick_min": float(np.min(wick_k)),
            "k_wick_mean": float(np.mean(wick_k)),
            "k_wick_max": float(np.max(wick_k)),
        }
        rows.append(row)
        step_count += 1

        if step_count % max(int(record_every), 1) == 0:
            time_history.append(float(current_t))
            temperature_history.append(np.array(T_2d, copy=True))
            outer_wall_temperature_history.append(np.array(T_2d[-1, :], copy=True))
            wick_k_history.append(np.array(wick_k, copy=True))
            q_eva_input_distribution_history.append(q_eva_in_dist)
            q_aba_loss_distribution_history.append(q_aba_loss_dist)
            q_con_rejection_distribution_history.append(q_con_rejection_dist)
            q_con_boundary_inflow_distribution_history.append(q_con_boundary_inflow_dist)

        if step_count == 1 or step_count % max(int(print_every), 1) == 0:
            print(
                f"t = {current_t:8.2f} s | "
                f"Qeva = {row['Q_eva_in']:10.3f} W | "
                f"Qcon = {row['Q_con_rejection']:10.3f} W | "
                f"Tmin = {row['T_hp_min']:8.3f} K | "
                f"Tmax = {row['T_hp_max']:8.3f} K | "
                f"k_wick_max = {row['k_wick_max']:10.3f} W/mK"
            )

    elapsed = time.time() - start_cpu
    final = rows[-1]

    print("-" * 78)
    print(f"Completed {step_count} steps in {elapsed:.2f} s")
    print(f"Final Q_eva_in       = {final['Q_eva_in']:.6f} W")
    print(f"Final Q_aba_loss     = {final['Q_aba_loss']:.6f} W")
    print(f"Final Q_con_reject   = {final['Q_con_rejection']:.6f} W")
    print(f"Final T_hp_min       = {final['T_hp_min']:.6f} K")
    print(f"Final T_hp_mean      = {final['T_hp_mean']:.6f} K")
    print(f"Final T_hp_max       = {final['T_hp_max']:.6f} K")
    print(f"Final k_wick_mean    = {final['k_wick_mean']:.6f} W/mK")
    print(f"Final k_wick_max     = {final['k_wick_max']:.6f} W/mK")

    if csv_path:
        write_csv(csv_path, rows)
        print(f"CSV written to: {csv_path}")

    if npz_path:
        write_npz(
            npz_path,
            {
                "time": np.asarray(time_history, dtype=float),
                "temperature_history": np.asarray(temperature_history, dtype=float),
                "outer_wall_temperature_history": np.asarray(outer_wall_temperature_history, dtype=float),
                "wick_k_history": np.asarray(wick_k_history, dtype=float),
                "q_eva_input_distribution_history": np.asarray(q_eva_input_distribution_history, dtype=float),
                "q_aba_loss_distribution_history": np.asarray(q_aba_loss_distribution_history, dtype=float),
                "q_con_rejection_distribution_history": np.asarray(q_con_rejection_distribution_history, dtype=float),
                "q_con_boundary_inflow_distribution_history": np.asarray(q_con_boundary_inflow_distribution_history, dtype=float),
                "x_centers": np.asarray(hp.mesh.x_centers, dtype=float),
                "y_centers": np.asarray(hp.mesh.y_centers, dtype=float),
                "initial_temp": np.array([initial_temp], dtype=float),
                "input_power": np.array([input_power], dtype=float),
                "condenser_h": np.array([condenser_h], dtype=float),
                "condenser_temp": np.array([condenser_temp], dtype=float),
                "dt": np.array([dt], dtype=float),
                "record_every": np.array([record_every], dtype=int),
                "axial_refine_factor": np.array([axial_refine_factor], dtype=float),
                "n_eva": np.array([hp.n_eva], dtype=int),
                "n_aba": np.array([hp.n_aba], dtype=int),
                "n_con": np.array([hp.n_con], dtype=int),
            }
        )
        print(f"NPZ written to: {npz_path}")

    return rows


def parse_args():
    parser = argparse.ArgumentParser(
        description="Single HP startup with evaporator power and condenser convection."
    )
    parser.add_argument("--t-end", type=float, default=3600.0)
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--initial-temp", type=float, default=293.15)
    parser.add_argument("--input-power", type=float, default=1800.0)
    parser.add_argument("--condenser-h", type=float, default=147.67)
    parser.add_argument("--condenser-temp", type=float, default=293.15)
    parser.add_argument("--axial-refine-factor", type=float, default=1.0)
    parser.add_argument("--n-eva", type=int, default=None)
    parser.add_argument("--n-aba", type=int, default=None)
    parser.add_argument("--n-con", type=int, default=None)
    parser.add_argument("--print-every", type=int, default=100)
    parser.add_argument("--record-every", type=int, default=1)
    parser.add_argument(
        "--csv-path",
        default=os.path.join(CURRENT_DIR, "single_hp_power_convection_startup_3600s.csv"),
    )
    parser.add_argument(
        "--npz-path",
        default=os.path.join(CURRENT_DIR, "single_hp_power_convection_startup_3600s.npz"),
    )
    parser.add_argument("--no-csv", action="store_true")
    parser.add_argument("--no-npz", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_case(
        t_end=args.t_end,
        dt=args.dt,
        initial_temp=args.initial_temp,
        input_power=args.input_power,
        condenser_h=args.condenser_h,
        condenser_temp=args.condenser_temp,
        axial_refine_factor=args.axial_refine_factor,
        n_eva=args.n_eva,
        n_aba=args.n_aba,
        n_con=args.n_con,
        print_every=args.print_every,
        record_every=args.record_every,
        csv_path=None if args.no_csv else args.csv_path,
        npz_path=None if args.no_npz else args.npz_path,
    )
