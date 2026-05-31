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

from Solvers.SystemManager import SystemManager
from testModule.test_single_hp_fin_energy_conservation import (
    build_hp_radiator,
    clear_boundary_conditions,
    create_dummy_fluid_network,
)


def prepare_cold_start_hp(initial_temp: float, eva_hot_side_temp: float, eva_r_ext: float = 1.0e-8):
    """
    低温启动单热管装配。
    说明：
    - 热管本体和辐射边界完全沿用现有单热管测试算例。
    - 蒸发段热侧仍采用当前测试算例的“定温热端 + 极小外热阻”口径，
      用 864 K 热侧代表高温流体换热环境。
    """
    hp_radiator = build_hp_radiator()

    eva_boundary = hp_radiator.hp.boundaries["outer_eva"]
    clear_boundary_conditions(eva_boundary)
    eva_boundary.add_resistance_condition(T_ext=eva_hot_side_temp, R_ext=eva_r_ext)

    hp_radiator.hp.T.fill(initial_temp)
    hp_radiator.hp.current_time = 0.0
    hp_radiator.last_fin_temperature.fill(initial_temp)
    hp_radiator.last_fin_radiation_distribution.fill(0.0)
    hp_radiator.last_fin_absorption_distribution.fill(0.0)
    hp_radiator.last_fin_net_from_root_distribution.fill(0.0)
    hp_radiator.last_fin_conductance_distribution.fill(0.0)
    hp_radiator.last_fin_effective_temperature_distribution.fill(hp_radiator.T_space)
    hp_radiator.last_fin_equivalent_resistance_distribution.fill(1.0e15)

    hp_radiator.hp.initialize_state()
    return hp_radiator


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
             initial_temp: float = 300.0,
             eva_hot_side_temp: float = 864.0,
             print_every: int = 20,
             csv_path: str = None,
             npz_path: str = None,
             record_every: int = 1):
    hp_radiator = prepare_cold_start_hp(
        initial_temp=initial_temp,
        eva_hot_side_temp=eva_hot_side_temp,
    )
    dummy_net = create_dummy_fluid_network()

    sys_manager = SystemManager(fluid_network=dummy_net, start_time=0.0)
    sys_manager.add_component(hp_radiator)
    sys_manager.initialize_system(dt_init=dt)

    rows = []
    time_history = []
    temperature_history = []
    outer_wall_temperature_history = []
    wick_k_history = []
    q_eva_in_distribution_history = []
    q_aba_rad_distribution_history = []
    q_con_gross_distribution_history = []
    q_con_bare_distribution_history = []
    q_con_fin_rad_distribution_history = []
    q_con_fin_net_distribution_history = []
    q_con_fin_abs_distribution_history = []
    step_count = 0
    start_cpu = time.time()

    print("=" * 72)
    print("Single HP cold-start test")
    print(f"initial_temp = {initial_temp:.3f} K")
    print(f"evaporator hot-side temperature = {eva_hot_side_temp:.3f} K")
    print(f"t_end = {t_end:.3f} s, dt = {dt:.3f} s")
    print("=" * 72)

    while sys_manager.global_time + 0.5 * dt <= t_end:
        sys_manager.step(dt=dt, inner_iter=1)

        current_t = sys_manager.global_time
        q_aba_dist, q_con_total_dist = hp_radiator.get_heat_rejection_distribution()
        breakdown = hp_radiator.get_heat_exchange_breakdown()
        T_2d = hp_radiator.get_temperature_distribution()
        k_2d = hp_radiator.hp.k_node.reshape(hp_radiator.hp.shape_nodes)
        wick_k = k_2d[:hp_radiator.hp.n_wick, :]

        Q_eva_in = float(np.sum(hp_radiator.hp.boundaries["outer_eva"].current_flux))
        Q_aba_rad = float(np.sum(q_aba_dist))
        Q_con_gross = float(np.sum(q_con_total_dist))
        Q_con_bare = float(np.sum(breakdown["bare_radiation"]))
        Q_con_fin_rad = float(np.sum(breakdown["fin_radiation"]))
        Q_con_fin_net = float(np.sum(breakdown["fin_net_from_root"]))
        Q_con_net = float(np.sum(breakdown["net_rejection"]))

        row = {
            "step": step_count + 1,
            "t": float(current_t),
            "Q_eva_in": Q_eva_in,
            "Q_aba_rad": Q_aba_rad,
            "Q_con_gross": Q_con_gross,
            "Q_con_bare": Q_con_bare,
            "Q_con_fin_rad": Q_con_fin_rad,
            "Q_con_fin_net": Q_con_fin_net,
            "Q_con_net": Q_con_net,
            "T_hp_min": float(np.min(hp_radiator.hp.T)),
            "T_hp_mean": float(np.mean(hp_radiator.hp.T)),
            "T_hp_max": float(np.max(hp_radiator.hp.T)),
            "T_eva_wall_mean": float(np.mean(T_2d[-1, hp_radiator.hp._slice_eva])),
            "T_con_wall_mean": float(np.mean(T_2d[-1, hp_radiator.hp._slice_con])),
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
            q_eva_in_distribution_history.append(
                np.array(hp_radiator.hp.boundaries["outer_eva"].current_flux, copy=True)
            )
            q_aba_rad_distribution_history.append(np.array(q_aba_dist, copy=True))
            q_con_gross_distribution_history.append(np.array(q_con_total_dist, copy=True))
            q_con_bare_distribution_history.append(np.array(breakdown["bare_radiation"], copy=True))
            q_con_fin_rad_distribution_history.append(np.array(breakdown["fin_radiation"], copy=True))
            q_con_fin_net_distribution_history.append(np.array(breakdown["fin_net_from_root"], copy=True))
            q_con_fin_abs_distribution_history.append(np.array(breakdown["fin_absorption"], copy=True))

        if step_count == 1 or step_count % max(print_every, 1) == 0:
            print(
                f"t = {current_t:7.3f} s | "
                f"Qeva = {Q_eva_in:10.3f} W | "
                f"Qcon = {Q_con_gross:10.3f} W | "
                f"Tmin = {row['T_hp_min']:8.3f} K | "
                f"Tmax = {row['T_hp_max']:8.3f} K | "
                f"k_wick_max = {row['k_wick_max']:10.3f} W/mK"
            )

    elapsed = time.time() - start_cpu
    final = rows[-1]

    print("-" * 72)
    print(f"Completed {step_count} steps in {elapsed:.2f} s")
    print(f"Final Q_eva_in      = {final['Q_eva_in']:.6f} W")
    print(f"Final Q_aba_rad     = {final['Q_aba_rad']:.6f} W")
    print(f"Final Q_con_gross   = {final['Q_con_gross']:.6f} W")
    print(f"Final T_hp_min      = {final['T_hp_min']:.6f} K")
    print(f"Final T_hp_mean     = {final['T_hp_mean']:.6f} K")
    print(f"Final T_hp_max      = {final['T_hp_max']:.6f} K")
    print(f"Final k_wick_mean   = {final['k_wick_mean']:.6f} W/mK")
    print(f"Final k_wick_max    = {final['k_wick_max']:.6f} W/mK")

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
                "q_eva_in_distribution_history": np.asarray(q_eva_in_distribution_history, dtype=float),
                "q_aba_rad_distribution_history": np.asarray(q_aba_rad_distribution_history, dtype=float),
                "q_con_gross_distribution_history": np.asarray(q_con_gross_distribution_history, dtype=float),
                "q_con_bare_distribution_history": np.asarray(q_con_bare_distribution_history, dtype=float),
                "q_con_fin_rad_distribution_history": np.asarray(q_con_fin_rad_distribution_history, dtype=float),
                "q_con_fin_net_distribution_history": np.asarray(q_con_fin_net_distribution_history, dtype=float),
                "q_con_fin_abs_distribution_history": np.asarray(q_con_fin_abs_distribution_history, dtype=float),
                "x_centers": np.asarray(hp_radiator.hp.mesh.x_centers, dtype=float),
                "y_centers": np.asarray(hp_radiator.hp.mesh.y_centers, dtype=float),
                "initial_temp": np.array([initial_temp], dtype=float),
                "eva_hot_side_temp": np.array([eva_hot_side_temp], dtype=float),
                "dt": np.array([dt], dtype=float),
                "record_every": np.array([record_every], dtype=int),
            }
        )
        print(f"NPZ written to: {npz_path}")

    return rows


def parse_args():
    parser = argparse.ArgumentParser(description="Single heat pipe cold-start test.")
    parser.add_argument("--t-end", type=float, default=20.0)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--initial-temp", type=float, default=300.0)
    parser.add_argument("--eva-hot-side-temp", type=float, default=864.0)
    parser.add_argument("--print-every", type=int, default=20)
    parser.add_argument(
        "--csv-path",
        default=os.path.join(CURRENT_DIR, "single_hp_cold_start_history.csv"),
    )
    parser.add_argument(
        "--npz-path",
        default=os.path.join(CURRENT_DIR, "single_hp_cold_start_history.npz"),
    )
    parser.add_argument("--record-every", type=int, default=1)
    parser.add_argument("--no-csv", action="store_true")
    parser.add_argument("--no-npz", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_case(
        t_end=args.t_end,
        dt=args.dt,
        initial_temp=args.initial_temp,
        eva_hot_side_temp=args.eva_hot_side_temp,
        print_every=args.print_every,
        csv_path=None if args.no_csv else args.csv_path,
        npz_path=None if args.no_npz else args.npz_path,
        record_every=args.record_every,
    )
