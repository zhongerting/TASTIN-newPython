import argparse
import os
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_full_ringhp as model


DEFAULT_CASE_NAME = "collector_ring_full_ringhp_1p3kgs_inlet_init_steady"


def collect_history_row(m, dt):
    inlet_junction = m["inlet_junction"]
    outlet_junction = m["outlet_junction"]
    hot_leg_to_ring = m["hot_leg_to_ring"]
    ring_to_manifold = m["ring_to_manifold"]
    manifolds = m["manifolds"]
    ring_channel = m["ring_channel"]
    ring_closure = m["ring_closure"]
    inlet_buffer_channel = m["inlet_buffer_channel"]
    outlet_buffer_channel = m["outlet_buffer_channel"]

    t_outlet_list = [channel.volumes[-1].T for channel in manifolds]
    row = {
        "time": float(m["sys_mgr"].global_time),
        "dt": float(dt),
        "W_in_total": float(inlet_junction.W),
        "W_out_total": float(outlet_junction.W),
        "W_ring_in_total": float(sum(j.W for j in hot_leg_to_ring)),
        "W_ring_out_total": float(sum(j.W for j in ring_to_manifold)),
        "W_ring_closure": float(ring_closure.W),
        "T_out_avg": float(np.mean(t_outlet_list)),
        "T_inlet_buffer_out": float(inlet_buffer_channel.volumes[-1].T),
        "T_outlet_buffer_out": float(outlet_buffer_channel.volumes[-1].T),
    }
    for idx, temp in enumerate(t_outlet_list, start=1):
        row[f"T_outlet_{idx}"] = float(temp)
    for idx in range(3):
        row[f"T_ring_inlet_node_{idx + 1}"] = float(
            ring_channel.volumes[model.INLET_NODE_INDICES[idx]].T
        )
        row[f"T_ring_outlet_node_{idx + 1}"] = float(
            ring_channel.volumes[model.OUTLET_NODE_INDICES[idx]].T
        )
    return row


def steady_slope(history, window):
    if len(history) < 2:
        return None
    end_time = history[-1]["time"]
    start_time = max(history[0]["time"], end_time - window)
    tail = [row for row in history if row["time"] >= start_time]
    if len(tail) < 2:
        return None
    dt = tail[-1]["time"] - tail[0]["time"]
    if dt <= 0.0:
        return None
    return (tail[-1]["T_out_avg"] - tail[0]["T_out_avg"]) / dt


def run_case(args):
    old_values = {
        "W_TOTAL": model.W_TOTAL,
        "W_BRANCH_TOTAL": model.W_BRANCH_TOTAL,
        "W_INLET_LEG_INIT": model.W_INLET_LEG_INIT,
        "T_INIT": model.T_INIT,
        "HP_INITIAL_TEMP": model.HP_INITIAL_TEMP,
    }
    try:
        init_temp = model.T_INLET if args.init_temp is None else float(args.init_temp)
        hp_init_temp = init_temp if args.hp_init_temp is None else float(args.hp_init_temp)
        model.W_TOTAL = float(args.total_flow)
        model.W_BRANCH_TOTAL = model.W_TOTAL / 3.0
        model.W_INLET_LEG_INIT = model.W_BRANCH_TOTAL
        model.T_INIT = init_temp
        model.HP_INITIAL_TEMP = hp_init_temp

        m = model.build_model()
        sys_mgr = m["sys_mgr"]
        model.print_pre_run_summary(m, args.case_name)
        if args.restart_from is not None:
            sys_mgr.load_global_state(args.restart_from)
            m["inlet_boundary"].set_boundary_state(P=model.P_OUTLET + 5000.0, T=model.T_INLET)
            model.sync_boundary_to_network(m["network"], m["inlet_boundary"])
            m["outlet_boundary"].set_boundary_state(P=model.P_OUTLET)
            model.sync_boundary_to_network(m["network"], m["outlet_boundary"])
            m["inlet_junction"].set_flow_rate(model.W_TOTAL)
            print(f"Restart loaded from: {args.restart_from}")
            print(f"Restart time: {sys_mgr.global_time:.6f} s")
        else:
            sys_mgr.initialize_system()
            print("System initialized from configured initial condition.")
        print(f"  T_INIT                  = {init_temp:.6f} K")
        print(f"  HP_INITIAL_TEMP         = {hp_init_temp:.6f} K")
        print(f"  W_TOTAL target           = {model.W_TOTAL:.6f} kg/s")
        print(f"  max time                 = {args.max_time:.6f} s")
        print(f"  steady window / tol      = {args.steady_window:.6f} s / {args.steady_tol:.6e} K/s")
        print(
            f"  startup dt control       = t<{args.startup_time:.6f}s, "
            f"max_dt={args.startup_max_dt:.6e}s, safety={args.startup_safety_factor:.3f}"
        )

        history = []
        next_print_time = model.next_event_time(sys_mgr.global_time, args.print_every)
        next_restart_save_time = model.next_event_time(sys_mgr.global_time, args.restart_save_every)
        converged = False

        while sys_mgr.global_time < args.max_time:
            if sys_mgr.global_time < args.startup_time:
                max_dt = min(args.max_dt, args.startup_max_dt)
                safety_factor = min(args.safety_factor, args.startup_safety_factor)
            else:
                max_dt = args.max_dt
                safety_factor = args.safety_factor

            dt = sys_mgr.compute_adaptive_dt(
                min_dt=args.min_dt,
                max_dt=max_dt,
                safety_factor=safety_factor,
            )
            if next_restart_save_time is not None:
                if sys_mgr.global_time < next_restart_save_time < sys_mgr.global_time + dt:
                    dt = next_restart_save_time - sys_mgr.global_time
            if next_print_time is not None:
                if sys_mgr.global_time < next_print_time < sys_mgr.global_time + dt:
                    dt = min(dt, next_print_time - sys_mgr.global_time)
            dt = min(dt, args.max_time - sys_mgr.global_time)

            sys_mgr.step(dt=dt, inner_iter=args.inner_iter)
            row = collect_history_row(m, dt)
            history.append(row)
            slope = steady_slope(history, args.steady_window)

            should_print = row["time"] >= args.max_time
            if next_print_time is not None and row["time"] >= next_print_time - 1.0e-12:
                should_print = True
            if should_print:
                slope_text = "n/a" if slope is None else f"{slope:.6e} K/s"
                print(
                    f"t = {row['time']:9.3f} s | "
                    f"T_out_avg = {row['T_out_avg']:.3f} K | "
                    f"W_in_total = {row['W_in_total']:.4f} kg/s | "
                    f"W_ring_in_total = {row['W_ring_in_total']:.4f} kg/s | "
                    f"dTdt_window = {slope_text}"
                )
                while next_print_time is not None and row["time"] >= next_print_time - 1.0e-12:
                    next_print_time += args.print_every

            if (
                args.restart_save_path is not None
                and args.restart_save_every > 0.0
                and next_restart_save_time is not None
                and row["time"] >= next_restart_save_time - 1.0e-12
            ):
                checkpoint_path = model.restart_checkpoint_path(
                    args.restart_save_path,
                    next_restart_save_time,
                )
                sys_mgr.save_global_state(checkpoint_path)
                print(f"Restart saved at t={row['time']:.3f} s: {checkpoint_path}")
                while next_restart_save_time is not None and row["time"] >= next_restart_save_time - 1.0e-12:
                    next_restart_save_time += args.restart_save_every

            if (
                slope is not None
                and row["time"] >= args.min_steady_time
                and abs(slope) <= args.steady_tol
            ):
                converged = True
                print(
                    f"Steady criterion met at t={row['time']:.6f} s: "
                    f"|dT_out_avg/dt|={abs(slope):.6e} K/s"
                )
                break

        model.write_history_csv(args.csv_path, history)
        if args.restart_save_path is not None:
            sys_mgr.save_global_state(args.restart_save_path)
            print(f"Final restart saved: {args.restart_save_path}")

        final_slope = steady_slope(history, args.steady_window)
        final = history[-1]
        print("=" * 70)
        print(f"Case completed: {args.case_name}")
        print(f"  converged      : {converged}")
        print(f"  final time     : {final['time']:.9f} s")
        print(f"  T_out_avg      : {final['T_out_avg']:.9f} K")
        print(f"  W_in_total     : {final['W_in_total']:.9f} kg/s")
        print(f"  W_ring_in_total: {final['W_ring_in_total']:.9f} kg/s")
        if final_slope is not None:
            print(f"  dT_out_avg/dt  : {final_slope:.9e} K/s")
        print("=" * 70)

        return m, history
    finally:
        for name, value in old_values.items():
            setattr(model, name, value)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run the 1.3 kg/s full-ring collector-ring + RingHP case from "
            "an all-inlet-temperature initial condition until steady or max time."
        )
    )
    parser.add_argument("--total-flow", type=float, default=1.3)
    parser.add_argument("--init-temp", type=float, default=None)
    parser.add_argument("--hp-init-temp", type=float, default=None)
    parser.add_argument("--max-time", type=float, default=100.0)
    parser.add_argument("--min-steady-time", type=float, default=10.0)
    parser.add_argument("--steady-window", type=float, default=10.0)
    parser.add_argument("--steady-tol", type=float, default=1.0e-3)
    parser.add_argument("--case-name", default=DEFAULT_CASE_NAME)
    parser.add_argument("--min-dt", type=float, default=2.0e-4)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--safety-factor", type=float, default=1.0)
    parser.add_argument("--startup-time", type=float, default=0.1)
    parser.add_argument("--startup-max-dt", type=float, default=5.0e-3)
    parser.add_argument("--startup-safety-factor", type=float, default=0.1)
    parser.add_argument("--inner-iter", type=int, default=2)
    parser.add_argument("--print-every", type=float, default=5.0)
    parser.add_argument("--restart-save-every", type=float, default=20.0)
    parser.add_argument("--restart-from", default=None)
    parser.add_argument(
        "--csv-path",
        default=os.path.join(current_dir, f"{DEFAULT_CASE_NAME}_history.csv"),
    )
    parser.add_argument(
        "--restart-save-path",
        default=os.path.join(current_dir, f"{DEFAULT_CASE_NAME}_restart.npz"),
    )
    return parser.parse_args()


def main():
    run_case(parse_args())


if __name__ == "__main__":
    main()
