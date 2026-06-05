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


DEFAULT_RESTART_CANDIDATES = [
    "collector_ring_full_ringhp_buffered_half_ringflow_500s_resume_from200s_restart_t0500s.npz",
    "collector_ring_full_ringhp_buffered_half_ringflow_500s_resume_from200s_restart.npz",
    "collector_ring_full_ringhp_buffered_half_ringflow_200s_resume_from60s_restart_t0200s.npz",
    "collector_ring_full_ringhp_buffered_half_ringflow_200s_resume_from60s_restart.npz",
]


def existing_default_restart():
    for name in DEFAULT_RESTART_CANDIDATES:
        path = os.path.join(current_dir, name)
        if os.path.exists(path):
            return path
    return os.path.join(current_dir, DEFAULT_RESTART_CANDIDATES[0])


def restart_time(restart_path):
    with np.load(restart_path, allow_pickle=True) as data:
        if "System/global_time" not in data:
            raise KeyError(f"Restart has no System/global_time: {restart_path}")
        return float(data["System/global_time"][0])


def summarize_steady(history, window, tolerance):
    if not history:
        raise ValueError("No history rows were produced.")

    end_time = float(history[-1]["time"])
    start_time = max(float(history[0]["time"]), end_time - window)
    tail = [row for row in history if float(row["time"]) >= start_time]
    if len(tail) < 2:
        tail = history

    t0 = float(tail[0]["time"])
    t1 = float(tail[-1]["time"])
    temp0 = float(tail[0]["T_out_avg"])
    temp1 = float(tail[-1]["T_out_avg"])
    dt = max(t1 - t0, 1.0e-30)
    slope = (temp1 - temp0) / dt
    flow_in = float(history[-1]["W_in_total"])
    flow_out = float(history[-1]["W_out_total"])
    ring_in = float(history[-1]["W_ring_in_total"])
    ring_out = float(history[-1]["W_ring_out_total"])

    print("-" * 78)
    print("Steady Debug Summary")
    print(f"  time range used       : {t0:.6f} s -> {t1:.6f} s")
    print(f"  T_out_avg final       : {temp1:.9f} K")
    print(f"  dT_out_avg/dt         : {slope:.9e} K/s")
    print(f"  abs slope tolerance   : {tolerance:.9e} K/s")
    print(f"  steady_by_T_out_avg   : {abs(slope) <= tolerance}")
    print(f"  W_in_total final      : {flow_in:.9f} kg/s")
    print(f"  W_out_total final     : {flow_out:.9f} kg/s")
    print(f"  W_in - W_out final    : {flow_in - flow_out:.9e} kg/s")
    print(f"  W_ring_in final       : {ring_in:.9f} kg/s")
    print(f"  W_ring_out final      : {ring_out:.9f} kg/s")
    print(f"  W_ring_in - out final : {ring_in - ring_out:.9e} kg/s")
    print(f"  W_ring_closure final  : {float(history[-1]['W_ring_closure']):.9f} kg/s")
    print("-" * 78)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Resume the current full-ring collector-ring + RingHP model and "
            "print a short steady-state diagnostic."
        )
    )
    parser.add_argument(
        "--restart-from",
        default=existing_default_restart(),
        help="Restart .npz path. Defaults to the latest known full_ringhp checkpoint.",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=1.0,
        help="Seconds to continue from the restart time when --target-time is omitted.",
    )
    parser.add_argument(
        "--target-time",
        type=float,
        default=None,
        help="Absolute target time in seconds. Overrides --duration.",
    )
    parser.add_argument("--case-name", default="collector_ring_full_ringhp_steady_debug")
    parser.add_argument("--print-every", type=float, default=1.0)
    parser.add_argument("--min-dt", type=float, default=1.0e-3)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--safety-factor", type=float, default=1.0)
    parser.add_argument("--inner-iter", type=int, default=2)
    parser.add_argument("--steady-window", type=float, default=1.0)
    parser.add_argument("--steady-tol", type=float, default=1.0e-3)
    parser.add_argument("--restart-save-every", type=float, default=0.0)
    parser.add_argument("--csv-path", default=None)
    parser.add_argument("--restart-save-path", default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    restart_path = os.path.abspath(args.restart_from)
    if not os.path.exists(restart_path):
        raise FileNotFoundError(restart_path)

    start_time = restart_time(restart_path)
    if args.target_time is None:
        target_time = start_time + args.duration
    else:
        target_time = args.target_time
    if target_time <= start_time:
        raise ValueError(
            f"Target time must be greater than restart time: {target_time} <= {start_time}"
        )

    csv_path = args.csv_path
    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{args.case_name}_history.csv")

    restart_save_path = args.restart_save_path
    if restart_save_path is None:
        restart_save_path = os.path.join(current_dir, f"{args.case_name}_restart.npz")

    _, history = model.run_case(
        case_name=args.case_name,
        t_end=target_time,
        min_dt=args.min_dt,
        max_dt=args.max_dt,
        safety_factor=args.safety_factor,
        inner_iter=args.inner_iter,
        print_every_time=args.print_every,
        csv_path=csv_path,
        restart_from=restart_path,
        restart_save_path=restart_save_path,
        restart_save_every=args.restart_save_every,
    )
    summarize_steady(history, args.steady_window, args.steady_tol)


if __name__ == "__main__":
    main()
