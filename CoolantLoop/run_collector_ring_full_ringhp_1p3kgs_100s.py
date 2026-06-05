import argparse
import os
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_full_ringhp as model


DEFAULT_CASE_NAME = "collector_ring_full_ringhp_buffered_1p3kgs_100s"
DEFAULT_TOTAL_FLOW = 1.3


def run_with_total_flow(args):
    old_w_total = model.W_TOTAL
    old_w_branch_total = model.W_BRANCH_TOTAL
    old_w_inlet_leg_init = model.W_INLET_LEG_INIT
    try:
        model.W_TOTAL = float(args.total_flow)
        model.W_BRANCH_TOTAL = model.W_TOTAL / 3.0
        model.W_INLET_LEG_INIT = model.W_BRANCH_TOTAL
        print("=" * 78)
        print("Full-ring collector-ring + RingHP flow-variant case")
        print(f"  total inlet flow target : {model.W_TOTAL:.9f} kg/s")
        print(f"  branch nominal flow     : {model.W_BRANCH_TOTAL:.9f} kg/s")
        print(f"  target time             : {args.t_end:.6f} s")
        print("=" * 78)

        return model.run_case(
            case_name=args.case_name,
            t_end=args.t_end,
            min_dt=args.min_dt,
            max_dt=args.max_dt,
            safety_factor=args.safety_factor,
            inner_iter=args.inner_iter,
            print_every_time=args.print_every,
            csv_path=args.csv_path,
            restart_from=args.restart_from,
            restart_save_path=args.restart_save_path,
            restart_save_every=args.restart_save_every,
        )
    finally:
        model.W_TOTAL = old_w_total
        model.W_BRANCH_TOTAL = old_w_branch_total
        model.W_INLET_LEG_INIT = old_w_inlet_leg_init


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run a separate full-ring collector-ring + RingHP case at 1.3 kg/s."
    )
    parser.add_argument("--total-flow", type=float, default=DEFAULT_TOTAL_FLOW)
    parser.add_argument("--t-end", type=float, default=100.0)
    parser.add_argument("--case-name", default=DEFAULT_CASE_NAME)
    parser.add_argument("--min-dt", type=float, default=2.0e-4)
    parser.add_argument("--max-dt", type=float, default=5.0e-3)
    parser.add_argument("--safety-factor", type=float, default=0.1)
    parser.add_argument("--inner-iter", type=int, default=2)
    parser.add_argument("--print-every", type=float, default=10.0)
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
    run_with_total_flow(parse_args())


if __name__ == "__main__":
    main()
