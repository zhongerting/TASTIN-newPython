import argparse
import os
import sys
from typing import Optional

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from test_core_assemble_v7_caseA import run_test_v7_case_a_heated
from test_core_assemble_v7_caseA_faststeady import (
    FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE,
    FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE_SCOPE,
    FASTSTEADY_TOTAL_POWER_W,
    compute_faststeady_energy_audit,
)


DEFAULT_RESTART_IN = "test_core_assemble_v7_caseA_faststeady_restart_t8800.npz"
DEFAULT_RESTART_OUT = "test_core_assemble_v7_caseA_faststeady_restart_t18800.npz"
DEFAULT_TARGET_TIME_S = 18800.0


def _optional_restart_path(path: str, no_final_restart: bool) -> Optional[str]:
    if no_final_restart:
        return None
    cleaned = path.strip()
    return cleaned or None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Continue V7 CaseA fast-steady transient from an existing restart."
    )
    parser.add_argument("--restart-in", default=DEFAULT_RESTART_IN)
    parser.add_argument("--restart-out", default=DEFAULT_RESTART_OUT)
    parser.add_argument("--target-time", type=float, default=DEFAULT_TARGET_TIME_S)
    parser.add_argument("--duration", type=float, default=None)
    parser.add_argument("--safety-factor", type=float, default=280.0)
    parser.add_argument("--max-dt", type=float, default=5.0)
    parser.add_argument(
        "--solid-heat-capacity-scale",
        type=float,
        default=FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE,
    )
    parser.add_argument(
        "--solid-heat-capacity-scale-scope",
        default=FASTSTEADY_SOLID_HEAT_CAPACITY_SCALE_SCOPE,
        choices=("all", "global_outer", "tfe_only"),
    )
    parser.add_argument("--global-outer-heat-capacity-scale", type=float, default=None)
    parser.add_argument("--save-interval", type=float, default=0.0)
    parser.add_argument("--keep-only-latest-restart", action="store_true")
    parser.add_argument("--no-final-restart", action="store_true")
    parser.add_argument("--disable-time-dt-cap", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    final_restart_file = _optional_restart_path(
        args.restart_out,
        args.no_final_restart,
    )
    result = run_test_v7_case_a_heated(
        run_duration_s=0.0 if args.duration is None else args.duration,
        target_time_s=args.target_time,
        total_power_w=FASTSTEADY_TOTAL_POWER_W,
        max_dt=args.max_dt,
        safety_factor=args.safety_factor,
        restart_file=args.restart_in,
        save_interval=args.save_interval,
        final_restart_file=final_restart_file,
        keep_only_latest_restart=args.keep_only_latest_restart,
        reset_design_flow_after_restart=True,
        solid_heat_capacity_scale=args.solid_heat_capacity_scale,
        solid_heat_capacity_scale_scope=args.solid_heat_capacity_scale_scope,
        global_outer_heat_capacity_scale=args.global_outer_heat_capacity_scale,
        use_time_dependent_dt_cap=not args.disable_time_dt_cap,
    )
    audit = compute_faststeady_energy_audit(result)
    print("FASTSTEADY_CONTINUE_ENERGY_AUDIT")
    print(audit)


if __name__ == "__main__":
    main()
