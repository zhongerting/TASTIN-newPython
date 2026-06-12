import argparse

import audit_v7_caseA_global_storage_continue as audit
from run_v8_caseA_common import (
    DEFAULT_SOLID_ODE_METHOD,
    build_loaded_case,
    parse_solid_ode_method,
    parse_v8_multipliers,
)


DEFAULT_RESTART = "testModule/v8_caseA_migrated/v8_caseA_migrated_latest_restart.npz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Continue V8 CaseA and audit radial storage energy.")
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", default="testModule/v8_caseA_global_storage_continue")
    parser.add_argument("--case-prefix", default="v8_caseA_global_storage_continue")
    parser.add_argument("--duration", type=float, default=100.0)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--inner-iter", type=int, default=1)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument(
        "--solid-ode-method",
        type=parse_solid_ode_method,
        default=DEFAULT_SOLID_ODE_METHOD,
    )
    parser.add_argument(
        "--ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=False),
        default=parse_v8_multipliers("1,6,12,15,3"),
    )
    parser.add_argument(
        "--tec-ring-multipliers",
        type=lambda text: parse_v8_multipliers(text, allow_zero=True),
        default=parse_v8_multipliers("1,6,12,15,0", allow_zero=True),
    )
    return parser.parse_args()


def main() -> None:
    audit.build_loaded_case = build_loaded_case
    audit.parse_args = parse_args
    audit.main()


if __name__ == "__main__":
    main()
