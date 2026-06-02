import argparse

import audit_v7_caseA_fluid_energy_network as audit
from run_v8_caseA_common import build_loaded_case, parse_v8_multipliers


DEFAULT_RESTART = "testModule/v8_caseA_long_energy/v8_caseA_long_energy_latest_restart.npz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit V8 CaseA fluid matrix energy closure.")
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--previous-restart-in", default=None)
    parser.add_argument("--output-dir", default="testModule/v8_caseA_fluid_energy_audit")
    parser.add_argument("--case-prefix", default="v8_caseA_fluid_energy_audit")
    parser.add_argument("--target-voltage", type=float, default=27.2)
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
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--dt-probe", type=float, default=None)
    parser.add_argument("--dt-probe-sweep", default="1e-4,1e-2,0.8,10.0")
    parser.add_argument("--run-regression-tests", action="store_true")
    return parser.parse_args()


def main() -> None:
    audit.build_loaded_case = build_loaded_case
    audit.parse_args = parse_args
    audit.main()


if __name__ == "__main__":
    main()
