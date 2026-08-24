"""Run the V14 accident with one half of the radiator sectors failed."""

from pathlib import Path
import sys

REPO_ROOT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "testModule").is_dir()
)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.heatpipe_failure_accident import (  # noqa: E402
    run_cli,
)


CASE_DIR = Path(__file__).resolve().parent
CASE_NAME = "V14_210kW_heatpipe_half_radiator_failure"
FAILURE_MODE = "upper_lower_A1_A3_half_radiator_transfer_0pct"

# A1-A3 are one contiguous half of the six-sector radiator ring.  Applying
# the same map to both physical rings disables 6 of 12 sectors and 171 of the
# nominal 340 heat pipes, which is the closest exact geometric half in this
# topology while preserving all nominal hydraulic paths.
FAILURE_MAP = {
    "upper": {
        0: {0: 0.0, 1: 0.0, 2: 0.0},
        1: {0: 0.0, 1: 0.0, 2: 0.0},
        2: {0: 0.0, 1: 0.0, 2: 0.0},
    },
    "lower": {
        0: {0: 0.0, 1: 0.0, 2: 0.0},
        1: {0: 0.0, 1: 0.0, 2: 0.0},
        2: {0: 0.0, 1: 0.0, 2: 0.0},
    },
}
DEFAULT_OUTPUT_DIR = CASE_DIR / "runs" / "from_t019865s"


if __name__ == "__main__":
    raise SystemExit(run_cli(
        case_name=CASE_NAME,
        failure_mode=FAILURE_MODE,
        failure_map=FAILURE_MAP,
        default_output_dir=DEFAULT_OUTPUT_DIR,
    ))
