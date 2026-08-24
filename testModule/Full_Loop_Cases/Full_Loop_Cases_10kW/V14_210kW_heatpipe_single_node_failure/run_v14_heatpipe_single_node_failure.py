"""Run the V14 A5 corresponding upper/lower single-node failure."""

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
CASE_NAME = "V14_210kW_heatpipe_single_node_failure"
FAILURE_MODE = "corresponding_upper_lower_node_transfer_0pct"
FAILURE_MAP = {"upper": {4: {2: 0.0}}, "lower": {4: {2: 0.0}}}
DEFAULT_OUTPUT_DIR = CASE_DIR / "runs" / "from_t019865s"


if __name__ == "__main__":
    raise SystemExit(run_cli(
        case_name=CASE_NAME,
        failure_mode=FAILURE_MODE,
        failure_map=FAILURE_MAP,
        default_output_dir=DEFAULT_OUTPUT_DIR,
    ))
