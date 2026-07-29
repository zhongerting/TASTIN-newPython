"""Run the fixed-power-until-trip TEC open-circuit accident."""

from pathlib import Path
import sys

REPO_ROOT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "testModule").is_dir()
)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.tec_open_circuit_accident import run_cli


CASE_NAME = "V14_210kW_TEC_open_circuit_accident_fixed_power"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "runs" / "final_from_t019865s"


if __name__ == "__main__":
    raise SystemExit(run_cli("fixed_power", CASE_NAME, DEFAULT_OUTPUT_DIR))
