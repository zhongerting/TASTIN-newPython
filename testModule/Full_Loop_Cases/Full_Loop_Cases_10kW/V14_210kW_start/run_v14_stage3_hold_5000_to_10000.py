"""Continue the accepted V14 Stage 3 fixed-current hold from 5000 to 10000 s."""

from pathlib import Path

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_start import (
    run_v14_startup_stages as startup,
)


CASE_DIR = Path(__file__).resolve().parent
RESTART_IN = (
    CASE_DIR
    / "startup_5000s_fixed_r_switchrestart_20260806"
    / "stage_3_cesium_tec_to_5000s"
    / "final_restart.npz"
)
OUTPUT_DIR = (
    CASE_DIR
    / "startup_10000s_fixed_i_continuation_20260806"
    / "stage_3_fixed_i_hold_5000_to_10000s"
)


def main() -> int:
    startup.STAGES[3] = (5000.0, 10000.0, OUTPUT_DIR.name)
    startup.run_stage(
        3,
        restart_in=RESTART_IN,
        output_dir=OUTPUT_DIR,
        max_dt_s=0.05,
        record_interval_s=1.0,
        checkpoint_interval_s=50.0,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
