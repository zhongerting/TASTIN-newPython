import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from test_core_assemble_v7_caseA import run_test_v7_case_a_heated


if __name__ == "__main__":
    result = run_test_v7_case_a_heated(
        target_time_s=6000.0,
        total_power_w=115000.0,
        max_dt=1.0,
        safety_factor=20.0,
        restart_file="test_core_assemble_v7_caseA_heated_restart_t4200.npz",
        save_interval=20.0,
        final_restart_file="test_core_assemble_v7_caseA_heated_restart_t6000.npz",
        keep_only_latest_restart=True,
    )
    print("FINAL_SUMMARY")
    print(result["final_summary"])
