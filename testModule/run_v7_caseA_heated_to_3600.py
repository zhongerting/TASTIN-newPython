import os
import re
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from test_core_assemble_v7_caseA import run_test_v7_case_a_heated


def _latest_restart():
    pattern = re.compile(r"test_core_assemble_v7_caseA_heated_restart_t(\d+)\.npz$")
    candidates = []
    for name in os.listdir(root_dir):
        match = pattern.match(name)
        if match:
            candidates.append((int(match.group(1)), name))
    if not candidates:
        return None
    candidates.sort()
    return candidates[-1][1]


if __name__ == "__main__":
    restart_file = _latest_restart()
    result = run_test_v7_case_a_heated(
        target_time_s=3600.0,
        total_power_w=115000.0,
        max_dt=1.0,
        safety_factor=20.0,
        restart_file=restart_file,
        save_interval=300.0,
        final_restart_file="test_core_assemble_v7_caseA_heated_restart_t3600.npz",
        keep_only_latest_restart=True,
    )
    print("BACKGROUND_FINAL_SUMMARY")
    print(result["final_summary"])
