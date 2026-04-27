import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from model_collector_ring_6segment import run_case


if __name__ == "__main__":
    case_name = "collector_ring_6segment_buffered_160s"
    run_case(
        case_name=case_name,
        t_end=160.0,
        print_every_time=1.0,
        csv_path=os.path.join(current_dir, f"{case_name}_history.csv"),
        restart_save_path=os.path.join(current_dir, f"{case_name}_restart.npz"),
        restart_save_every=40.0,
        restart_from=None,
        profiler_summary_path=os.path.join(
            current_dir, f"{case_name}_profiler_summary.csv"
        ),
        profiler_snapshot_path=os.path.join(
            current_dir, f"{case_name}_profiler_snapshots.csv"
        ),
        profiler_report_path=os.path.join(
            current_dir, f"{case_name}_profiler_report.txt"
        ),
    )
