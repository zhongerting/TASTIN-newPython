import os
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_full_ringhp as model


def main():
    case_name = "collector_ring_full_ringhp_buffered_half_ringflow_500s_resume_from200s"
    model.run_case(
        case_name=case_name,
        t_end=500.0,
        print_every_time=1.0,
        csv_path=os.path.join(
            current_dir,
            f"{case_name}_history.csv",
        ),
        restart_from=os.path.join(
            current_dir,
            "collector_ring_full_ringhp_buffered_half_ringflow_200s_resume_from60s_restart_t0200s.npz",
        ),
        restart_save_path=os.path.join(
            current_dir,
            f"{case_name}_restart.npz",
        ),
        restart_save_every=100.0,
    )


if __name__ == "__main__":
    main()
