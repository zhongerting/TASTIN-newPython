import os
import sys
import time

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import CoolantLoop.test_coolant_loop_v5 as ring_case


def main():
    t_end = float(os.environ.get("LONG_CASE_T_END", "1800.0"))
    dt = float(os.environ.get("LONG_CASE_DT", "0.25"))
    inner_iter = int(os.environ.get("LONG_CASE_INNER_ITER", "1"))
    print_every = int(os.environ.get("LONG_CASE_PRINT_EVERY", "400"))
    case_name = os.environ.get("LONG_CASE_NAME", "ringhp_v5_fin_tangent_1800s_dt025")
    csv_path = os.environ.get(
        "LONG_CASE_CSV",
        os.path.join("HeatPipe优化重构", f"{case_name}_history.csv"),
    )

    start = time.perf_counter()
    history = ring_case.run_case(
        case_name=case_name,
        t_end=t_end,
        min_dt=dt,
        max_dt=dt,
        safety_factor=1.0,
        inner_iter=inner_iter,
        print_every=print_every,
        csv_path=csv_path,
        restart_from=None,
        restart_save_path=None,
        restart_save_every=0.0,
    )
    elapsed = time.perf_counter() - start
    print(f"LONG_CASE_ROWS={len(history)}")
    print(f"LONG_CASE_WALL_TIME_S={elapsed:.3f}")


if __name__ == "__main__":
    main()
