import os
import sys
import time

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from profiler import TEASAProfiler
from Components.HPwithFin import HPwithFin
from Components.basicComponents.HeatPipe2D import HeatPipe2D
from Components.RingHP import RingHP, SingleVolumeProxy


def env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    return float(raw)


def env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    return int(raw)


def env_str(name: str, default: str) -> str:
    raw = os.getenv(name)
    if raw is None or raw == "":
        return default
    return raw


def install_profile_hook(cls, method_name: str) -> None:
    method = getattr(cls, method_name)
    if getattr(method, "_teasa_profile_installed", False):
        return

    base_method = getattr(method, "__wrapped__", method)
    wrapped = TEASAProfiler.profile(base_method)
    wrapped._teasa_profile_installed = True
    setattr(cls, method_name, wrapped)


def install_profile_hooks() -> None:
    targets = [
        (RingHP, "pre_step"),
        (HPwithFin, "pre_step"),
        (HPwithFin, "_solve_fin_quasi_steady"),
        (HeatPipe2D, "_update_properties"),
        (HeatPipe2D, "_update_boundaries_state"),
        (HeatPipe2D, "_compute_fluxes"),
        (SingleVolumeProxy, "add_coupling_source_distribution"),
    ]
    for cls, method_name in targets:
        install_profile_hook(cls, method_name)


def print_focus_stats(limit: int = 12) -> None:
    prefixes = (
        "RingHP.",
        "HPwithFin.",
        "HeatPipe2D.",
        "SingleVolumeProxy.",
    )
    stats = [
        (name, data)
        for name, data in TEASAProfiler.stats.items()
        if name.startswith(prefixes)
    ]
    stats.sort(key=lambda item: item[1]["time"], reverse=True)

    print("\n" + "=" * 78)
    print(f"{'HeatPipe Function':<48} | {'Calls':<8} | {'Total Time (s)':<14}")
    print("-" * 78)
    for name, data in stats[:limit]:
        print(f"{name:<48} | {data['count']:<8} | {data['time']:<14.6f}")
    print("=" * 78)


def run_benchmark():
    install_profile_hooks()

    import CoolantLoop.test_coolant_loop_v5 as ring_case

    dt = env_float("BENCH_DT", 0.05)
    n_steps = env_int("BENCH_N_STEPS", 36000)
    inner_iter = env_int("BENCH_INNER_ITER", 1)
    print_every = env_int("BENCH_PRINT_EVERY", 50)
    top_n = env_int("BENCH_TOP", 12)
    case_name = env_str("BENCH_CASE_NAME", "ringhp_v5_baseline")

    csv_path = os.path.join(CURRENT_DIR, f"{case_name}_history.csv")
    total_time = dt * n_steps

    print("=" * 78)
    print("HeatPipe benchmark based on CoolantLoop/test_coolant_loop_v5.py")
    print(f"dt = {dt:.6f} s, n_steps = {n_steps}, inner_iter = {inner_iter}")
    print(f"target end time = {total_time:.6f} s")
    print("=" * 78)

    TEASAProfiler.stats.clear()

    wall_start = time.perf_counter()
    history = ring_case.run_case(
        case_name=case_name,
        t_end=total_time,
        min_dt=dt,
        max_dt=dt,
        safety_factor=2.0,
        inner_iter=inner_iter,
        print_every=print_every,
        csv_path=csv_path,
        restart_from=None,
        restart_save_path=None,
        restart_save_every=0.0,
    )
    wall_time = time.perf_counter() - wall_start

    print(f"Benchmark wall time: {wall_time:.6f} s")
    print(f"CSV saved to: {csv_path}")
    print_focus_stats(limit=top_n)
    TEASAProfiler.report()

    return {
        "history": history,
        "wall_time": wall_time,
        "csv_path": csv_path,
        "stats": dict(TEASAProfiler.stats),
    }


if __name__ == "__main__":
    run_benchmark()
