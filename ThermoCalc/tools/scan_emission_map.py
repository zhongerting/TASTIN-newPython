from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "emission_scan_outputs"


FLOAT_FIELDS = (
    "J",
    "Vd",
    "delta_V",
    "phiE",
    "phiC",
    "obstructed_residual",
    "transition_residual",
    "saturation_residual",
)

INT_FIELDS = (
    "regime",
    "iteration_count",
    "obstructed_iterations",
    "transition_iterations",
    "saturation_iterations",
    "error_code",
)

BOOL_FIELDS = (
    "converged",
    "finite_flag",
    "done",
)


def import_te_solver() -> Any:
    candidates = []
    env_dir = os.environ.get("THERMOCALC_TE_SOLVER_DIR")
    if env_dir:
        candidates.append(Path(env_dir))
    candidates.extend((ROOT / "build_cp312" / "Release", ROOT / "build" / "Release", ROOT))

    for candidate in reversed(candidates):
        if candidate.exists():
            sys.path.insert(0, str(candidate))

    import te_solver  # type: ignore

    if not hasattr(te_solver, "calc_emission_point"):
        searched = ", ".join(str(path) for path in candidates)
        raise RuntimeError(
            "te_solver.calc_emission_point is missing. Build the updated "
            f"ThermoCalc extension first. Searched: {searched}"
        )
    return te_solver


def axis_from_args(start: float, stop: float, count: int) -> np.ndarray:
    if count <= 0:
        raise ValueError("Axis point count must be positive.")
    return np.linspace(start, stop, count, dtype=np.float64)


def cesium_pressure_from_tcs(tcs: float) -> float:
    return 2.45e8 / math.sqrt(tcs) * math.exp(-8910.0 / tcs)


def tcs_from_cesium_pressure(target_pressure: float) -> float:
    if target_pressure <= 0.0:
        raise ValueError("Cesium pressure must be positive.")
    lo = 300.0
    hi = 1200.0
    while cesium_pressure_from_tcs(lo) > target_pressure:
        lo *= 0.8
    while cesium_pressure_from_tcs(hi) < target_pressure:
        hi *= 1.2
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if cesium_pressure_from_tcs(mid) < target_pressure:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def pressure_axis_from_args(start: float, stop: float, count: int, spacing: str) -> np.ndarray:
    if count <= 0:
        raise ValueError("Pressure point count must be positive.")
    if start <= 0.0 or stop <= 0.0:
        raise ValueError("Pressure limits must be positive.")
    if spacing == "log":
        return np.geomspace(start, stop, count, dtype=np.float64)
    if spacing == "linear":
        return np.linspace(start, stop, count, dtype=np.float64)
    raise ValueError("Pressure spacing must be 'log' or 'linear'.")


def new_arrays(shape: tuple[int, ...]) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for name in FLOAT_FIELDS:
        arrays[name] = np.full(shape, np.nan, dtype=np.float64)
    for name in INT_FIELDS:
        arrays[name] = np.full(shape, -1, dtype=np.int32)
    for name in BOOL_FIELDS:
        arrays[name] = np.zeros(shape, dtype=bool)
    return arrays


def axes_match(data: dict[str, np.ndarray], axes: dict[str, np.ndarray]) -> bool:
    return all(np.array_equal(data[name], values) for name, values in axes.items())


def load_or_create(path: Path, axes: dict[str, np.ndarray], resume: bool) -> dict[str, np.ndarray]:
    shape = tuple(len(axes[name]) for name in ("TE_axis", "TC_axis", "Vo_axis", "Tcs_axis"))
    if resume and path.exists():
        with np.load(path, allow_pickle=False) as loaded:
            loaded_dict = {name: loaded[name] for name in loaded.files}
        if not axes_match(loaded_dict, axes):
            raise ValueError("Checkpoint axes do not match the requested scan grid.")
        return loaded_dict

    arrays = new_arrays(shape)
    arrays.update(axes)
    return arrays


def save_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    np.savez_compressed(tmp_path, **arrays)
    if tmp_path.exists():
        tmp_path.replace(path)
    else:
        tmp_npz = Path(str(tmp_path) + ".npz")
        tmp_npz.replace(path)


def summarize(arrays: dict[str, np.ndarray], elapsed_s: float | None = None) -> dict[str, Any]:
    done = arrays["done"]
    total = int(done.size)
    done_count = int(np.count_nonzero(done))
    finite_done = done & arrays["finite_flag"]
    converged_done = done & arrays["converged"]
    failed_done = done & ~(arrays["finite_flag"] & arrays["converged"])

    summary: dict[str, Any] = {
        "total_points": total,
        "done_points": done_count,
        "remaining_points": total - done_count,
        "finite_points": int(np.count_nonzero(finite_done)),
        "converged_points": int(np.count_nonzero(converged_done)),
        "failed_or_nonfinite_points": int(np.count_nonzero(failed_done)),
        "finite_rate_done": float(np.count_nonzero(finite_done) / done_count) if done_count else 0.0,
        "converged_rate_done": float(np.count_nonzero(converged_done) / done_count) if done_count else 0.0,
    }

    if elapsed_s is not None:
        summary["elapsed_s"] = float(elapsed_s)
        summary["points_per_s"] = float(done_count / elapsed_s) if elapsed_s > 0 else 0.0

    if done_count:
        iterations = arrays["iteration_count"][done]
        summary["iteration_min"] = int(np.min(iterations))
        summary["iteration_mean"] = float(np.mean(iterations))
        summary["iteration_max"] = int(np.max(iterations))

        regimes = arrays["regime"][done]
        for regime in (-1, 0, 1, 2):
            summary[f"regime_{regime}_count"] = int(np.count_nonzero(regimes == regime))

        j_values = arrays["J"][finite_done]
        if j_values.size:
            summary["J_min"] = float(np.nanmin(j_values))
            summary["J_mean"] = float(np.nanmean(j_values))
            summary["J_max"] = float(np.nanmax(j_values))

    return summary


def write_summary_files(output_dir: Path, summary: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "emission_scan_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    with (output_dir / "emission_scan_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(("key", "value"))
        for key, value in summary.items():
            writer.writerow((key, value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan local thermionic-emission phase space.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--checkpoint-interval", type=int, default=1000)
    parser.add_argument("--progress-interval", type=int, default=1000)
    parser.add_argument("--max-runtime-s", type=float, default=0.0)
    parser.add_argument("--d-gap", type=float, default=0.5)

    parser.add_argument("--te-min", type=float, default=800.0)
    parser.add_argument("--te-max", type=float, default=2400.0)
    parser.add_argument("--te-count", type=int, default=25)
    parser.add_argument("--tc-min", type=float, default=600.0)
    parser.add_argument("--tc-max", type=float, default=1000.0)
    parser.add_argument("--tc-count", type=int, default=15)
    parser.add_argument("--vo-min", type=float, default=0.0)
    parser.add_argument("--vo-max", type=float, default=2.5)
    parser.add_argument("--vo-count", type=int, default=25)
    parser.add_argument("--tcs-min", type=float, default=500.0)
    parser.add_argument("--tcs-max", type=float, default=800.0)
    parser.add_argument("--tcs-count", type=int, default=20)
    parser.add_argument("--pcs-min", type=float, default=None)
    parser.add_argument("--pcs-max", type=float, default=None)
    parser.add_argument("--pcs-count", type=int, default=None)
    parser.add_argument("--pcs-spacing", choices=("log", "linear"), default="log")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    te_solver = import_te_solver()

    axes = {
        "TE_axis": axis_from_args(args.te_min, args.te_max, args.te_count),
        "TC_axis": axis_from_args(args.tc_min, args.tc_max, args.tc_count),
        "Vo_axis": axis_from_args(args.vo_min, args.vo_max, args.vo_count),
    }
    if args.pcs_min is not None or args.pcs_max is not None or args.pcs_count is not None:
        if args.pcs_min is None or args.pcs_max is None:
            raise ValueError("--pcs-min and --pcs-max must be provided together.")
        pcs_count = args.pcs_count if args.pcs_count is not None else args.tcs_count
        pcs_axis = pressure_axis_from_args(args.pcs_min, args.pcs_max, pcs_count, args.pcs_spacing)
        axes["Tcs_axis"] = np.array([tcs_from_cesium_pressure(float(p)) for p in pcs_axis], dtype=np.float64)
        axes["Pcs_axis"] = pcs_axis
    else:
        axes["Tcs_axis"] = axis_from_args(args.tcs_min, args.tcs_max, args.tcs_count)
    checkpoint_path = args.output_dir / "emission_scan_checkpoint.npz"
    final_path = args.output_dir / "emission_scan_map.npz"
    arrays = load_or_create(checkpoint_path, axes, args.resume)

    shape = arrays["done"].shape
    total = int(arrays["done"].size)
    already_done = int(np.count_nonzero(arrays["done"]))
    start_time = time.perf_counter()
    completed_this_run = 0

    print(f"Scanning {total} points; already done: {already_done}; d_gap={args.d_gap}")

    try:
        for idx in np.ndindex(shape):
            if arrays["done"][idx]:
                continue

            te = float(arrays["TE_axis"][idx[0]])
            tc = float(arrays["TC_axis"][idx[1]])
            vo = float(arrays["Vo_axis"][idx[2]])
            tcs = float(arrays["Tcs_axis"][idx[3]])

            try:
                result = te_solver.calc_emission_point(te, tc, vo, tcs, args.d_gap)
                for name in FLOAT_FIELDS:
                    arrays[name][idx] = float(result[name])
                for name in INT_FIELDS:
                    if name != "error_code":
                        arrays[name][idx] = int(result[name])
                arrays["error_code"][idx] = 0
                arrays["converged"][idx] = bool(result["converged"])
                arrays["finite_flag"][idx] = bool(result["finite_flag"])
            except Exception as exc:  # noqa: BLE001 - scan must keep going.
                arrays["error_code"][idx] = 1
                arrays["regime"][idx] = -1
                arrays["converged"][idx] = False
                arrays["finite_flag"][idx] = False
                if completed_this_run < 10:
                    print(f"Point failed at {idx}: {exc}")

            arrays["done"][idx] = True
            completed_this_run += 1

            done_count = already_done + completed_this_run
            if args.progress_interval > 0 and completed_this_run % args.progress_interval == 0:
                elapsed = time.perf_counter() - start_time
                rate = completed_this_run / elapsed if elapsed > 0 else 0.0
                print(f"Progress: {done_count}/{total} ({rate:.2f} points/s)")

            if args.checkpoint_interval > 0 and completed_this_run % args.checkpoint_interval == 0:
                save_npz(checkpoint_path, arrays)

            if args.max_runtime_s > 0 and (time.perf_counter() - start_time) >= args.max_runtime_s:
                print(f"Stopping after max runtime {args.max_runtime_s:.1f} s; checkpoint saved.")
                break

    finally:
        save_npz(checkpoint_path, arrays)

    elapsed = time.perf_counter() - start_time
    save_npz(final_path, arrays)
    summary = summarize(arrays, elapsed)
    write_summary_files(args.output_dir, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
