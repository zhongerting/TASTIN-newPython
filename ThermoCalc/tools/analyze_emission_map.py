from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "emission_scan_outputs" / "emission_scan_map.npz"


def load_scan(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as loaded:
        return {name: loaded[name] for name in loaded.files}


def summarize(arrays: dict[str, np.ndarray]) -> dict[str, Any]:
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


def write_summary(output_dir: Path, summary: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "emission_scan_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    with (output_dir / "emission_scan_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(("key", "value"))
        for key, value in summary.items():
            writer.writerow((key, value))


def projection_rows(
    arrays: dict[str, np.ndarray],
    axis_a: str,
    axis_b: str,
) -> list[dict[str, Any]]:
    axes = ("TE_axis", "TC_axis", "Vo_axis", "Tcs_axis")
    axis_index = {name: i for i, name in enumerate(axes)}
    ia = axis_index[axis_a]
    ib = axis_index[axis_b]
    reduce_axes = tuple(i for i in range(4) if i not in (ia, ib))

    done = arrays["done"]
    finite = done & arrays["finite_flag"]
    converged = done & arrays["converged"]
    failed = done & ~(arrays["finite_flag"] & arrays["converged"])

    done_count = np.sum(done, axis=reduce_axes)
    finite_count = np.sum(finite, axis=reduce_axes)
    converged_count = np.sum(converged, axis=reduce_axes)
    failed_count = np.sum(failed, axis=reduce_axes)
    iter_sum = np.sum(np.where(done, arrays["iteration_count"], 0), axis=reduce_axes)
    j_sum = np.sum(np.where(finite, arrays["J"], 0.0), axis=reduce_axes)

    rows: list[dict[str, Any]] = []
    axis_a_values = arrays[axis_a]
    axis_b_values = arrays[axis_b]
    for i, a_value in enumerate(axis_a_values):
        for j, b_value in enumerate(axis_b_values):
            idx = (i, j)
            n_done = int(done_count[idx])
            n_finite = int(finite_count[idx])
            rows.append(
                {
                    axis_a.replace("_axis", ""): float(a_value),
                    axis_b.replace("_axis", ""): float(b_value),
                    "done_count": n_done,
                    "finite_rate": float(n_finite / n_done) if n_done else 0.0,
                    "converged_rate": float(converged_count[idx] / n_done) if n_done else 0.0,
                    "failed_or_nonfinite_rate": float(failed_count[idx] / n_done) if n_done else 0.0,
                    "iteration_mean": float(iter_sum[idx] / n_done) if n_done else 0.0,
                    "J_mean_finite": float(j_sum[idx] / n_finite) if n_finite else float("nan"),
                }
            )
    return rows


def write_projection(output_dir: Path, name: str, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path = output_dir / f"projection_{name}.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze a ThermoCalc emission scan NPZ.")
    parser.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir or args.input.parent
    arrays = load_scan(args.input)

    summary = summarize(arrays)
    write_summary(output_dir, summary)
    write_projection(output_dir, "TE_TC", projection_rows(arrays, "TE_axis", "TC_axis"))
    write_projection(output_dir, "Vo_Tcs", projection_rows(arrays, "Vo_axis", "Tcs_axis"))
    write_projection(output_dir, "TE_Vo", projection_rows(arrays, "TE_axis", "Vo_axis"))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
