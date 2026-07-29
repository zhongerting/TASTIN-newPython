"""Check fixed-I candidate histories against the 20% acceptance window."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def _number(value):
    return float(value)


def _true(value):
    return value is True or str(value).lower() == "true"


def _relative_drift(times, values):
    mean_t = sum(times) / len(times)
    mean_v = sum(values) / len(values)
    denominator = sum((time - mean_t) ** 2 for time in times)
    slope = (
        sum((time - mean_t) * (value - mean_v) for time, value in zip(times, values))
        / denominator
        if denominator else 0.0
    )
    return abs(slope) * (max(times) - min(times)) / mean_v


def analyze_rows(rows, *, window_s=300.0):
    if not rows:
        raise ValueError("history is empty")
    end_time = max(_number(row["elapsed_s"]) for row in rows)
    window = [
        row for row in rows
        if _number(row["elapsed_s"]) >= end_time - float(window_s) - 1.0e-9
    ]
    times = [_number(row["elapsed_s"]) for row in window]
    powers = [_number(row["tec_main_electric_power_W"]) for row in window]
    finite = all(math.isfinite(value) for value in times + powers)
    converged = all(
        _true(row["tec_main_converged"]) and _true(row["fluid_converged"])
        for row in window
    )
    temperature_fields = sorted(
        key for key in window[0]
        if key.endswith("_T_K") or ("temperature" in key and key.endswith("_K"))
    )
    variations = {}
    for field in temperature_fields:
        values = [_number(row[field]) for row in window]
        mean = sum(values) / len(values)
        variations[field] = (max(values) - min(values)) / (2.0 * mean)
        finite = finite and all(math.isfinite(value) for value in values)
    mean_power = sum(powers) / len(powers)
    power_variation = (max(powers) - min(powers)) / (2.0 * mean_power)
    power_drift = _relative_drift(times, powers)
    full_window = max(times) - min(times) >= float(window_s) - 1.0e-6
    accepted = (
        full_window
        and finite
        and converged
        and min(powers) >= 2000.0
        and max(powers) <= 2200.0
        and power_variation <= 0.03
        and power_drift <= 0.01
        and all(value <= 0.03 for value in variations.values())
    )
    return {
        "accepted": accepted,
        "window_start_s": min(times),
        "window_end_s": max(times),
        "record_count": len(window),
        "full_window": full_window,
        "electric_power_min_W": min(powers),
        "electric_power_mean_W": mean_power,
        "electric_power_max_W": max(powers),
        "electric_power_half_range_fraction": power_variation,
        "electric_power_drift_fraction": power_drift,
        "converged": converged,
        "finite": finite,
        "temperature_half_range_fraction": variations,
    }


def analyze_history(path, *, window_s=300.0):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        result = analyze_rows(list(csv.DictReader(stream)), window_s=window_s)
    return {"history": str(path), **result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histories", nargs="+", type=Path)
    parser.add_argument("--window", type=float, default=300.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = [analyze_history(path, window_s=args.window) for path in args.histories]
    text = json.dumps(result, ensure_ascii=False, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
