from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np

from benke_thermal_network import BENKE_AVERAGE_THERMOCOUPLE_INDICES

SLEEVE_DIGITIZED_CSV = "benke_sleeve_thermocouple_12pt_digitized.csv"
WATER_BALANCE_DIGITIZED_CSV = "benke_water_balance_digitized.csv"

REFERENCE_RANGES = {
    "active_zone_power_w": (3003.0 - 1.0, 3003.0 + 1.0),
    "regulated_he_effective_k_w_m_k": (0.073, 0.087),
    "water_h_w_m2_k": (528.0, 1012.0),
}


def _range_check(value: float, low: float, high: float) -> dict[str, Any]:
    finite = math.isfinite(float(value))
    return {
        "value": float(value),
        "min": float(low),
        "max": float(high),
        "passed": bool(finite and low <= float(value) <= high),
    }


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _mean_for_indices(values_by_index: dict[int, float], indices: tuple[int, ...]) -> tuple[float, list[int]]:
    values = []
    missing = []
    for idx in indices:
        value = values_by_index.get(idx, float("nan"))
        if math.isfinite(value):
            values.append(float(value))
        else:
            missing.append(idx)
    if not values:
        return float("nan"), missing
    return float(np.mean(np.asarray(values, dtype=float))), missing


def _compare_sleeve_thermocouples(result, data_dir: Path) -> dict[str, Any]:
    path = data_dir / SLEEVE_DIGITIZED_CSV
    if not path.exists():
        return {
            "status": "missing",
            "expected_file": str(path),
            "message": "Digitized Benke 12-point sleeve thermocouple data is not available yet.",
        }
    rows = _read_csv_rows(path)
    required = {"thermocouple_index", "sleeve_outer_k"}
    if not rows or not required.issubset(rows[0].keys()):
        raise ValueError(f"{path} must contain columns: {sorted(required)}")
    exp_by_index = {int(row["thermocouple_index"]): float(row["sleeve_outer_k"]) for row in rows}
    calc = np.asarray(result.sleeve_thermocouple_temperature_k, dtype=float)
    calc_by_index = {idx + 1: float(value) for idx, value in enumerate(calc)}
    exp = []
    calc_matched = []
    ignored_indices = []
    for idx in range(1, calc.size + 1):
        if idx in exp_by_index:
            exp_value = exp_by_index[idx]
            calc_value = calc_by_index[idx]
            if math.isfinite(exp_value) and math.isfinite(calc_value):
                exp.append(exp_value)
                calc_matched.append(calc_value)
            else:
                ignored_indices.append(idx)
    if not exp:
        raise ValueError(f"{path} has no finite matching thermocouple values in 1..{calc.size}")
    exp_arr = np.asarray(exp, dtype=float)
    calc_arr = np.asarray(calc_matched, dtype=float)
    diff = calc_arr - exp_arr

    exp_benke_avg, missing_exp_avg = _mean_for_indices(exp_by_index, BENKE_AVERAGE_THERMOCOUPLE_INDICES)
    calc_benke_avg, missing_calc_avg = _mean_for_indices(calc_by_index, BENKE_AVERAGE_THERMOCOUPLE_INDICES)
    excluded_from_benke_average = [idx for idx in range(1, calc.size + 1) if idx not in set(BENKE_AVERAGE_THERMOCOUPLE_INDICES)]

    return {
        "status": "compared",
        "source_file": str(path),
        "point_count": int(exp_arr.size),
        "expected_point_count": int(calc.size),
        "ignored_indices": ignored_indices,
        "mae_k": float(np.mean(np.abs(diff))),
        "rmse_k": float(np.sqrt(np.mean(diff**2))),
        "max_abs_error_k": float(np.max(np.abs(diff))),
        "mean_error_k": float(np.mean(diff)),
        "benke_average_point_count": int(len(BENKE_AVERAGE_THERMOCOUPLE_INDICES) - len(set(missing_exp_avg + missing_calc_avg))),
        "benke_average_indices": list(BENKE_AVERAGE_THERMOCOUPLE_INDICES),
        "excluded_from_benke_average_indices": excluded_from_benke_average,
        "experimental_benke_average_k": exp_benke_avg,
        "calculated_benke_average_k": calc_benke_avg,
        "benke_average_error_k": calc_benke_avg - exp_benke_avg,
        "benke_average_abs_error_k": abs(calc_benke_avg - exp_benke_avg),
        "missing_benke_average_exp_indices": missing_exp_avg,
        "missing_benke_average_calc_indices": missing_calc_avg,
    }


def _compare_water_balance(result, data_dir: Path) -> dict[str, Any]:
    path = data_dir / WATER_BALANCE_DIGITIZED_CSV
    if not path.exists():
        return {
            "status": "missing",
            "expected_file": str(path),
            "message": "Digitized Benke water outlet or water heat-balance data is not available yet.",
        }
    rows = _read_csv_rows(path)
    if not rows:
        raise ValueError(f"{path} contains no rows")
    row = rows[0]
    comparison: dict[str, Any] = {"status": "compared", "source_file": str(path)}
    if "water_outlet_k" in row and row["water_outlet_k"] != "":
        exp = float(row["water_outlet_k"])
        calc = float(result.water_bulk_outlet_k)
        comparison["water_outlet_error_k"] = calc - exp
        comparison["water_outlet_abs_error_k"] = abs(calc - exp)
    if "water_delta_t_k" in row and row["water_delta_t_k"] != "":
        exp = float(row["water_delta_t_k"])
        calc = float(result.water_bulk_outlet_k - result.config.water_inlet_temperature_k)
        comparison["water_delta_t_error_k"] = calc - exp
        comparison["water_delta_t_abs_error_k"] = abs(calc - exp)
    if len(comparison) == 2:
        raise ValueError(f"{path} must contain water_outlet_k or water_delta_t_k")
    return comparison


def evaluate_benke_validation(result, data_dir: Path | str) -> dict[str, Any]:
    data_dir = Path(data_dir)
    range_checks = {
        "active_zone_power_w": _range_check(
            result.active_zone_power_w,
            *REFERENCE_RANGES["active_zone_power_w"],
        ),
        "regulated_he_effective_k_w_m_k": _range_check(
            result.config.regulated_he_effective_k_w_m_k,
            *REFERENCE_RANGES["regulated_he_effective_k_w_m_k"],
        ),
        "water_h_w_m2_k": _range_check(
            result.config.water_h_w_m2_k,
            *REFERENCE_RANGES["water_h_w_m2_k"],
        ),
        "energy_balance_error_w": {
            "value": float(result.energy_balance_error_w),
            "tolerance_abs_w": 1.0e-6,
            "passed": bool(abs(float(result.energy_balance_error_w)) <= 1.0e-6),
        },
    }
    sleeve = _compare_sleeve_thermocouples(result, data_dir)
    water = _compare_water_balance(result, data_dir)
    sleeve_compared = sleeve["status"] == "compared"
    water_compared = water["status"] == "compared"
    compared = sleeve_compared or water_compared
    all_range_passed = all(item.get("passed", False) for item in range_checks.values())
    if sleeve_compared and water_compared:
        status = "complete_with_digitized_data"
    elif compared:
        status = "quantitative_partial_with_digitized_data"
    else:
        status = "partial_missing_digitized_data"
    return {
        "status": status,
        "range_check_status": "passed" if all_range_passed else "failed",
        "range_checks": range_checks,
        "sleeve_thermocouple_comparison": sleeve,
        "water_balance_comparison": water,
        "data_dir": str(data_dir),
        "missing_data_note": None
        if compared
        else "Only literature range checks were possible. Add digitized Benke CSV files to enable quantitative MAE/RMSE validation.",
    }
