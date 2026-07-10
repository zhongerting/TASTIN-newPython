from __future__ import annotations

import argparse
import csv
import json
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Sequence

from benke_parameter_scan import DEFAULT_REGULATED_HE_K_VALUES, DEFAULT_WATER_H_VALUES
from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalCase, BenkeThermalNetworkConfig, solve_benke_thermal_network
from benke_validation import evaluate_benke_validation

CASE_DIR = Path(__file__).resolve().parent
DEFAULT_COOLANT_HEAT_FRACTION_VALUES = (1.0,)


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _require_sleeve_comparison(validation: dict) -> dict:
    comparison = validation["sleeve_thermocouple_comparison"]
    if comparison["status"] != "compared":
        raise FileNotFoundError(comparison["expected_file"])
    return comparison


def calibrate_to_sleeve_thermocouples(
    experimental_data_dir: Path | str,
    regulated_he_k_values: Iterable[float] = DEFAULT_REGULATED_HE_K_VALUES,
    water_h_values: Iterable[float] = DEFAULT_WATER_H_VALUES,
    base_config: BenkeThermalNetworkConfig | None = None,
) -> tuple[dict, list[dict]]:
    best, rows = calibrate_adjustable_thermal_parameters(
        experimental_data_dir,
        case=BENKE_TYPICAL_CASE,
        base_config=BenkeThermalNetworkConfig() if base_config is None else base_config,
        regulated_he_k_values=regulated_he_k_values,
        water_h_values=water_h_values,
        coolant_heat_fraction_values=((BenkeThermalNetworkConfig() if base_config is None else base_config).coolant_heat_fraction,),
    )
    compact_rows = [
        {
            "regulated_he_effective_k_w_m_k": row["regulated_he_effective_k_w_m_k"],
            "water_h_w_m2_k": row["water_h_w_m2_k"],
            "sleeve_mae_k": row["sleeve_mae_k"],
            "sleeve_rmse_k": row["sleeve_rmse_k"],
            "sleeve_max_abs_error_k": row["sleeve_max_abs_error_k"],
            "sleeve_mean_error_k": row["sleeve_mean_error_k"],
            "range_check_status": row["range_check_status"],
        }
        for row in rows
    ]
    compact_best = min(compact_rows, key=lambda row: (row["sleeve_rmse_k"], row["sleeve_mae_k"]))
    return compact_best, compact_rows


def calibrate_adjustable_thermal_parameters(
    experimental_data_dir: Path | str,
    *,
    case: BenkeThermalCase = BENKE_TYPICAL_CASE,
    base_config: BenkeThermalNetworkConfig | None = None,
    regulated_he_k_values: Iterable[float] = DEFAULT_REGULATED_HE_K_VALUES,
    water_h_values: Iterable[float] = DEFAULT_WATER_H_VALUES,
    coolant_heat_fraction_values: Iterable[float] = DEFAULT_COOLANT_HEAT_FRACTION_VALUES,
) -> tuple[dict, list[dict]]:
    experimental_data_dir = Path(experimental_data_dir)
    base_config = BenkeThermalNetworkConfig() if base_config is None else base_config
    rows: list[dict] = []
    for he_k in regulated_he_k_values:
        for water_h in water_h_values:
            for coolant_fraction in coolant_heat_fraction_values:
                config = replace(
                    base_config,
                    regulated_he_effective_k_w_m_k=float(he_k),
                    water_h_w_m2_k=float(water_h),
                    coolant_heat_fraction=float(coolant_fraction),
                )
                result = solve_benke_thermal_network(case, config)
                validation = evaluate_benke_validation(result, experimental_data_dir)
                sleeve = _require_sleeve_comparison(validation)
                water = validation["water_balance_comparison"]
                rows.append(
                    {
                        "regulated_he_effective_k_w_m_k": float(he_k),
                        "water_h_w_m2_k": float(water_h),
                        "coolant_heat_fraction": float(coolant_fraction),
                        "active_zone_power_w": float(result.active_zone_power_w),
                        "coolant_heat_input_w": float(result.heat_source_w.sum()),
                        "water_outlet_k": float(result.water_bulk_outlet_k),
                        "water_delta_t_k": float(result.water_bulk_outlet_k - config.water_inlet_temperature_k),
                        "sleeve_mae_k": float(sleeve["mae_k"]),
                        "sleeve_rmse_k": float(sleeve["rmse_k"]),
                        "sleeve_max_abs_error_k": float(sleeve["max_abs_error_k"]),
                        "sleeve_mean_error_k": float(sleeve["mean_error_k"]),
                        "sleeve_point_count": int(sleeve["point_count"]),
                        "water_outlet_abs_error_k": water.get("water_outlet_abs_error_k"),
                        "water_delta_t_abs_error_k": water.get("water_delta_t_abs_error_k"),
                        "range_check_status": validation["range_check_status"],
                    }
                )
    if not rows:
        raise ValueError("Calibration grid is empty.")
    best = min(rows, key=lambda row: (row["sleeve_rmse_k"], row["sleeve_mae_k"]))
    return best, rows


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Calibrate Benke thermal-network parameters against digitized sleeve thermocouple data.")
    parser.add_argument("--experimental-data-dir", type=Path, default=CASE_DIR / "experimental_data")
    parser.add_argument("--regulated-he-k-values", nargs="*", type=float, default=list(DEFAULT_REGULATED_HE_K_VALUES))
    parser.add_argument("--water-h-values", nargs="*", type=float, default=list(DEFAULT_WATER_H_VALUES))
    parser.add_argument("--coolant-heat-fraction-values", nargs="*", type=float, default=list(DEFAULT_COOLANT_HEAT_FRACTION_VALUES))
    parser.add_argument("--tisa-power-w", type=float, default=BENKE_TYPICAL_CASE.tisa_power_w)
    parser.add_argument("--regulated-he-pressure-torr", type=float, default=BENKE_TYPICAL_CASE.regulated_he_pressure_torr)
    parser.add_argument("--active-length-m", type=float, default=BenkeThermalNetworkConfig.active_length_m)
    parser.add_argument("--tisa-heated-length-m", type=float, default=BenkeThermalNetworkConfig.tisa_heated_length_m)
    parser.add_argument("--water-inlet-temperature-k", type=float, default=BenkeThermalNetworkConfig.water_inlet_temperature_k)
    parser.add_argument("--water-mass-flow-kg-s", type=float, default=BenkeThermalNetworkConfig.water_mass_flow_kg_s)
    parser.add_argument("--output-dir", type=Path, default=CASE_DIR / "calibration_results")
    args = parser.parse_args(argv)

    case = BenkeThermalCase(
        name="benke_calibration_case",
        tisa_power_w=args.tisa_power_w,
        regulated_he_pressure_torr=args.regulated_he_pressure_torr,
    )
    base_config = BenkeThermalNetworkConfig(
        active_length_m=args.active_length_m,
        tisa_heated_length_m=args.tisa_heated_length_m,
        water_inlet_temperature_k=args.water_inlet_temperature_k,
        water_mass_flow_kg_s=args.water_mass_flow_kg_s,
    )
    best, rows = calibrate_adjustable_thermal_parameters(
        args.experimental_data_dir,
        case=case,
        base_config=base_config,
        regulated_he_k_values=args.regulated_he_k_values,
        water_h_values=args.water_h_values,
        coolant_heat_fraction_values=args.coolant_heat_fraction_values,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(args.output_dir / "benke_adjustable_parameter_calibration_grid.csv", rows)
    with (args.output_dir / "benke_adjustable_parameter_calibration_best.json").open("w", encoding="utf-8") as handle:
        json.dump(best, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps(best, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
