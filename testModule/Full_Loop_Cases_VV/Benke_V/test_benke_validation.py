from __future__ import annotations

import csv
import sys
import shutil
import unittest
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class BenkeValidationComparisonTests(unittest.TestCase):
    def test_missing_digitized_data_reports_partial_status(self) -> None:
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network
        from benke_validation import evaluate_benke_validation

        result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig())
        data_dir = CASE_DIR / "_tmp_missing_data_test"
        if data_dir.exists():
            shutil.rmtree(data_dir)
        data_dir.mkdir()
        try:
            report = evaluate_benke_validation(result, data_dir)
        finally:
            shutil.rmtree(data_dir)

        self.assertEqual(report["status"], "partial_missing_digitized_data")
        self.assertTrue(report["range_checks"]["active_zone_power_w"]["passed"])
        self.assertTrue(report["range_checks"]["regulated_he_effective_k_w_m_k"]["passed"])
        self.assertEqual(report["sleeve_thermocouple_comparison"]["status"], "missing")

    def test_digitized_sleeve_data_computes_error_metrics_but_remains_partial_without_water_data(self) -> None:
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network
        from benke_validation import evaluate_benke_validation

        result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig())
        data_dir = CASE_DIR / "_tmp_digitized_sleeve_data_test"
        if data_dir.exists():
            shutil.rmtree(data_dir)
        data_dir.mkdir()
        try:
            csv_path = data_dir / "benke_sleeve_thermocouple_12pt_digitized.csv"
            with csv_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["thermocouple_index", "sleeve_outer_k"])
                writer.writeheader()
                for idx, value in enumerate(result.sleeve_thermocouple_temperature_k, start=1):
                    writer.writerow({"thermocouple_index": idx, "sleeve_outer_k": float(value) + 2.0})

            report = evaluate_benke_validation(result, data_dir)
        finally:
            shutil.rmtree(data_dir)

        comparison = report["sleeve_thermocouple_comparison"]
        self.assertEqual(report["status"], "quantitative_partial_with_digitized_data")
        self.assertEqual(comparison["status"], "compared")
        self.assertEqual(comparison["point_count"], 11)
        self.assertEqual(comparison["ignored_indices"], [9])
        self.assertAlmostEqual(comparison["mae_k"], 2.0)
        self.assertAlmostEqual(comparison["rmse_k"], 2.0)

    def test_digitized_sleeve_and_water_data_reports_complete_status(self) -> None:
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network
        from benke_validation import evaluate_benke_validation

        result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig())
        data_dir = CASE_DIR / "_tmp_digitized_complete_data_test"
        if data_dir.exists():
            shutil.rmtree(data_dir)
        data_dir.mkdir()
        try:
            sleeve_path = data_dir / "benke_sleeve_thermocouple_12pt_digitized.csv"
            with sleeve_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["thermocouple_index", "sleeve_outer_k"])
                writer.writeheader()
                for idx, value in enumerate(result.sleeve_thermocouple_temperature_k, start=1):
                    writer.writerow({"thermocouple_index": idx, "sleeve_outer_k": float(value)})

            water_path = data_dir / "benke_water_balance_digitized.csv"
            with water_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["water_outlet_k", "water_delta_t_k"])
                writer.writeheader()
                writer.writerow(
                    {
                        "water_outlet_k": float(result.water_bulk_outlet_k) + 1.5,
                        "water_delta_t_k": float(result.water_bulk_outlet_k - result.config.water_inlet_temperature_k) - 0.5,
                    }
                )

            report = evaluate_benke_validation(result, data_dir)
        finally:
            shutil.rmtree(data_dir)

        self.assertEqual(report["status"], "complete_with_digitized_data")
        self.assertEqual(report["sleeve_thermocouple_comparison"]["status"], "compared")
        self.assertEqual(report["water_balance_comparison"]["status"], "compared")
        self.assertAlmostEqual(report["water_balance_comparison"]["water_outlet_abs_error_k"], 1.5)
        self.assertAlmostEqual(report["water_balance_comparison"]["water_delta_t_abs_error_k"], 0.5)



    def test_nan_sleeve_thermocouple_values_are_ignored(self) -> None:
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network
        from benke_validation import evaluate_benke_validation

        result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig())
        data_dir = CASE_DIR / "_tmp_digitized_nan_sleeve_data_test"
        if data_dir.exists():
            shutil.rmtree(data_dir)
        data_dir.mkdir()
        try:
            csv_path = data_dir / "benke_sleeve_thermocouple_12pt_digitized.csv"
            with csv_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["thermocouple_index", "sleeve_outer_k"])
                writer.writeheader()
                for idx, value in enumerate(result.sleeve_thermocouple_temperature_k, start=1):
                    measured = "NaN" if idx == 9 else float(value) + 1.0
                    writer.writerow({"thermocouple_index": idx, "sleeve_outer_k": measured})

            report = evaluate_benke_validation(result, data_dir)
        finally:
            shutil.rmtree(data_dir)

        comparison = report["sleeve_thermocouple_comparison"]
        self.assertEqual(report["status"], "quantitative_partial_with_digitized_data")
        self.assertEqual(comparison["point_count"], 11)
        self.assertEqual(comparison["expected_point_count"], 12)
        self.assertEqual(comparison["ignored_indices"], [9])
        self.assertAlmostEqual(comparison["mae_k"], 1.0)

    def test_benke_average_statistics_exclude_end_and_failed_thermocouples(self) -> None:
        from benke_thermal_network import BenkeThermalCase, BenkeThermalNetworkConfig, solve_benke_thermal_network
        from benke_validation import evaluate_benke_validation

        case = BenkeThermalCase(name="benke_user_data_row", tisa_power_w=3412.0, regulated_he_pressure_torr=10.0)
        result = solve_benke_thermal_network(
            case,
            BenkeThermalNetworkConfig(active_length_m=0.410, water_inlet_temperature_k=289.71, water_mass_flow_kg_s=0.043518),
        )
        report = evaluate_benke_validation(result, CASE_DIR / "experimental_data")
        sleeve = report["sleeve_thermocouple_comparison"]

        self.assertEqual(sleeve["benke_average_point_count"], 9)
        self.assertEqual(sleeve["excluded_from_benke_average_indices"], [1, 9, 12])
        self.assertAlmostEqual(sleeve["experimental_benke_average_k"], 726.8288888888889)
        self.assertIn("calculated_benke_average_k", sleeve)
        self.assertIn("benke_average_error_k", sleeve)
if __name__ == "__main__":
    unittest.main()
