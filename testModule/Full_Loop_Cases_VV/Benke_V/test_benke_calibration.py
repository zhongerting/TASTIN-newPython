from __future__ import annotations

import csv
import sys
import unittest
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class BenkeCalibrationTests(unittest.TestCase):
    def test_calibration_recovers_known_he_k_and_water_h_from_synthetic_sleeve_data(self) -> None:
        from benke_calibration import calibrate_to_sleeve_thermocouples
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network

        target_config = BenkeThermalNetworkConfig(
            regulated_he_effective_k_w_m_k=0.08,
            water_h_w_m2_k=800.0,
        )
        target = solve_benke_thermal_network(BENKE_TYPICAL_CASE, target_config)
        data_dir = CASE_DIR / "_tmp_calibration_test"
        data_dir.mkdir(exist_ok=True)
        try:
            csv_path = data_dir / "benke_sleeve_thermocouple_12pt_digitized.csv"
            with csv_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["thermocouple_index", "sleeve_outer_k"])
                writer.writeheader()
                for idx, value in enumerate(target.sleeve_thermocouple_temperature_k, start=1):
                    writer.writerow({"thermocouple_index": idx, "sleeve_outer_k": float(value)})

            best, rows = calibrate_to_sleeve_thermocouples(
                data_dir,
                regulated_he_k_values=[0.073, 0.08, 0.087],
                water_h_values=[528.0, 800.0, 1012.0],
            )
        finally:
            if csv_path.exists():
                csv_path.unlink()
            data_dir.rmdir()

        self.assertEqual(len(rows), 9)
        self.assertAlmostEqual(best["regulated_he_effective_k_w_m_k"], 0.08)
        self.assertAlmostEqual(best["water_h_w_m2_k"], 800.0)
        self.assertAlmostEqual(best["sleeve_mae_k"], 0.0)
        self.assertAlmostEqual(best["sleeve_rmse_k"], 0.0)



    def test_adjustable_parameter_scan_reduces_real_benke_sleeve_error(self) -> None:
        from benke_calibration import calibrate_adjustable_thermal_parameters
        from benke_thermal_network import BenkeThermalCase, BenkeThermalNetworkConfig, solve_benke_thermal_network
        from benke_validation import evaluate_benke_validation

        case = BenkeThermalCase(name="benke_user_data_row", tisa_power_w=3412.0, regulated_he_pressure_torr=10.0)
        base_config = BenkeThermalNetworkConfig(
            water_inlet_temperature_k=289.71,
            water_mass_flow_kg_s=0.043518,
            water_h_w_m2_k=800.0,
            regulated_he_effective_k_w_m_k=0.08,
        )
        baseline = evaluate_benke_validation(
            solve_benke_thermal_network(case, base_config),
            CASE_DIR / "experimental_data",
        )["sleeve_thermocouple_comparison"]

        best, rows = calibrate_adjustable_thermal_parameters(
            CASE_DIR / "experimental_data",
            case=case,
            base_config=base_config,
            regulated_he_k_values=[0.08, 0.10, 0.15],
            water_h_values=[800.0, 1500.0, 5000.0],
            coolant_heat_fraction_values=[0.94, 1.0],
        )

        self.assertEqual(len(rows), 18)
        self.assertLess(best["sleeve_rmse_k"], baseline["rmse_k"])
        self.assertLess(best["sleeve_mae_k"], baseline["mae_k"])
        self.assertIn("coolant_heat_fraction", best)
if __name__ == "__main__":
    unittest.main()
