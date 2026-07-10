from __future__ import annotations

import sys
import unittest
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class FakeThermoCalcModel:
    def __init__(self, n_elements: int, n_nodes: int) -> None:
        self.n_elements = n_elements
        self.n_nodes = n_nodes
        self.target_voltage = None
        self.i_guess = None
        self._input_data = type("InputData", (), {})()

    def setup_circuit_mode(self, mode_str: str, target_value: float, I_guess: float = 150.0) -> None:
        if mode_str != "fixed_u":
            raise AssertionError(f"unexpected mode {mode_str}")
        self.target_voltage = float(target_value)
        self.i_guess = float(I_guess)

    def set_temperatures(self, _temitter, _tcollector) -> None:
        return None

    def set_tcs(self, _tcs) -> None:
        return None

    def calculate(self, verbose: bool = False) -> float:
        return 1.25

    def get_global_results(self):
        voltage = float(self.target_voltage)
        current = max(0.0, 10.0 - 5.0 * voltage)
        return {
            "Iout": current,
            "Uout": voltage,
            "converged": True,
            "iteration_count": 3,
            "zero_emission_skipped": False,
        }




class CapturingThermoCalcModel(FakeThermoCalcModel):
    captured_emitter_mean = None
    captured_tcs_mean = None

    def set_temperatures(self, temitter, tcollector) -> None:
        import numpy as np

        CapturingThermoCalcModel.captured_emitter_mean = float(np.mean(temitter))

    def set_tcs(self, tcs) -> None:
        import numpy as np

        CapturingThermoCalcModel.captured_tcs_mean = float(np.mean(tcs))
class TargetMismatchThermoCalcModel(FakeThermoCalcModel):
    def get_global_results(self):
        results = super().get_global_results()
        if abs(float(self.target_voltage) - 2.0) < 1.0e-12:
            results["Uout"] = 0.5
        return results

class VenableValidationRunnerTests(unittest.TestCase):
    def test_runner_repo_root_points_to_tastin_python_package_root(self) -> None:
        from venable_validation_runner import REPO_ROOT

        self.assertTrue((REPO_ROOT / "ThermoCalc").is_dir())
        self.assertTrue((REPO_ROOT / "testModule").is_dir())
    def test_voltage_scan_selects_maximum_power_and_computes_error(self) -> None:
        from venable_table71_data import TABLE71_CASES
        from venable_validation_runner import VoltageScanConfig, scan_case_fixed_voltage

        result = scan_case_fixed_voltage(
            TABLE71_CASES[0],
            VoltageScanConfig(start_v=0.0, stop_v=2.0, points=5),
            model_factory=FakeThermoCalcModel,
        )

        self.assertEqual(result.case_id, "venable_table71_qaz_0892w")
        self.assertAlmostEqual(result.u_max_calc_v, 1.0)
        self.assertAlmostEqual(result.i_at_max_calc_a, 5.0)
        self.assertAlmostEqual(result.p_out_calc_w, 5.0)
        self.assertAlmostEqual(result.eta_calc_percent, 100.0 * 5.0 / 892.0)
        self.assertAlmostEqual(result.p_out_rel_error, (5.0 - 10.23) / 10.23)
        self.assertTrue(result.converged_all)
        self.assertEqual(result.scan_points, 5)

    def test_scan_accepts_custom_thermal_closure(self) -> None:
        from venable_single_tfe_model import VenableThermalClosure
        from venable_table71_data import TABLE71_CASES
        from venable_validation_runner import VoltageScanConfig, scan_case_fixed_voltage

        CapturingThermoCalcModel.captured_emitter_mean = None
        closure = VenableThermalClosure(
            emitter_mean_min_k=1200.0,
            emitter_mean_max_k=1200.0,
            collector_mean_min_k=700.0,
            collector_mean_max_k=700.0,
            axial_shape_amplitude=0.0,
        )
        scan_case_fixed_voltage(
            TABLE71_CASES[0],
            VoltageScanConfig(start_v=0.1, stop_v=0.2, points=2),
            model_factory=CapturingThermoCalcModel,
            thermal_closure=closure,
        )

        self.assertAlmostEqual(CapturingThermoCalcModel.captured_emitter_mean, 1200.0)

    def test_scan_accepts_thermal_network_closure(self) -> None:
        from venable_single_tfe_model import THERMAL_MODEL_THERMAL_NETWORK_V1, VenableThermalClosure
        from venable_table71_data import TABLE71_CASES
        from venable_validation_runner import VoltageScanConfig, scan_case_fixed_voltage

        CapturingThermoCalcModel.captured_emitter_mean = None
        closure = VenableThermalClosure(
            thermal_model_mode=THERMAL_MODEL_THERMAL_NETWORK_V1,
            cooling_water_mass_flow_kg_s=0.03,
        )
        scan_case_fixed_voltage(
            TABLE71_CASES[-1],
            VoltageScanConfig(start_v=0.1, stop_v=0.2, points=2),
            model_factory=CapturingThermoCalcModel,
            thermal_closure=closure,
        )

        self.assertIsNotNone(CapturingThermoCalcModel.captured_emitter_mean)
        self.assertGreater(CapturingThermoCalcModel.captured_emitter_mean, 500.0)

    def test_scan_accepts_pressure_formula_tcs_mode(self) -> None:
        from venable_single_tfe_model import TCS_MODE_PRESSURE_FORMULA, tcs_from_cesium_pressure
        from venable_table71_data import TABLE71_CASES
        from venable_validation_runner import VoltageScanConfig, scan_case_fixed_voltage

        CapturingThermoCalcModel.captured_tcs_mean = None
        scan_case_fixed_voltage(
            TABLE71_CASES[0],
            VoltageScanConfig(start_v=0.1, stop_v=0.2, points=2),
            model_factory=CapturingThermoCalcModel,
            tcs_mode=TCS_MODE_PRESSURE_FORMULA,
        )

        self.assertAlmostEqual(
            CapturingThermoCalcModel.captured_tcs_mean,
            tcs_from_cesium_pressure(TABLE71_CASES[0].pcs_torr),
        )
    def test_voltage_scan_marks_target_voltage_mismatch_invalid(self) -> None:
        from venable_table71_data import TABLE71_CASES
        from venable_validation_runner import VoltageScanConfig, scan_case_fixed_voltage

        result = scan_case_fixed_voltage(
            TABLE71_CASES[0],
            VoltageScanConfig(start_v=1.0, stop_v=2.0, points=2),
            model_factory=TargetMismatchThermoCalcModel,
        )

        self.assertFalse(result.converged_all)
        self.assertFalse(result.scan_records[-1]["target_matched"])
        self.assertTrue(result.scan_records[0]["target_matched"])
    def test_validation_summary_reports_all_cases_and_error_metrics(self) -> None:
        from venable_table71_data import TABLE71_CASES
        from venable_validation_runner import (
            VoltageScanConfig,
            run_validation_cases,
            summarize_validation_results,
        )

        results = run_validation_cases(
            TABLE71_CASES[:2],
            VoltageScanConfig(start_v=0.0, stop_v=2.0, points=5),
            model_factory=FakeThermoCalcModel,
        )
        summary = summarize_validation_results(results)

        self.assertEqual(len(results), 2)
        self.assertEqual(summary["case_count"], 2)
        self.assertGreater(summary["mean_abs_p_out_rel_error"], 0.0)
        self.assertIn("max_abs_p_out_rel_error", summary)


if __name__ == "__main__":
    unittest.main()
