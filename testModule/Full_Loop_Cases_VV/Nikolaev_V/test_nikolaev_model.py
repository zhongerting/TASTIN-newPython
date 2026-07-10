from __future__ import annotations

import sys
import unittest
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class NikolaevModelTests(unittest.TestCase):
    def test_source_data_contains_table2_operating_points(self) -> None:
        from nikolaev_source_data import SOURCE_DOI, TABLE2_OPERATING_POINTS

        self.assertEqual(SOURCE_DOI, "10.1063/1.47120")
        self.assertEqual(len(TABLE2_OPERATING_POINTS), 3)
        self.assertAlmostEqual(TABLE2_OPERATING_POINTS[0].voltage_v, 0.7)
        self.assertAlmostEqual(TABLE2_OPERATING_POINTS[-1].emitter_temperature_k, 1910.0)

    def test_operating_model_matches_table2_anchor_values(self) -> None:
        from nikolaev_single_tfe_model import calculate_operating_point

        low = calculate_operating_point(0.7)
        mid = calculate_operating_point(0.8)
        high = calculate_operating_point(0.9)

        self.assertAlmostEqual(low.current_a, 300.0 / 0.7)
        self.assertAlmostEqual(mid.current_a, 375.0)
        self.assertAlmostEqual(high.current_a, 300.0 / 0.9)
        self.assertAlmostEqual(low.thermal_power_kw, 4.2)
        self.assertAlmostEqual(mid.thermal_power_kw, 4.1)
        self.assertAlmostEqual(high.thermal_power_kw, 4.1)
        self.assertAlmostEqual(low.emitter_temperature_k, 1880.0)
        self.assertAlmostEqual(mid.emitter_temperature_k, 1890.0)
        self.assertAlmostEqual(high.emitter_temperature_k, 1910.0)

    def test_table3_and_table4_grids_close_at_published_points(self) -> None:
        from nikolaev_single_tfe_model import capillary_limit_diameter_mm, fuel_max_temperature_k
        from nikolaev_source_data import TABLE3_FUEL_TEMPERATURES, TABLE4_CAPILLARY_LIMITS

        for point in TABLE3_FUEL_TEMPERATURES:
            self.assertAlmostEqual(
                fuel_max_temperature_k(point.free_volume_percent, point.radial_factor, TABLE3_FUEL_TEMPERATURES),
                point.max_fuel_temperature_k,
            )
        for point in TABLE4_CAPILLARY_LIMITS:
            self.assertAlmostEqual(
                capillary_limit_diameter_mm(point.free_volume_percent, point.radial_factor, TABLE4_CAPILLARY_LIMITS),
                point.max_capillary_diameter_mm,
            )

    def test_table3_interpolates_between_published_points(self) -> None:
        from nikolaev_single_tfe_model import fuel_max_temperature_k
        from nikolaev_source_data import TABLE3_FUEL_TEMPERATURES

        interpolated = fuel_max_temperature_k(25.0, 1.075, TABLE3_FUEL_TEMPERATURES)

        self.assertGreater(interpolated, 2070.0)
        self.assertLess(interpolated, 2220.0)


if __name__ == "__main__":
    unittest.main()

class NikolaevThermoCalcPathTests(unittest.TestCase):
    def test_build_thermocalc_case_uses_fixed_literature_geometry_and_heat_input(self) -> None:
        from nikolaev_thermocalc_model import NikolaevThermalNetworkConfig, build_thermocalc_case
        from nikolaev_source_data import TABLE2_OPERATING_POINTS

        case = build_thermocalc_case(TABLE2_OPERATING_POINTS[1], NikolaevThermalNetworkConfig(n_nodes=12))

        self.assertAlmostEqual(case.geometry.interelectrode_gap_m, 0.0005)
        self.assertAlmostEqual(case.geometry.heater_length_m, 0.350)
        self.assertAlmostEqual(case.geometry.emitter_outer_diameter_m, 0.0173)
        self.assertEqual(case.arrays.temitter_k.shape, (1, 12))
        self.assertEqual(case.arrays.tcollector_k.shape, (1, 12))
        self.assertGreater(float(case.arrays.temitter_k.mean()), float(case.arrays.tcollector_k.mean()))
        self.assertAlmostEqual(float(case.arrays.tcollector_k.mean()), 870.0, delta=30.0)

    def test_fixed_voltage_scan_takes_current_from_model_factory_not_table_power_formula(self) -> None:
        from nikolaev_thermocalc_runner import ThermoCalcVoltageConfig, scan_table2_point
        from nikolaev_thermocalc_model import NikolaevThermalNetworkConfig
        from nikolaev_source_data import TABLE2_OPERATING_POINTS

        class FakeThermoCalc:
            def __init__(self, n_elements: int, n_nodes: int) -> None:
                self._input_data = type("InputData", (), {})()
                self.target_voltage = None

            def setup_circuit_mode(self, mode_str: str, target_value: float, I_guess: float = 150.0) -> None:
                self.target_voltage = float(target_value)

            def set_temperatures(self, _te, _tc) -> None:
                return None

            def set_tcs(self, _tcs) -> None:
                return None

            def calculate(self, verbose: bool = False) -> float:
                return 0.0

            def get_global_results(self):
                return {"Iout": 123.0, "Uout": self.target_voltage, "converged": True, "iteration_count": 2}

        result = scan_table2_point(
            TABLE2_OPERATING_POINTS[1],
            thermal_config=NikolaevThermalNetworkConfig(n_nodes=8),
            voltage_config=ThermoCalcVoltageConfig(i_guess_a=300.0),
            model_factory=FakeThermoCalc,
        )

        self.assertAlmostEqual(result.current_calc_a, 123.0)
        self.assertNotAlmostEqual(result.current_calc_a, 300.0 / TABLE2_OPERATING_POINTS[1].voltage_v)
        self.assertTrue(result.converged)

class NikolaevThermoElectricClosedLoopTests(unittest.TestCase):
    def test_closed_loop_iterates_electronic_heat_terms_back_to_temperatures(self) -> None:
        import numpy as np
        from nikolaev_source_data import TABLE2_OPERATING_POINTS
        from nikolaev_thermocalc_model import NikolaevThermalNetworkConfig
        from nikolaev_thermoelectric_closed_loop import (
            ClosedLoopConfig,
            solve_closed_loop_point,
        )

        class FakeThermoCalc:
            calls = 0

            def __init__(self, n_elements: int, n_nodes: int) -> None:
                self._input_data = type("InputData", (), {})()
                self.n_nodes = n_nodes
                self.te_history = []
                self.target_voltage = None

            def setup_circuit_mode(self, mode_str: str, target_value: float, I_guess: float = 150.0) -> None:
                self.target_voltage = float(target_value)

            def set_temperatures(self, te, tc) -> None:
                self.te_history.append(float(te.mean()))

            def set_tcs(self, _tcs) -> None:
                return None

            def calculate(self, verbose: bool = False) -> float:
                FakeThermoCalc.calls += 1
                return 0.0

            def get_global_results(self):
                return {"Iout": 100.0 + FakeThermoCalc.calls, "Uout": self.target_voltage, "converged": True, "iteration_count": 2}

            def get_tec_results(self, idx: int):
                n = self.n_nodes
                return {
                    "J": np.full(n, 0.005),
                    "phiE": np.full(n, 4.9),
                    "TE": np.full(n, 1900.0),
                    "UE": np.full(n, 0.9),
                    "UC": np.zeros(n),
                    "joulePowerE": np.full(n, 1.5),
                    "joulePowerC": np.full(n, 0.25),
                }

        result = solve_closed_loop_point(
            TABLE2_OPERATING_POINTS[1],
            thermal_config=NikolaevThermalNetworkConfig(n_nodes=6),
            closed_loop_config=ClosedLoopConfig(max_iterations=4, relaxation=0.5, temperature_tolerance_k=0.0),
            model_factory=FakeThermoCalc,
        )

        self.assertEqual(result.outer_iterations, 4)
        self.assertGreater(FakeThermoCalc.calls, 1)
        self.assertGreater(abs(result.emitter_temperature_mean_k - result.initial_emitter_temperature_mean_k), 1.0)
        self.assertLess(result.electron_cooling_power_w, 0.0)
        self.assertGreater(result.collector_electron_heating_power_w, 0.0)
        self.assertAlmostEqual(result.joule_power_emitter_w, 9.0)
        self.assertAlmostEqual(result.joule_power_collector_w, 1.5)

    def test_closed_loop_uses_authoritative_joule_power_arrays_not_voltage_gradient(self) -> None:
        import numpy as np
        from nikolaev_thermoelectric_closed_loop import extract_thermoelectric_feedback
        from nikolaev_thermocalc_model import NikolaevTfeGeometry

        n = 5
        geometry = NikolaevTfeGeometry(n_nodes=n)
        tec_result = {
            "J": np.full(n, 0.01),
            "phiE": np.full(n, 4.9),
            "TE": np.full(n, 1900.0),
            "UE": np.linspace(0.0, 100.0, n),
            "UC": np.linspace(10.0, -10.0, n),
            "joulePowerE": np.full(n, 2.0),
            "joulePowerC": np.full(n, 0.5),
        }

        feedback = extract_thermoelectric_feedback(tec_result, geometry)

        self.assertAlmostEqual(float(feedback.joule_power_emitter_w.sum()), 10.0)
        self.assertAlmostEqual(float(feedback.joule_power_collector_w.sum()), 2.5)
        self.assertGreater(float(np.ptp(tec_result["UE"])), 50.0)

class NikolaevPhysicalThermalHydraulicLoopTests(unittest.TestCase):
    def test_physical_loop_tracks_coolant_flow_and_energy_gain(self) -> None:
        import numpy as np
        from nikolaev_source_data import TABLE2_OPERATING_POINTS
        from nikolaev_physical_tfe_loop import PhysicalLoopConfig, solve_physical_tfe_point

        class FakeThermoCalc:
            calls = 0

            def __init__(self, n_elements: int, n_nodes: int) -> None:
                self._input_data = type("InputData", (), {})()
                self.n_nodes = n_nodes
                self.target_voltage = None

            def setup_circuit_mode(self, mode_str: str, target_value: float, I_guess: float = 150.0) -> None:
                self.target_voltage = float(target_value)

            def set_temperatures(self, _te, _tc) -> None:
                return None

            def set_tcs(self, _tcs) -> None:
                return None

            def calculate(self, verbose: bool = False) -> float:
                FakeThermoCalc.calls += 1
                return 0.0

            def get_global_results(self):
                return {"Iout": 250.0, "Uout": self.target_voltage, "converged": True, "iteration_count": 2}

            def get_tec_results(self, idx: int):
                n = self.n_nodes
                return {
                    "J": np.full(n, 0.006),
                    "phiE": np.full(n, 4.9),
                    "TE": np.full(n, 1880.0),
                    "UE": np.full(n, 0.8),
                    "UC": np.zeros(n),
                    "joulePowerE": np.full(n, 0.8),
                    "joulePowerC": np.full(n, 0.2),
                }

        result = solve_physical_tfe_point(
            TABLE2_OPERATING_POINTS[1],
            config=PhysicalLoopConfig(n_nodes=8, max_iterations=5, relaxation=0.5, temperature_tolerance_k=0.0, coolant_mass_flow_kg_s=0.02),
            model_factory=FakeThermoCalc,
        )

        self.assertGreater(result.outer_iterations, 1)
        self.assertGreater(result.coolant_outlet_temperature_k, result.coolant_inlet_temperature_k)
        self.assertGreater(result.coolant_heat_gain_w, 0.0)
        self.assertAlmostEqual(result.current_calc_a, 250.0)
        self.assertLess(abs(result.thermal_balance_residual_w), 1.0e-6)
        self.assertLess(result.collector_temperature_mean_k, result.emitter_temperature_mean_k)

    def test_lower_coolant_mass_flow_raises_outlet_temperature(self) -> None:
        from nikolaev_source_data import TABLE2_OPERATING_POINTS
        from nikolaev_physical_tfe_loop import PhysicalLoopConfig, solve_thermal_hydraulic_update
        from nikolaev_thermoelectric_closed_loop import ThermoElectricFeedback
        import numpy as np

        point = TABLE2_OPERATING_POINTS[1]
        n = 10
        heat_source = np.full(n, point.thermal_power_kw * 1000.0 / n)
        feedback = ThermoElectricFeedback(
            electron_emitter_power_w=np.full(n, -10.0),
            electron_collector_power_w=np.full(n, 7.0),
            electron_emitter_flux_w_m2=np.zeros(n),
            electron_collector_flux_w_m2=np.zeros(n),
            joule_power_emitter_w=np.full(n, 0.5),
            joule_power_collector_w=np.full(n, 0.2),
        )
        fast = solve_thermal_hydraulic_update(heat_source, feedback, PhysicalLoopConfig(n_nodes=n, coolant_mass_flow_kg_s=0.04))
        slow = solve_thermal_hydraulic_update(heat_source, feedback, PhysicalLoopConfig(n_nodes=n, coolant_mass_flow_kg_s=0.01))

        self.assertGreater(slow.coolant_outlet_temperature_k, fast.coolant_outlet_temperature_k)
        self.assertAlmostEqual(slow.coolant_heat_gain_w, fast.coolant_heat_gain_w)


class NikolaevAxialConductionTests(unittest.TestCase):
    def test_physical_axial_conduction_uses_material_k_and_geometry_area(self) -> None:
        import numpy as np
        from Materials.Solids.MoNb import MoNb
        from Materials.Solids.Molybdenum import Molybdenum
        from nikolaev_physical_tfe_loop import axial_face_conductance_w_per_k
        from nikolaev_thermocalc_model import NikolaevTfeGeometry

        geometry = NikolaevTfeGeometry(n_nodes=5)
        temperature = np.array([1800.0, 1900.0, 2000.0, 1900.0, 1800.0], dtype=float)

        emitter_g = axial_face_conductance_w_per_k(temperature, geometry.emitter_cross_area_m2, geometry.active_length_m / geometry.n_nodes, MoNb())
        collector_g = axial_face_conductance_w_per_k(temperature, geometry.collector_cross_area_m2, geometry.active_length_m / geometry.n_nodes, Molybdenum())

        self.assertEqual(emitter_g.shape, (4,))
        self.assertEqual(collector_g.shape, (4,))
        self.assertAlmostEqual(float(emitter_g[1]), float(MoNb().conductivity(1950.0)) * geometry.emitter_cross_area_m2 / 0.1)
        self.assertAlmostEqual(float(collector_g[1]), float(Molybdenum().conductivity(1950.0)) * geometry.collector_cross_area_m2 / 0.1)

    def test_physical_update_uses_real_axial_conduction_not_temperature_smoothing(self) -> None:
        import numpy as np
        from nikolaev_physical_tfe_loop import PhysicalLoopConfig, solve_thermal_hydraulic_update, zero_feedback

        heat_source = np.array([0.0, 2000.0, 0.0, 2000.0, 0.0], dtype=float)
        no_axial = solve_thermal_hydraulic_update(
            heat_source,
            zero_feedback(5),
            PhysicalLoopConfig(n_nodes=5, axial_conduction_enabled=False, max_temperature_k=10000.0),
        )
        with_axial = solve_thermal_hydraulic_update(
            heat_source,
            zero_feedback(5),
            PhysicalLoopConfig(n_nodes=5, axial_conduction_enabled=True, max_temperature_k=10000.0),
        )

        self.assertLess(float(np.ptp(with_axial.emitter_temperature_k)), float(np.ptp(no_axial.emitter_temperature_k)))
        self.assertLess(float(np.ptp(with_axial.collector_temperature_k)), float(np.ptp(no_axial.collector_temperature_k)))
        self.assertAlmostEqual(with_axial.coolant_heat_gain_w, float(np.sum(with_axial.heat_to_coolant_w)))
        self.assertAlmostEqual(with_axial.coolant_heat_gain_w, float(np.sum(heat_source)), delta=1.0e-5)
