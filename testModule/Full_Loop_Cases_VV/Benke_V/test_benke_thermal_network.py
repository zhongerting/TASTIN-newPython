from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class BenkeThermalNetworkTests(unittest.TestCase):
    def test_tisa_power_is_reduced_to_active_zone_power(self) -> None:
        from benke_thermal_network import BENKE_ACTIVE_ZONE_FRACTION, active_zone_power_w

        self.assertAlmostEqual(BENKE_ACTIVE_ZONE_FRACTION, 0.88)
        self.assertAlmostEqual(active_zone_power_w(3412.5), 3003.0, places=9)

    def test_centered_tisa_heat_source_conserves_active_zone_power(self) -> None:
        from benke_thermal_network import build_centered_tisa_heat_source_w

        heat = build_centered_tisa_heat_source_w(3003.0, n_nodes=60, active_length_m=0.375, heated_length_m=0.300)

        self.assertEqual(heat.shape, (60,))
        self.assertAlmostEqual(float(heat.sum()), 3003.0, places=9)
        self.assertGreater(float(heat[len(heat) // 2]), float(heat[0]))
        self.assertAlmostEqual(float(heat[0]), float(heat[-1]), places=12)

    def test_typical_case_solves_water_energy_balance_and_temperature_ordering(self) -> None:
        from benke_thermal_network import BenkeThermalNetworkConfig, BENKE_TYPICAL_CASE, solve_benke_thermal_network

        result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig())

        self.assertAlmostEqual(result.active_zone_power_w, 3003.0, delta=0.05)
        self.assertLess(abs(result.energy_balance_error_w), 1.0e-8)
        self.assertTrue(np.all(result.collector_inner_temperature_k > result.sleeve_outer_temperature_k))
        self.assertTrue(np.all(result.sleeve_outer_temperature_k > result.water_bulk_temperature_k))
        self.assertGreater(result.water_bulk_outlet_k, result.config.water_inlet_temperature_k)
        self.assertEqual(result.sleeve_thermocouple_temperature_k.shape, (12,))

    def test_lower_regulated_he_conductivity_raises_collector_temperature(self) -> None:
        from benke_thermal_network import BenkeThermalNetworkConfig, BENKE_TYPICAL_CASE, solve_benke_thermal_network

        base = solve_benke_thermal_network(
            BENKE_TYPICAL_CASE,
            BenkeThermalNetworkConfig(regulated_he_effective_k_w_m_k=0.08),
        )
        lower_k = solve_benke_thermal_network(
            BENKE_TYPICAL_CASE,
            BenkeThermalNetworkConfig(regulated_he_effective_k_w_m_k=0.073),
        )

        self.assertGreater(
            float(lower_k.collector_inner_temperature_k.mean()),
            float(base.collector_inner_temperature_k.mean()),
        )
        self.assertGreater(
            float(lower_k.sleeve_outer_temperature_k.mean()),
            float(base.sleeve_outer_temperature_k.mean()),
        )



    def test_coolant_heat_fraction_reduces_water_temperature_rise_without_changing_active_power(self) -> None:
        from benke_thermal_network import BenkeThermalNetworkConfig, BENKE_TYPICAL_CASE, solve_benke_thermal_network

        base = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig())
        reduced = solve_benke_thermal_network(
            BENKE_TYPICAL_CASE,
            BenkeThermalNetworkConfig(coolant_heat_fraction=0.9),
        )

        self.assertAlmostEqual(base.active_zone_power_w, reduced.active_zone_power_w)
        self.assertAlmostEqual(float(reduced.heat_source_w.sum()), 0.9 * reduced.active_zone_power_w)
        self.assertLess(reduced.water_bulk_outlet_k, base.water_bulk_outlet_k)
        self.assertLess(float(reduced.sleeve_outer_temperature_k.mean()), float(base.sleeve_outer_temperature_k.mean()))

    def test_default_axial_validation_domain_matches_benke_thermocouple_span(self) -> None:
        from benke_thermal_network import BenkeThermalNetworkConfig

        config = BenkeThermalNetworkConfig()

        self.assertAlmostEqual(config.active_length_m, 0.410)
        self.assertAlmostEqual(config.tisa_heated_length_m, 0.300)
    def test_benke_thermocouple_positions_and_average_indices_are_explicit(self) -> None:
        from benke_thermal_network import (
            BENKE_AVERAGE_THERMOCOUPLE_INDICES,
            BENKE_THERMOCOUPLE_POSITIONS_MM,
        )

        self.assertEqual(BENKE_THERMOCOUPLE_POSITIONS_MM, (-205.0, -163.0, -108.0, -55.0, -55.0, 0.0, 0.0, 55.0, None, 108.0, 163.0, 205.0))
        self.assertEqual(BENKE_AVERAGE_THERMOCOUPLE_INDICES, (2, 3, 4, 5, 6, 7, 8, 10, 11))

    def test_thermocouple_sampling_uses_internal_sleeve_radius_and_marks_benke_average_points(self) -> None:
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig, solve_benke_thermal_network

        result = solve_benke_thermal_network(BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig(active_length_m=0.410))

        self.assertEqual(result.sleeve_thermocouple_temperature_k.shape, (12,))
        self.assertTrue(np.isnan(result.sleeve_thermocouple_temperature_k[8]))
        self.assertEqual(result.sleeve_thermocouple_z_mm[0], -205.0)
        self.assertEqual(result.sleeve_thermocouple_z_mm[8], None)
        self.assertTrue(result.sleeve_thermocouple_included_in_benke_average[1])
        self.assertFalse(result.sleeve_thermocouple_included_in_benke_average[0])
        self.assertFalse(result.sleeve_thermocouple_included_in_benke_average[8])
        self.assertFalse(result.sleeve_thermocouple_included_in_benke_average[11])
        self.assertLess(result.sleeve_thermocouple_radius_m, result.config.sleeve_outer_radius_m)
        self.assertGreater(result.sleeve_thermocouple_radius_m, result.config.sleeve_outer_radius_m - 0.0031)
if __name__ == "__main__":
    unittest.main()
