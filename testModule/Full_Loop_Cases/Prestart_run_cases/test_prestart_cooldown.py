import unittest

import numpy as np

from Solvers.Couplers import FluidSolidCouple
from testModule.Full_Loop_Cases.Prestart_run_cases.run_prestart_cooldown import (
    OPERATING_PARAMETERS,
    _build_case,
    _radiator_rejection_w,
    _set_local_implicit_coupling,
    _update_external_heat,
    _update_shield_heat,
    aggregate_flux_to_six,
)


class TestPrestartCooldown(unittest.TestCase):
    def test_flux_aggregation(self):
        np.testing.assert_allclose(aggregate_flux_to_six(np.arange(18.0)), [1, 4, 7, 10, 13, 16])
        np.testing.assert_allclose(aggregate_flux_to_six(np.arange(78.0)), [6, 19, 32, 45, 58, 71])

    def test_v14_shield_interface(self):
        build = _build_case("v14")
        flux = _update_shield_heat(build)
        shield = build["radiator_thermal_shield"]
        self.assertEqual(len(shield.radiator_units), 18)
        self.assertEqual(shield.solar_heat_flux_w_m2.shape, (18,))
        self.assertEqual(flux.shape, (6,))
        shield.pre_step(0.0, 0.0)
        self.assertTrue(np.isfinite(shield.last_effective_background_mean_k))

    def test_v15_shield_mount(self):
        build = _build_case("v15")
        _update_shield_heat(build)
        shield = build["radiator_thermal_shield"]
        self.assertEqual(len(shield.radiator_units), 78)
        self.assertEqual(shield.solar_heat_flux_w_m2.shape, (78,))
        self.assertFalse(build["external_heat_enabled"])

    def test_direct_radiator_external_heat_mounts_without_shield(self):
        for case, expected_units in (("v14", 18), ("v15", 78)):
            build = _build_case(case, external_heat_target="radiator")
            flux = _update_external_heat(build)
            units = build.get("radiator_units")
            if units is None:
                units = [hp for ring in build["ring_hps"] for hp in ring.hp_units]
            self.assertTrue(build["external_heat_enabled"])
            self.assertNotIn("radiator_thermal_shield", build)
            self.assertEqual(len(units), expected_units)
            self.assertEqual(build["external_heat_period_s"], 6552.0)
            self.assertEqual(build["radiator_config"].external_heat_scale_factor, 1.0)
            self.assertEqual(build["diagnostic_external_heat_source"].scale_factor, 1.0)
            self.assertEqual(flux.shape, (6,))
            self.assertTrue(np.all(np.isfinite(flux)))

    def test_shield_external_heat_uses_outer_surface_absorptivity(self):
        for case in ("v14", "v15"):
            build = _build_case(case, external_heat_target="shield")
            self.assertEqual(build["diagnostic_external_heat_source"].scale_factor, 0.1)

    def test_v14_rejection_uses_double_ring_and_includes_ring_wall(self):
        build = _build_case("v14", external_heat_target="radiator")
        symmetric = build["radiator_config"].symmetric_ring_multiplier
        heatpipes = symmetric * sum(
            ring.get_total_heat_rejection_scaled()
            for ring in build["ring_hps"]
        )
        ring_wall = symmetric * sum(
            max(0.0, -float(np.sum(solid.boundaries["right"].current_flux)))
            for solid in build["ring_solids"]
        )
        self.assertAlmostEqual(_radiator_rejection_w(build), heatpipes + ring_wall)

    def test_initial_temperature_is_configurable(self):
        build = _build_case("v15", initial_temperature_k=310.0)
        np.testing.assert_allclose(build["system"].fluid_solver.T_vec, 310.0, atol=2.0e-3)
        for solid in build["system"].solid_components.values():
            np.testing.assert_allclose(solid.T, 310.0)

    def test_current_operating_parameters_are_retained(self):
        for case in ("v14", "v15"):
            build = _build_case(case, external_heat_target="radiator", initial_temperature_k=300.0)
            expected = OPERATING_PARAMETERS[case]
            self.assertAlmostEqual(build["pump_total_head_pa"], expected["pump_head_pa"])
            self.assertAlmostEqual(build["prestart_emissivity"], expected["emissivity"])

    def test_prestart_cases_use_local_implicit_fluid_solid_coupling(self):
        for case in ("v14", "v15"):
            build = _build_case(case)
            _set_local_implicit_coupling(build["system"])
            schemes = [
                coupler.coupling_time_scheme
                for coupler in build["system"].couplers
                if isinstance(coupler, FluidSolidCouple)
                and coupler.solid_node_capacitance is not None
            ]
            self.assertTrue(schemes)
            self.assertEqual(set(schemes), {"local_implicit"})


if __name__ == "__main__":
    unittest.main()
