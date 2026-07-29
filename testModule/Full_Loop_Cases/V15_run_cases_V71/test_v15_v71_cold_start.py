import unittest

import numpy as np

from .run_v15_v71_cold_start import ColdStartConfig, build_cold_start_case


class V15V71ColdStartTests(unittest.TestCase):
    def test_build_uses_requested_cold_start_contract(self):
        config = ColdStartConfig(duration_s=1.0, dt_s=0.1)
        build = build_cold_start_case(config, initialize_hydraulics=False)

        self.assertFalse(build["core"].enable_tec_coupled)
        self.assertAlmostEqual(float(np.sum(build["axial_power_profile"])), 1.0)
        self.assertTrue(np.allclose(build["system"].fluid_solver.T_vec, 723.0))
        for solid in build["system"].solid_components.values():
            self.assertTrue(np.allclose(solid.T, 723.0))
            self.assertEqual(solid.ode_method, "implicit_euler")
        for radiator in build["radiator_units"]:
            self.assertEqual(radiator.T_space, 4.0)


if __name__ == "__main__":
    unittest.main()
