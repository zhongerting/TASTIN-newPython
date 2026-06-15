import unittest

import numpy as np

from Solvers.Couplers import FluidSolidCouple


class LocalImplicitHeatExchangeTests(unittest.TestCase):
    def test_wall_hotter_reduces_delta_without_sign_reversal(self):
        q_to_fluid, delta_new, c_eff = FluidSolidCouple.compute_local_implicit_exchange(
            delta_old=np.array([100.0]),
            lambda_vals=np.array([10.0]),
            C_solid=np.array([20.0]),
            C_fluid=np.array([5.0]),
            dt=100.0,
        )

        self.assertGreater(q_to_fluid[0], 0.0)
        self.assertGreater(delta_new[0], 0.0)
        self.assertLess(delta_new[0], 100.0)
        self.assertAlmostEqual(c_eff[0], 4.0)

    def test_fluid_hotter_is_symmetric(self):
        q_to_fluid, delta_new, _ = FluidSolidCouple.compute_local_implicit_exchange(
            delta_old=np.array([-80.0]),
            lambda_vals=np.array([7.0]),
            C_solid=np.array([30.0]),
            C_fluid=np.array([10.0]),
            dt=50.0,
        )

        self.assertLess(q_to_fluid[0], 0.0)
        self.assertLess(delta_new[0], 0.0)
        self.assertLess(abs(delta_new[0]), 80.0)

    def test_small_dt_approaches_explicit_heat_rate(self):
        delta_old = np.array([12.0])
        lambda_vals = np.array([3.5])
        q_to_fluid, _, _ = FluidSolidCouple.compute_local_implicit_exchange(
            delta_old=delta_old,
            lambda_vals=lambda_vals,
            C_solid=np.array([1000.0]),
            C_fluid=np.array([900.0]),
            dt=1.0e-6,
        )

        self.assertAlmostEqual(q_to_fluid[0], (lambda_vals * delta_old)[0], places=4)

    def test_rejects_non_positive_dt(self):
        with self.assertRaises(ValueError):
            FluidSolidCouple.compute_local_implicit_exchange(
                delta_old=np.array([1.0]),
                lambda_vals=np.array([1.0]),
                C_solid=np.array([1.0]),
                C_fluid=np.array([1.0]),
                dt=0.0,
            )


if __name__ == "__main__":
    unittest.main()
