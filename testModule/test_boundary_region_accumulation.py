import unittest
from unittest.mock import patch

import numpy as np

from Solvers.HeatConduction.Boundary import BoundaryRegion


class BoundaryRegionAccumulationTests(unittest.TestCase):
    def test_zero_external_resistance_behaves_as_dirichlet_boundary(self):
        boundary = BoundaryRegion(shape=(1,), area_array=np.array([1.0]))
        boundary.clear_conditions()
        boundary.add_resistance_condition(T_ext=500.0, R_ext=0.0)
        boundary.update_internal_state(
            T_node=np.array([300.0]),
            R_int=np.array([2.0]),
        )

        flux = boundary.compute_net_flux_for_solver()

        self.assertAlmostEqual(float(flux[0]), 100.0)
        self.assertTrue(bool(boundary._resistance_mask[0]))
        self.assertAlmostEqual(float(boundary.T_eff[0]), 500.0)
        self.assertLess(float(boundary.R_eff[0]), 1.0e-10)

    def test_accumulate_bc_tracks_flux_conditions_in_flux_cache(self):
        boundary = BoundaryRegion(shape=(1,), area_array=np.array([1.0]))
        boundary.clear_conditions()
        boundary.add_flux_condition(np.array([12.5]))

        self.assertAlmostEqual(float(boundary.Q_sum_flux[0]), 12.5)
        self.assertAlmostEqual(float(boundary.J_sum[0]), 0.0)

    def test_zero_external_resistance_does_not_feed_inf_to_nan_to_num(self):
        original_nan_to_num = np.nan_to_num
        nonfinite_inputs = []

        def monitored_nan_to_num(x, *args, **kwargs):
            arr = np.asarray(x)
            if np.issubdtype(arr.dtype, np.number) and np.any(~np.isfinite(arr)):
                nonfinite_inputs.append(arr.copy())
            return original_nan_to_num(x, *args, **kwargs)

        boundary = BoundaryRegion(shape=(1,), area_array=np.array([1.0]))
        boundary.clear_conditions()

        with patch("numpy.nan_to_num", monitored_nan_to_num):
            boundary.add_resistance_condition(T_ext=500.0, R_ext=0.0)
            boundary.update_internal_state(
                T_node=np.array([300.0]),
                R_int=np.array([2.0]),
            )
            boundary.compute_net_flux_for_solver()

        self.assertEqual(nonfinite_inputs, [])


if __name__ == "__main__":
    unittest.main()
