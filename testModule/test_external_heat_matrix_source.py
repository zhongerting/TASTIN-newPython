import unittest

import numpy as np

from Components.ExternalHeatSources import ExternalHeatFluxBC, OrbitalMatrixHeatSource


class OrbitalMatrixHeatSourceTests(unittest.TestCase):
    def test_n6_sum_returns_first_sample(self):
        source = OrbitalMatrixHeatSource(
            shape=(6,),
            matrix_key="is58p5_w0_8p12_N6_sum",
        )

        expected = np.array([
            383.00606911,
            217.43880587,
            33.37829324,
            703.08817091,
            1392.08891942,
            887.12881210,
        ])

        np.testing.assert_allclose(source.get_heat_flux(1.0), expected, rtol=0.0, atol=1e-8)

    def test_n6_sum_interpolates_between_samples(self):
        source = OrbitalMatrixHeatSource(
            shape=(6,),
            matrix_key="is58p5_w0_8p12_N6_sum",
        )

        row_1 = np.array([
            383.00606911,
            217.43880587,
            33.37829324,
            703.08817091,
            1392.08891942,
            887.12881210,
        ])
        row_2 = np.array([
            382.95033209,
            217.40939586,
            33.37719445,
            701.00633623,
            1390.14738600,
            885.48148633,
        ])

        np.testing.assert_allclose(source.get_heat_flux(1.5), 0.5 * (row_1 + row_2), rtol=0.0, atol=1e-8)

    def test_periodic_time_wraps_to_first_sample(self):
        source = OrbitalMatrixHeatSource(
            shape=(6,),
            matrix_key="is58p5_w0_8p12_N6_sum",
        )

        np.testing.assert_allclose(source.get_heat_flux(361.0), source.get_heat_flux(1.0), rtol=0.0, atol=1e-8)


    def test_n18_sum_matrix_is_available_from_csv(self):
        source = OrbitalMatrixHeatSource(
            shape=(18,),
            matrix_key="is58p5_w0_8p12_N18_sum",
        )

        self.assertEqual(source.get_heat_flux(1.0).shape, (18,))
        self.assertAlmostEqual(float(source.get_heat_flux(1.0)[0]), 383.00606911, places=8)
        self.assertAlmostEqual(float(source.get_heat_flux(1.0)[17]), 365.98554992, places=8)
    def test_shape_mismatch_raises(self):
        with self.assertRaisesRegex(ValueError, "column count"):
            OrbitalMatrixHeatSource(
                shape=(5,),
                matrix_key="is58p5_w0_8p12_N6_sum",
            )

    def test_external_heat_flux_bc_keeps_area_conversion_contract(self):
        source = OrbitalMatrixHeatSource(
            shape=(6,),
            matrix_key="is58p5_w0_8p12_N6_sum",
        )
        bc = ExternalHeatFluxBC(
            heat_source=source,
            area_array=np.full(6, 2.0),
        )

        bc.current_time = 1.0
        bc.update_state(np.full(6, 300.0))

        np.testing.assert_allclose(bc.q_flux, 2.0 * source.get_heat_flux(1.0), rtol=0.0, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
