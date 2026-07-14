import unittest

import numpy as np

from Components.circumferential_mapping import (
    build_uniform_circumferential_mapping,
    map_circumferential_intensive,
)


class CircumferentialMappingTests(unittest.TestCase):
    def test_aligned_18_to_12_has_expected_overlap_and_column_conservation(self):
        mapping = build_uniform_circumferential_mapping(18, 12)

        self.assertEqual(mapping.shape, (12, 18))
        np.testing.assert_allclose(mapping.sum(axis=0), 1.0, atol=1.0e-14)
        np.testing.assert_allclose(mapping[0, :3], [1.0, 0.5, 0.0])
        np.testing.assert_allclose(mapping[1, :3], [0.0, 0.5, 1.0])

    def test_offset_mapping_wraps_across_zero_degrees(self):
        mapping = build_uniform_circumferential_mapping(
            4,
            4,
            target_offset_deg=45.0,
        )

        np.testing.assert_allclose(mapping[0], [0.5, 0.5, 0.0, 0.0])
        np.testing.assert_allclose(mapping[3], [0.5, 0.0, 0.0, 0.5])
        np.testing.assert_allclose(mapping.sum(axis=0), 1.0, atol=1.0e-14)

    def test_extensive_mapping_preserves_total(self):
        mapping = build_uniform_circumferential_mapping(18, 12)
        source_power = np.arange(1.0, 19.0)
        target_power = mapping @ source_power

        self.assertAlmostEqual(float(np.sum(target_power)), float(np.sum(source_power)))

    def test_v14_multipliers_weight_t4_by_represented_heat_pipe_count(self):
        mapping = build_uniform_circumferential_mapping(18, 12)
        multipliers = np.array([
            5, 6, 6,
            5, 5, 6,
            5, 6, 6,
            5, 5, 6,
            5, 6, 6,
            5, 6, 6,
        ], dtype=float)
        source_t4 = np.arange(1.0, 19.0) ** 4

        target_t4 = map_circumferential_intensive(
            source_t4,
            mapping,
            source_weights=multipliers,
        )

        expected_0 = (5.0 * source_t4[0] + 3.0 * source_t4[1]) / 8.0
        expected_1 = (3.0 * source_t4[1] + 6.0 * source_t4[2]) / 9.0
        self.assertAlmostEqual(target_t4[0], expected_0)
        self.assertAlmostEqual(target_t4[1], expected_1)

    def test_constant_intensive_field_is_preserved(self):
        mapping = build_uniform_circumferential_mapping(18, 12)
        result = map_circumferential_intensive(
            np.full(18, 1234.5),
            mapping,
            source_weights=np.arange(1.0, 19.0),
        )

        np.testing.assert_allclose(result, 1234.5)

    def test_batched_values_use_final_axis_as_source_segments(self):
        mapping = build_uniform_circumferential_mapping(18, 12)
        source = np.vstack([np.arange(18.0), np.arange(18.0) + 100.0])

        result = map_circumferential_intensive(source, mapping)

        self.assertEqual(result.shape, (2, 12))
        np.testing.assert_allclose(result[1] - result[0], 100.0)

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            build_uniform_circumferential_mapping(0, 12)
        with self.assertRaises(ValueError):
            build_uniform_circumferential_mapping(18, 12, source_offset_deg=np.nan)

        mapping = build_uniform_circumferential_mapping(18, 12)
        with self.assertRaises(ValueError):
            map_circumferential_intensive(np.ones(17), mapping)
        with self.assertRaises(ValueError):
            map_circumferential_intensive(np.full(18, np.nan), mapping)
        with self.assertRaises(ValueError):
            map_circumferential_intensive(
                np.ones(18),
                mapping,
                source_weights=-np.ones(18),
            )
        with self.assertRaises(ValueError):
            map_circumferential_intensive(
                np.ones(18),
                mapping,
                source_weights=np.zeros(18),
            )


if __name__ == "__main__":
    unittest.main()
