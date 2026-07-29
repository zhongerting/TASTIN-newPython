import unittest
from pathlib import Path

from .run_fast_steady import FastSteadyRunConfig, ScaledHeatCapacityMaterial


class FastSteadyTests(unittest.TestCase):
    def test_scaled_material_changes_only_heat_capacity(self):
        class Material:
            def heat_capacity(self, temperature):
                return 500.0 + temperature

            def conductivity(self, temperature):
                return 20.0 + 0.01 * temperature

            def density(self, temperature):
                return 8000.0 - 0.1 * temperature

        base = Material()
        scaled = ScaledHeatCapacityMaterial(base, 0.01)

        self.assertAlmostEqual(
            scaled.heat_capacity(800.0),
            0.01 * base.heat_capacity(800.0),
        )
        self.assertEqual(scaled.conductivity(800.0), base.conductivity(800.0))
        self.assertEqual(scaled.density(800.0), base.density(800.0))

    def test_accelerated_coupling_controls_are_configurable(self):
        config = FastSteadyRunConfig(
            restart_in=Path('baseline.npz'),
            inner_iter=10,
            interface_relaxation=0.1,
        )

        self.assertEqual(config.inner_iter, 10)
        self.assertEqual(config.interface_relaxation, 0.1)


if __name__ == "__main__":
    unittest.main()
