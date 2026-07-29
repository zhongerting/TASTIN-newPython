import unittest

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_external_heat_2orbits import (
    run_v14_210kw_fixed_power_external_heat_2orbits as runner,
)


class FixedPowerExternalHeatTwoOrbitTests(unittest.TestCase):
    def test_duration_is_exactly_two_orbits(self):
        self.assertAlmostEqual(
            runner.TOTAL_DURATION_S,
            2.0 * 5668.144369,
        )

    def test_default_restart_exists(self):
        self.assertTrue(runner.DEFAULT_RESTART.is_file())


if __name__ == "__main__":
    unittest.main()
