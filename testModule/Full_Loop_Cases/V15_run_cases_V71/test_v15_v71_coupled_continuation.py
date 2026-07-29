import unittest

from .run_v15_v71_cold_start import ColdStartConfig, build_cold_start_case
from .run_v15_v71_coupled_continuation import ContinuationConfig, configure_tec


class V15V71CoupledContinuationTests(unittest.TestCase):
    def test_configures_lookup_fixed_voltage_and_update_interval(self):
        config = ContinuationConfig(output_dir=".", dt_s=0.05, thermo_update_interval_s=0.8)
        build = build_cold_start_case(ColdStartConfig(dt_s=0.05), initialize_hydraulics=False)
        configure_tec(build, config)

        core = build["core"]
        self.assertTrue(core.enable_tec_coupled)
        self.assertTrue(core.tec_lookup_enabled)
        self.assertEqual(core.tec_circuit_mode, "fixed_u")
        self.assertAlmostEqual(core.tec_circuit_groups["main"].thermo_calc._input_data.target_val, 27.2)
        self.assertAlmostEqual(core.thermo_update_interval, 0.8)


if __name__ == "__main__":
    unittest.main()
