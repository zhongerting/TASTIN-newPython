"""Focused tests for selectable low-power ramp shapes."""

from __future__ import annotations

import unittest

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import (
    run_v14_low_power_fixed_i as runner,
)


class RampShapeTests(unittest.TestCase):
    def test_default_cubic_and_quintic_weights(self):
        self.assertEqual(runner.LowPowerRunConfig().ramp_shape, "cubic")
        self.assertAlmostEqual(
            runner.smooth_ramp(12.5, 10.0, 10.0, 0.0, 1.0), 0.15625)
        self.assertAlmostEqual(
            runner.smooth_ramp(
                12.5, 10.0, 10.0, 0.0, 1.0, shape="quintic"),
            0.103515625,
        )

    def test_shape_validation_rejects_unknown_value(self):
        runner.validate_run_config(runner.LowPowerRunConfig(ramp_shape="quintic"))
        with self.assertRaises(ValueError):
            runner.validate_run_config(runner.LowPowerRunConfig(ramp_shape="linear"))

    def test_control_setpoints_pass_shape_and_keep_exact_endpoints(self):
        midpoint = runner.control_setpoints_at_step_end(
            step_end_elapsed_s=12.5, hold_before_ramp_s=10.0,
            ramp_duration_s=10.0, ramp_shape="quintic",
            initial_power_w=210000.0, final_power_w=110000.0,
            initial_flow_kg_s=2.46, final_flow_kg_s=2.0,
        )
        self.assertAlmostEqual(midpoint[0], 199648.4375)
        endpoint = runner.control_setpoints_at_step_end(
            step_end_elapsed_s=20.0, hold_before_ramp_s=10.0,
            ramp_duration_s=10.0, ramp_shape="quintic",
            initial_power_w=210000.0, final_power_w=110000.0,
            initial_flow_kg_s=2.46, final_flow_kg_s=2.0,
        )
        self.assertEqual(endpoint, (110000.0, 2.0))

    def test_manifest_persists_shape_and_old_manifest_keeps_caller_value(self):
        config = runner.LowPowerRunConfig(ramp_shape="quintic")
        manifest = runner.make_manifest(
            config=config, source_config={}, candidate_start_time_s=1.0,
            initial_current_a=2.0, initial_electric_power_w=3.0,
            initial_outlet_k=4.0, initial_power_w=5.0,
            initial_flow_kg_s=6.0,
        )
        self.assertEqual(manifest["trajectory"]["ramp_shape"], "quintic")
        restored = runner.restore_trajectory_config(
            runner.LowPowerRunConfig(), manifest)
        self.assertEqual(restored.ramp_shape, "quintic")
        del manifest["trajectory"]["ramp_shape"]
        restored = runner.restore_trajectory_config(
            runner.LowPowerRunConfig(ramp_shape="quintic"), manifest)
        self.assertEqual(restored.ramp_shape, "quintic")
        self.assertEqual(
            runner.restore_trajectory_config(
                runner.LowPowerRunConfig(), manifest).ramp_shape,
            "cubic",
        )

    def test_cli_accepts_ramp_shape(self):
        args = runner._parser().parse_args(["--ramp-shape", "quintic"])
        self.assertEqual(args.ramp_shape, "quintic")


if __name__ == "__main__":
    unittest.main()
