"""Focused tests for the independent low-power TEC update interval."""

from __future__ import annotations

import math
import unittest
from unittest import mock

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import (
    run_v14_low_power_fixed_i as runner,
)


def _manifest(config: runner.LowPowerRunConfig) -> dict:
    return runner.make_manifest(
        config=config,
        source_config={"wire_resistance_scale": 0.335},
        candidate_start_time_s=19265.0,
        initial_current_a=213.4,
        initial_electric_power_w=10812.0,
        initial_outlet_k=852.8,
        initial_power_w=210000.0,
        initial_flow_kg_s=2.46,
    )


class TecUpdateIntervalTests(unittest.TestCase):
    def test_default_is_half_second_and_validation_requires_finite_positive(self):
        self.assertEqual(runner.LowPowerRunConfig().tec_update_interval_s, 0.5)
        runner.validate_run_config(runner.LowPowerRunConfig(tec_update_interval_s=1.0))
        for value in (0.0, -0.1, math.nan, math.inf):
            with self.subTest(value=value), self.assertRaises(ValueError):
                runner.validate_run_config(
                    runner.LowPowerRunConfig(tec_update_interval_s=value))

    def test_manifest_persists_and_resume_restores_interval(self):
        manifest = _manifest(
            runner.LowPowerRunConfig(tec_update_interval_s=1.25))
        self.assertEqual(manifest["trajectory"]["tec_update_interval_s"], 1.25)
        restored = runner.restore_trajectory_config(
            runner.LowPowerRunConfig(tec_update_interval_s=0.75), manifest)
        self.assertEqual(restored.tec_update_interval_s, 1.25)

    def test_old_manifest_keeps_requested_or_default_interval(self):
        manifest = _manifest(runner.LowPowerRunConfig())
        del manifest["trajectory"]["tec_update_interval_s"]
        configured = runner.restore_trajectory_config(
            runner.LowPowerRunConfig(tec_update_interval_s=0.75), manifest)
        defaulted = runner.restore_trajectory_config(
            runner.LowPowerRunConfig(), manifest)
        self.assertEqual(configured.tec_update_interval_s, 0.75)
        self.assertEqual(defaulted.tec_update_interval_s, 0.5)

    def test_cli_accepts_independent_interval(self):
        args = runner._parser().parse_args(["--tec-update-interval", "1.5"])
        self.assertEqual(args.tec_update_interval, 1.5)

    def test_resume_sets_interval_after_immediate_refresh(self):
        events = []

        class Core:
            thermo_update_interval = None

            def setup_tec_circuit(self, *args, **kwargs):
                events.append("fixed_i")

        core = Core()
        config = runner.LowPowerRunConfig(
            dt_s=0.05, tec_update_interval_s=1.25)
        manifest = _manifest(config)
        metrics = {
            "tec_main_converged": True,
            "tec_main_current_A": 213.4,
            "tec_main_voltage_V": 20.0,
        }
        with mock.patch.object(runner, "_apply_fixed_core_power"), \
             mock.patch.object(runner, "set_total_flow_target"), \
             mock.patch.object(runner, "_apply_wire_resistance"), \
             mock.patch.object(
                 runner, "_refresh_tec",
                 side_effect=lambda *_: events.append("refresh")), \
             mock.patch.object(runner, "collect_metrics", return_value=metrics):
            runner.prepare_candidate_state(
                {"core": core}, {"wire_resistance_scale": 0.335}, config,
                current_time_s=19265.0, resume_manifest=manifest)
        self.assertIn("refresh", events)
        self.assertEqual(core.thermo_update_interval, 1.25)


if __name__ == "__main__":
    unittest.main()
