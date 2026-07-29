"""Smoke tests for the V14_10kW 210 kW debug runner."""

from __future__ import annotations

import json
import inspect
import tempfile
import unittest
from pathlib import Path

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    _load_restart_with_coupling_dt,
    DebugRunConfig,
    build_debug_case,
    run_debug_case,
)


class V14DebugRunnerSmokeTest(unittest.TestCase):
    def test_build_debug_case_can_skip_fixed_power_handoff(self):
        parameters = inspect.signature(build_debug_case).parameters
        self.assertIn('apply_fixed_power', parameters)
        parameter = parameters['apply_fixed_power']
        self.assertTrue(parameter.default)

    def test_restart_sync_receives_saved_positive_dt(self):
        class System:
            def __init__(self):
                self.received_dt = None

            def _run_couplers(self, dt=None):
                self.received_dt = dt

            def load_global_state(self, unused_path):
                self._last_dt = 0.05
                self._run_couplers()

        system = System()
        _load_restart_with_coupling_dt(
            system,
            'restart.npz',
            fallback_dt=0.01,
        )

        self.assertEqual(system.received_dt, 0.05)
        self.assertNotIn('_run_couplers', system.__dict__)

    def test_runner_writes_restart_and_key_metrics(self):
        with tempfile.TemporaryDirectory(dir=r"E:\tmp") as tmp:
            out_dir = Path(tmp) / "v14_debug_smoke"
            result = run_debug_case(DebugRunConfig(
                output_dir=out_dir,
                stage_durations_s=(0.02,),
                dt_s=0.02,
                record_interval_s=0.02,
                inner_iter=1,
                fluid_max_iter=20,
            ))
            self.assertTrue((out_dir / "stage_01_restart.npz").exists())
            self.assertTrue((out_dir / "stage_01_summary.json").exists())
            self.assertTrue((out_dir / "history.csv").exists())
            self.assertTrue((out_dir / "latest_state.json").exists())

            config = json.loads((out_dir / "run_config.json").read_text(encoding="utf-8"))
            self.assertFalse(config["tec_electrical_calculation_enabled"])
            self.assertFalse(config["external_heat_enabled"])
            self.assertEqual(config["solid_ode_method"], "implicit_euler")
            self.assertAlmostEqual(config["power_w"], 210000.0)
            self.assertAlmostEqual(config["target_flow_kg_s"], 2.46)
            self.assertAlmostEqual(config["initial_temperature_k"], 754.15)
            self.assertAlmostEqual(config["space_temperature_k"], 4.0)

            latest = result["latest_metrics"]
            for key in (
                "core_inlet_T_K",
                "core_outlet_T_K",
                "radiator_heat_rejection_W",
                "upper_ring_heatpipe_rejection_W",
                "lower_ring_heatpipe_rejection_W",
                "pump_required_head_total_Pa",
            ):
                self.assertIn(key, latest)
            self.assertAlmostEqual(latest["pump_a_flow_kg_s"], 2.46, places=6)
            self.assertAlmostEqual(latest["pump_b_flow_kg_s"], 2.46, places=6)


if __name__ == "__main__":
    unittest.main()
