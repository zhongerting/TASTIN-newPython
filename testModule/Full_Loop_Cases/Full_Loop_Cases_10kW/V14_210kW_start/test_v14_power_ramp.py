from types import SimpleNamespace
import unittest

import numpy as np

from .run_v14_power_ramp import _update_threshold_events
from ..run_v14_shield_radiator_startup import build_case


class _Shield:
    def __init__(self):
        self.active = True
        self.pre_step_calls = 0

    def set_active(self, active):
        self.active = bool(active)

    def pre_step(self, _dt, _time):
        self.pre_step_calls += 1


class _Pump:
    def __init__(self):
        self.target_W = 0.615

    def set_flow_rate(self, flow):
        self.target_W = float(flow)


class TestPowerRampThresholds(unittest.TestCase):
    def test_threshold_events_latch_once(self):
        system = SimpleNamespace(
            global_time=10.0,
            fluid_solver=SimpleNamespace(T_vec=np.array([300.0, 400.0])),
        )
        shield = _Shield()
        pump = _Pump()
        build = {
            "system": system,
            "core_outlet_connector": SimpleNamespace(T=400.0),
            "radiator_thermal_shield": shield,
            "pump_a": pump,
        }
        state = {
            "shield_jettisoned": False,
            "shield_threshold_k": 373.0,
            "shield_jettison_time_s": None,
            "flow_increased": False,
            "outlet_threshold_k": 500.0,
            "flow_increase_time_s": None,
            "high_flow_kg_s": 1.23,
        }


        _update_threshold_events(build, state)
        self.assertTrue(shield.active)
        self.assertEqual(pump.target_W, 0.615)

        system.global_time = 20.0
        system.fluid_solver.T_vec[:] = 373.0
        build["core_outlet_connector"].T = 500.0
        _update_threshold_events(build, state)
        self.assertFalse(shield.active)
        self.assertEqual(shield.pre_step_calls, 1)
        self.assertEqual(pump.target_W, 1.23)
        self.assertEqual(state["shield_jettison_time_s"], 20.0)
        self.assertEqual(state["flow_increase_time_s"], 20.0)

    def test_startup_uses_helium_tec_gap(self):
        build = build_case()

        for tfe in build["tfes"].values():
            config = tfe.gap_configs["tec_gap"]
            couple = tfe.couplers["tec_couple"]
            self.assertEqual(type(config.material).__name__, "Helium")
            self.assertEqual(config.h_eq, 5678.0)
            self.assertAlmostEqual(couple.k_gas / couple.gap, 5678.0)


if __name__ == "__main__":
    unittest.main()
