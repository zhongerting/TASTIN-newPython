import unittest

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.tec_open_circuit_accident import (
    HALF_ORBIT_S,
    ORBIT_PERIOD_S,
    disable_all_tec,
    record_interval_s,
    required_end_elapsed_s,
)


class _TFE:
    def __init__(self):
        self.cleared = False

    def clear_tec_sources(self):
        self.cleared = True


class _System:
    global_time = 10.0

    def _refresh_solid_boundary_cache(self, **kwargs):
        self.refreshed = kwargs


class _Core:
    enable_tec_coupled = True


class TecOpenCircuitAccidentTests(unittest.TestCase):
    def test_record_schedule_and_combined_end_condition(self):
        self.assertEqual(record_interval_s(0.0), 0.1)
        self.assertEqual(record_interval_s(20.0), 1.0)
        self.assertEqual(record_interval_s(19.9999999985), 1.0)
        self.assertEqual(record_interval_s(100.0), 10.0)
        self.assertEqual(record_interval_s(99.9999999985), 10.0)
        self.assertEqual(
            required_end_elapsed_s(ORBIT_PERIOD_S, 4000.0, HALF_ORBIT_S),
            4000.0 + HALF_ORBIT_S,
        )
        self.assertEqual(
            required_end_elapsed_s(ORBIT_PERIOD_S, 10.0, HALF_ORBIT_S),
            ORBIT_PERIOD_S,
        )

    def test_disable_all_tec_clears_every_tfe(self):
        tfes = {"a": _TFE(), "b": _TFE()}
        system = _System()
        core = _Core()
        disable_all_tec({"tfes": tfes, "system": system, "core": core})
        self.assertTrue(all(tfe.cleared for tfe in tfes.values()))
        self.assertFalse(core.enable_tec_coupled)
        self.assertTrue(system.refreshed["update_flux"])


if __name__ == "__main__":
    unittest.main()
