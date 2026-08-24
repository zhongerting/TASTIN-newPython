from types import SimpleNamespace
import csv
import io
from pathlib import Path
import unittest

from unittest.mock import MagicMock, patch

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fast_shutdown import (
    run_v14_210kw_fast_shutdown as runner,
)


class FastShutdownTests(unittest.TestCase):
    def test_default_restart_is_completed_two_orbit_state(self):
        self.assertTrue(runner.DEFAULT_RESTART.is_file())
        self.assertEqual(runner.DEFAULT_RESTART.name, "stage_01_restart.npz")

    def test_minus_two_dollars_uses_delayed_neutron_fraction(self):
        core = SimpleNamespace(point_reactor=SimpleNamespace(beta_total=0.0079321))
        self.assertAlmostEqual(runner.external_reactivity(core, -2.0), -0.0158642)

    def test_positive_scram_reactivity_is_rejected(self):
        with self.assertRaises(ValueError):
            runner._validate(runner.FastShutdownRunConfig(
                scram_reactivity_dollars=0.1,
            ))

    def test_tec_open_uses_global_iout_field(self):
        tfe = SimpleNamespace(clear_tec_sources=MagicMock())
        thermo_calc = SimpleNamespace(
            get_global_results=lambda: {"Iout": 0.0, "Uout": 50.65},
        )
        core = SimpleNamespace(enable_tec_coupled=True, thermo_calc=thermo_calc)
        state = {
            "tec_open_circuit_active": False,
            "tec_open_circuit_time_s": float("nan"),
            "tec_open_circuit_trigger_current_A": float("nan"),
        }
        opened = runner.maybe_open_tec(
            {"core": core, "tfes": {"Ring1": tfe}}, state, 0.01, 127.0,
        )
        self.assertTrue(opened)
        self.assertFalse(core.enable_tec_coupled)
        self.assertEqual(state["tec_open_circuit_time_s"], 127.0)
        tfe.clear_tec_sources.assert_called_once_with()

    def test_zero_current_history_restores_open_tec(self):
        tfe = SimpleNamespace(clear_tec_sources=MagicMock())
        core = SimpleNamespace(enable_tec_coupled=True)
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=[
            "shutdown_elapsed_s", "tec_main_current_A",
        ])
        writer.writeheader()
        writer.writerows([
            {"shutdown_elapsed_s": 126.0, "tec_main_current_A": 0.02},
            {"shutdown_elapsed_s": 127.0, "tec_main_current_A": 0.0},
            {"shutdown_elapsed_s": 5668.0, "tec_main_current_A": 0.0},
        ])
        stream.seek(0)
        with patch.object(Path, "is_file", return_value=True), \
                patch.object(Path, "open", return_value=stream):
            restored, elapsed = runner._restore_open_tec_continuation(
                {"core": core, "tfes": {"Ring1": tfe}}, Path("run_config.json"), 0.01,
            )
        self.assertTrue(restored)
        self.assertEqual(elapsed, 127.0)
        self.assertFalse(core.enable_tec_coupled)
        tfe.clear_tec_sources.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
