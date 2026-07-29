import csv
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1 import (
    run_v14_210kw_fixed_power_loca_1 as runner,
)


class FixedPowerLocaTests(unittest.TestCase):
    def test_default_restart_exists(self):
        self.assertTrue(runner.DEFAULT_RESTART.is_file())

    def test_default_accident_parameters(self):
        self.assertEqual(runner.DEFAULT_EMISSIVITY, 0.8)
        self.assertEqual(runner.DEFAULT_RECORD_INTERVAL_S, 0.2)

    def test_classification(self):
        self.assertEqual(runner.classify_fluid_name("Chan_Ring1_Vol_01"), "core")
        self.assertEqual(runner.classify_fluid_name("Upper_A1_RingHP_Header_01"), "collector_ring")
        self.assertEqual(runner.classify_fluid_name("RadiatorInletHeader_Vol_01"), "ordinary_pipe")
        self.assertEqual(runner.classify_solid_name("Upper_A1_RingWall"), "pipe_wall")
        self.assertEqual(runner.classify_solid_name("Upper_A1_HPwithFin_01"), "heat_pipe")

    def test_invalid_emissivity_is_rejected(self):
        with self.assertRaises(ValueError):
            runner._validate_config(runner.LocARunConfig(radiation_emissivity=0.0))

    def test_failure_boundary_uses_existing_coolant_only(self):
        system = SimpleNamespace(
            solid_components={
                "Ring1_Collector": SimpleNamespace(T=np.array([1499.0, 1500.0])),
            },
            fluid_solver=SimpleNamespace(T_vec=np.array([2000.0])),
        )
        reason, metrics = runner.evaluate_failure_reason(
            system, runner.LocARunConfig(), coolant_present=False,
        )
        self.assertEqual(reason, "collector_temperature_limit")
        self.assertTrue(np.isnan(metrics["coolant_max_T_K"]))

    def test_staged_recording_and_minus_two_dollar_scram(self):
        self.assertEqual(
            [runner.staged_record_interval(t) for t in (0.0, 20.0, 100.0, 400.0, 600.0)],
            [0.5, 2.0, 5.0, 10.0, 20.0],
        )
        core = SimpleNamespace(point_reactor=SimpleNamespace(beta_total=0.0079321))
        config = runner.LocARunConfig(
            enable_reactivity_feedback=True, scram_time_s=5.0,
            scram_reactivity_dollars=-2.0,
        )
        self.assertEqual(runner.external_reactivity_for_elapsed(core, config, 4.95), (0.0, 0.0))
        rho, dollars = runner.external_reactivity_for_elapsed(core, config, 5.0)
        self.assertAlmostEqual(rho, -0.0158642)
        self.assertEqual(dollars, -2.0)

    def test_post_scram_low_current_disables_and_clears_tec(self):
        tfe = SimpleNamespace(cleared=False)
        tfe.clear_tec_sources = lambda: setattr(tfe, "cleared", True)
        core = SimpleNamespace(enable_tec_coupled=True)
        build = {"core": core, "tfes": {"Ring1": tfe}}
        accident = {
            "tec_open_circuit_active": False,
            "tec_open_circuit_time_s": np.nan,
            "tec_open_circuit_trigger_current_A": np.nan,
        }
        config = runner.LocARunConfig(
            enable_reactivity_feedback=True, scram_time_s=5.0,
            tec_open_circuit_current_threshold_a=0.01,
        )
        with patch.object(runner, "_tec_main_metrics", return_value={"tec_main_current_A": 0.005}):
            self.assertTrue(
                runner.maybe_transition_tec_to_open_circuit(build, accident, config, 6.0)
            )
        self.assertFalse(core.enable_tec_coupled)
        self.assertTrue(tfe.cleared)
        self.assertTrue(accident["tec_open_circuit_active"])

    def test_postprocessing_histories_are_written_separately(self):
        payload = {
            "time_s": np.array([10.0]), "accident_elapsed_s": np.array([1.0]),
            "coolant_present": np.array([0]), "fluid_volume_names": np.array(["v"]),
            "fluid_volume_categories": np.array(["core"]),
            "fluid_temperature_K": np.array([np.nan]), "fluid_pressure_Pa": np.array([np.nan]),
            "fluid_enthalpy_J_kg": np.array([np.nan]),
            "fluid_reference_temperature_K": np.array([700.0]),
            "fluid_reference_pressure_Pa": np.array([1.0e5]),
            "fluid_reference_enthalpy_J_kg": np.array([1.0]),
            "fluid_junction_names": np.array(["j"]),
            "fluid_junction_categories": np.array(["ordinary_pipe"]),
            "fluid_mass_flow_kg_s": np.array([0.0]), "fluid_velocity_m_s": np.array([0.0]),
            "fluid_reference_mass_flow_kg_s": np.array([2.0]),
            "fluid_reference_velocity_m_s": np.array([3.0]),
            "solid_names": np.array(["s"]), "solid_categories": np.array(["pipe_wall"]),
            "solid_shapes": np.array(["[2]"]), "solid_offsets": np.array([0, 2]),
            "solid_temperature_K": np.array([800.0, 801.0]),
            "tfe_names": np.array(["Ring1"]), "tfe_multipliers": np.array([6]),
            "tec_main_current_A": np.array([4.0]), "tec_main_voltage_V": np.array([5.0]),
            "tec_main_electric_power_W": np.array([20.0]), "tec_main_converged": np.array([1]),
            "feedback_fuel": np.array([1.0]), "feedback_electrode": np.array([2.0]),
            "feedback_moderator": np.array([3.0]), "feedback_reflector": np.array([4.0]),
            "feedback_total_absolute": np.array([10.0]),
            "feedback_total_change_from_accident": np.array([0.5]),
            "point_kinetics_enabled": np.array([1]), "core_power_W": np.array([210000.0]),
            "fission_power_W": np.array([197000.0]), "decay_power_W": np.array([13000.0]),
            "external_reactivity": np.array([0.0]),
            "external_reactivity_dollars": np.array([0.0]),
            "effective_temperature_feedback": np.array([0.0]),
            "total_reactivity": np.array([0.0]),
            "tec_open_circuit_active": np.array([0]),
            "tec_open_circuit_time_s": np.array([np.nan]),
        }
        for key in (
            "tec_current_density_A_m2", "tec_emitter_potential_V",
            "tec_collector_potential_V", "tec_emitter_collector_voltage_drop_V",
            "tec_electron_cooling_flux_W_m2", "tec_electron_heating_flux_W_m2",
            "tec_electron_cooling_power_W", "tec_electron_heating_power_W",
            "tec_emitter_joule_power_axial_W", "tec_collector_joule_power_axial_W",
        ):
            payload[key] = np.ones((1, 2))
        with tempfile.TemporaryDirectory(dir=runner.CASE_DIR / "runs") as directory:
            out_dir = Path(directory)
            runner.append_postprocessing_histories(out_dir, payload)
            expected_rows = {
                "history_coolant.csv": 2,
                "history_solids.csv": 2,
                "history_electrical.csv": 2,
                "history_reactivity.csv": 1,
            }
            for name, count in expected_rows.items():
                with (out_dir / name).open(encoding="utf-8", newline="") as stream:
                    self.assertEqual(len(list(csv.DictReader(stream))), count)


if __name__ == "__main__":
    unittest.main()
