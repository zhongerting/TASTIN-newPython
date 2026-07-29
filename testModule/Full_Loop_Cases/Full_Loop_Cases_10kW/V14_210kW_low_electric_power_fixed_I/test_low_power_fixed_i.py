"""Focused tests for the V14 fixed-current low-power runner."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import shutil
import unittest
from unittest import mock
import uuid

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import (
    run_v14_low_power_fixed_i as runner,
)


class _Pump:
    def __init__(self) -> None:
        self.target_W = 2.46

    def set_flow_rate(self, value: float) -> None:
        self.target_W = float(value)


class _UncontrolledPump:
    pass


class LowPowerControlTests(unittest.TestCase):
    def test_manifest_persists_explicit_current_and_external_heat_switch(self):
        config = runner.LowPowerRunConfig(
            fixed_current_a=213.469, external_heat_enabled=False)
        manifest = runner.make_manifest(
            config=config, source_config={"wire_resistance_scale": 0.335},
            candidate_start_time_s=10.0, initial_current_a=213.469,
            initial_electric_power_w=2000.0, initial_outlet_k=730.0,
            initial_power_w=122000.0, initial_flow_kg_s=2.46,
        )
        self.assertEqual(manifest["trajectory"]["fixed_current_a"], 213.469)
        self.assertFalse(manifest["trajectory"]["external_heat_enabled"])

    def test_low_power_restart_loads_original_source_config_from_manifest(self):
        root = runner.CASE_DIR / f"_test_source_manifest_{uuid.uuid4().hex}"
        try:
            root.mkdir()
            restart = root / "final_restart.npz"
            restart.write_bytes(b"restart")
            (root / "run_config.json").write_text(
                json.dumps({"final_power_w": 118500.0}), encoding="utf-8")
            (root / "run_manifest.json").write_text(json.dumps({
                "source_build_config": {
                    "power_w": 210000.0,
                    "wire_resistance_scale": 0.335,
                },
            }), encoding="utf-8")
            self.assertEqual(
                runner.load_source_config(restart)["wire_resistance_scale"], 0.335)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_smooth_ramp_has_exact_endpoints(self):
        self.assertEqual(runner.smooth_ramp(0.0, 10.0, 20.0, 210000.0, 150000.0), 210000.0)
        self.assertEqual(runner.smooth_ramp(30.0, 10.0, 20.0, 210000.0, 150000.0), 150000.0)
        midpoint = runner.smooth_ramp(20.0, 10.0, 20.0, 210000.0, 150000.0)
        self.assertAlmostEqual(midpoint, 180000.0)

    def test_flow_target_updates_only_controlled_series_pump(self):
        build = {'pump_a': _Pump(), 'pump_b': _UncontrolledPump()}
        runner.set_total_flow_target(build, 1.75)
        self.assertEqual(build['pump_a'].target_W, 1.75)
        self.assertFalse(hasattr(build['pump_b'], 'target_W'))

    def test_flow_target_requires_a_controlled_pump(self):
        with self.assertRaises(RuntimeError):
            runner.set_total_flow_target(
                {'pump_a': _UncontrolledPump(), 'pump_b': _UncontrolledPump()},
                1.75,
            )

    def test_fixed_i_requires_convergence_and_current_match(self):
        runner.validate_fixed_i_result({
            "tec_main_converged": True,
            "tec_main_current_A": 213.0,
            "tec_main_voltage_V": 49.0,
        }, target_current_a=213.0)
        with self.assertRaises(RuntimeError):
            runner.validate_fixed_i_result({
                "tec_main_converged": False,
                "tec_main_current_A": 0.0,
                "tec_main_voltage_V": 60.0,
            }, target_current_a=213.0)
        with self.assertRaises(RuntimeError):
            runner.validate_fixed_i_result({
                "tec_main_converged": True,
                "tec_main_current_A": 210.0,
                "tec_main_voltage_V": 49.0,
            }, target_current_a=213.0)

    def test_step_end_setpoints_apply_exact_ramp_endpoint(self):
        q_set, w_set = runner.control_setpoints_at_step_end(
            step_end_elapsed_s=110.0,
            hold_before_ramp_s=10.0,
            ramp_duration_s=100.0,
            initial_power_w=210000.0,
            final_power_w=125000.0,
            initial_flow_kg_s=2.46,
            final_flow_kg_s=1.90,
        )
        self.assertEqual(q_set, 125000.0)
        self.assertEqual(w_set, 1.90)

    def test_staged_recording_uses_loca_schedule(self):
        expected = {0.0: 0.5, 20.0: 0.5, 20.1: 2.0, 100.1: 5.0,
                    400.1: 10.0, 600.1: 20.0}
        for elapsed, interval in expected.items():
            self.assertEqual(runner.record_interval_for_elapsed(elapsed, staged=True), interval)
        self.assertEqual(runner.record_interval_for_elapsed(999.0, staged=False, default_s=7.0), 7.0)

    def test_history_row_is_visible_immediately(self):
        tmp = runner.CASE_DIR / f"_test_history_{uuid.uuid4().hex}"
        try:
            tmp.mkdir()
            path = tmp / "history.csv"
            runner.append_history_row(path, {"time_s": 1.0, "value": 2.0}, ("time_s", "value"))
            with path.open(newline="", encoding="utf-8") as stream:
                self.assertEqual(list(csv.DictReader(stream)), [{"time_s": "1.0", "value": "2.0"}])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    def test_checkpoint_save_replaces_complete_file(self):
        class _System:
            def save_global_state(self, path: str) -> None:
                Path(path).write_bytes(b"new checkpoint")

        tmp = runner.CASE_DIR / f"_test_checkpoint_{uuid.uuid4().hex}"
        try:
            tmp.mkdir()
            path = tmp / "latest_restart.npz"
            path.write_bytes(b"old checkpoint")
            runner.save_checkpoint_atomic(_System(), path)
            self.assertEqual(path.read_bytes(), b"new checkpoint")
            self.assertFalse((tmp / "latest_restart.tmp.npz").exists())
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    def test_hard_trip_rejects_temperature_solver_and_nonfinite_states(self):
        safe = {
            "fluid_converged": True, "tec_main_converged": True,
            "tec_main_current_A": 213.0, "tec_main_voltage_V": 20.0,
            "core_outlet_T_K": 850.0, "fuel_max_T_K": 2000.0,
            "collector_max_T_K": 1000.0, "emitter_max_T_K": 2200.0,
            "coolant_max_T_K": 900.0, "moderator_max_T_K": 800.0,
            "reflector_max_T_K": 850.0, "channel_wall_max_T_K": 900.0,
        }
        self.assertIsNone(runner.evaluate_hard_trip(safe, initial_outlet_k=850.0))
        for field, value, reason in (
            ("fluid_converged", False, "hydraulic_nonconvergence"),
            ("tec_main_converged", False, "tec_nonconvergence"),
            ("fuel_max_T_K", 2700.1, "fuel_temperature_limit"),
            ("core_outlet_T_K", 850.6, "core_outlet_temperature_limit"),
            ("tec_main_voltage_V", math.nan, "nonfinite_metric"),
        ):
            row = dict(safe)
            row[field] = value
            self.assertEqual(runner.evaluate_hard_trip(row, initial_outlet_k=850.0)["stop_reason"], reason)

    def test_outlet_soft_limit_can_be_deferred_during_fixed_i_handoff(self):
        metrics = {
            "fluid_converged": True, "tec_main_converged": True,
            "tec_main_current_A": 213.0, "tec_main_voltage_V": 20.0,
            "core_outlet_T_K": 854.0, "fuel_max_T_K": 2000.0,
            "collector_max_T_K": 1000.0, "emitter_max_T_K": 2200.0,
            "coolant_max_T_K": 900.0, "moderator_max_T_K": 800.0,
            "reflector_max_T_K": 850.0,
        }
        self.assertFalse(runner.outlet_limit_is_active(0.1))
        self.assertTrue(runner.outlet_limit_is_active(1.0))
        self.assertIsNone(runner.evaluate_hard_trip(
            metrics, initial_outlet_k=852.8, enforce_outlet_limit=False))
        self.assertEqual(runner.evaluate_hard_trip(
            metrics, initial_outlet_k=852.8, enforce_outlet_limit=True
        )["stop_reason"], "core_outlet_temperature_limit")

    def test_interval_due_accepts_floating_point_endpoint(self):
        self.assertTrue(runner.interval_due(0.0999999999985, 0.0, 0.1))

    def test_record_due_accepts_floating_point_endpoint(self):
        self.assertTrue(runner.record_due(
            elapsed_s=0.0999999999985, last_record_s=0.0,
            interval_s=0.1, duration_s=0.1, stopped=False,
        ))

    def test_failed_result_maps_to_nonzero_exit(self):
        self.assertEqual(runner.exit_code_for_result({"stop_reason": "completed"}), 0)
        self.assertNotEqual(runner.exit_code_for_result({"stop_reason": "tec_nonconvergence"}), 0)

    def test_resume_preparation_directly_restores_saved_fixed_i(self):
        class _Core:
            def __init__(self):
                self.calls = []

            def setup_tec_circuit(self, *args, **kwargs):
                self.calls.append((args, kwargs))

        core = _Core()
        manifest = {
            "candidate_start_time_s": 1000.0,
            "baseline": {"I0_A": 213.4, "Pe0_W": 10812.0, "Tout0_K": 852.8,
                         "Q0_W": 210000.0, "W0_kg_s": 2.46},
        }
        fixed_i = {"tec_main_converged": True, "tec_main_current_A": 213.4,
                   "tec_main_voltage_V": 20.0}
        with mock.patch.object(runner, "_apply_fixed_core_power"), \
             mock.patch.object(runner, "set_total_flow_target"), \
             mock.patch.object(runner, "_apply_wire_resistance"), \
             mock.patch.object(runner, "_refresh_tec"), \
             mock.patch.object(runner, "collect_metrics", return_value=fixed_i) as collect:
            state = runner.prepare_candidate_state(
                {"core": core}, {"wire_resistance_scale": 0.335},
                runner.LowPowerRunConfig(), current_time_s=1120.0,
                resume_manifest=manifest,
            )
        self.assertEqual(core.calls[0][0][:2], ("fixed_i", 213.4))
        self.assertEqual(collect.call_count, 1)
        self.assertEqual(state["elapsed_s"], 120.0)
        self.assertEqual(state["candidate_start_time_s"], 1000.0)
    def test_config_rejects_invalid_time_and_control_values(self):
        for change in (
            {"dt_s": 0.0}, {"duration_s": 0.0}, {"record_interval_s": -1.0},
            {"hold_before_ramp_s": -1.0}, {"ramp_duration_s": -1.0},
            {"checkpoint_interval_s": -1.0}, {"final_power_w": 0.0},
            {"final_flow_kg_s": math.nan}, {"fluid_max_iter": 0},
        ):
            values = dict(runner.LowPowerRunConfig().__dict__)
            values.update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                runner.validate_run_config(runner.LowPowerRunConfig(**values))

    def test_manifest_persists_source_baseline_and_trajectory(self):
        config = runner.LowPowerRunConfig(
            hold_before_ramp_s=20.0, ramp_duration_s=1000.0,
            final_power_w=120000.0, final_flow_kg_s=2.1,
        )
        manifest = runner.make_manifest(
            config=config, source_config={"wire_resistance_scale": 0.335},
            candidate_start_time_s=19264.65, initial_current_a=213.4,
            initial_electric_power_w=10812.0, initial_outlet_k=852.8,
            initial_power_w=210000.0, initial_flow_kg_s=2.46,
        )
        self.assertEqual(manifest["source_build_config"]["wire_resistance_scale"], 0.335)
        self.assertEqual(manifest["baseline"], {
            "I0_A": 213.4, "Pe0_W": 10812.0, "Tout0_K": 852.8,
            "Q0_W": 210000.0, "W0_kg_s": 2.46,
        })
        self.assertEqual(manifest["candidate_start_time_s"], 19264.65)
        self.assertEqual(manifest["trajectory"]["ramp_duration_s"], 1000.0)
    def test_resume_restores_saved_trajectory(self):
        manifest = {"trajectory": {
            "hold_before_ramp_s": 20.0, "ramp_duration_s": 1000.0,
            "final_power_w": 120000.0, "final_flow_kg_s": 2.1,
            "duration_s": 1600.0, "dt_s": 0.2,
        }}
        restored = runner.restore_trajectory_config(
            runner.LowPowerRunConfig(), manifest)
        self.assertEqual(restored.duration_s, 1600.0)
        self.assertEqual(restored.ramp_duration_s, 1000.0)
        self.assertEqual(restored.final_power_w, 120000.0)
        self.assertEqual(restored.resume_from, runner.LowPowerRunConfig().resume_from)

    def test_resume_build_does_not_reapply_source_fixed_power(self):
        debug = object()
        with mock.patch.object(runner, "build_debug_case", return_value={}) as build:
            runner.build_case_for_run(debug, resuming=True)
        build.assert_called_once_with(debug, apply_fixed_power=False)

    def test_resume_applies_current_qw_before_tec_refresh_and_labels_first_row(self):
        events = []

        class _Core:
            def setup_tec_circuit(self, *args, **kwargs):
                events.append("fixed_i")

        manifest = {
            "candidate_start_time_s": 1000.0,
            "baseline": {"I0_A": 213.4, "Pe0_W": 10812.0, "Tout0_K": 852.8,
                         "Q0_W": 210000.0, "W0_kg_s": 2.46},
        }
        config = runner.LowPowerRunConfig(
            hold_before_ramp_s=20.0, ramp_duration_s=1000.0,
            final_power_w=120000.0, final_flow_kg_s=2.1,
        )
        fixed_i = {"tec_main_converged": True, "tec_main_current_A": 213.4,
                   "tec_main_voltage_V": 20.0}
        with mock.patch.object(runner, "_apply_fixed_core_power",
                               side_effect=lambda build, q: events.append(("q", q))), \
             mock.patch.object(runner, "set_total_flow_target",
                               side_effect=lambda build, w: events.append(("w", w))), \
             mock.patch.object(runner, "_apply_wire_resistance"), \
             mock.patch.object(runner, "_refresh_tec",
                               side_effect=lambda core, t: events.append("refresh")), \
             mock.patch.object(runner, "collect_metrics", return_value=fixed_i):
            state = runner.prepare_candidate_state(
                {"core": _Core()}, {"wire_resistance_scale": 0.335}, config,
                current_time_s=1120.0, resume_manifest=manifest,
            )
        expected_q, expected_w = runner.control_setpoints_at_step_end(
            step_end_elapsed_s=120.0, hold_before_ramp_s=20.0,
            ramp_duration_s=1000.0, initial_power_w=210000.0,
            final_power_w=120000.0, initial_flow_kg_s=2.46,
            final_flow_kg_s=2.1,
        )
        self.assertEqual(events[0], ("q", expected_q))
        self.assertEqual(events[1], ("w", expected_w))
        self.assertLess(events.index(("w", expected_w)), events.index("refresh"))
        self.assertEqual(state["resume_power_setpoint_W"], expected_q)
        self.assertEqual(state["resume_flow_setpoint_kg_s"], expected_w)

    def test_resume_rewinds_history_to_checkpoint_before_rewriting_that_time(self):
        root = runner.CASE_DIR / f"_test_resume_rewind_{uuid.uuid4().hex}"
        try:
            root.mkdir()
            path = root / "history.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=("elapsed_s", "value"))
                writer.writeheader()
                for elapsed in (0.0, 0.1, 0.15, 0.2):
                    writer.writerow({"elapsed_s": elapsed, "value": elapsed})
            runner.rewind_history_for_resume(path, checkpoint_elapsed_s=0.1)
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([float(row["elapsed_s"]) for row in rows], [0.0])
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def test_resume_history_requires_exact_existing_header(self):
        path = runner.CASE_DIR / "_test_resume_history_header.csv"
        try:
            path.write_text("a,b\n1,2\n", encoding="utf-8")
            self.assertEqual(runner.resolve_history_fields(path, ("a", "b"), resuming=True),
                             ("a", "b"))
            with self.assertRaises(RuntimeError):
                runner.resolve_history_fields(path, ("a", "c"), resuming=True)
        finally:
            if path.exists():
                path.unlink()

    def test_cladding_temperatures_are_monitored_without_an_invented_limit(self):
        metrics = {
            "fluid_converged": True, "tec_main_converged": True,
            "tec_main_current_A": 213.0, "tec_main_voltage_V": 20.0,
            "core_outlet_T_K": 850.0, "fuel_max_T_K": 2000.0,
            "collector_max_T_K": 1000.0, "emitter_max_T_K": 2200.0,
            "coolant_max_T_K": 900.0, "moderator_max_T_K": 800.0,
            "reflector_max_T_K": 850.0,
            "inner_clad_max_T_K": 1400.0, "outer_clad_max_T_K": 1500.0,
        }
        self.assertIsNone(runner.evaluate_hard_trip(metrics, initial_outlet_k=850.0))

    def test_final_setpoint_validation_checks_controlled_pump(self):
        core = type('Core', (), {'last_total_core_power': 120000.0})()
        build = {'pump_a': _Pump(), 'pump_b': _UncontrolledPump()}
        build['pump_a'].target_W = 1.9
        runner.validate_final_setpoints(build, core, 120000.0, 1.9)
        build['pump_a'].target_W = 2.0
        with self.assertRaises(RuntimeError):
            runner.validate_final_setpoints(build, core, 120000.0, 1.9)
if __name__ == "__main__":
    unittest.main()
