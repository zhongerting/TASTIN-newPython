import json
import shutil
from pathlib import Path
from types import SimpleNamespace
import unittest
import uuid
from unittest import mock

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import run_v14_low_power_fixed_i as runner


class FailureArtifactTest(unittest.TestCase):
    def test_exception_before_manifest_is_not_resumable(self):
        out_dir = runner.CASE_DIR / f"_test_failure_{uuid.uuid4().hex}"
        system = mock.Mock()
        try:
            try:
                raise RuntimeError("step exploded")
            except RuntimeError as exc:
                runner.write_failure_artifacts(out_dir, exc, system=system)
            payload = json.loads((out_dir / "failure.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["exception_type"], "RuntimeError")
            self.assertIn("step exploded", payload["traceback"])
            self.assertFalse(payload["restart_resumable"])
            self.assertFalse((out_dir / "emergency_restart.npz").exists())
            system.save_global_state.assert_not_called()
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)

    def test_exception_after_manifest_writes_resumable_restart(self):
        out_dir = runner.CASE_DIR / f"_test_failure_{uuid.uuid4().hex}"

        class _System:
            def save_global_state(self, path: str) -> None:
                Path(path).write_bytes(b"restart")

        try:
            out_dir.mkdir(parents=True)
            (out_dir / "run_manifest.json").write_text("{}", encoding="utf-8")
            try:
                raise RuntimeError("step exploded")
            except RuntimeError as exc:
                runner.write_failure_artifacts(out_dir, exc, system=_System())
            payload = json.loads((out_dir / "failure.json").read_text(encoding="utf-8"))
            self.assertTrue(payload["restart_resumable"])
            self.assertTrue((out_dir / "emergency_restart.npz").is_file())
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)

    def test_initial_preflight_trip_returns_controlled_result(self):
        out_dir = runner.CASE_DIR / f"_test_preflight_{uuid.uuid4().hex}"
        resume_path = out_dir / "resume.npz"
        manifest = {
            "source_build_config": {},
            "candidate_start_time_s": 100.0,
            "baseline": {
                "I0_A": 200.0, "Pe0_W": 10000.0, "Tout0_K": 850.0,
                "Q0_W": 210000.0, "W0_kg_s": 2.46,
            },
            "trajectory": {
                "hold_before_ramp_s": 0.0, "ramp_duration_s": 1000.0,
                "final_power_w": 120000.0, "final_flow_kg_s": 2.2,
                "duration_s": 10.0, "dt_s": 0.1,
            },
        }

        class _System:
            global_time = 100.0
            fluid_solver = SimpleNamespace(
                T_vec=[850.0], P_vec=[1.0], h_vec=[1.0], rho_vec=[1.0], W_vec=[1.0])
            solid_components = {}

            def save_global_state(self, path: str) -> None:
                Path(path).write_bytes(b"restart")

        fixed_i = {
            "tec_main_electric_power_W": 10000.0,
            "core_outlet_T_K": 850.0,
        }
        state = {
            "fixed_i_metrics": fixed_i,
            "initial_current_A": 200.0,
            "initial_electric_power_W": 10000.0,
            "initial_outlet_T_K": 850.0,
            "initial_power_W": 210000.0,
            "initial_flow_kg_s": 2.46,
            "candidate_start_time_s": 100.0,
            "elapsed_s": 0.0,
            "resume_power_setpoint_W": 210000.0,
            "resume_flow_setpoint_kg_s": 2.46,
        }
        trip = {"stop_reason": "coolant_temperature_limit", "actual": 1060.0}
        try:
            out_dir.mkdir(parents=True)
            resume_path.write_bytes(b"source")
            (out_dir / "run_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8")
            (out_dir / "history_control.csv").write_text(
                "elapsed_s\n0.0\n", encoding="utf-8")
            config = runner.LowPowerRunConfig(resume_from=resume_path, output_dir=out_dir)
            with mock.patch.object(runner, "_load_debug_config",
                                   return_value=(SimpleNamespace(
                                       inner_iter=1, external_heat_period_s=5668.144369,
                                       external_heat_time_origin_s=0.0), {})), \
                 mock.patch.object(runner, "build_case_for_run",
                                   return_value={"system": _System(), "core": object()}), \
                 mock.patch.object(runner, "prepare_candidate_state", return_value=state), \
                 mock.patch.object(runner, "_collect_safety_metrics", return_value={}), \
                 mock.patch.object(runner, "resolve_history_fields",
                                   side_effect=lambda path, fields, resuming: tuple(fields)), \
                 mock.patch.object(runner, "append_history_row"), \
                 mock.patch.object(runner, "evaluate_hard_trip", return_value=trip):
                result = runner.run_low_power_case(config)
            self.assertEqual(result["stop_reason"], "coolant_temperature_limit")
            self.assertEqual(result["trip"], trip)
            self.assertTrue((out_dir / "limit_trip.json").is_file())
            self.assertTrue((out_dir / "emergency_restart.npz").is_file())
            self.assertTrue((out_dir / "run_summary.json").is_file())
            self.assertFalse((out_dir / "failure.json").exists())
        finally:
            shutil.rmtree(out_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()