import json
import tempfile
import unittest
from pathlib import Path

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import (
    supervise_low_power_candidates as supervisor,
)


class SupervisorTest(unittest.TestCase):
    def candidate(self, root: Path) -> supervisor.Candidate:
        return supervisor.Candidate(
            name="q060_w100",
            output_dir=root / "q060_w100",
            restart_in=root / "initial.npz",
            duration_s=500.0,
            dt_s=0.05,
            tec_update_interval_s=0.5,
            record_interval_s=2.0,
            hold_before_ramp_s=20.0,
            ramp_duration_s=100.0,
            ramp_shape="quintic",
            final_power_w=126000.0,
            final_flow_kg_s=2.46,
            checkpoint_interval_s=20.0,
            staged_recording=True,
        )

    def test_resume_command_preserves_all_candidate_recording_parameters(self):
        with tempfile.TemporaryDirectory() as raw:
            candidate = self.candidate(Path(raw))
            resume = candidate.output_dir / "latest_restart.npz"
            command = supervisor.build_command(candidate, resume_from=resume)
            joined = " ".join(map(str, command))
            self.assertEqual(command[0], str(supervisor.PYTHON_EXE))
            self.assertEqual(command[1:3], ["-m", supervisor.RUNNER_MODULE])
            for expected in (
                "--duration 500.0", "--dt 0.05",
                "--tec-update-interval 0.5", "--record-interval 2.0",
                "--hold-before-ramp 20.0", "--ramp-duration 100.0",
                "--ramp-shape quintic",
                "--final-power 126000.0", "--final-flow 2.46",
                "--checkpoint-interval 20.0", "--staged-recording",
                f"--resume-from {resume}",
            ):
                self.assertIn(expected, joined)

    def test_candidate_defaults_to_cubic_ramp(self):
        candidate = self.candidate(Path("unused"))
        values = dict(candidate.__dict__)
        values.pop("ramp_shape")
        self.assertEqual(supervisor.Candidate(**values).ramp_shape, "cubic")

    def test_exception_restart_prefers_resumable_emergency_then_latest(self):
        with tempfile.TemporaryDirectory() as raw:
            out = Path(raw)
            history = out / "history_control.csv"
            latest = out / "latest_restart.npz"
            emergency = out / "emergency_restart.npz"
            history.write_text("elapsed_s\n1.0\n", encoding="utf-8")
            latest.write_bytes(b"latest")
            emergency.write_bytes(b"emergency")
            (out / "failure.json").write_text(
                json.dumps({"restart_resumable": True}), encoding="utf-8")
            self.assertEqual(
                supervisor.select_restart(out, exceptional_exit=True), emergency)
            (out / "failure.json").write_text(
                json.dumps({"restart_resumable": False}), encoding="utf-8")
            self.assertEqual(
                supervisor.select_restart(out, exceptional_exit=True), latest)
            history.unlink()
            self.assertIsNone(
                supervisor.select_restart(out, exceptional_exit=True))

    def test_stall_uses_initialization_grace_until_history_exists(self):
        with tempfile.TemporaryDirectory() as raw:
            out = Path(raw)
            self.assertFalse(supervisor.is_stalled(
                out, started_at_s=100.0, now_s=999.0,
                initialization_grace_s=900.0, stall_timeout_s=1800.0))
            self.assertTrue(supervisor.is_stalled(
                out, started_at_s=100.0, now_s=1000.0,
                initialization_grace_s=900.0, stall_timeout_s=1800.0))
            history = out / "history_control.csv"
            history.write_text("elapsed_s\n0.0\n", encoding="utf-8")
            self.assertFalse(supervisor.is_stalled(
                out, started_at_s=100.0, now_s=1899.0,
                initialization_grace_s=900.0, stall_timeout_s=1800.0,
                history_updated_at_s=100.0))
            self.assertTrue(supervisor.is_stalled(
                out, started_at_s=100.0, now_s=1900.0,
                initialization_grace_s=900.0, stall_timeout_s=1800.0,
                history_updated_at_s=100.0))

    def test_retry_stops_after_two_retries(self):
        self.assertTrue(supervisor.can_retry(0, max_retries=2))
        self.assertTrue(supervisor.can_retry(1, max_retries=2))
        self.assertFalse(supervisor.can_retry(2, max_retries=2))

    def test_completed_and_limit_trip_are_terminal(self):
        with tempfile.TemporaryDirectory() as raw:
            out = Path(raw)
            (out / "run_summary.json").write_text(
                json.dumps({"stop_reason": "completed"}), encoding="utf-8")
            self.assertTrue(supervisor.is_terminal_result(out))
            (out / "run_summary.json").unlink()
            (out / "limit_trip.json").write_text("{}", encoding="utf-8")
            self.assertTrue(supervisor.is_terminal_result(out))


if __name__ == "__main__":
    unittest.main()
