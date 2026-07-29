import shutil
import unittest
import uuid
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import (
    extract_phase_baseline as baseline,
)


CASE_DIR = Path(__file__).resolve().parent


@contextmanager
def temporary_case_directory():
    path = CASE_DIR / f"_tmp_phase_baseline_{uuid.uuid4().hex}"
    path.mkdir()
    try:
        yield path
    finally:
        shutil.rmtree(path)


class PhaseBaselineTest(unittest.TestCase):
    def _write_checkpoint(self, path: Path, time_s: float, extra=None) -> None:
        payload = {"System/global_time": np.array([time_s])}
        for level in ("Upper", "Lower"):
            for segment in range(1, 7):
                stem = f"Solid_{level}_A{segment}_path"
                payload[f"{stem}_Solid/T"] = np.full(3, 700.0 + segment)
                for node in range(3):
                    temperature = np.empty((3, 14))
                    temperature[:, 0] = np.array([800.0, 900.0, 1000.0])
                    temperature[:, 1] = np.array([810.0, 910.0, 1010.0])
                    temperature[:, 2:] = np.array([820.0, 920.0, 1020.0])[:, None]
                    payload[f"{stem}_RingHP_HP_node{node}_HP_inner/T"] = temperature.ravel()
        payload.update(extra or {})
        np.savez_compressed(path, **payload)

    def test_orbital_phase_uses_explicit_origin(self):
        self.assertAlmostEqual(
            baseline.orbital_phase(13864.2 + 5668.144369 + 12.5, 13864.2),
            12.5,
            places=8,
        )

    def test_incompatible_checkpoint_schema_is_rejected(self):
        with temporary_case_directory() as tmp:
            root = Path(tmp)
            self._write_checkpoint(root / "checkpoint_t000001s.npz", 1.0)
            self._write_checkpoint(
                root / "checkpoint_t000002s.npz", 2.0,
                {"unexpected": np.ones(2)},
            )
            with self.assertRaisesRegex(ValueError, "schema differs"):
                baseline.extract_directory(root, phase_origin_s=0.0)

    def test_extracts_topology_weighted_sections_and_marks_fins_absent(self):
        with temporary_case_directory() as tmp:
            root = Path(tmp)
            self._write_checkpoint(root / "checkpoint_t000010s.npz", 10.0)
            self._write_checkpoint(root / "stage_01_restart.npz", 20.0)
            result = baseline.extract_directory(root, phase_origin_s=0.0)

        self.assertEqual(result["checkpoint_count"], 2)
        self.assertEqual(result["checkpoints"][1]["file"], "stage_01_restart.npz")
        row = result["checkpoints"][0]
        radial_weights = np.pi * np.diff(np.array([0.0075, 0.0080, 0.0085, 0.0090]) ** 2)
        expected_eva = np.average([800.0, 900.0, 1000.0], weights=radial_weights)
        self.assertAlmostEqual(
            row["heat_pipe_temperature_k"]["all"]["evaporator"]["volume_weighted_mean"],
            expected_eva,
        )
        self.assertEqual(
            row["collector_ring_wall_temperature_k"]["all"]["cell_count"], 36)
        self.assertEqual(result["fin_temperature_k"]["available"], False)
        self.assertIn("not stored", result["fin_temperature_k"]["reason"])


if __name__ == "__main__":
    unittest.main()
