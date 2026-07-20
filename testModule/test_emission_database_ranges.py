from pathlib import Path
import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
THERMOCALC_DIR = REPO_ROOT / "ThermoCalc"
THERMOCALC_TOOLS_DIR = THERMOCALC_DIR / "tools"
for path in (THERMOCALC_DIR, THERMOCALC_TOOLS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from tools import emission_database
from tools.emission_database import REGIONS, points_for_axes, region_axes
from ThermoCalc import ThermoCalcWrapper


class EmissionDatabaseRangeTests(unittest.TestCase):
    def test_accident_region_extends_collector_temperature_to_1500_k(self):
        axes = region_axes(REGIONS["accident"])

        np.testing.assert_allclose(axes["TC_axis"], np.arange(500.0, 1501.0, 10.0))
        self.assertEqual(points_for_axes(axes), 9_693_576)

    def test_full_database_unique_point_count_after_tc_extension(self):
        total = sum(points_for_axes(region_axes(spec)) for spec in REGIONS.values())

        self.assertEqual(total, 22_576_428)

    def test_pyd_info_honors_explicit_solver_directory(self):
        pyd_path = THERMOCALC_DIR / "te_solver.cp312-win_amd64.pyd"
        with (
            mock.patch.dict(os.environ, {"THERMOCALC_TE_SOLVER_DIR": str(THERMOCALC_DIR)}),
            mock.patch.object(emission_database, "ROOT", REPO_ROOT / "missing"),
        ):
            info = emission_database.pyd_info()

        self.assertEqual(info["path"], str(pyd_path.relative_to(REPO_ROOT)))
        self.assertEqual(info["sha256"], emission_database.file_sha256(pyd_path))

    def test_raw_chunk_iterator_excludes_optimized_sidecars(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            chunk_dir = Path(temp_dir) / "chunks" / "core"
            chunk_dir.mkdir(parents=True)
            (chunk_dir / "core_0000.npz").touch()
            (chunk_dir / "core_0000.optimized.npz").touch()

            paths = list(emission_database.iter_raw_chunk_paths(Path(temp_dir)))

        self.assertEqual([path.name for path in paths], ["core_0000.npz"])

    def test_default_runtime_database_prefers_tc1500_extension(self):
        def manifest_exists(path):
            return path.name == "runtime_dense_manifest.json" and path.parent.name == "pcs_0p02_5torr_tc1500"

        with mock.patch.object(Path, "exists", autospec=True, side_effect=manifest_exists):
            selected = ThermoCalcWrapper._find_default_lookup_database()

        self.assertIsNotNone(selected)
        self.assertEqual(selected.name, "pcs_0p02_5torr_tc1500")


if __name__ == "__main__":
    unittest.main()
