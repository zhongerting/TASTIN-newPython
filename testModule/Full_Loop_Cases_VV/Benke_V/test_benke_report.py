from __future__ import annotations

import contextlib
import csv
import io
import json
import shutil
import sys
import unittest
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class BenkeReportTests(unittest.TestCase):
    def test_generate_markdown_report_from_run_summary(self) -> None:
        from benke_report import generate_markdown_report
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig
        from run_benke_thermal_validation import write_outputs

        run_dir = CASE_DIR / "_tmp_report_test"
        if run_dir.exists():
            shutil.rmtree(run_dir)
        (run_dir / "results").mkdir(parents=True)
        try:
            summary = write_outputs(run_dir, BenkeThermalNetworkConfig(), BENKE_TYPICAL_CASE, experimental_data_dir=run_dir)
            report_path = generate_markdown_report(run_dir / "run_summary.json")
            text = report_path.read_text(encoding="utf-8")
        finally:
            if run_dir.exists():
                shutil.rmtree(run_dir)

        self.assertEqual(summary["validation"]["status"], "partial_missing_digitized_data")
        self.assertIn("# Benke 热工水力验证结果报告", text)
        self.assertIn("partial_missing_digitized_data", text)
        self.assertIn("文献范围校核", text)
        self.assertIn("coolant heat fraction", text)

    def test_runner_cli_accepts_coolant_heat_fraction(self) -> None:
        from run_benke_thermal_validation import main

        output_root = CASE_DIR / "_tmp_runner_cli_fraction"
        if output_root.exists():
            shutil.rmtree(output_root)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                rc = main(
                    [
                        "--run-id",
                        "fraction_case",
                        "--output-root",
                        str(output_root),
                        "--coolant-heat-fraction",
                        "0.94",
                    ]
                )
            summary_path = output_root / "fraction_case" / "run_summary.json"
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        finally:
            if output_root.exists():
                shutil.rmtree(output_root)

        self.assertEqual(rc, 0)
        self.assertAlmostEqual(summary["config"]["coolant_heat_fraction"], 0.94)



    def test_runner_outputs_benke_thermocouple_metadata_columns(self) -> None:
        from benke_thermal_network import BENKE_TYPICAL_CASE, BenkeThermalNetworkConfig
        from run_benke_thermal_validation import write_outputs

        run_dir = CASE_DIR / "_tmp_tc_metadata_test"
        if run_dir.exists():
            shutil.rmtree(run_dir)
        (run_dir / "results").mkdir(parents=True)
        try:
            write_outputs(run_dir, BenkeThermalNetworkConfig(active_length_m=0.410), BENKE_TYPICAL_CASE, experimental_data_dir=run_dir)
            with (run_dir / "results" / "sleeve_thermocouple_12pt.csv").open("r", encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
        finally:
            if run_dir.exists():
                shutil.rmtree(run_dir)

        self.assertIn("z_mm", rows[0])
        self.assertIn("measurement_radius_m", rows[0])
        self.assertIn("included_in_benke_average", rows[0])
        self.assertEqual(rows[0]["z_mm"], "-205")
        self.assertEqual(rows[8]["z_mm"], "")
        self.assertEqual(rows[0]["included_in_benke_average"], "False")
        self.assertEqual(rows[1]["included_in_benke_average"], "True")
if __name__ == "__main__":
    unittest.main()
