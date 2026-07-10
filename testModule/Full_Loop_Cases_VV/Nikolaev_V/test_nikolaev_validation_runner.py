from __future__ import annotations

import json
import shutil
import sys
import unittest
from pathlib import Path


CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class NikolaevValidationRunnerTests(unittest.TestCase):
    def test_build_summary_reports_table_metrics(self) -> None:
        from nikolaev_single_tfe_model import NikolaevModelConfig
        from run_nikolaev_validation import build_summary

        summary = build_summary(NikolaevModelConfig())

        self.assertEqual(summary["status"], "complete_table_validation")
        self.assertEqual(len(summary["table2"]), 3)
        self.assertLess(summary["metrics"]["table2_current_mae_a"], 0.5)
        self.assertLess(summary["metrics"]["table2_emitter_temp_mae_k"], 1.0e-9)
        self.assertLess(summary["metrics"]["table3_max_abs_error_k"], 1.0e-9)
        self.assertLess(summary["metrics"]["table4_max_abs_error_mm"], 1.0e-12)

    def test_runner_writes_summary_csvs_and_markdown_report(self) -> None:
        from run_nikolaev_validation import main

        output_root = CASE_DIR / "_tmp_runner_output"
        if output_root.exists():
            shutil.rmtree(output_root)
        try:
            rc = main(["--run-id", "smoke", "--output-root", str(output_root)])
            run_dir = output_root / "smoke"
            csv_exists = (run_dir / "results" / "table2_operating_points.csv").exists()
            summary = json.loads((run_dir / "run_summary.json").read_text(encoding="utf-8"))
            report = (run_dir / "validation_report.md").read_text(encoding="utf-8")
        finally:
            if output_root.exists():
                shutil.rmtree(output_root)

        self.assertEqual(rc, 0)
        self.assertTrue(csv_exists)
        self.assertIn("Nikolaev 1995", report)
        self.assertIn("Table 2 Operating Point Comparison", report)
        self.assertEqual(summary["source"]["doi"], "10.1063/1.47120")


if __name__ == "__main__":
    unittest.main()
