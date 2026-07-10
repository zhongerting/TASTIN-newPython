from __future__ import annotations

import sys
import unittest
from pathlib import Path

CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class BenkeParameterScanTests(unittest.TestCase):
    def test_parameter_envelope_scan_returns_all_grid_points(self) -> None:
        from benke_parameter_scan import scan_benke_parameter_envelope

        rows, summary = scan_benke_parameter_envelope(
            regulated_he_k_values=[0.073, 0.08, 0.087],
            water_h_values=[528.0, 800.0, 1012.0],
        )

        self.assertEqual(len(rows), 9)
        self.assertEqual(summary["grid_point_count"], 9)
        self.assertLess(summary["sleeve_outer_mean_k_min"], summary["sleeve_outer_mean_k_max"])
        self.assertLess(summary["collector_inner_mean_k_min"], summary["collector_inner_mean_k_max"])
        self.assertTrue(all(row["range_check_status"] == "passed" for row in rows))

    def test_lower_he_k_or_lower_water_h_increases_temperatures(self) -> None:
        from benke_parameter_scan import scan_benke_parameter_envelope

        rows, _summary = scan_benke_parameter_envelope(
            regulated_he_k_values=[0.073, 0.087],
            water_h_values=[528.0, 1012.0],
        )
        by_key = {(row["regulated_he_effective_k_w_m_k"], row["water_h_w_m2_k"]): row for row in rows}

        coldest = by_key[(0.087, 1012.0)]
        hottest = by_key[(0.073, 528.0)]
        self.assertGreater(hottest["sleeve_outer_mean_k"], coldest["sleeve_outer_mean_k"])
        self.assertGreater(hottest["collector_inner_mean_k"], coldest["collector_inner_mean_k"])


if __name__ == "__main__":
    unittest.main()
