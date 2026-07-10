import json
import math
import unittest
from pathlib import Path

from .run_v15_v71_flow_path_smoke import run_smoke_case


class V15V71FlowPathSmokeTests(unittest.TestCase):
    def test_hydraulic_only_smoke_writes_finite_results(self):
        output_path = Path(__file__).with_name("v15_v71_flow_path_smoke_result.json")
        if output_path.exists():
            output_path.unlink()

        result = run_smoke_case(output_path=output_path)

        self.assertTrue(output_path.exists())
        saved = json.loads(output_path.read_text(encoding="utf-8"))
        self.assertEqual(saved["case"], "V15_V71_flow_path_hydraulic_only_smoke")
        self.assertTrue(saved["hydraulic_init_converged"])
        self.assertTrue(saved["hydraulic_step_completed"])
        self.assertEqual(saved["case_version"], "v15_v71_center0p30_uniform_pipefin_full_loop")
        self.assertEqual(saved["axial_power_profile_name"], "center_0p30m_uniform")
        self.assertAlmostEqual(saved["axial_power_profile_sum"], 1.0, places=12)
        self.assertGreater(saved["axial_power_profile_heated_node_count"], 0)
        self.assertEqual(saved["pressure_reference_volumes"], ["CoreInletConnector"])
        self.assertEqual(saved["fixed_pressure_boundary_volumes"], [])
        self.assertEqual(saved["nonfinite_flow_junctions"], [])
        self.assertEqual(saved["nonfinite_pressure_volumes"], [])
        self.assertEqual(saved["radiator_tube_count"], 78)
        self.assertEqual(saved["cold_return_branch_count"], 3)
        self.assertGreater(saved["max_abs_flow_kg_s_after_step"], 0.0)
        self.assertTrue(math.isfinite(saved["max_abs_flow_kg_s_after_step"]))
        self.assertAlmostEqual(
            saved["pump_total_head_pa"],
            result["pump_total_head_pa"],
            places=9,
        )


if __name__ == "__main__":
    unittest.main()
