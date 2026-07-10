import json
import math
import unittest
from pathlib import Path

from .run_v14_hydraulic_stability import run_stability_case


class V14HydraulicStabilityTests(unittest.TestCase):
    def test_hydraulic_stability_writes_bounded_flow_drift(self):
        output_path = Path(__file__).with_name("v14_hydraulic_stability_result.json")
        if output_path.exists():
            output_path.unlink()

        result = run_stability_case(output_path=output_path, n_steps=80, last_window=20)

        self.assertTrue(output_path.exists())
        saved = json.loads(output_path.read_text(encoding="utf-8"))
        self.assertEqual(saved["case"], "V14_10kW_hydraulic_stability")
        self.assertTrue(saved["hydraulic_init_converged"])
        self.assertTrue(saved["stable_flow_reached"])
        self.assertEqual(saved["nonfinite_flow_junctions"], [])
        self.assertEqual(saved["nonfinite_pressure_volumes"], [])
        self.assertLess(saved["last_window_max_abs_pump_flow_change_kg_s"], 1.0e-6)
        self.assertTrue(math.isfinite(saved["final_pump_a_flow_kg_s"]))
        self.assertAlmostEqual(saved["final_pump_a_flow_kg_s"], result["final_pump_a_flow_kg_s"], places=12)


if __name__ == "__main__":
    unittest.main()
