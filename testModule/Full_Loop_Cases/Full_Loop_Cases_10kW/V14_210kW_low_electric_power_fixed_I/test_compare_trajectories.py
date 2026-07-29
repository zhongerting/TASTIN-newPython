import csv
import json
import shutil
import unittest
import uuid
from pathlib import Path

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I import (
    compare_trajectories,
)


CASE_DIR = Path(__file__).resolve().parent


class CompareTrajectoriesTest(unittest.TestCase):
    def setUp(self):
        self.run_dir = CASE_DIR / f"_tmp_compare_trajectories_{uuid.uuid4().hex}"
        self.run_dir.mkdir()

    def tearDown(self):
        shutil.rmtree(self.run_dir)

    def test_reports_endpoint_and_smoothness_metrics_using_elapsed_time(self):
        manifest = {
            "trajectory": {
                "hold_before_ramp_s": 20.0,
                "ramp_duration_s": 80.0,
                "final_power_w": 100.0,
                "final_flow_kg_s": 2.0,
            }
        }
        (self.run_dir / "run_manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )

        fields = [
            "elapsed_s",
            "electric_power_ratio",
            "thermal_power_setpoint_W",
            "flow_setpoint_kg_s",
            "core_outlet_T_K",
            "hp_evaporator_temperature_mean_K",
            "hp_condenser_temperature_mean_K",
            "collector_ring_wall_temperature_mean_K",
            "radiator_fin_temperature_mean_K",
        ]
        with (self.run_dir / "history_control.csv").open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for elapsed_s in range(0, 201, 20):
                if elapsed_s <= 100:
                    ratio = 1.0 - 0.0061 * elapsed_s
                    power = 200.0 - elapsed_s
                    flow = 3.0 - 0.01 * elapsed_s
                else:
                    ratio = 0.39 + 0.0001 * (elapsed_s - 100)
                    power = 100.0
                    flow = 2.0
                writer.writerow(
                    {
                        "elapsed_s": elapsed_s,
                        "electric_power_ratio": ratio,
                        "thermal_power_setpoint_W": power,
                        "flow_setpoint_kg_s": flow,
                        "core_outlet_T_K": 850.0 + 0.01 * elapsed_s,
                        "hp_evaporator_temperature_mean_K": 750.0 - 0.02 * elapsed_s,
                        "hp_condenser_temperature_mean_K": 745.0 - 0.01 * elapsed_s,
                        "collector_ring_wall_temperature_mean_K": 787.0 - 0.03 * elapsed_s,
                        "radiator_fin_temperature_mean_K": 731.0 - 0.04 * elapsed_s,
                    }
                )

        result = compare_trajectories.analyze_history(
            self.run_dir / "history_control.csv"
        )

        self.assertEqual(result["freeze_elapsed_s"], 100.0)
        self.assertTrue(result["endpoint_power_compliant"])
        self.assertTrue(result["frozen_power_band_compliant"])
        self.assertTrue(result["setpoints_frozen"])
        self.assertAlmostEqual(result["power_ratio_end"], 0.4)
        self.assertAlmostEqual(result["power_ratio_final_100s_slope_per_s"], 0.0001)
        self.assertAlmostEqual(result["power_ratio_tv_over_net_change"], 0.62 / 0.6)
        self.assertAlmostEqual(result["power_ratio_peak_60s_slope_per_s"], 0.0061)
        self.assertAlmostEqual(result["core_outlet_T_max_K"], 852.0)
        self.assertAlmostEqual(result["core_outlet_T_final_100s_slope_K_per_s"], 0.01)
        self.assertAlmostEqual(result["hp_evaporator_temperature_mean_end_K"], 746.0)
        self.assertAlmostEqual(result["collector_ring_wall_temperature_mean_delta_K"], -6.0)
        self.assertAlmostEqual(result["radiator_fin_temperature_mean_final_100s_slope_K_per_s"], -0.04)

    def test_rejects_history_missing_a_required_column(self):
        (self.run_dir / "run_manifest.json").write_text(
            json.dumps(
                {"trajectory": {"hold_before_ramp_s": 0.0, "ramp_duration_s": 0.0}}
            ),
            encoding="utf-8",
        )
        (self.run_dir / "history_control.csv").write_text(
            "elapsed_s,electric_power_ratio\n0,1\n", encoding="utf-8"
        )

        with self.assertRaisesRegex(ValueError, "missing required columns"):
            compare_trajectories.analyze_history(self.run_dir / "history_control.csv")


if __name__ == "__main__":
    unittest.main()
