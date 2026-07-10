import argparse
import csv
import tempfile
import unittest
from pathlib import Path

from testModule import run_single_hp_cold_start_800k as run800


class TestSingleHeatPipeColdStartRunScript(unittest.TestCase):
    def test_default_geometry_uses_refined_axial_wall_nodes(self):
        hp = run800.build_single_heat_pipe(
            heater_temperature_k=800.0,
            heater_ha_w_per_k=0.3,
            condenser_emissivity=0.03,
            condenser_background_k=4.0,
            initial_temperature_k=300.0,
            n_eva=15,
            n_con=135,
        )

        self.assertEqual(hp.n_eva, 15)
        self.assertEqual(hp.n_con, 135)
        self.assertEqual(hp.shape_nodes[1], 150)

        wall_columns = run800.outer_wall_temperature_columns(hp)
        self.assertEqual(len(wall_columns), hp.shape_nodes[1])
        self.assertTrue(all(name.startswith("wall_outer_z") for name in wall_columns))

    def test_short_run_writes_outer_wall_axial_temperature_history(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = argparse.Namespace(
                duration_s=0.2,
                dt_s=0.1,
                record_interval_s=0.1,
                heater_temperature_k=800.0,
                heater_ha_w_per_k=0.3,
                condenser_emissivity=0.03,
                condenser_background_k=4.0,
                initial_temperature_k=300.0,
                n_eva=15,
                n_con=135,
                output_dir=tmpdir,
            )

            summary = run800.run_case(args)
            history_csv = Path(summary["history_csv"])
            self.assertTrue(history_csv.exists())

            with history_csv.open(newline="", encoding="utf-8") as f:
                header = next(csv.reader(f))

            wall_columns = [name for name in header if name.startswith("wall_outer_z")]
            self.assertEqual(len(wall_columns), 150)
            self.assertIn("outer_wall_plot_png", summary)
            self.assertTrue(Path(summary["outer_wall_plot_png"]).exists())


if __name__ == "__main__":
    unittest.main()
