import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from testModule.Full_Loop_Cases_10kW import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    build_v14_case_a_system,
)


class FullLoop10kWCoreGeometryTests(unittest.TestCase):
    def test_core_defaults_describe_five_ring_series_tec_layout(self):
        cfg = FullLoopCoreConfig(main_tec_enabled=False)

        self.assertEqual(tuple(cfg.representative_names), ("Center", "Ring1", "Ring2", "Ring3", "Ring4"))
        self.assertEqual(tuple(cfg.ring_multipliers), (1, 6, 9, 18, 24))
        self.assertEqual(tuple(cfg.tec_ring_multipliers), (1, 6, 9, 18, 24))
        self.assertEqual(sum(cfg.ring_multipliers), 58)
        self.assertEqual(sum(cfg.tec_ring_multipliers), 58)
        self.assertEqual(cfg.physical_ring_count, 5)
        self.assertAlmostEqual(cfg.main_tec_target_value, 50.5)
        self.assertEqual(
            cfg.representative_ring_mapping,
            {"Center": 0, "Ring1": 1, "Ring2": 2, "Ring3": 3, "Ring4": 4},
        )

    def test_v14_10kw_tec_lookup_config_reaches_thermocalc(self):
        created = []

        class FakeThermoCalcModel:
            def __init__(self, n_elements, n_nodes, lookup_db=None, enable_lookup=None, lookup_regions=None):
                created.append({
                    "n_elements": n_elements,
                    "n_nodes": n_nodes,
                    "lookup_db": lookup_db,
                    "enable_lookup": enable_lookup,
                    "lookup_regions": lookup_regions,
                })
                self._input_data = SimpleNamespace()

            def setup_circuit_mode(self, mode_str, target_value, I_guess=150.0):
                self.circuit_mode = (mode_str, target_value, I_guess)

        with patch("Components.ReactorCore.ThermoCalcModel", FakeThermoCalcModel):
            build_v14_case_a_system(
                core_config=FullLoopCoreConfig(
                    tec_lookup_enabled=True,
                    tec_lookup_db="lookup-db",
                    tec_lookup_regions=("hot", "cold"),
                ),
                flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
                pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
            )

        self.assertTrue(created)
        self.assertEqual(created[0]["lookup_db"], "lookup-db")
        self.assertIs(created[0]["enable_lookup"], True)
        self.assertEqual(created[0]["lookup_regions"], ("hot", "cold"))
    def test_v14_10kw_core_geometry_builds_hydraulic_only(self):
        build = build_v14_case_a_system(
            core_config=FullLoopCoreConfig(main_tec_enabled=False),
            flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
            pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
        )
        core = build["core"]

        self.assertEqual(set(build["fluid_channels"]), {"Center", "Ring1", "Ring2", "Ring3", "Ring4"})
        self.assertEqual(sum(build["ring_multipliers"].values()), 58)
        self.assertEqual(sum(build["tec_ring_multipliers"].values()), 58)
        self.assertEqual(len(core.mod_rings), 5)

        expected_edges = np.array([0.0, 21.0e-3, 53.25e-3, 86.75e-3, 120.5e-3, 164.0e-3])
        for ring, r_in, r_out in zip(core.mod_rings, expected_edges[:-1], expected_edges[1:]):
            self.assertAlmostEqual(float(ring.mesh.x_faces[0]), float(r_in))
            self.assertAlmostEqual(float(ring.mesh.x_faces[-1]), float(r_out))

        self.assertAlmostEqual(float(core.barrel.mesh.x_faces[0]), 164.0e-3)
        self.assertAlmostEqual(float(core.barrel.mesh.x_faces[-1]), 166.0e-3)
        self.assertAlmostEqual(float(core.reflector.mesh.x_faces[0]), 166.0e-3)
        self.assertAlmostEqual(float(core.reflector.mesh.x_faces[-1]), 261.0e-3)

        net = build["system"].fluid_solver
        self.assertTrue(net.initialize_hydraulics(dt=0.01, tol=1.0e-4, max_iter=1000))
        net.step_hydraulic(1.0e-4)
        self.assertTrue(math.isfinite(float(np.max(np.abs(net.W_vec)))))


if __name__ == "__main__":
    unittest.main()
