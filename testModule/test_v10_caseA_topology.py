import os
import sys
import unittest

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Solvers.Hydrodynamics.BoundaryVolume import InletJunction
from Solvers.Hydrodynamics.Components import MacroFlowJunction
from test_core_assemble_v10_caseA import (
    V10_CASE_VERSION,
    build_v10_case_a_system,
    v10_flow_diagnostics,
)


class TestV10CaseATopology(unittest.TestCase):
    def setUp(self):
        self.build = build_v10_case_a_system(pipe_n_nodes=1, external_pipe_n_nodes=1)
        self.net = self.build["system"].fluid_solver

    def test_case_identity_and_open_boundaries(self):
        self.assertEqual(self.build["case_version"], V10_CASE_VERSION)
        fixed_pressure = [
            vol.name for vol in self.net.volumes_obj
            if bool(getattr(vol, "is_pressure_boundary", False))
        ]
        self.assertEqual(fixed_pressure, ["V10_OutletBoundary_FixedPressure"])
        self.assertIsInstance(self.build["j_inlet"], InletJunction)
        self.assertAlmostEqual(self.build["j_inlet"].target_W, 1.3)
        self.assertFalse(any(bool(getattr(j, "is_pump_junction", False)) for j in self.net.junctions_obj))

    def test_no_v9_temporary_hot_outlet_merge(self):
        volume_names = {getattr(vol, "name", "") for vol in self.net.volumes_obj}
        junction_names = {getattr(junc, "name", "") for junc in self.net.junctions_obj}
        self.assertNotIn("V9_HotOutletMerge", volume_names)
        self.assertNotIn("J_HotOutletMerge_to_OutletBoundary", junction_names)
        self.assertNotIn("V9_OutletBoundary_FixedPressure", volume_names)

    def test_hot_outlets_are_explicit_and_ring_interface_is_scaled(self):
        self.assertEqual(len(self.build["hot_outlet_branches"]), 3)
        self.assertEqual(len(self.build["hot_outlet_to_ring_junctions"]), 3)
        for junc in self.build["hot_outlet_to_ring_junctions"]:
            self.assertIsInstance(junc, MacroFlowJunction)
            self.assertEqual(junc.multiplier, 2.0)
        self.assertEqual(set(self.build["inlet_mix_nodes"]), {"I1", "I2", "I3"})
        self.assertEqual(set(self.build["outlet_mix_nodes"]), {"O1", "O2", "O3"})

    def test_collector_ring_and_cold_return_geometry_exist_once(self):
        self.assertEqual(len(self.build["ring"]["sectors"]), 6)
        self.assertEqual(len(self.build["manifolds"]), 3)
        self.assertAlmostEqual(self.build["manifolds"][0].total_length, 0.40911)
        self.assertAlmostEqual(self.build["manifolds"][1].total_length, 1.41912)
        self.assertAlmostEqual(self.build["manifolds"][2].total_length, 1.41912)
        self.assertIsInstance(self.build["j_cold_23_in"], MacroFlowJunction)
        self.assertIsInstance(self.build["j_cold_23_out"], MacroFlowJunction)
        self.assertEqual(self.build["j_cold_23_in"].multiplier, 2.0)
        self.assertEqual(self.build["j_cold_23_out"].multiplier, 2.0)

    def test_initial_flow_diagnostics_close_at_design_values(self):
        diag = v10_flow_diagnostics(self.build)
        self.assertEqual(diag["fixed_pressure_boundary_count"], 1.0)
        self.assertEqual(diag["fixed_pressure_boundaries"], ["V10_OutletBoundary_FixedPressure"])
        self.assertFalse(diag["has_hot_outlet_merge"])
        self.assertAlmostEqual(diag["inlet_total_flow_kg_s"], 1.3)
        self.assertAlmostEqual(sum(diag["hot_outlet_macro_flows_kg_s"]), 1.3)
        self.assertAlmostEqual(diag["single_ring_in_total_kg_s"], 0.65)
        self.assertAlmostEqual(diag["single_ring_out_total_kg_s"], 0.65)


if __name__ == "__main__":
    unittest.main()
