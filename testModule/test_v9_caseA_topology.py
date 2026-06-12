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
from test_core_assemble_v9_caseA import (
    V9_CASE_VERSION,
    build_v9_case_a_system,
    v9_flow_diagnostics,
)


class TestV9CaseATopology(unittest.TestCase):
    def setUp(self):
        self.build = build_v9_case_a_system(pipe_n_nodes=1, external_pipe_n_nodes=1)
        self.net = self.build["system"].fluid_solver

    def test_case_identity_and_v8_core_are_preserved(self):
        self.assertEqual(self.build["case_version"], V9_CASE_VERSION)
        self.assertEqual(
            list(self.build["tfes"]),
            ["Center", "Ring1", "Ring2", "Ring3_TEC", "Ring3_Open"],
        )
        self.assertEqual(sum(self.build["ring_multipliers"].values()), 37)
        self.assertEqual(sum(self.build["tec_ring_multipliers"].values()), 34)

    def test_open_loop_boundary_conditions_are_not_overconstrained(self):
        fixed_pressure = [
            vol for vol in self.net.volumes_obj
            if bool(getattr(vol, "is_pressure_boundary", False))
        ]
        self.assertEqual([vol.name for vol in fixed_pressure], ["V9_OutletBoundary_FixedPressure"])
        self.assertIsInstance(self.build["j_inlet"], InletJunction)
        self.assertAlmostEqual(self.build["j_inlet"].target_W, 1.3)
        self.assertFalse(any(bool(getattr(j, "is_pump_junction", False)) for j in self.net.junctions_obj))

    def test_asymmetric_radiator_outlet_branches_are_not_merged(self):
        b38 = self.build["radiator_outlet_branch_38"]
        b4450 = self.build["radiator_outlet_branch_44_50_rep"]
        self.assertAlmostEqual(b38.total_length, 0.40911)
        self.assertAlmostEqual(b4450.total_length, 1.41912)
        self.assertIsInstance(self.build["j_rad_44_50_in"], MacroFlowJunction)
        self.assertEqual(self.build["j_rad_44_50_in"].multiplier, 2.0)

    def test_cold_return_and_hot_outlet_simplifications_match_design(self):
        self.assertIsInstance(self.build["j_cold_23_in"], MacroFlowJunction)
        self.assertEqual(self.build["j_cold_23_in"].multiplier, 2.0)
        self.assertEqual(len(self.build["hot_outlet_branches"]), 3)
        self.assertTrue(all(j.multiplier == 2.0 for j in [self.build["j_cold_23_in"]]))
        self.assertTrue(all(type(j).__name__ == "FlowJunction" for j in self.build["hot_outlet_entry_junctions"]))

    def test_initial_flow_diagnostics_close_at_design_values(self):
        diag = v9_flow_diagnostics(self.build)
        self.assertEqual(diag["fixed_pressure_boundary_count"], 1.0)
        self.assertFalse(diag["has_pump_junction"])
        self.assertAlmostEqual(diag["inlet_total_flow_kg_s"], 1.3)
        self.assertAlmostEqual(diag["radiator_branch_38_flow_kg_s"], 1.3 / 3.0)
        self.assertAlmostEqual(diag["radiator_branch_44_50_macro_flow_kg_s"], 2.0 * 1.3 / 3.0)
        self.assertAlmostEqual(diag["cold_return_branch_2_3_macro_flow_kg_s"], 2.0 * 1.3 / 3.0)
        self.assertAlmostEqual(diag["tfe_total_macro_flow_kg_s"], 1.3)


if __name__ == "__main__":
    unittest.main()
