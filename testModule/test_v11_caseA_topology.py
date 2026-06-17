import os
import sys
import unittest

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Solvers.Hydrodynamics.Components import MacroFlowJunction, PumpJunction
from test_core_assemble_v11_caseA import (
    V11_CASE_VERSION,
    V11_DEFAULT_PUMP_TOTAL_HEAD_PA,
    build_v11_case_a_system,
    v11_flow_diagnostics,
)


class TestV11CaseATopology(unittest.TestCase):
    def setUp(self):
        self.build = build_v11_case_a_system(pipe_n_nodes=1, external_pipe_n_nodes=1)
        self.net = self.build["system"].fluid_solver

    def test_case_identity_and_closed_pressure_reference(self):
        self.assertEqual(self.build["case_version"], V11_CASE_VERSION)
        fixed_pressure = [
            vol.name for vol in self.net.volumes_obj
            if bool(getattr(vol, "is_pressure_boundary", False))
        ]
        pressure_reference = [
            vol.name for vol in self.net.volumes_obj
            if bool(getattr(vol, "is_pressure_reference", False))
        ]
        self.assertEqual(fixed_pressure, [])
        self.assertEqual(pressure_reference, ["CoreInletConnector"])

    def test_open_boundaries_are_removed(self):
        volume_names = {getattr(vol, "name", "") for vol in self.net.volumes_obj}
        junction_names = {getattr(junc, "name", "") for junc in self.net.junctions_obj}
        self.assertNotIn("V10_InletBoundary_FixedFlow", volume_names)
        self.assertNotIn("V10_OutletBoundary_FixedPressure", volume_names)
        self.assertNotIn("J_V10_InletBoundary_to_CoreInletConnector", junction_names)
        self.assertNotIn("J_ColdReturnOutletMerge_to_OutletBoundary", junction_names)

    def test_two_equal_series_pumps_between_outer_header_and_distributor(self):
        self.assertIsInstance(self.build["pump_a"], PumpJunction)
        self.assertIsInstance(self.build["pump_b"], PumpJunction)
        self.assertIs(self.build["pump_a"].from_vol, self.build["radiator_outer_header_52"].volumes[-1])
        self.assertIs(self.build["pump_a"].to_vol, self.build["pump_mid_node"])
        self.assertIs(self.build["pump_b"].from_vol, self.build["pump_mid_node"])
        self.assertIs(self.build["pump_b"].to_vol, self.build["pump_outlet_distributor"])
        self.assertAlmostEqual(self.build["pump_a"].delta_p, V11_DEFAULT_PUMP_TOTAL_HEAD_PA / 2.0)
        self.assertAlmostEqual(self.build["pump_b"].delta_p, V11_DEFAULT_PUMP_TOTAL_HEAD_PA / 2.0)

    def test_pump_outlet_distributor_feeds_cold_return_branches(self):
        self.assertEqual(self.build["pump_outlet_distributor"].name, "V11_PumpOutletDistributor_51")
        self.assertIs(self.build["j_cold_1_in"].from_vol, self.build["pump_outlet_distributor"])
        self.assertIs(self.build["j_cold_23_in"].from_vol, self.build["pump_outlet_distributor"])
        self.assertIsInstance(self.build["j_cold_23_in"], MacroFlowJunction)
        self.assertEqual(self.build["j_cold_23_in"].multiplier, 2.0)
        self.assertIs(self.build["j_cold_merge_to_core_inlet"].to_vol, self.build["core_inlet_connector"])

    def test_collector_ring_interface_is_unchanged_from_v10(self):
        self.assertEqual(len(self.build["hot_outlet_branches"]), 3)
        self.assertEqual(len(self.build["hot_outlet_to_ring_junctions"]), 3)
        for junc in self.build["hot_outlet_to_ring_junctions"]:
            self.assertIsInstance(junc, MacroFlowJunction)
            self.assertEqual(junc.multiplier, 2.0)
        self.assertEqual(len(self.build["ring"]["sectors"]), 6)
        self.assertEqual(len(self.build["manifolds"]), 3)

    def test_initial_flow_diagnostics_close_at_design_values(self):
        diag = v11_flow_diagnostics(self.build)
        self.assertEqual(diag["fixed_pressure_boundary_count"], 0.0)
        self.assertEqual(diag["fixed_pressure_boundaries"], [])
        self.assertEqual(diag["pressure_reference_count"], 1.0)
        self.assertEqual(diag["pressure_reference"], ["CoreInletConnector"])
        self.assertTrue(diag["has_pump_junction"])
        self.assertEqual(diag["pump_junction_count"], 2.0)
        self.assertAlmostEqual(diag["pump_total_head_pa"], V11_DEFAULT_PUMP_TOTAL_HEAD_PA)
        self.assertAlmostEqual(diag["pump_mean_flow_kg_s"], 1.3)
        self.assertAlmostEqual(diag["closed_loop_flow_kg_s"], 1.3)
        self.assertAlmostEqual(sum(diag["hot_outlet_macro_flows_kg_s"]), 1.3)
        self.assertAlmostEqual(diag["single_ring_in_total_kg_s"], 0.65)
        self.assertAlmostEqual(diag["single_ring_out_total_kg_s"], 0.65)


if __name__ == "__main__":
    unittest.main()
