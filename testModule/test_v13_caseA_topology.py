import os
import sys
import unittest

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from test_core_assemble_v13_caseA import (
    V13_CASE_VERSION,
    V13_DEFAULT_PUMP_TOTAL_HEAD_PA,
    build_v13_case_a_system,
    v13_basic_diagnostics,
)


class V13CaseATopologyTest(unittest.TestCase):
    def test_closed_loop_pump_and_pipefin_radiator_path(self):
        build = build_v13_case_a_system(
            pipe_n_nodes=2,
            n_tubes=8,
            n_axial=2,
            n_fin_width=4,
            total_inlet_flow_kg_s=1.3,
        )
        diag = v13_basic_diagnostics(build)
        net = build["system"].fluid_solver
        names = {getattr(vol, "name", "") for vol in net.volumes_obj}
        junction_names = {getattr(junc, "name", "") for junc in net.junctions_obj}

        self.assertEqual(build["case_version"], V13_CASE_VERSION)
        self.assertEqual(diag["fixed_pressure_boundary_count"], 0.0)
        self.assertEqual(diag["pressure_reference_count"], 1.0)
        self.assertIn("V12_CoreInletConnector", diag["pressure_reference"])
        self.assertNotIn("V12_InletBoundary_FixedFlow", names)
        self.assertNotIn("V12_OutletBoundary_FixedPressure", names)
        self.assertIn("V13_PumpMidNode", names)
        self.assertIn("V13_PumpOutletNode", names)
        self.assertNotIn("J_V12_InletBoundary_to_Pipe11", junction_names)
        self.assertNotIn("J_Pipe09_to_OutletBoundary", junction_names)
        self.assertIn("J_Pipe09_to_V13_PumpA", junction_names)
        self.assertIn("J_V13_PumpA_to_PumpB", junction_names)
        self.assertIn("J_V13_PumpOutlet_to_Pipe11", junction_names)
        self.assertEqual(diag["pump_junction_count"], 2.0)
        self.assertAlmostEqual(diag["pump_total_head_pa"], V13_DEFAULT_PUMP_TOTAL_HEAD_PA)
        self.assertEqual(diag["radiator_tube_count"], 8)
        self.assertTrue(diag["tec_coupled_enabled"])
        self.assertFalse(any(name.startswith("V10_") for name in names))


if __name__ == "__main__":
    unittest.main()
