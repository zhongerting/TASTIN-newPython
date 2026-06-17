import os
import sys
import unittest

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from test_core_assemble_v12_caseA import (
    V12_CASE_VERSION,
    build_v12_case_a_system,
    v12_basic_diagnostics,
)


class V12CaseATopologyTest(unittest.TestCase):
    def test_open_loop_boundaries_and_radiator_path(self):
        build = build_v12_case_a_system(
            pipe_n_nodes=2,
            n_tubes=8,
            n_axial=2,
            n_fin_width=4,
            total_inlet_flow_kg_s=1.3,
        )
        diag = v12_basic_diagnostics(build)
        net = build["system"].fluid_solver
        names = {getattr(vol, "name", "") for vol in net.volumes_obj}
        junction_names = {getattr(junc, "name", "") for junc in net.junctions_obj}

        self.assertEqual(build["case_version"], V12_CASE_VERSION)
        self.assertEqual(diag["fixed_pressure_boundary_count"], 1.0)
        self.assertIn("V12_InletBoundary_FixedFlow", names)
        self.assertIn("V12_OutletBoundary_FixedPressure", names)
        self.assertIn("V12_CoreInletConnector", names)
        self.assertIn("V12_CoreOutletConnector", names)
        self.assertIn("V12_RadiatorInletSplit", names)
        self.assertIn("V12_RadiatorOutletMix", names)
        self.assertIn("J_V12_InletBoundary_to_Pipe11", junction_names)
        self.assertIn("J_CoreOutletConnector_to_Pipe05", junction_names)
        self.assertIn("J_Pipe05_to_RadiatorInletSplit", junction_names)
        self.assertIn("J_Pipe09_to_OutletBoundary", junction_names)
        self.assertEqual(diag["radiator_tube_count"], 8)
        self.assertFalse(diag["tec_coupled_enabled"])
        self.assertFalse(any(name.startswith("V10_") for name in names))


if __name__ == "__main__":
    unittest.main()
