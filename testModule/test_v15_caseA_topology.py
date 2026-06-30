import ast
import unittest
from pathlib import Path

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V15PipeFinRadiatorConfig,
    build_v15_case_a_system,
    v15_basic_diagnostics,
)


BANNED_IMPORTS = {
    "testModule.test_core_assemble_v8_caseA",
    "testModule.test_core_assemble_v9_caseA",
    "testModule.test_core_assemble_v10_caseA",
    "testModule.test_core_assemble_v11_caseA",
    "testModule.test_core_assemble_v12_caseA",
    "testModule.test_core_assemble_v13_caseA",
    "CoolantLoop.model_collector_ring_6segment_v9_interface",
}
BANNED_NAME_FRAGMENTS = (
    "V12_",
    "Pipe05",
    "Pipe06",
    "Pipe07",
    "Pipe08",
    "Pipe09",
    "Pipe11",
)


class V15CaseATopologyTests(unittest.TestCase):
    def build_case(self):
        return build_v15_case_a_system(
            core_config=FullLoopCoreConfig(main_tec_enabled=False),
            flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
            pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
            radiator_config=V15PipeFinRadiatorConfig(),
        )

    def test_v15_builds_full_pipefin_radiator_topology(self):
        build = self.build_case()
        self.assertEqual(build["case_version"], "v15_pipefin_radiator_full_loop")

        net = build["system"].fluid_solver
        volume_names = {vol.name for vol in net.volumes_obj}
        junction_names = {junc.name for junc in net.junctions_obj}
        component_names = {component.name for component in build["system"].components}

        self.assertEqual(len(build["radiator_upper_headers"]), 78)
        self.assertEqual(len(build["radiator_lower_headers"]), 78)
        self.assertEqual(len(build["radiator_tube_channels"]), 78)
        self.assertEqual(len(build["radiator_units"]), 78)
        self.assertEqual(len(build["cold_return_branches"]), 3)

        for idx in (1, 40, 78):
            self.assertIn(f"RadiatorUpperHeader_{idx:02d}_Vol_01", volume_names)
            self.assertIn(f"RadiatorLowerHeader_{idx:02d}_Vol_01", volume_names)
            self.assertIn(f"RadiatorTubeFluid_{idx:02d}_Vol_01", volume_names)
            self.assertIn(f"RadiatorTube_{idx:02d}", component_names)
            self.assertIn(f"J_RadiatorUpperRing_{idx:02d}_to_{(idx % 78) + 1:02d}", junction_names)
            self.assertIn(f"J_RadiatorLowerRing_{idx:02d}_to_{(idx % 78) + 1:02d}", junction_names)
            self.assertIn(f"J_RadiatorUpper_to_Tube_{idx:02d}", junction_names)
            self.assertIn(f"J_RadiatorTube_{idx:02d}_to_Lower", junction_names)

        for name in (
            "RadiatorInletDistributor",
            "RadiatorInnerHeader",
            "RadiatorOuterHeader_Vol_01",
            "PumpOutletDistributor",
            "ColdReturnBranch_1_Vol_01",
            "ColdReturnBranch_2_Vol_01",
            "ColdReturnBranch_3_Vol_01",
        ):
            self.assertIn(name, volume_names)

        self.assertIn("J_PumpA", junction_names)
        self.assertIn("J_PumpB", junction_names)
        self.assertIn("PumpMidNode", volume_names)
        self.assertIn("PumpOutletNode", volume_names)
        self.assertIn("J_PumpOutletNode_to_PumpOutletDistributor", junction_names)
        self.assertNotIn("J_PumpOutletNode_to_CoreInletSegment", junction_names)
        self.assertNotIn("CoreInletSegment_Vol_01", volume_names)

        pressure_references = [
            vol.name for vol in net.volumes_obj if getattr(vol, "is_pressure_reference", False)
        ]
        pressure_boundaries = [
            vol.name for vol in net.volumes_obj if getattr(vol, "is_pressure_boundary", False)
        ]
        self.assertEqual(pressure_references, ["CoreInletConnector"])
        self.assertEqual(pressure_boundaries, [])

        all_names = volume_names | junction_names | component_names
        bad_names = sorted(
            name for name in all_names for fragment in BANNED_NAME_FRAGMENTS if fragment in name
        )
        self.assertEqual(bad_names, [])

        diagnostics = v15_basic_diagnostics(build)
        self.assertEqual(diagnostics["radiator_tube_count"], 78)
        self.assertEqual(diagnostics["cold_return_branch_count"], 3)
        self.assertAlmostEqual(diagnostics["pump_total_head_pa"], 6466.56)

    def test_full_loop_cases_do_not_import_legacy_case_builders(self):
        root = Path(__file__).resolve().parent / "Full_Loop_Cases"
        violations = []
        for py_file in root.glob("*.py"):
            tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name in BANNED_IMPORTS:
                            violations.append((py_file.name, alias.name))
                elif isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                    if not node.level and module in BANNED_IMPORTS:
                        violations.append((py_file.name, module))
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
