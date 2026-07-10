import ast
import os
import sys
import unittest


current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)


class FullLoopCommonBuilderTest(unittest.TestCase):
    def test_common_builder_constructs_closed_loop_with_configurable_tec_and_pumps(self):
        from testModule.Full_Loop_Cases import (
            FullLoopCoreConfig,
            FullLoopFlowConfig,
            FullLoopPumpConfig,
            build_full_loop_common_base,
            full_loop_common_diagnostics,
        )

        build = build_full_loop_common_base(
            core_config=FullLoopCoreConfig(
                inlet_temperature_k=727.0,
                reference_pressure_pa=205000.0,
                tec_ring_multipliers=(1, 6, 12, 15, 3),
                main_tec_enabled=False,
            ),
            flow_config=FullLoopFlowConfig(
                total_flow_kg_s=1.3,
                radiator_header_area_m2=4.2e-4,
                radiator_header_dh_m=0.021,
                radiator_inlet_header_length_m=0.31,
                radiator_outlet_header_length_m=0.27,
                radiator_header_n_nodes=3,
            ),
            pump_config=FullLoopPumpConfig(
                pump_total_head_pa=8123.0,
                pump_flow_control=False,
            ),
            close_with_placeholder_bridge=True,
        )

        diag = full_loop_common_diagnostics(build)
        net = build["system"].fluid_solver
        volume_names = {getattr(vol, "name", "") for vol in net.volumes_obj}
        junction_names = {getattr(junc, "name", "") for junc in net.junctions_obj}

        self.assertEqual(diag["case_version"], "full_loop_common_base")
        self.assertEqual(diag["fixed_pressure_boundary_count"], 0)
        self.assertEqual(diag["pressure_reference"], ["CoreInletConnector"])
        self.assertEqual(build["tec_ring_multipliers"]["Ring3_Open"], 3)
        self.assertEqual(set(build["fluid_channels"]), {"Center", "Ring1", "Ring2", "Ring3_TEC", "Ring3_Open"})
        self.assertIn("RadiatorInletHeader_Vol_01", volume_names)
        self.assertIn("RadiatorOutletHeader_Vol_01", volume_names)
        self.assertIn("J_RadiatorPlaceholderBridge", junction_names)
        self.assertIn("CoreInletSegment_Vol_01", volume_names)
        self.assertIn("J_PumpOutletNode_to_CoreInletSegment", junction_names)
        self.assertIn("J_CoreInletSegment_to_CoreInletConnector", junction_names)
        self.assertNotIn("J_PumpOutletNode_to_CoreInletConnector", junction_names)
        self.assertIs(build["j_pump_outlet_to_core_inlet_segment"].from_vol, build["pump_outlet_node"])
        self.assertIs(build["j_core_inlet_segment_to_core_inlet"].to_vol, build["core_inlet_connector"])
        self.assertAlmostEqual(build["radiator_inlet_header"].area, 4.2e-4)
        self.assertAlmostEqual(build["radiator_inlet_header"].d_h, 0.021)
        self.assertAlmostEqual(build["radiator_inlet_header"].total_length, 0.31)
        self.assertEqual(build["radiator_inlet_header"].n_nodes, 3)
        self.assertAlmostEqual(diag["pump_total_head_pa"], 8123.0)
        self.assertAlmostEqual(build["pump_a"].delta_p, 8123.0 / 2.0)
        self.assertAlmostEqual(build["pump_b"].delta_p, 8123.0 / 2.0)

    def test_common_package_does_not_import_legacy_case_builders(self):
        package_dir = os.path.join(root_dir, "testModule", "Full_Loop_Cases")
        banned = {
            "test_core_assemble_v7_caseA",
            "test_core_assemble_v8_caseA",
            "test_core_assemble_v11_caseA",
            "test_core_assemble_v12_caseA",
            "test_core_assemble_v13_caseA",
        }
        offenders = []
        for filename in os.listdir(package_dir):
            if not filename.endswith(".py"):
                continue
            path = os.path.join(package_dir, filename)
            with open(path, "r", encoding="utf-8") as f:
                tree = ast.parse(f.read(), filename=path)
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        root = alias.name.split(".")[0]
                        if root in banned:
                            offenders.append((filename, alias.name))
                elif isinstance(node, ast.ImportFrom):
                    module = (node.module or "").split(".")[0]
                    if module in banned:
                        offenders.append((filename, node.module))

        self.assertEqual(offenders, [])

    def test_flow_controlled_pumps_receive_target_flow(self):
        from testModule.Full_Loop_Cases import (
            FullLoopCoreConfig,
            FullLoopFlowConfig,
            FullLoopPumpConfig,
            build_full_loop_common_base,
        )

        build = build_full_loop_common_base(
            core_config=FullLoopCoreConfig(main_tec_enabled=False),
            flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
            pump_config=FullLoopPumpConfig(
                pump_total_head_pa=9000.0,
                pump_flow_control=True,
                target_flow_kg_s=1.15,
            ),
            close_with_placeholder_bridge=True,
        )

        self.assertTrue(hasattr(build["pump_a"], "target_W"))
        self.assertTrue(hasattr(build["pump_b"], "target_W"))
        self.assertAlmostEqual(build["pump_a"].target_W, 1.15)
        self.assertAlmostEqual(build["pump_b"].target_W, 1.15)

if __name__ == "__main__":
    unittest.main()

