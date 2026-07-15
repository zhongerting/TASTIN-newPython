import ast
import math
import unittest
from pathlib import Path

import numpy as np

from Components.ExternalHeatSources import ExternalHeatFluxBC
from Solvers.Hydrodynamics.Components import MacroFlowJunction

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    build_v14_case_a_system,
    v14_basic_diagnostics,
)


BANNED_MODULES = {
    "testModule.test_core_assemble_v7_caseA",
    "testModule.test_core_assemble_v8_caseA",
    "testModule.test_core_assemble_v10_caseA",
    "testModule.test_core_assemble_v11_caseA",
    "testModule.test_core_assemble_v12_caseA",
    "testModule.test_core_assemble_v13_caseA",
    "CoolantLoop.model_collector_ring_6segment_v9_interface",
}


class V14CaseATopologyTests(unittest.TestCase):
    def build_case(self):
        return build_v14_case_a_system(
            core_config=FullLoopCoreConfig(main_tec_enabled=False),
            flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
            pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
            radiator_config=V14HeatPipeRadiatorConfig(
                hot_branch_n_nodes=1,
                manifold_node_counts=(1, 1, 1),
                hp_n_con=2,
                n_fin_height=3,
            ),
        )

    def test_v14_builds_closed_heatpipe_radiator_loop_without_old_case_imports(self):
        build = self.build_case()
        self.assertEqual(build["case_version"], "v14_heatpipe_radiator_full_loop")

        volume_names = {vol.name for vol in build["system"].fluid_solver.volumes_obj}
        junction_names = {junc.name for junc in build["system"].fluid_solver.junctions_obj}

        for name in ("HotOutletBranch_1_Vol_01", "HotOutletBranch_2_Vol_01", "HotOutletBranch_3_Vol_01"):
            self.assertIn(name, volume_names)
        for name in ("InletMix_I1", "InletMix_I2", "InletMix_I3", "OutletMix_O1", "OutletMix_O2", "OutletMix_O3"):
            self.assertIn(name, volume_names)
        for name in ("Manifold_1_Vol_01", "Manifold_2_Vol_01", "Manifold_3_Vol_01"):
            self.assertIn(name, volume_names)
        for name in (
            "A1_I1_to_O1_Channel_Vol_01",
            "A2_O1_to_I2_Channel_Vol_01",
            "A3_I2_to_O2_Channel_Vol_01",
            "A4_O2_to_I3_Channel_Vol_01",
            "A5_I3_to_O3_Channel_Vol_01",
            "A6_O3_to_I1_Channel_Vol_01",
        ):
            self.assertIn(name, volume_names)

        self.assertIn("J_PumpOutletNode_to_CoreInletSegment", junction_names)
        self.assertIn("J_CoreInletSegment_to_CoreInletConnector", junction_names)
        self.assertNotIn("J_PumpOutletNode_to_CoreInletConnector", junction_names)

        pressure_references = [vol.name for vol in build["system"].fluid_solver.volumes_obj if getattr(vol, "is_pressure_reference", False)]
        pressure_boundaries = [vol.name for vol in build["system"].fluid_solver.volumes_obj if getattr(vol, "is_pressure_boundary", False)]
        self.assertEqual(pressure_references, ["CoreInletConnector"])
        self.assertEqual(pressure_boundaries, [])

        self.assertEqual(len(build["hot_outlet_branches"]), 3)
        self.assertEqual(len(build["ring_sectors"]), 6)
        self.assertEqual(len(build["ring_hps"]), 6)
        self.assertEqual(len(build["ring_solids"]), 6)
        self.assertEqual(len(build["manifolds"]), 3)
        self.assertFalse(build["external_heat_enabled"])
        self.assertTrue(all(not np.any(ring_hp._hp_external_heat_enabled) for ring_hp in build["ring_hps"]))
        for ring_hp in build['ring_hps']:
            for hp in ring_hp.hp_units:
                external_bcs = [
                    condition
                    for condition in hp.hp.boundaries['outer_con'].conditions
                    if isinstance(condition, ExternalHeatFluxBC)
                ]
                self.assertEqual(len(external_bcs), 0)

        self.assertEqual(len(build["hot_outlet_to_ring_junctions"]), 3)
        self.assertTrue(all(isinstance(j, MacroFlowJunction) for j in build["hot_outlet_to_ring_junctions"]))
        self.assertTrue(all(math.isclose(j.multiplier, 2.0) for j in build["hot_outlet_to_ring_junctions"]))
        self.assertEqual(len(build["manifold_to_outlet_header_junctions"]), 3)
        self.assertTrue(all(isinstance(j, MacroFlowJunction) for j in build["manifold_to_outlet_header_junctions"]))
        self.assertTrue(all(math.isclose(j.multiplier, 2.0) for j in build["manifold_to_outlet_header_junctions"]))

        diagnostics = v14_basic_diagnostics(build)
        self.assertAlmostEqual(diagnostics["pump_total_head_pa"], 6466.56)
        self.assertAlmostEqual(diagnostics["single_ring_in_total_kg_s"], 0.65)
        self.assertAlmostEqual(diagnostics["single_ring_out_total_kg_s"], 0.65)

    def test_v14_keeps_original_heat_pipe_view_factors_explicitly(self):
        build = self.build_case()
        hp = build["ring_hps"][0].hp_units[0]

        self.assertAlmostEqual(hp.radiation_outer_up_view_factor, 1.0)
        self.assertAlmostEqual(hp.radiation_outer_down_view_factor, 1.0)
        self.assertAlmostEqual(hp.radiation_inner_up_view_factor, 0.0)
        self.assertAlmostEqual(hp.radiation_inner_down_view_factor, 0.3)
        self.assertAlmostEqual(hp.radiation_outer_view_factor, 1.0)
        self.assertAlmostEqual(hp.radiation_inner_view_factor, 0.3)
        self.assertAlmostEqual(hp.radiation_face_average_factor, 0.65)
        self.assertAlmostEqual(hp.effective_emissivity, 0.75 * 0.65)
        self.assertAlmostEqual(hp.effective_fin_emissivity, 0.75 * 0.65)

    def test_v14_external_heat_uses_one_sided_absorption_and_ring_multipliers_once(self):
        config = V14HeatPipeRadiatorConfig(
            external_heat_enabled=True,
            hot_branch_n_nodes=1,
            manifold_node_counts=(1, 1, 1),
            hp_n_con=2,
            n_fin_height=3,
        )
        build = build_v14_case_a_system(
            core_config=FullLoopCoreConfig(main_tec_enabled=False),
            flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
            pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
            radiator_config=config,
        )
        ring_hp = build["ring_hps"][0]
        hp = ring_hp.hp_units[0]
        absorption_factor = config.external_heat_absorption_efficiency * config.hp_emissivity
        expected_wall_area = hp.tube_bare_area * 0.5 * absorption_factor
        expected_fin_area = hp.fin_illuminated_area_array * absorption_factor

        np.testing.assert_allclose(hp.wall_external_absorption_area, expected_wall_area)
        np.testing.assert_allclose(hp.fin_external_absorption_area, expected_fin_area)
        external_bcs = [
            condition
            for condition in hp.hp.boundaries['outer_con'].conditions
            if isinstance(condition, ExternalHeatFluxBC)
        ]
        self.assertEqual(len(external_bcs), 1)
        np.testing.assert_allclose(external_bcs[0].area_array, expected_wall_area)
        self.assertAlmostEqual(hp.fin_external_area_scale, absorption_factor)

        wall_abs, fin_abs, total_abs = hp.get_external_heat_absorption_distribution(0.0)
        raw_q = hp.external_heat_accounting_source.get_heat_flux(0.0)
        np.testing.assert_allclose(wall_abs, raw_q * expected_wall_area)
        np.testing.assert_allclose(fin_abs, raw_q * expected_fin_area)
        np.testing.assert_allclose(total_abs, wall_abs + fin_abs)
        hp.pre_step(dt=0.01, current_time=0.0)
        np.testing.assert_allclose(
            hp.last_fin_absorption_distribution,
            raw_q * expected_fin_area,
        )

        local_total = ring_hp.get_total_external_heat_absorption(0.0)
        scaled_total = ring_hp.get_total_external_heat_absorption_scaled(0.0)
        expected_scaled = 0.0
        hp_position = 0
        for node_index, present in enumerate(ring_hp._hp_presence_mask):
            if present:
                expected_scaled += ring_hp._hp_multipliers[node_index] * float(
                    np.sum(ring_hp.hp_units[hp_position].get_external_heat_absorption_distribution(0.0)[2])
                )
                hp_position += 1
        self.assertGreater(local_total, 0.0)
        self.assertAlmostEqual(scaled_total, expected_scaled)
    def test_full_loop_cases_do_not_import_legacy_case_builders(self):
        root = Path(__file__).resolve().parent / "Full_Loop_Cases"
        violations = []
        for py_file in root.glob("*.py"):
            tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
            for node in ast.walk(tree):
                module = None
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        module = alias.name
                        if module in BANNED_MODULES:
                            violations.append((py_file.name, module))
                elif isinstance(node, ast.ImportFrom):
                    module = node.module or ""
                    if node.level:
                        continue
                    if module in BANNED_MODULES:
                        violations.append((py_file.name, module))
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
