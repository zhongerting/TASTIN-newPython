import unittest

from Solvers.Hydrodynamics.Components import MacroFlowJunction

from testModule.Full_Loop_Cases_10kW import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    build_v14_case_a_system,
    v14_basic_diagnostics,
)


class V14_10kWRadiatorTopologyTests(unittest.TestCase):
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

    def test_v14_10kw_uses_explicit_upper_lower_heatpipe_rings(self):
        build = self.build_case()
        cfg = build["radiator_config"]
        volume_names = {vol.name for vol in build["system"].fluid_solver.volumes_obj}
        junctions = build["system"].fluid_solver.junctions_obj
        junction_names = {junc.name for junc in junctions}

        self.assertEqual(cfg.hp_l_eva_m, 0.0605)
        self.assertEqual(cfg.hp_l_aba_m, 0.0415)
        self.assertEqual(cfg.hp_l_con_m, 0.47)
        self.assertEqual(build["upper_heatpipe_count"], 154)
        self.assertEqual(build["lower_heatpipe_count"], 186)
        self.assertEqual(len(build["upper_ring_sectors"]), 6)
        self.assertEqual(len(build["lower_ring_sectors"]), 6)
        self.assertEqual(len(build["ring_sectors"]), 12)
        self.assertEqual(len(build["ring_hps"]), 12)

        for prefix in ("Upper", "Lower"):
            for inlet_node in ("I1", "I2", "I3"):
                self.assertIn(f"{prefix}_InletMix_{inlet_node}", volume_names)
            for sector in ("A1_I1_to_O1", "A2_O1_to_I2", "A3_I2_to_O2", "A4_O2_to_I3", "A5_I3_to_O3", "A6_O3_to_I1"):
                self.assertIn(f"{prefix}_{sector}_Channel_Vol_01", volume_names)

        for outlet_node in ("O1", "O2", "O3"):
            self.assertIn(f"OutletMix_{outlet_node}", volume_names)
            self.assertNotIn(f"Upper_OutletMix_{outlet_node}", volume_names)
            self.assertNotIn(f"Lower_OutletMix_{outlet_node}", volume_names)

        for idx in (1, 2, 3):
            self.assertIn(f"J_HotOutletBranch_{idx}_to_Upper_InletMix_I{idx}", junction_names)
            self.assertIn(f"J_HotOutletBranch_{idx}_to_Lower_InletMix_I{idx}", junction_names)
            self.assertIn(f"J_OutletMix_O{idx}_to_Manifold_{idx}", junction_names)
            self.assertIn(f"J_Manifold_{idx}_to_RadiatorOutletHeader", junction_names)

        self.assertFalse(any(isinstance(j, MacroFlowJunction) and getattr(j, "multiplier", 1.0) == 2.0 for j in junctions))

        diag = v14_basic_diagnostics(build)
        self.assertAlmostEqual(diag["upper_ring_in_total_kg_s"], 1.3 * 154.0 / 340.0)
        self.assertAlmostEqual(diag["lower_ring_in_total_kg_s"], 1.3 * 186.0 / 340.0)



    def test_v14_10kw_loads_n18_external_heat_on_upper_and_lower_rings(self):
        build = self.build_case()
        ring_hps = build["ring_hps"]

        self.assertEqual(build["radiator_external_heat_matrix_key"], "is58p5_w0_8p12_N18_sum")
        self.assertEqual(len(ring_hps), 12)
        self.assertTrue(all(ring_hp._hp_external_heat_enabled.all() for ring_hp in ring_hps))
        self.assertGreater(ring_hps[0].get_total_external_heat_absorption_scaled(1.0), 0.0)
        self.assertAlmostEqual(
            ring_hps[0].get_total_external_heat_absorption_scaled(1.0),
            ring_hps[6].get_total_external_heat_absorption_scaled(1.0),
            places=9,
        )

if __name__ == "__main__":
    unittest.main()
