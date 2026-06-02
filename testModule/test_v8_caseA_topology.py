import os
import sys
import unittest

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from run_v8_caseA_common import passive_tec_source_totals
from test_core_assemble_v8_caseA import build_v8_case_a_system


class TestV8CaseATopology(unittest.TestCase):
    def setUp(self):
        self.build = build_v8_case_a_system(
            pipe_n_nodes=1,
            solid_heat_capacity_scale=1.0,
            solid_heat_capacity_scale_scope="global_outer",
        )

    def test_outer_ring_is_split_without_adding_physical_moderator_ring(self):
        core = self.build["core"]
        self.assertEqual(
            list(self.build["tfes"]),
            ["Center", "Ring1", "Ring2", "Ring3_TEC", "Ring3_Open"],
        )
        self.assertEqual(len(core.mod_rings), 4)
        self.assertEqual(core.get_ring_member_names(3), ["Ring3_TEC", "Ring3_Open"])
        self.assertEqual(sum(self.build["ring_multipliers"].values()), 37)
        self.assertEqual(sum(self.build["tec_ring_multipliers"].values()), 34)

    def test_outer_ring_split_preserves_normalized_power_factors(self):
        core = self.build["core"]
        self.assertAlmostEqual(core.power_factor_weighted_sum, 1.0, places=12)
        self.assertEqual(
            core.tfe_power_factors["Ring3_TEC"],
            core.tfe_power_factors["Ring3_Open"],
        )

    def test_passive_outer_tfe_clears_only_active_tec_sources(self):
        tfe = self.build["tfes"]["Ring3_Open"]
        tfe.update_joule_power_sources(
            Q_emitter_axial=np.ones(tfe.mesh.n_axial),
            Q_collector_axial=np.ones(tfe.mesh.n_axial),
        )
        tfe.update_plasma_flux(
            q_e_flux=-np.ones(tfe.mesh.n_axial),
            q_c_flux=np.ones(tfe.mesh.n_axial),
        )
        tfe.clear_tec_sources()
        self.assertTrue(all(value == 0.0 for value in passive_tec_source_totals(self.build).values()))
        self.assertGreater(float(tfe.couplers["tec_couple"].k_gas), 0.0)


if __name__ == "__main__":
    unittest.main()
