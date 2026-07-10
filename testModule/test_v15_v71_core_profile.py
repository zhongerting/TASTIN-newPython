import unittest

import numpy as np

from testModule.Full_Loop_Cases import (
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V15PipeFinRadiatorConfig,
    build_v15_v71_case_a_system,
    build_center_uniform_axial_power_profile,
)


class V15V71CoreProfileTests(unittest.TestCase):
    def test_center_uniform_profile_uses_center_0p30m_overlap(self):
        profile = build_center_uniform_axial_power_profile(6, 25, 6, 0.30)

        self.assertEqual(len(profile), 37)
        self.assertAlmostEqual(float(np.sum(profile)), 1.0)
        self.assertEqual(float(np.sum(profile[:6])), 0.0)
        self.assertEqual(float(np.sum(profile[-6:])), 0.0)
        heated = np.flatnonzero(profile > 0.0)
        self.assertEqual(heated.tolist(), list(range(int(heated[0]), int(heated[-1]) + 1)))
        self.assertAlmostEqual(
            float(np.sum(profile * (np.arange(37) + 0.5))),
            18.5,
            places=12,
        )

    def test_v71_builder_applies_center_uniform_profile_to_all_tfes(self):
        build = build_v15_v71_case_a_system(
            core_config=FullLoopCoreConfig(main_tec_enabled=False),
            flow_config=FullLoopFlowConfig(total_flow_kg_s=1.3),
            pump_config=FullLoopPumpConfig(pump_total_head_pa=6466.56),
            radiator_config=V15PipeFinRadiatorConfig(tube_emissivity=0.0, fin_emissivity=0.0),
        )
        expected = build_center_uniform_axial_power_profile(6, 25, 6, 0.30)

        self.assertEqual(build["case_version"], "v15_v71_center0p30_uniform_pipefin_full_loop")
        self.assertEqual(build["axial_power_profile_name"], "center_0p30m_uniform")
        for tfe in build["tfes"].values():
            np.testing.assert_allclose(tfe.axial_power_profile, expected)
            np.testing.assert_allclose(tfe.solids["pellet"].power_allocation_weights.reshape(5, 37).sum(axis=0), expected)


if __name__ == "__main__":
    unittest.main()
