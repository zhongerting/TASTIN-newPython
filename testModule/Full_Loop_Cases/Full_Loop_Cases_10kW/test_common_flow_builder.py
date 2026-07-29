import unittest

from Materials.Fluids.SodiumPotassium78 import SodiumPotassium78
from Solvers.Hydrodynamics.Components import PumpJunction

from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig
from .common_flow_builder import FlowControlledPumpJunction, build_common_flow_objects


class CommonFlowBuilderTests(unittest.TestCase):
    def test_two_series_pumps_use_one_flow_controller_and_keep_total_head(self):
        total_head = 6466.56
        build = build_common_flow_objects(
            FullLoopCoreConfig(),
            FullLoopFlowConfig(total_flow_kg_s=2.46),
            FullLoopPumpConfig(
                pump_total_head_pa=total_head,
                pump_flow_control=True,
                target_flow_kg_s=2.46,
            ),
            material=SodiumPotassium78(),
            close_with_placeholder_bridge=True,
        )

        self.assertIsInstance(build['pump_a'], FlowControlledPumpJunction)
        self.assertIsInstance(build['pump_b'], PumpJunction)
        self.assertNotIsInstance(build['pump_b'], FlowControlledPumpJunction)
        self.assertFalse(hasattr(build['pump_b'], 'target_W'))
        self.assertAlmostEqual(
            build['pump_a'].compute_pump_head(0.0)
            + build['pump_b'].compute_pump_head(0.0),
            total_head,
        )


if __name__ == '__main__':
    unittest.main()
