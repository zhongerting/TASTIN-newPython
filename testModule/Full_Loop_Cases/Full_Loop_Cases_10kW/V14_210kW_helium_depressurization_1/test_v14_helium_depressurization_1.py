from types import SimpleNamespace
import unittest

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_fixed_power_LOCA_1 import (
    run_v14_210kw_fixed_power_loca_1 as loca,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization_1 import (
    run_v14_helium_depressurization_1 as runner,
)


class HeliumDepressurization1Tests(unittest.TestCase):
    def test_defaults_match_requested_case(self):
        config = runner.HeliumRunConfig()
        self.assertTrue(config.restart_in.is_file())
        self.assertEqual(config.duration_s, 2000.0)
        self.assertEqual(config.collector_failure_temperature_k, 1500.0)
        self.assertEqual(config.emitter_failure_temperature_k, 3000.0)
        self.assertEqual(config.coolant_failure_temperature_k, 1058.0)
        self.assertEqual(config.moderator_failure_temperature_k, 930.0)
        self.assertEqual(config.reflector_failure_temperature_k, 1000.0)

    def test_live_coolant_history_uses_current_state(self):
        net = SimpleNamespace(
            volumes_obj=[SimpleNamespace(name="v")], junctions_obj=[],
            T_vec=np.array([800.0]), P_vec=np.array([2.0]),
            h_vec=np.array([3.0]), W_vec=np.array([]),
        )
        payload = loca._fluid_payload(
            SimpleNamespace(fluid_solver=net),
            {"T": np.array([700.0]), "P": np.array([1.0]),
             "h": np.array([1.5]), "W": np.array([])},
            coolant_present=True,
        )
        self.assertEqual(payload["fluid_temperature_K"][0], 800.0)
        self.assertEqual(payload["fluid_reference_temperature_K"][0], 700.0)


if __name__ == "__main__":
    unittest.main()
