import os
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from ThermoCalc import ThermoCalcWrapper as wrapper


class ThermoCalcLookupConfigTests(unittest.TestCase):
    def test_explicit_lookup_false_overrides_environment(self):
        calls = []
        enabled = []
        def fake_load(*args, **kwargs):
            calls.append((args, kwargs))
            return 7

        fake_solver = SimpleNamespace(
            InputData=lambda: SimpleNamespace(),
            CalculationMode=SimpleNamespace(FixedVoltage="fixed_u"),
            set_emission_lookup_enabled=lambda value: enabled.append(bool(value)),
        )
        env = {"THERMOCALC_LOOKUP_DB": "env-db", "THERMOCALC_ENABLE_LOOKUP": "1"}
        with patch.dict(os.environ, env):
            with patch.object(wrapper, "te_solver", fake_solver), patch.object(wrapper, "load_emission_lookup_database", fake_load):
                model = wrapper.ThermoCalcModel(1, 1, lookup_db="case-db", enable_lookup=False)

        self.assertEqual(calls, [])
        self.assertEqual(enabled, [False])
        self.assertEqual(model.lookup_db, "case-db")
        self.assertFalse(model.lookup_enabled)

    def test_explicit_lookup_true_loads_case_database_and_regions(self):
        calls = []

        def fake_load(*args, **kwargs):
            calls.append((args, kwargs))
            return 3

        fake_solver = SimpleNamespace(InputData=lambda: SimpleNamespace(), CalculationMode=SimpleNamespace(FixedVoltage="fixed_u"))
        with patch.object(wrapper, "te_solver", fake_solver), patch.object(wrapper, "load_emission_lookup_database", fake_load):
            model = wrapper.ThermoCalcModel(
                1,
                1,
                lookup_db="case-db",
                enable_lookup=True,
                lookup_regions=("hot", "cold"),
            )

        self.assertEqual(calls, [(('case-db',), {"enable": True, "regions": ("hot", "cold")})])
        self.assertEqual(model.lookup_loaded_blocks, 3)
        self.assertTrue(model.lookup_enabled)


if __name__ == "__main__":
    unittest.main()
