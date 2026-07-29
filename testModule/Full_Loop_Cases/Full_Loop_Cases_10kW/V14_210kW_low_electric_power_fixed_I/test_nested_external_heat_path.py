"""Regression test for the nested V14 case package data path."""

import unittest

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW import v14_heatpipe_radiator


class NestedExternalHeatPathTest(unittest.TestCase):
    def test_external_heat_csv_exists(self):
        self.assertTrue(v14_heatpipe_radiator.EXTERNAL_HEAT_CSV.is_file())


if __name__ == "__main__":
    unittest.main()
