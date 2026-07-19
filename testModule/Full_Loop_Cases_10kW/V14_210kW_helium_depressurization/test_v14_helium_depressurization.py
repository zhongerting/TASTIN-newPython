import unittest

from testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization import (
    run_v14_helium_depressurization as runner,
)


class Gap:
    def __init__(self):
        self.gap = 5.0e-5
        self.k_gas = 5678.0 * self.gap
        self.eps1 = 0.60
        self.eps2 = 0.80


class TFE:
    def __init__(self):
        self.couplers = {'collector_iclad_gap': Gap()}


class HeliumGapTests(unittest.TestCase):
    def make_build(self):
        names = runner.REPRESENTATIVE_NAMES
        return {
            'tfes': {name: TFE() for name in names},
            'ring_multipliers': dict(zip(names, runner.EXPECTED_MULTIPLIERS)),
        }

    def test_collects_exactly_five_expected_gaps(self):
        gaps = runner.collect_helium_gaps(self.make_build())
        self.assertEqual(tuple(gaps), runner.REPRESENTATIVE_NAMES)
        self.assertEqual(
            tuple(multiplier for _, multiplier in gaps.values()),
            runner.EXPECTED_MULTIPLIERS,
        )

    def test_instantaneous_loss_only_clears_gas_conduction(self):
        gaps = runner.collect_helium_gaps(self.make_build())
        runner.set_helium_h_eq(gaps, 0.0)
        for gap, _ in gaps.values():
            self.assertEqual(gap.k_gas, 0.0)
            self.assertEqual(gap.eps1, 0.60)
            self.assertEqual(gap.eps2, 0.80)

    def test_source_accident_marker_is_required_and_boolean(self):
        self.assertFalse(runner.read_source_accident_state(
            {'helium_accident_active': False}
        ))
        self.assertTrue(runner.read_source_accident_state(
            {'helium_accident_active': True}
        ))
        with self.assertRaises(ValueError):
            runner.read_source_accident_state({})
        with self.assertRaises(ValueError):
            runner.read_source_accident_state({'helium_accident_active': 'false'})


if __name__ == '__main__':
    unittest.main()
