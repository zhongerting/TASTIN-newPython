import unittest

import numpy as np

from testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization import (
    run_v14_helium_depressurization as runner,
)


class Gap:
    def __init__(self):
        self.gap = 5.0e-5
        self.k_gas = 5678.0 * self.gap
        self.eps1 = 0.60
        self.eps2 = 0.80
        self.R_gap_total = np.array([1.0, 2.0])
        self.bc1 = type('Boundary', (), {
            'current_flux': np.array([-1.0, -2.0]),
        })()


class TFE:
    def __init__(self):
        self.couplers = {'collector_iclad_gap': Gap()}
        self.solids = {
            'collector': Solid([1000.0, 1010.0]),
            'inner_clad': Solid([900.0, 910.0]),
        }


class Mesh:
    class Geom:
        node_centers_y = np.array([0.1, 0.2])

    geom_data = Geom()


class Solid:
    def __init__(self, values):
        self.T = np.asarray(values, dtype=float)
        self.mesh = Mesh()


class CoreForLimits:
    def __init__(self):
        self.tfes = {
            name: type('TFEForLimit', (), {'solids': {
                'inner_clad': Solid([900.0, 1000.0]),
                'outer_clad': Solid([890.0, 910.0]),
                'pellet': Solid([2600.0, 2690.0]),
                'collector': Solid(
                    [1000.0, 1024.0]
                    if name == 'Ring3' else [1000.0, 1010.0]
                ),
                'moderator': Solid([800.0, 850.0]),
            }})()
            for name in runner.REPRESENTATIVE_NAMES
        }
        self.mod_rings = [Solid([870.0, 900.0])]
        self.reflector = Solid([900.0, 950.0])


class HeliumGapTests(unittest.TestCase):
    def make_build(self):
        names = runner.REPRESENTATIVE_NAMES
        return {
            'tfes': {name: TFE() for name in names},
            'ring_multipliers': dict(zip(names, runner.EXPECTED_MULTIPLIERS)),
            'system': type('System', (), {'global_time': 15.0})(),
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

    def test_collector_limit_reports_ring_and_axial_position(self):
        peaks = runner.collect_temperature_peaks(CoreForLimits())
        trip = runner.find_limit_trip(peaks)
        self.assertEqual(trip['component'], 'collector')
        self.assertEqual(trip['representative'], 'Ring3')
        self.assertEqual(trip['limit_k'], 1023.0)
        self.assertEqual(trip['actual_k'], 1024.0)
        self.assertEqual(trip['axial_position_m'], 0.2)

    def test_nonfinite_temperature_trips_before_numeric_limits(self):
        core = CoreForLimits()
        core.tfes['Center'].solids['collector'].T[0] = np.nan
        trip = runner.find_limit_trip(runner.collect_temperature_peaks(core))
        self.assertEqual(trip['component'], 'nonfinite_temperature')
        self.assertEqual(trip['source_component'], 'collector')
        self.assertEqual(trip['representative'], 'Center')

    def test_helium_metrics_scale_representative_gap_heat(self):
        build = self.make_build()
        gaps = runner.collect_helium_gaps(build)
        runner.set_helium_h_eq(gaps, 0.0)
        row = runner.collect_helium_metrics(
            build,
            gaps,
            accident_time_s=10.0,
            active=True,
        )
        self.assertEqual(row['accident_elapsed_s'], 5.0)
        self.assertEqual(row['helium_h_eq_W_m2K'], 0.0)
        self.assertEqual(row['helium_gap_heat_out_scaled_W'], 174.0)
        self.assertEqual(row['helium_gap_R_total_min_K_W'], 1.0)
        self.assertEqual(row['helium_gap_R_total_max_K_W'], 2.0)


if __name__ == '__main__':
    unittest.main()
