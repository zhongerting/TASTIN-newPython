import unittest
from pathlib import Path

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
        self.bound1 = type('Boundary', (), {
            'get_coupling_surface_snapshot': lambda self: (
                np.array([1001.0, 1004.0]), np.ones(2)
            ),
        })()
        self.bound2 = type('Boundary', (), {
            'get_coupling_surface_snapshot': lambda self: (
                np.array([1000.0, 1000.0]), np.ones(2)
            ),
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
    def __init__(self, violation='collector'):
        self.tfes = {
            name: type('TFEForLimit', (), {'solids': {
                'inner_clad': Solid([900.0, 1000.0]),
                'outer_clad': Solid([890.0, 910.0]),
                'pellet': Solid([2600.0, 2690.0]),
                'collector': Solid([1000.0, 1010.0]),
                'moderator': Solid([800.0, 850.0]),
            }})()
            for name in runner.REPRESENTATIVE_NAMES
        }
        self.mod_rings = [Solid([870.0, 900.0])]
        self.reflector = Solid([900.0, 950.0])
        if violation == 'collector':
            self.tfes['Ring3'].solids['collector'].T[1] = 1024.0
        elif violation == 'channel_wall':
            self.tfes['Center'].solids['inner_clad'].T[1] = 1059.0
        elif violation == 'pellet':
            self.tfes['Center'].solids['pellet'].T[1] = 2701.0
        elif violation == 'moderator':
            self.mod_rings[0].T[1] = 931.0
        elif violation == 'reflector':
            self.reflector.T[1] = 1001.0
        elif violation is not None:
            raise ValueError(violation)


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

    def test_each_temperature_limit_class_trips(self):
        for component in (
                'channel_wall', 'pellet', 'collector',
                'moderator', 'reflector'):
            with self.subTest(component=component):
                core = CoreForLimits(violation=component)
                trip = runner.find_limit_trip(
                    runner.collect_temperature_peaks(core)
                )
                self.assertEqual(trip['component'], component)

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

    def test_default_run_config_matches_approved_case(self):
        config = runner.HeliumAccidentRunConfig(restart_in=Path('steady.npz'))
        self.assertEqual(config.duration_s, 100.0)
        self.assertEqual(config.dt_s, 0.05)
        self.assertEqual(config.tec_update_interval_s, 0.05)
        self.assertEqual(config.record_interval_s, 0.1)
        self.assertEqual(config.checkpoint_interval_s, 10.0)
        self.assertEqual(config.wall_limit_k, 1058.0)
        self.assertEqual(config.pellet_limit_k, 2700.0)
        self.assertEqual(config.collector_limit_k, 1023.0)
        self.assertEqual(config.moderator_limit_k, 930.0)
        self.assertEqual(config.reflector_limit_k, 1000.0)

    def test_accident_tec_update_interval_is_applied_to_core(self):
        core = type('Core', (), {'thermo_update_interval': 0.8})()
        scheduler_threshold = runner.set_tec_update_interval(core, 0.05)
        elapsed_at_real_timestamp = 13864.25 - 13864.20
        self.assertLess(scheduler_threshold, 0.05)
        self.assertEqual(core.thermo_update_interval, scheduler_threshold)
        self.assertGreaterEqual(
            elapsed_at_real_timestamp,
            core.thermo_update_interval,
        )
        with self.assertRaises(ValueError):
            runner.set_tec_update_interval(core, 0.0)

    def test_nonfinite_trip_scans_all_numeric_diagnostics(self):
        row = {
            'core_total_power_W': 210000.0,
            'fission_power_W': 197000.0,
            'decay_power_W': 13000.0,
            'effective_temperature_feedback': 0.0,
            'total_reactivity': 0.0,
            'min_fluid_T_K': np.nan,
        }
        trip = runner._nonfinite_metric_trip(row)
        self.assertEqual(trip['component'], 'nonfinite_metric')
        self.assertEqual(trip['field'], 'min_fluid_T_K')

    def test_state_trip_fails_closed_on_solver_nonconvergence(self):
        config = runner.HeliumAccidentRunConfig(restart_in=Path('steady.npz'))
        row = {
            'core_total_power_W': 210000.0,
            'min_fluid_T_K': 700.0,
            'fluid_converged': False,
            'tec_main_converged': True,
        }
        reason, trip = runner.evaluate_state_trip(
            row,
            peaks=[],
            config=config,
            baseline_power_w=210000.0,
            require_tec_convergence=True,
        )
        self.assertEqual(reason, 'hydraulic_nonconvergence')
        self.assertEqual(trip['component'], 'fluid_solver')

        row['fluid_converged'] = True
        row['tec_main_converged'] = False
        reason, trip = runner.evaluate_state_trip(
            row,
            peaks=[],
            config=config,
            baseline_power_w=210000.0,
            require_tec_convergence=True,
        )
        self.assertEqual(reason, 'tec_nonconvergence')
        self.assertEqual(trip['component'], 'tec_main')

    def test_raw_state_scanner_catches_partial_nan(self):
        fluid = type('Fluid', (), {
            'T_vec': np.array([700.0, 710.0]),
            'P_vec': np.array([1.0e5, 1.0e5]),
            'h_vec': np.array([1.0, 2.0]),
            'rho_vec': np.array([800.0, 800.0]),
            'W_vec': np.array([2.46, 2.46]),
        })()
        bad_solid = Solid([900.0, np.nan])
        system = type('System', (), {
            'fluid_solver': fluid,
            'solid_components': {'bad_solid': bad_solid},
        })()
        trip = runner.find_nonfinite_model_state({'system': system})
        self.assertEqual(trip['component'], 'nonfinite_model_state')
        self.assertEqual(trip['field'], 'solid:bad_solid.T')
        self.assertEqual(trip['flat_index'], 1)

        bad_solid.T[:] = 900.0
        fluid.W_vec[1] = np.inf
        trip = runner.find_nonfinite_model_state({'system': system})
        self.assertEqual(trip['field'], 'fluid.W_vec')
        self.assertEqual(trip['flat_index'], 1)

    def test_refresh_tec_now_updates_and_applies_main_group(self):
        calc = type('Calc', (), {
            'calculate': lambda self, verbose=False: setattr(
                self, 'calculate_called', True
            ),
        })()
        group = type('Group', (), {
            'name': 'main',
            'thermo_calc': calc,
            'last_update_time': -999.0,
        })()
        core = type('Core', (), {
            'enable_tec_coupled': True,
            '_last_thermo_update_time': -999.0,
            'iter_tec_circuit_groups': lambda self: [group],
            '_sync_tec_group_temperatures': lambda self, item: setattr(
                self, 'synced_group', item
            ),
            '_apply_tec_group_results': lambda self, item: setattr(
                self, 'applied_group', item
            ),
        })()
        runner.refresh_tec_now(core, 13864.2)
        self.assertTrue(calc.calculate_called)
        self.assertIs(core.synced_group, group)
        self.assertIs(core.applied_group, group)
        self.assertEqual(group.last_update_time, 13864.2)
        self.assertEqual(core._last_thermo_update_time, 13864.2)

    def test_accident_restart_is_reapplied_without_retrigger(self):
        gaps = runner.collect_helium_gaps(self.make_build())
        event = runner.restore_or_trigger_accident(
            gaps,
            source_config={
                'helium_accident_active': True,
                'helium_accident_time_absolute_s': 13864.2,
            },
            current_time_s=13874.2,
        )
        self.assertFalse(event['triggered_now'])
        self.assertEqual(event['accident_time_absolute_s'], 13864.2)
        self.assertTrue(all(gap.k_gas == 0.0 for gap, _ in gaps.values()))


if __name__ == '__main__':
    unittest.main()
