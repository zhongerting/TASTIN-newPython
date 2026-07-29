import json
from pathlib import Path
import tempfile
import unittest

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control import (
    run_v14_210kw_reactivity_control as runner_module,
)


class V14ReactivityControlRunnerTests(unittest.TestCase):
    def test_runner_module_exists(self):
        runner = Path(__file__).with_name('run_v14_210kw_reactivity_control.py')
        self.assertTrue(runner.exists())

    def test_loads_model_settings_from_restart_directory(self):
        self.assertTrue(hasattr(runner_module, 'ReactivityControlRunConfig'))
        self.assertTrue(hasattr(runner_module, 'load_baseline_debug_config'))
        runtime_type = runner_module.ReactivityControlRunConfig
        load_config = runner_module.load_baseline_debug_config
        with tempfile.TemporaryDirectory(dir=r'E:\tmp') as tmp:
            root = Path(tmp)
            restart = root / 'steady_restart.npz'
            restart.write_bytes(b'placeholder')
            source = {
                'checkpoint_interval_s': 20.0,
                'external_heat_enabled': False,
                'fluid_max_iter': 77,
                'hp_up_view_factor': 0.0,
                'initial_temperature_k': 754.15,
                'inner_iter': 1,
                'lower_hp_down_view_factor': 0.4,
                'point_kinetics_enabled': False,
                'power_w': 210000.0,
                'radiator_emissivity': 0.7475,
                'space_temperature_k': 4.0,
                'target_flow_kg_s': 2.46,
                'tec_current_guess_a': 206.0,
                'tec_electrical_calculation_enabled': True,
                'tec_lookup_db': 'lookup-db',
                'tec_lookup_enabled': True,
                'tec_lookup_regions': ['core', 'accident'],
                'tec_voltage_v': 50.65,
                'upper_hp_down_view_factor': 0.3,
                'wire_resistance_scale': 0.335,
            }
            (root / 'run_config.json').write_text(json.dumps(source), encoding='utf-8')
            runtime = runtime_type(
                restart_in=restart,
                output_dir=root / 'output',
            )

            debug, loaded = load_config(runtime)

            self.assertEqual(loaded, source)
            self.assertEqual(debug.restart_in, restart)
            self.assertEqual(debug.output_dir, runtime.output_dir)
            self.assertAlmostEqual(debug.power_w, 210000.0)
            self.assertAlmostEqual(debug.tec_voltage_v, 50.65)
            self.assertAlmostEqual(debug.wire_resistance_scale, 0.335)
            self.assertTrue(debug.tec_electrical_enabled)
            self.assertEqual(tuple(debug.tec_lookup_regions), ('core', 'accident'))

    def test_fixed_power_restart_initializes_point_kinetics_once(self):
        self.assertTrue(hasattr(runner_module, 'prepare_reactivity_control'))

        class Core:
            has_point_reactor = False
            last_total_core_power = 210000.0
            initialized_power = None
            control_drum_reactivity_model = type('Drum', (), {'enabled': False})()

            def initialize_point_reactor(self, total_power_initial):
                self.initialized_power = total_power_initial
                self.has_point_reactor = True

            def get_effective_reactivity_feedback(self):
                return 0.0

        core = Core()
        handoff = runner_module.prepare_reactivity_control(
            core,
            source_point_kinetics_enabled=False,
            expected_power_w=210000.0,
        )

        self.assertEqual(handoff, 'fixed_power_handoff')
        self.assertAlmostEqual(core.initialized_power, 210000.0)

    def test_reactivity_restart_preserves_existing_point_kinetics(self):
        self.assertTrue(hasattr(runner_module, 'prepare_reactivity_control'))

        class Core:
            has_point_reactor = True
            last_total_core_power = 209500.0
            initialize_calls = 0
            control_drum_reactivity_model = type('Drum', (), {'enabled': False})()

            def initialize_point_reactor(self, total_power_initial):
                self.initialize_calls += 1

        core = Core()
        handoff = runner_module.prepare_reactivity_control(
            core,
            source_point_kinetics_enabled=True,
            expected_power_w=210000.0,
        )

        self.assertEqual(handoff, 'reactivity_continuation')
        self.assertEqual(core.initialize_calls, 0)

    def test_diagnostics_use_only_relative_temperature_feedback(self):
        self.assertTrue(hasattr(runner_module, 'collect_reactivity_diagnostics'))

        class Feedback:
            fuel = -1.0e-4
            electrode = 2.0e-5
            moderator = -3.0e-5
            reflector = 1.0e-5
            total = -1.0e-4

        class Reference:
            total = -1.2e-4

        class PointReactor:
            total_power = 209500.0
            fission_power = 205000.0
            decay_power = 4500.0

        class Core:
            point_reactor = PointReactor()
            feedback_reference_result = Reference()
            control_drum_reactivity_model = type('Drum', (), {'enabled': False})()

            def compute_reactivity_feedback(self):
                return Feedback()

            def get_control_drum_reactivity(self):
                return 0.0

        diagnostics = runner_module.collect_reactivity_diagnostics(
            Core(),
            handoff_type='fixed_power_handoff',
            initial_power_w=210000.0,
        )

        self.assertEqual(diagnostics['external_reactivity'], 0.0)
        self.assertEqual(diagnostics['control_drum_reactivity'], 0.0)
        self.assertAlmostEqual(diagnostics['effective_temperature_feedback'], 2.0e-5)
        self.assertAlmostEqual(diagnostics['total_reactivity'], 2.0e-5)
        self.assertAlmostEqual(diagnostics['power_relative_change'], -500.0 / 210000.0)

    def test_sub_nanosecond_end_residual_does_not_create_tiny_step(self):
        self.assertTrue(hasattr(runner_module, '_next_step_dt'))
        self.assertIsNone(runner_module._next_step_dt(
            current_time=9173.849999998616,
            end_time=9173.849999998762,
            requested_dt=0.05,
        ))
        self.assertAlmostEqual(runner_module._next_step_dt(
            current_time=9173.8,
            end_time=9173.85,
            requested_dt=0.05,
        ), 0.05)


if __name__ == '__main__':
    unittest.main()
