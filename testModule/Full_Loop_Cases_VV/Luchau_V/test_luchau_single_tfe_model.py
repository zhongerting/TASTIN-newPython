import math
import unittest
from types import SimpleNamespace

import numpy as np

from testModule.Full_Loop_Cases_VV.Luchau_V.luchau_single_tfe_model import (
    LUCHAU_AXIAL_NODE_ALLOCATION,
    LuchauSingleTFEConfig,
    LuchauSingleTFEThermoCalcCoupler,
    build_center_heater_profile,
    build_luchau_single_tfe,
    build_node_lengths,
    configure_luchau_thermocalc,
    tcs_from_cesium_pressure,
)


class LuchauSingleTFEModelTests(unittest.TestCase):
    def test_center_heater_profile_is_normalized_over_center_30cm(self):
        cfg = LuchauSingleTFEConfig(thermal_power_w=3000.0, target_voltage_v=1.0)

        profile = build_center_heater_profile(cfg)
        node_lengths = build_node_lengths()

        self.assertEqual(profile.shape, (sum(LUCHAU_AXIAL_NODE_ALLOCATION),))
        self.assertAlmostEqual(float(np.sum(profile)), 1.0, places=12)
        self.assertGreater(int(np.count_nonzero(profile)), 0)

        implied_overlap_lengths = profile * cfg.heater_length_m
        self.assertAlmostEqual(float(np.sum(implied_overlap_lengths)), cfg.heater_length_m, places=12)
        self.assertLess(float(np.min(implied_overlap_lengths[profile > 0.0])), float(np.max(node_lengths)))

    def test_build_single_tfe_uses_requested_flow_temperature_and_power(self):
        cfg = LuchauSingleTFEConfig(thermal_power_w=3000.0, target_voltage_v=1.0)

        build = build_luchau_single_tfe(cfg)

        self.assertEqual(build["tfe"].mesh.n_axial, 37)
        self.assertAlmostEqual(float(build["tfe"].geom.height), 0.507, places=12)
        self.assertAlmostEqual(float(build["inlet"].T), 727.0, places=3)
        self.assertAlmostEqual(float(build["inlet_junction"].W), 1.3 / 37.0, places=12)
        self.assertAlmostEqual(float(np.sum(build["tfe"].solids["pellet"].Q_source)), 3000.0, places=9)

        pellet_q = build["tfe"].solids["pellet"].Q_source.reshape(
            build["tfe"].mesh.n_r_pellet,
            build["tfe"].mesh.n_axial,
        )
        axial_q = np.sum(pellet_q, axis=0)
        self.assertTrue(np.all(axial_q[build["heater_profile"] == 0.0] == 0.0))

    def test_power_and_voltage_are_required_runtime_parameters(self):
        with self.assertRaisesRegex(ValueError, "thermal_power_w"):
            LuchauSingleTFEConfig(thermal_power_w=None, target_voltage_v=1.0)
        with self.assertRaisesRegex(ValueError, "target_voltage_v"):
            LuchauSingleTFEConfig(thermal_power_w=3000.0, target_voltage_v=None)

    def test_configure_thermocalc_uses_single_fixed_voltage_geometry(self):
        cfg = LuchauSingleTFEConfig(thermal_power_w=3000.0, target_voltage_v=1.2)
        build = build_luchau_single_tfe(cfg)
        fake = FakeThermoCalc()

        configure_luchau_thermocalc(fake, build, cfg)

        self.assertEqual(fake.mode, "fixed_u")
        self.assertAlmostEqual(fake.target_value, 1.2, places=12)
        self.assertEqual(fake.emitter_temperatures.shape, (1, 37))
        self.assertEqual(fake.collector_temperatures.shape, (1, 37))
        self.assertEqual(fake.tcs.shape, (1, 37))
        self.assertTrue(np.all(np.isfinite(fake.tcs)))
        self.assertAlmostEqual(float(fake.tcs[0, 0]), tcs_from_cesium_pressure(0.4), places=9)

        self.assertEqual(fake._input_data.dlE.shape, (1, 37))
        self.assertEqual(fake._input_data.sideAreaE.shape, (1, 37))
        self.assertAlmostEqual(float(fake._input_data.d_gap[0]), 0.5, places=12)

        geom = build["tfe"].geom
        expected_emitter_cross = math.pi * (geom.r_emitter_outer**2 - geom.r_fission_gas_outer**2)
        expected_collector_cross = math.pi * (geom.r_collector_outer**2 - geom.r_collector_inner**2)
        self.assertAlmostEqual(float(fake._input_data.crossAreaE[0]), expected_emitter_cross, places=15)
        self.assertAlmostEqual(float(fake._input_data.crossAreaC[0]), expected_collector_cross, places=15)

    def test_runner_requires_power_and_voltage_arguments(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import build_parser

        parser = build_parser()
        args = parser.parse_args([
            "--thermal-power-w",
            "3000",
            "--target-voltage-v",
            "1.2",
            "--skip-thermocalc-calc",
        ])
        self.assertAlmostEqual(args.thermal_power_w, 3000.0, places=12)
        self.assertAlmostEqual(args.target_voltage_v, 1.2, places=12)
        self.assertTrue(args.skip_thermocalc_calc)

        with self.assertRaises(SystemExit):
            parser.parse_args(["--target-voltage-v", "1.2"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--thermal-power-w", "3000"])

    def test_single_tfe_thermocalc_coupler_updates_temperatures_and_sources(self):
        cfg = LuchauSingleTFEConfig(thermal_power_w=3000.0, target_voltage_v=0.3)
        build = build_luchau_single_tfe(cfg)
        fake = FakeThermoCalc()
        configure_luchau_thermocalc(fake, build, cfg)
        fake.tec_result = {
            "UE": np.zeros(37),
            "UC": np.zeros(37),
            "rhoE": np.ones(37),
            "rhoC": np.ones(37),
            "joulePowerE": np.full(37, 0.1),
            "joulePowerC": np.full(37, 0.2),
            "J": np.full(37, 0.001),
            "phiE": np.full(37, 2.0),
            "TE": np.full(37, 1200.0),
        }
        coupler = LuchauSingleTFEThermoCalcCoupler(
            name="fake_coupler",
            tfe=build["tfe"],
            thermo_model=fake,
            alpha_tec=1.0,
        )

        coupler.pre_step(dt=0.1, current_time=0.0)

        self.assertEqual(fake.calculate_count, 1)
        self.assertEqual(fake.emitter_temperatures.shape, (1, 37))
        self.assertTrue(np.all(build["tfe"].electric_data.emitter_joule_heat > 0.0))
        self.assertTrue(np.any(build["tfe"].plasma_data.electron_cooling_flux < 0.0))
        self.assertIs(coupler.last_global_results, fake.global_results)

    def test_steady_runner_stops_when_temperature_rate_is_below_threshold(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_single_tfe import _advance_system_until_steady

        system = FakeSystem()
        build = {"tfe": SimpleNamespace(solids={"emitter": SimpleNamespace(T=np.array([1000.0005]))})}
        result = _advance_system_until_steady(
            system=system,
            build=build,
            duration_s=10.0,
            dt_s=0.1,
            steady_dtemp_k_s=1.0e-2,
            steady_window_steps=1,
        )

        self.assertTrue(result["steady_reached"])
        self.assertEqual(result["steps"], 1)
        self.assertAlmostEqual(result["final_time_s"], 0.1, places=12)
    def test_preheat_sweep_voltage_points_include_endpoint(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_preheat_then_sweep import _voltage_points

        points = _voltage_points(0.30, 0.40, 0.02)

        self.assertEqual(points, [0.30, 0.32, 0.34, 0.36, 0.38, 0.40])
    def test_preheat_sweep_defaults_to_corrected_dense_lookup_database(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_luchau_preheat_then_sweep import (
            DEFAULT_LOOKUP_DB,
            DEFAULT_LOOKUP_REGIONS,
        )

        self.assertTrue(DEFAULT_LOOKUP_DB.name == "pcs_0p02_5torr")
        self.assertTrue((DEFAULT_LOOKUP_DB / "runtime_dense_manifest.json").exists())
        self.assertEqual(DEFAULT_LOOKUP_REGIONS, ("core", "startup", "high_power", "accident"))
    def test_fit_scan_parse_float_list_accepts_scalar_values(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.fit_all_powers_cs_wire import _parse_float_list

        self.assertEqual(_parse_float_list(1.35, positive=True), [1.35])
        self.assertEqual(_parse_float_list("1.2,1.35", positive=True), [1.2, 1.35])

    def test_fit_scan_progress_interval_uses_at_least_one_step(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.fit_all_powers_cs_wire import _progress_interval_steps

        self.assertEqual(_progress_interval_steps(duration_s=200.0, dt_s=0.05, progress_interval_s=10.0), 200)
        self.assertEqual(_progress_interval_steps(duration_s=1.0, dt_s=0.05, progress_interval_s=0.0), 20)
        self.assertEqual(_progress_interval_steps(duration_s=1.0, dt_s=2.0, progress_interval_s=0.1), 1)
    def test_fit_scan_builds_manual_single_power_voltage_targets(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.fit_all_powers_cs_wire import _manual_targets_from_args

        rows = _manual_targets_from_args(
            power_kwt=3.26,
            voltage_list=[0.30, 0.40],
            target_current_list=[410.0, 300.0],
        )

        self.assertEqual([row["power_kwt"] for row in rows], [3.26, 3.26])
        self.assertEqual([row["voltage_v"] for row in rows], [0.30, 0.40])
        self.assertEqual([row["target_current_a"] for row in rows], [410.0, 300.0])
    def test_validation_scan_uses_power_indexed_directory_and_restart_name(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_validation_power_scan import (
            validation_directory,
            validation_restart_path,
        )

        base = validation_directory(5.6)
        restart = validation_restart_path(5.6, 1.0469)

        self.assertEqual(base.name, "Validation_6-5.6")
        self.assertEqual(restart.parent.name, "restarts")
        self.assertEqual(restart.name, "restart_5p6_1p0469.npz")

    def test_validation_scan_selects_first_or_continuation_schedule(self):
        from pathlib import Path
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_validation_power_scan import select_run_schedule

        first = select_run_schedule(Path("missing_restart.npz"))
        self.assertTrue(first["is_first_run"])
        self.assertEqual(first["duration_s"], 800.0)
        self.assertEqual(first["dt_s"], 0.05)

        existing = Path(__file__)
        continuation = select_run_schedule(existing)
        self.assertFalse(continuation["is_first_run"])
        self.assertEqual(continuation["duration_s"], 100.0)
        self.assertEqual(continuation["dt_s"], 0.1)

    def test_validation_scan_has_default_targets_for_all_powers(self):
        from testModule.Full_Loop_Cases_VV.Luchau_V.run_validation_power_scan import (
            DEFAULT_TARGET_SERIES_BY_POWER,
            _parse_float_list,
        )

        self.assertEqual(set(DEFAULT_TARGET_SERIES_BY_POWER), {3.26, 3.52, 4.01, 4.49, 5.0, 5.6})
        for power, series in DEFAULT_TARGET_SERIES_BY_POWER.items():
            voltages = _parse_float_list(series["voltage_list"])
            currents = _parse_float_list(series["target_current_list"])
            self.assertEqual(len(voltages), len(currents), power)
            self.assertGreater(len(voltages), 0, power)

        self.assertAlmostEqual(_parse_float_list(DEFAULT_TARGET_SERIES_BY_POWER[4.49]["voltage_list"])[0], 0.56202)
        self.assertAlmostEqual(_parse_float_list(DEFAULT_TARGET_SERIES_BY_POWER[3.26]["target_current_list"])[-1], 198.5581)
class FakeThermoCalc:
    def __init__(self):
        self._input_data = SimpleNamespace()
        self.mode = None
        self.target_value = None
        self.i_guess = None
        self.emitter_temperatures = None
        self.collector_temperatures = None
        self.tcs = None
        self.calculate_count = 0
        self.tec_result = None
        self.global_results = {"Iout": 1.0, "Uout": 0.3, "converged": True}

    def setup_circuit_mode(self, mode_str, target_value, I_guess=150.0):
        self.mode = mode_str
        self.target_value = float(target_value)
        self.i_guess = float(I_guess)

    def set_temperatures(self, emitter, collector):
        self.emitter_temperatures = np.asarray(emitter, dtype=float)
        self.collector_temperatures = np.asarray(collector, dtype=float)

    def set_tcs(self, tcs):
        self.tcs = np.asarray(tcs, dtype=float)

    def calculate(self, verbose=False):
        self.calculate_count += 1

    def get_tec_results(self, idx):
        return self.tec_result

    def get_global_results(self):
        return self.global_results


class FakeSystem:
    def __init__(self):
        self.global_time = 0.0
        self.step_count = 0

    def step(self, dt):
        self.global_time += float(dt)
        self.step_count += 1


if __name__ == "__main__":
    unittest.main()
