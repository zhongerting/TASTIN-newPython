import math
import os
import sys
from types import SimpleNamespace

import numpy as np


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)


from testModule.run_v13_start_case import (
    configure_startup_main_tec_circuit,
    effective_startup_main_tec_load_resistance,
    estimate_tec_wire_joule_loss_w,
    prepare_startup_step,
    recover_startup_tec_state_from_restart,
    maybe_switch_startup_tec_to_fixed_voltage,
    parse_args,
    startup_energy_residual_diagnostics,
    startup_tec_global_diagnostics,
)
from testModule.v13_startup_control import (
    MatrixColumnHeatSource,
    MatrixColumnLineHeatSource,
    V13StartupControlConfig,
    V13StartupController,
    apply_tec_gap_h_eq,
    reset_fluid_temperatures,
    reset_solid_temperatures,
    shield_qsss_from_matrix,
)

import testModule.run_v13_start_case as run_case_module


class FakeGapCoupler:
    def __init__(self, gap):
        self.gap = gap
        self.k_gas = -1.0


class FakeTFE:
    def __init__(self, gap):
        self.couplers = {"tec_couple": FakeGapCoupler(gap)}
class FakeFluidMaterial:
    def enthalpy(self, T, P):
        return 10.0 * float(T) + 1.0e-6 * float(P)

    def temperature_from_enthalpy(self, h, P):
        return (float(h) - 1.0e-6 * float(P)) / 10.0

    def density(self, T, P):
        return 900.0

    def viscosity(self, T, P):
        return 1.0e-3

    def heat_capacity(self, T, P):
        return 10.0

    def liquid_density_derivative_T(self, T):
        return 0.0


class FakeFluidVolume:
    def __init__(self, name, T, P):
        self.name = name
        self.T = T
        self.P = P
        self.material = FakeFluidMaterial()
        self.h = self.material.enthalpy(T, P)
        self.Q_wall = 1.0
        self.Q_vol = 2.0
        self.implicit_coeff = 3.0

    def update_properties(self, material):
        self.T = material.temperature_from_enthalpy(self.h, self.P)
        self.rho = material.density(self.T, self.P)
        self.mu = material.viscosity(self.T, self.P)


class FakeFluidNetwork:
    def __init__(self):
        self.volumes_obj = [FakeFluidVolume("a", 700.0, 1.0e5), FakeFluidVolume("b", 800.0, 2.0e5)]
        self.T_vec = np.zeros(2)
        self.h_vec = np.zeros(2)
        self.P_vec = np.zeros(2)
        self.properties_updated = False

    def _initialize_state_from_objects(self):
        for idx, vol in enumerate(self.volumes_obj):
            self.T_vec[idx] = vol.T
            self.h_vec[idx] = vol.h
            self.P_vec[idx] = vol.P

    def _update_fluid_properties(self):
        self.properties_updated = True


class FakeFluidSystem:
    def __init__(self):
        self.fluid_solver = FakeFluidNetwork()

class FakeBoundary:
    def __init__(self):
        self.T_surface = np.full((2,), 800.0)
        self.T_adj_node = np.full((2,), 790.0)
        self.current_flux = np.full((2,), 10.0)


class FakeSolid:
    def __init__(self):
        self.T = np.full(3, 900.0)
        self.dTdt = np.ones(3)
        self.current_time = 5.0
        self.last_trial_temperature_min = 900.0
        self.last_trial_temperature_max = 900.0
        self.last_trial_temperature_time = 5.0
        self.last_step_success = False
        self.last_step_failure_message = "old"
        self.boundaries = {"right": FakeBoundary()}
        self.updated_properties = False

    def _update_properties(self):
        self.updated_properties = True


class FakeSystem:
    def __init__(self):
        self.solid_components = {"fake": FakeSolid()}


class FakeStartupCore:
    def __init__(self):
        self.calls = []
        class FakeThermo:
            N_elem = 34

        self.thermo_calc = FakeThermo()

    def setup_tec_circuit(self, mode_str, target_value, I_guess=150.0, topology="series"):
        self.calls.append({
            "mode_str": mode_str,
            "target_value": float(target_value),
            "I_guess": float(I_guess),
            "topology": topology,
        })

def assert_close(actual, expected, tol=1.0e-9):
    assert math.isclose(float(actual), float(expected), rel_tol=tol, abs_tol=tol), (
        actual,
        expected,
    )


def test_titam_startup_schedule_matches_document_milestones():
    controller = V13StartupController(V13StartupControlConfig())

    before_critical = controller.evaluate(0.0, core_inlet_temperature_k=373.0, emitter_temperature_k=373.0)
    assert before_critical.phase == "SAFETY_DRUM_WITHDRAWAL"
    assert_close(before_critical.safety_drum_angle_deg, 0.0)
    assert_close(before_critical.control_drum_angle_deg, 0.0)
    assert before_critical.thermal_power_w == controller.config.source_power_w

    critical_time = controller.config.critical_time_s
    at_critical = controller.evaluate(critical_time, core_inlet_temperature_k=373.0, emitter_temperature_k=373.0)
    assert at_critical.phase == "INITIAL_SUPERCRITICAL_RAMP"
    assert_close(at_critical.safety_drum_angle_deg, 180.0)
    assert_close(at_critical.control_drum_angle_deg, 125.0)

    low_power = controller.evaluate(
        controller.config.low_power_hold_start_s,
        core_inlet_temperature_k=373.0,
        emitter_temperature_k=373.0,
    )
    assert low_power.phase == "LOW_POWER_HOLD"
    assert_close(low_power.control_drum_angle_deg, 145.0)
    assert_close(low_power.thermal_power_w, 5000.0)

    fast_end = controller.evaluate(
        controller.config.fast_ramp_end_s,
        core_inlet_temperature_k=373.0,
        emitter_temperature_k=373.0,
    )
    assert fast_end.phase == "SLOW_POWER_RAMP"
    assert_close(fast_end.thermal_power_w, 35000.0)

    slow_end = controller.evaluate(
        controller.config.slow_ramp_end_s,
        core_inlet_temperature_k=373.0,
        emitter_temperature_k=373.0,
    )
    assert slow_end.phase == "CRITICAL_POWER_HOLD"
    assert_close(slow_end.thermal_power_w, 110000.0)
    assert 87.0 <= slow_end.control_drum_angle_deg <= 90.0


def test_shield_jettison_is_temperature_triggered_and_latched():
    controller = V13StartupController(V13StartupControlConfig(shield_jettison_temperature_k=400.0))

    cold = controller.evaluate(10.0, core_inlet_temperature_k=399.9, emitter_temperature_k=373.0)
    assert cold.radiation_shield_active is True
    assert cold.shield_jettisoned is False

    hot = controller.evaluate(11.0, core_inlet_temperature_k=400.0, emitter_temperature_k=373.0)
    assert hot.radiation_shield_active is False
    assert hot.shield_jettisoned is True

    cooled = controller.evaluate(12.0, core_inlet_temperature_k=390.0, emitter_temperature_k=373.0)
    assert cooled.radiation_shield_active is False
    assert cooled.shield_jettisoned is True


def test_default_cesium_gap_h_eq_matches_v7_steady_value():
    config = V13StartupControlConfig()

    assert_close(config.cesium_gap_h_eq_w_m2_k, 29.0)


def test_v13_start_runner_default_cesium_gap_h_eq_matches_v7_steady_value():
    old_argv = sys.argv[:]
    try:
        sys.argv = ["run_v13_start_case.py"]
        args = parse_args()
    finally:
        sys.argv = old_argv

    assert_close(args.cesium_gap_h_eq_w_m2_k, 29.0)

def test_tfe_ignition_immediately_sets_cesium_gap_and_enables_fixed_r_tec():
    config = V13StartupControlConfig(
        tfe_start_after_critical_s=100.0,
        tfe_start_emitter_temperature_k=1040.0,
        tec_electrical_start_after_cesium_s=0.0,
        tec_electrical_start_cs_fraction=0.0,
        tec_electrical_start_emitter_temperature_k=0.0,
        cs_transition_tau_s=100.0,
        helium_gap_h_eq_w_m2_k=1200.0,
        cesium_gap_h_eq_w_m2_k=29.0,
    )
    controller = V13StartupController(config)
    ready_time = config.critical_time_s + config.tfe_start_after_critical_s

    before_ignition = controller.evaluate(ready_time - 1.0, core_inlet_temperature_k=725.0, emitter_temperature_k=900.0)
    assert before_ignition.cesium_conditioning_started is False
    assert before_ignition.tec_enabled is False
    assert_close(before_ignition.cs_fraction, 0.0)
    assert_close(before_ignition.tec_gap_h_eq_w_m2_k, 1200.0)

    ignition = controller.evaluate(ready_time, core_inlet_temperature_k=725.0, emitter_temperature_k=900.0)
    assert ignition.cesium_conditioning_started is True
    assert ignition.tec_enabled is True
    assert_close(ignition.cs_fraction, 1.0)
    assert_close(ignition.tec_gap_h_eq_w_m2_k, 29.0)

def test_apply_tec_gap_h_eq_updates_all_tfe_gap_couplers():
    tfes = {
        "Center": FakeTFE(gap=2.5e-4),
        "Ring1": FakeTFE(gap=3.0e-4),
    }

    updated = apply_tec_gap_h_eq(tfes, h_eq_w_m2_k=800.0)

    assert updated == 2
    assert_close(tfes["Center"].couplers["tec_couple"].k_gas, 800.0 * 2.5e-4)
    assert_close(tfes["Ring1"].couplers["tec_couple"].k_gas, 800.0 * 3.0e-4)
def test_reset_fluid_temperatures_updates_volume_and_network_state():
    system = FakeFluidSystem()

    count = reset_fluid_temperatures(system, 373.0)

    network = system.fluid_solver
    assert count == 2
    assert network.properties_updated is True
    np.testing.assert_allclose(network.T_vec, np.full(2, 373.0))
    for idx, vol in enumerate(network.volumes_obj):
        assert vol.T == 373.0
        assert vol.h == network.h_vec[idx]
        assert vol.Q_wall == 0.0
        assert vol.Q_vol == 0.0
        assert vol.implicit_coeff == 0.0

def test_reset_solid_temperatures_updates_solid_and_boundary_caches():
    system = FakeSystem()

    count = reset_solid_temperatures(system, 373.0, current_time_s=0.0)

    solid = system.solid_components["fake"]
    assert count == 1
    np.testing.assert_allclose(solid.T, np.full(3, 373.0))
    np.testing.assert_allclose(solid.dTdt, np.zeros(3))
    assert solid.current_time == 0.0
    assert solid.last_trial_temperature_min == 373.0
    assert solid.last_trial_temperature_max == 373.0
    assert solid.last_step_success is True
    assert solid.last_step_failure_message == ""
    assert solid.updated_properties is True
    boundary = solid.boundaries["right"]
    np.testing.assert_allclose(boundary.T_surface, np.full(2, 373.0))
    np.testing.assert_allclose(boundary.T_adj_node, np.full(2, 373.0))
    np.testing.assert_allclose(boundary.current_flux, np.zeros(2))

def test_matrix_column_heat_source_broadcasts_one_tube_flux_to_axial_nodes():
    source = MatrixColumnHeatSource(
        shape=(8,),
        matrix_key="is58p5_w0_8p12_N78_sum",
        column_index=0,
        scale_factor=0.25,
        periodic=True,
    )

    flux = source.get_heat_flux(123.4)

    assert flux.shape == (8,)
    assert np.all(np.isfinite(flux))
    np.testing.assert_allclose(flux, np.full(8, flux[0]))
    assert flux[0] >= 0.0

def test_matrix_column_line_heat_source_converts_w_per_m_to_equivalent_flux_density():
    area = np.array([0.2, 0.4, 0.5])
    length = np.array([0.1, 0.2, 0.25])
    source = MatrixColumnLineHeatSource(
        shape=(3,),
        matrix_key="is58p5_w0_8p12_N78_sum",
        column_index=0,
        node_lengths_m=length,
        area_array_m2=area,
        scale_factor=1.0,
        periodic=True,
    )

    density = source.get_heat_flux(123.4)
    q_line = MatrixColumnHeatSource(
        shape=(3,),
        matrix_key="is58p5_w0_8p12_N78_sum",
        column_index=0,
        scale_factor=1.0,
        periodic=True,
    ).get_heat_flux(123.4)[0]

    np.testing.assert_allclose(density * area, q_line * length)

def test_fixed_startup_power_overrides_ramp_schedule():
    controller = V13StartupController(V13StartupControlConfig(fixed_power_w=0.0))

    before_critical = controller.evaluate(10.0, core_inlet_temperature_k=373.0, emitter_temperature_k=373.0)
    after_ramp_time = controller.config.slow_ramp_end_s + 500.0
    after_ramp = controller.evaluate(after_ramp_time, core_inlet_temperature_k=373.0, emitter_temperature_k=373.0)

    assert before_critical.thermal_power_w == 0.0
    assert after_ramp.thermal_power_w == 0.0
    assert after_ramp.fission_power_w == 0.0

def test_shield_qsss_from_matrix_maps_six_side_partitions_to_eight_nodes():
    qsss = shield_qsss_from_matrix(
        "is58p5_w0_8p12_N6_sum",
        time_s=123.4,
        scale_factor=0.5,
    )

    assert qsss.shape == (8,)
    assert np.all(np.isfinite(qsss))
    assert np.all(qsss[:6] >= 0.0)
    np.testing.assert_allclose(qsss[6:], np.zeros(2))

def test_startup_main_tec_uses_fixed_resistance_then_switches_to_fixed_voltage():
    class Args:
        startup_main_tec_initial_mode = "fixed_r"
        startup_main_tec_load_resistance_ohm = 0.0044
        startup_main_tec_load_resistance_scope = "total"
        startup_main_tec_switch_voltage_v = 27.2
        startup_main_tec_i_guess_a = 150.0
        target_voltage = 27.2

    core = FakeStartupCore()
    build = {"core": core}

    configure_startup_main_tec_circuit(core, Args)

    assert core.calls[-1]["mode_str"] == "fixed_r"
    assert_close(core.calls[-1]["target_value"], 0.0044)
    assert build.get("startup_main_tec_switched_to_fixed_voltage") is None

    below = {"tec_main_voltage_v": 26.9, "tec_main_current_a": 6000.0}
    switched = maybe_switch_startup_tec_to_fixed_voltage(build, Args, below)
    assert switched is False
    assert len(core.calls) == 1

    above = {"tec_main_voltage_v": 27.21, "tec_main_current_a": 6180.0}
    switched = maybe_switch_startup_tec_to_fixed_voltage(build, Args, above)
    assert switched is True
    assert build["startup_main_tec_switched_to_fixed_voltage"] is True
    assert core.calls[-1]["mode_str"] == "fixed_u"
    assert_close(core.calls[-1]["target_value"], 27.2)
    assert_close(core.calls[-1]["I_guess"], 6180.0)

def test_startup_tec_global_diagnostics_reports_convergence_and_finiteness():
    class Thermo:
        def get_global_results(self):
            return {
                "mode": "fixed_r",
                "converged": True,
                "iteration_count": 4,
                "zero_emission_skipped": False,
                "zero_emission_reason": None,
                "Iout": 64.0,
                "Uout": 9.7,
                "Rload": 0.1496,
            }

    class Core:
        thermo_calc = Thermo()

    diag = startup_tec_global_diagnostics(Core())

    assert diag["tec_solver_mode"] == "fixed_r"
    assert diag["tec_solver_converged"] is True
    assert diag["tec_solver_iteration_count"] == 4
    assert diag["tec_solver_output_finite"] is True
    assert diag["tec_solver_zero_emission_skipped"] is False



def test_startup_energy_residual_diagnostics_tracks_storage_and_radiator_balance():
    record = {
        "core_heat_power_w": 110000.0,
        "coolant_enthalpy_rise_w": 102978.091,
        "tec_total_electric_power_w": 5692.628,
        "q_radiator_total_w": 102977.955,
    }

    diag = startup_energy_residual_diagnostics(record)

    assert_close(diag["core_heat_minus_coolant_enthalpy_minus_electric_w"], 1329.281)
    assert_close(diag["core_heat_minus_radiator_minus_electric_w"], 1329.417)
    assert_close(diag["radiator_minus_coolant_enthalpy_w"], -0.136)
    assert_close(diag["core_energy_storage_residual_rel"], 1329.281 / 110000.0)
    assert_close(diag["radiator_coolant_balance_rel"], -0.136 / 102978.091)


def test_startup_energy_residual_diagnostics_treats_missing_electric_as_zero():
    record = {
        "startup_thermal_power_w": 100.0,
        "coolant_enthalpy_rise_w": 80.0,
        "q_radiator_total_w": 79.5,
    }

    diag = startup_energy_residual_diagnostics(record)

    assert_close(diag["core_heat_minus_coolant_enthalpy_minus_electric_w"], 20.0)
    assert_close(diag["core_heat_minus_radiator_minus_electric_w"], 20.5)
    assert_close(diag["radiator_minus_coolant_enthalpy_w"], -0.5)

def test_startup_main_tec_load_resistance_scope_can_use_per_tec_value():
    class Args:
        startup_main_tec_load_resistance_ohm = 0.0044
        startup_main_tec_load_resistance_scope = "per_tec"

    core = FakeStartupCore()

    assert_close(effective_startup_main_tec_load_resistance(core, Args), 0.1496)

    class TotalArgs:
        startup_main_tec_load_resistance_ohm = 0.0044
        startup_main_tec_load_resistance_scope = "total"

    assert_close(effective_startup_main_tec_load_resistance(core, TotalArgs), 0.0044)


def test_startup_controller_can_seed_existing_cesium_conditioning_state():
    config = V13StartupControlConfig(
        cs_transition_tau_s=100.0,
        tec_electrical_start_after_cesium_s=100000.0,
        helium_gap_h_eq_w_m2_k=600.0,
        cesium_gap_h_eq_w_m2_k=250.0,
    )
    controller = V13StartupController(config)
    controller.seed_cesium_conditioning(absolute_time_s=1000.0, cs_fraction=0.98)

    seeded = controller.evaluate(1000.0, core_inlet_temperature_k=740.0, emitter_temperature_k=1180.0)

    assert seeded.cesium_conditioning_started is True
    assert seeded.tec_enabled is False
    assert_close(seeded.cs_fraction, 1.0)
    assert_close(seeded.tec_gap_h_eq_w_m2_k, 250.0)



def test_estimate_tec_wire_joule_loss_w_uses_current_and_wire_series():
    record = {"tec_total_current_a": 2000.0, "wire_resistance_ohm": np.array([0.001, 0.0012, 0.0006, 0.0002])}

    loss = estimate_tec_wire_joule_loss_w(record)

    expected_wire_resistance = float(np.sum(record["wire_resistance_ohm"]))
    expected = 34 * (2000.0 / 2.0) ** 2 * expected_wire_resistance
    assert_close(loss, expected)


def test_startup_energy_residual_diagnostics_includes_external_heat_and_wire_loss():
    record = {
        "core_heat_power_w": 110000.0,
        "coolant_enthalpy_rise_w": 102978.091,
        "tec_total_electric_power_w": 5692.628,
        "q_radiator_total_w": 102977.955,
        "radiator_tube_external_heat_w": 1200.0,
        "wire_resistance_ohm": np.array([0.001, 0.001]),
        "tec_total_current_a": 200.0,
    }
    diag = startup_energy_residual_diagnostics(record)
    wire_loss = estimate_tec_wire_joule_loss_w(record)
    assert_close(diag["tec_wire_joule_loss_w"], wire_loss)
    assert_close(diag["corrected_core_energy_residual_w"], 110000.0 - 102978.091 - 5692.628 - wire_loss)
    assert_close(diag["corrected_loop_energy_residual_w"], 110000.0 + 1200.0 - 102977.955 - 5692.628 - wire_loss)


def test_restart_loaded_fixed_u_tec_state_flags_skip_startup_setup():
    class RestartCore:
        enable_tec_coupled = True

        def __init__(self) -> None:
            self.calls = []
            self.tfes = {}

        def get_tec_circuit_global_results(self):
            return {"main": {"mode": "fixed_u"}}

        def setup_tec_circuit(self, *args, **kwargs) -> None:
            self.calls.append((args, kwargs))

        def set_thermo_update_time(self, _value) -> None:
            pass

        def update_neutronic_power(self, *args, **kwargs) -> None:
            pass

    class FakeController:
        def evaluate(self, *_args, **_kwargs):
            cmd = SimpleNamespace(
                tec_enabled=True,
                thermal_power_w=1.0,
                fission_power_w=0.0,
                decay_power_w=0.0,
                tec_gap_h_eq_w_m2_k=29.0,
            )
            return cmd

    class FakeConnector:
        def __init__(self, t):
            self.T = t

    class FakeSystem:
        def __init__(self, time):
            self.global_time = time

    args = SimpleNamespace(
        startup_main_tec_initial_mode="fixed_r",
        startup_main_tec_load_resistance_ohm=0.1,
        startup_main_tec_load_resistance_scope="total",
        startup_main_tec_switch_voltage_v=27.2,
        startup_main_tec_i_guess_a=150.0,
        target_voltage=27.2,
        thermo_update_interval=0.5,
    )

    core = RestartCore()
    build = {
        "core": core,
        "system": FakeSystem(12.0),
        "startup_controller": FakeController(),
        "core_inlet_connector": FakeConnector(710.0),
        "tec_has_been_enabled": False,
    }

    old_get_wire_resistance = run_case_module.get_wire_resistance
    run_case_module.get_wire_resistance = lambda _core: [0.001, 0.002, 0.0004, 0.0002]
    try:
        recover_startup_tec_state_from_restart(build, core)
        prepare_startup_step(build, args)
    finally:
        run_case_module.get_wire_resistance = old_get_wire_resistance

    assert build["tec_has_been_enabled"] is True
    assert build["startup_main_tec_switched_to_fixed_voltage"] is True
    assert build["wire_resistance_ohm"] == [0.001, 0.002, 0.0004, 0.0002]
    assert core.calls == []


if __name__ == "__main__":
    test_titam_startup_schedule_matches_document_milestones()
    test_shield_jettison_is_temperature_triggered_and_latched()
    test_default_cesium_gap_h_eq_matches_v7_steady_value()
    test_v13_start_runner_default_cesium_gap_h_eq_matches_v7_steady_value()
    test_tfe_ignition_immediately_sets_cesium_gap_and_enables_fixed_r_tec()
    test_apply_tec_gap_h_eq_updates_all_tfe_gap_couplers()
    test_reset_fluid_temperatures_updates_volume_and_network_state()
    test_reset_solid_temperatures_updates_solid_and_boundary_caches()
    test_matrix_column_heat_source_broadcasts_one_tube_flux_to_axial_nodes()
    test_matrix_column_line_heat_source_converts_w_per_m_to_equivalent_flux_density()
    test_fixed_startup_power_overrides_ramp_schedule()
    test_shield_qsss_from_matrix_maps_six_side_partitions_to_eight_nodes()
    test_startup_main_tec_uses_fixed_resistance_then_switches_to_fixed_voltage()
    test_startup_tec_global_diagnostics_reports_convergence_and_finiteness()
    test_estimate_tec_wire_joule_loss_w_uses_current_and_wire_series()
    test_startup_energy_residual_diagnostics_includes_external_heat_and_wire_loss()
    test_restart_loaded_fixed_u_tec_state_flags_skip_startup_setup()
    test_startup_energy_residual_diagnostics_tracks_storage_and_radiator_balance()
    test_startup_energy_residual_diagnostics_treats_missing_electric_as_zero()
    test_startup_main_tec_load_resistance_scope_can_use_per_tec_value()
    test_startup_controller_can_seed_existing_cesium_conditioning_state()
    print("V13 startup control checks passed.")
