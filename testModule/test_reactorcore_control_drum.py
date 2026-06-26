import math
import os
import sys

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from Components.ReactorCore import (
    ControlDrumReactivityModel,
    ReactivityFeedbackResult,
    ReactivityFeedbackTemperatureSummary,
    ReactorCore,
)


BETA_TOTAL = 0.0079321


class FakePointReactor:
    def __init__(self):
        self.beta_total = BETA_TOTAL
        self.fission_power = 1.0
        self.decay_power = 0.0
        self.total_power = 1.0
        self.last_reactivity_control = None
        self.last_reactivity_feedback = None

    def step(self, dt, reactivity_control, reactivity_feedback):
        self.last_reactivity_control = float(reactivity_control)
        self.last_reactivity_feedback = float(reactivity_feedback)
        return True


def make_minimal_core(model=None):
    core = ReactorCore.__new__(ReactorCore)
    core.point_reactor = FakePointReactor()
    core.feedback_reference_result = ReactivityFeedbackResult(
        temperatures=ReactivityFeedbackTemperatureSummary(),
        total=0.0,
    )
    core.last_reactivity_control = 0.0
    core.last_effective_reactivity_feedback = 0.0
    core.control_drum_reactivity_model = model or ControlDrumReactivityModel()
    core.compute_reactivity_feedback = lambda: ReactivityFeedbackResult(
        temperatures=ReactivityFeedbackTemperatureSummary(),
        total=0.0,
    )
    core.update_neutronic_power = lambda **kwargs: None
    return core


def assert_close(actual, expected, tol=1.0e-10):
    assert math.isclose(float(actual), float(expected), rel_tol=tol, abs_tol=tol), (
        actual,
        expected,
    )


def test_control_drum_polynomial_endpoints_and_cold_reference():
    model = ControlDrumReactivityModel()

    assert_close(model.polynomial_reactivity_dollars(0.0), -4.0)
    assert_close(model.polynomial_reactivity_dollars(180.0), 0.79592288, tol=1.0e-8)
    assert_close(model.delta_reactivity_dollars(), 0.0)

    model.theta_deg = 180.0
    assert_close(model.delta_reactivity_dollars(), 4.79592288, tol=1.0e-8)

    expected_cold_rho = (0.952 - 1.0) / 0.952
    expected_cold_dollars = expected_cold_rho / BETA_TOTAL
    assert_close(model.cold_reference_reactivity_dollars(BETA_TOTAL), expected_cold_dollars)


def test_control_drum_model_is_disabled_by_default():
    model = ControlDrumReactivityModel(theta_deg=180.0)
    assert model.enabled is False
    assert_close(model.reactivity(BETA_TOTAL), 0.0)

    core = make_minimal_core(model=model)
    core.advance_neutronics(dt=0.1, reactivity_control=1.23e-4)

    assert_close(core.point_reactor.last_reactivity_control, 1.23e-4)
    assert_close(core.point_reactor.last_reactivity_feedback, 0.0)
    assert_close(core.get_control_drum_reactivity_dollars(), 0.0)


def test_enabled_control_drum_adds_dimensional_reactivity_to_control():
    model = ControlDrumReactivityModel(enabled=True, theta_deg=180.0)
    core = make_minimal_core(model=model)

    external_control = 1.0e-5
    core.advance_neutronics(dt=0.1, reactivity_control=external_control)

    expected = external_control + model.reactivity(BETA_TOTAL)
    assert_close(core.point_reactor.last_reactivity_control, expected)
    assert_close(core.last_reactivity_control, expected)


def test_reactorcore_configure_and_diagnostics_api():
    core = make_minimal_core()

    model = core.configure_control_drum_reactivity(
        enabled=True,
        theta_deg=180.0,
        reference_theta_deg=0.0,
        cold_reference_keff=0.952,
    )

    assert model is core.control_drum_reactivity_model
    diagnostics = core.get_control_drum_diagnostics()
    assert diagnostics["control_drum_enabled"] is True
    assert_close(diagnostics["control_drum_theta_deg"], 180.0)
    assert_close(diagnostics["control_drum_delta_reactivity_dollars"], 4.79592288, tol=1.0e-8)
    assert_close(
        diagnostics["control_drum_reactivity"],
        diagnostics["control_drum_total_reactivity_dollars"] * BETA_TOTAL,
    )


def test_control_drum_model_state_roundtrip():
    model = ControlDrumReactivityModel(
        enabled=True,
        theta_deg=180.0,
        reference_theta_deg=0.0,
        cold_reference_keff=0.952,
    )

    restored = ControlDrumReactivityModel.from_state(model.to_state())

    assert restored.enabled is True
    assert_close(restored.theta_deg, 180.0)
    assert_close(restored.reference_theta_deg, 0.0)
    assert_close(restored.cold_reference_keff, 0.952)
    assert_close(restored.reactivity(BETA_TOTAL), model.reactivity(BETA_TOTAL))


def test_control_drum_angle_clamping():
    model = ControlDrumReactivityModel(theta_deg=181.0)
    assert_close(model.polynomial_reactivity_dollars(), model.polynomial_reactivity_dollars(180.0))

    model.set_angle(-5.0)
    assert_close(model.theta_deg, 0.0)


def test_old_restart_without_control_drum_keys_keeps_default_model():
    core = make_minimal_core()
    state = {
        "Core/enable_tec_coupled": [False],
        "Core/tec_topology": ["series"],
        "Core/tec_circuit_mode": ["fixed_u"],
        "Core/reserved_parallel_tec_enabled": [False],
        "Core/_last_thermo_update_time": [0.0],
        "Core/last_total_core_power": [1.0],
        "Core/last_effective_reactivity_feedback": [0.0],
        "Core/last_reactivity_control": [0.0],
        "Core/feedback_reference/total": [0.0],
    }
    core.enable_tec_coupled = True
    core.tec_topology = "series"
    core.tec_circuit_mode = "fixed_u"
    core.reserved_parallel_tec_enabled = False
    core._last_thermo_update_time = 0.0
    core.tfes = {}
    core.point_reactor = None

    core.load_state_dict(state, "Core")

    assert core.control_drum_reactivity_model.enabled is False
    assert_close(core.control_drum_reactivity_model.theta_deg, 0.0)


if __name__ == "__main__":
    test_control_drum_polynomial_endpoints_and_cold_reference()
    test_control_drum_model_is_disabled_by_default()
    test_enabled_control_drum_adds_dimensional_reactivity_to_control()
    test_reactorcore_configure_and_diagnostics_api()
    test_control_drum_model_state_roundtrip()
    test_control_drum_angle_clamping()
    test_old_restart_without_control_drum_keys_keeps_default_model()
    print("ReactorCore control drum checks passed.")
