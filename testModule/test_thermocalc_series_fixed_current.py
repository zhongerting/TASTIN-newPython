import os
import sys

import numpy as np


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel


N_ELEMENTS = 2
N_NODES = 12


def _make_model(mode: str, target: float, i_guess: float = 100.0, n_elements: int = N_ELEMENTS) -> ThermoCalcModel:
    model = ThermoCalcModel(n_elements=n_elements, n_nodes=N_NODES, enable_lookup=False)
    emitter = np.linspace(1500.0, 1750.0, N_NODES)[np.newaxis, :]
    collector = np.linspace(750.0, 900.0, N_NODES)[np.newaxis, :]
    model.set_temperatures(
        np.repeat(emitter, n_elements, axis=0),
        np.repeat(collector, n_elements, axis=0),
    )
    model.set_tcs(600.0)
    model.setup_circuit_mode(mode, target, I_guess=i_guess)
    return model


def _assert_finite_tec_state(model: ThermoCalcModel) -> None:
    for index in range(model.N_elem):
        result = model.get_tec_results(index)
        for field in ("J", "V", "UE", "UC", "joulePowerE", "joulePowerC"):
            assert np.all(np.isfinite(np.asarray(result[field], dtype=float))), field


def test_series_fixed_current_zero_target_returns_finite_open_circuit():
    for n_elements in (1, 2):
        model = _make_model("fixed_i", 0.0, i_guess=0.0, n_elements=n_elements)
        model.calculate(verbose=False)

        result = model.get_global_results()
        assert result["Iout"] == 0.0
        assert np.isfinite(result["Uout"])
        assert result["Uout"] > 0.0
        assert result["converged"] is True
        assert result["iteration_count"] == 1
        _assert_finite_tec_state(model)


def test_series_fixed_current_matches_fixed_voltage_operating_point():
    voltage_model = _make_model("fixed_u", 1.6)
    voltage_model.calculate(verbose=False)
    voltage_result = voltage_model.get_global_results()
    assert voltage_result["converged"] is True

    target_current = float(voltage_result["Iout"])
    current_model = _make_model("fixed_i", target_current, i_guess=target_current)
    current_model.calculate(verbose=False)
    current_result = current_model.get_global_results()

    assert current_result["converged"] is True
    assert np.isclose(current_result["Iout"], target_current, rtol=0.0, atol=1.0e-9)
    assert np.isclose(current_result["Uout"], voltage_result["Uout"], rtol=0.05, atol=0.05)
    assert current_result["Uout"] > 0.0
    _assert_finite_tec_state(current_model)


def test_series_fixed_current_rejects_non_generating_target_and_opens_circuit():
    voltage_model = _make_model("fixed_u", 1.6)
    voltage_model.calculate(verbose=False)
    reference_current = max(float(voltage_model.get_global_results()["Iout"]), 1.0)

    model = _make_model("fixed_i", reference_current * 50.0, i_guess=reference_current)
    model.calculate(verbose=False)
    result = model.get_global_results()

    assert result["Iout"] == 0.0
    assert np.isfinite(result["Uout"])
    assert result["Uout"] > 0.0
    assert result["converged"] is False
    assert result["iteration_count"] == 1
    _assert_finite_tec_state(model)


def test_series_fixed_current_double_failure_returns_finite_zero_output():
    voltage_model = _make_model("fixed_u", 1.6)
    voltage_model.calculate(verbose=False)
    reference_current = max(float(voltage_model.get_global_results()["Iout"]), 1.0)

    model = _make_model("fixed_i", reference_current * 1000.0, i_guess=reference_current)
    model.calculate(verbose=False)
    result = model.get_global_results()

    assert result["Iout"] == 0.0
    assert result["Uout"] == 0.0
    assert result["converged"] is False
    assert result["iteration_count"] == 1
    _assert_finite_tec_state(model)

def test_series_fixed_resistance_matches_fixed_voltage_operating_point():
    voltage_model = _make_model("fixed_u", 1.6)
    voltage_model.calculate(verbose=False)
    voltage_result = voltage_model.get_global_results()
    assert voltage_result["converged"] is True

    target_current = float(voltage_result["Iout"])
    resistance = float(voltage_result["Uout"]) / target_current
    resistance_model = _make_model("fixed_r", resistance, i_guess=target_current)
    resistance_model.build()
    resistance_model._circuit.Iout = 0.0
    resistance_model._circuit.Uout = 0.0
    resistance_model.calculate(verbose=False)
    resistance_result = resistance_model.get_global_results()

    assert resistance_result["converged"] is True
    assert resistance_result["Iout"] > 0.0
    assert resistance_result["Uout"] > 0.0
    assert abs(resistance_result["Uout"] - resistance_result["Iout"] * resistance) <= 1.0e-3
    assert np.isclose(resistance_result["Iout"], target_current, rtol=0.02, atol=0.1)
    _assert_finite_tec_state(resistance_model)

if __name__ == "__main__":
    test_series_fixed_current_zero_target_returns_finite_open_circuit()
    test_series_fixed_current_matches_fixed_voltage_operating_point()
    test_series_fixed_current_rejects_non_generating_target_and_opens_circuit()
    test_series_fixed_current_double_failure_returns_finite_zero_output()
    test_series_fixed_resistance_matches_fixed_voltage_operating_point()
    print("ThermoCalc series fixed-current checks passed.")
