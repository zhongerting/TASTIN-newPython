import os
import sys

import numpy as np


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel


N_NODES = 5


def _expect_value_error(callback, message):
    try:
        callback()
    except ValueError:
        return
    raise AssertionError(message)


def _make_model(*, nonuniform: bool):
    model = ThermoCalcModel(n_elements=1, n_nodes=N_NODES)
    model.setup_circuit_mode("fixed_u", 0.8, I_guess=20.0)
    model.set_temperatures(
        np.linspace(1500.0, 1700.0, N_NODES)[np.newaxis, :],
        np.linspace(800.0, 900.0, N_NODES)[np.newaxis, :],
    )
    model.set_tcs(610.0)
    if nonuniform:
        model._input_data.dlE = np.array([[0.03, 0.05, 0.08, 0.10, 0.12]])
        model._input_data.dlC = np.array([[0.03, 0.05, 0.08, 0.10, 0.12]])
        model._input_data.sideAreaE = np.array([[0.0003, 0.0005, 0.0008, 0.0010, 0.0012]])
        model._input_data.sideAreaC = np.array([[0.0004, 0.0006, 0.0009, 0.0011, 0.0013]])
        model._input_data.resistanceWire = np.array([[1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6]])
    model.calculate(verbose=False)
    return model


def _reconstruct_fvm_joule(potential, resistivity, lengths, cross_area, left_terminal, right_terminal):
    conductance = cross_area / (0.5 * (resistivity[:-1] + resistivity[1:]) * lengths[:-1])
    node_power = np.zeros_like(potential)
    face_power = conductance * np.diff(potential) ** 2
    node_power[:-1] += 0.5 * face_power
    node_power[1:] += 0.5 * face_power
    node_power[0] += cross_area / (resistivity[0] * lengths[0] * 0.5) * (potential[0] - left_terminal) ** 2
    node_power[-1] += (
        cross_area
        / (resistivity[-1] * lengths[-1] * 0.5)
        * (potential[-1] - right_terminal) ** 2
    )
    return node_power


def _assert_fvm_joule_outputs(model, results):
    for electrode in ("E", "C"):
        values = np.asarray(results[f"joulePower{electrode}"], dtype=float)
        assert values.shape == (N_NODES,)
        assert np.all(np.isfinite(values))
        assert np.all(values >= 0.0)
        expected = _reconstruct_fvm_joule(
            np.asarray(results[f"U{electrode}"], dtype=float),
            np.asarray(results[f"rho{electrode}"], dtype=float),
            np.asarray(getattr(model._input_data, f"dl{electrode}")[0], dtype=float),
            float(getattr(model._input_data, f"crossArea{electrode}")[0]),
            float(results[f"terminalPointU{electrode}1"]),
            float(results[f"terminalPointU{electrode}2"]),
        )
        assert np.allclose(values, expected, rtol=1.0e-12, atol=1.0e-12)


def test_input_shape_validation():
    model = ThermoCalcModel(n_elements=1, n_nodes=N_NODES)
    model._input_data.sideAreaE = np.ones(1)
    _expect_value_error(
        model.build,
        "create_circuit() accepted a one-dimensional sideAreaE array.",
    )


def test_fixed_current_is_explicitly_rejected():
    model = ThermoCalcModel(n_elements=1, n_nodes=N_NODES)
    _expect_value_error(
        lambda: model.setup_circuit_mode("fixed_i", 10.0),
        "fixed_I unexpectedly reached an unbound C++ enum.",
    )


def test_uniform_and_nonuniform_runtime_interfaces():
    uniform = _make_model(nonuniform=False)
    uniform_results = uniform.get_tec_results(0)
    for field in ("phiE", "phiC", "Vd"):
        values = np.asarray(uniform_results[field], dtype=float)
        assert values.shape == (N_NODES,)
        assert np.all(np.isfinite(values))
    _assert_fvm_joule_outputs(uniform, uniform_results)

    updated_tcs = np.linspace(620.0, 640.0, N_NODES)[np.newaxis, :]
    uniform.set_tcs(updated_tcs)
    assert np.allclose(np.asarray(uniform._circuit.TECs[0].Tcs), updated_tcs[0])

    updated_te = np.linspace(1520.0, 1720.0, N_NODES)[np.newaxis, :]
    updated_tc = np.linspace(810.0, 910.0, N_NODES)[np.newaxis, :]
    uniform.set_temperatures(updated_te, updated_tc)
    assert np.allclose(np.asarray(uniform._circuit.TECs[0].Temitter), updated_te[0])
    assert np.allclose(np.asarray(uniform._circuit.TECs[0].Tcollector), updated_tc[0])

    nonuniform = _make_model(nonuniform=True)
    nonuniform_results = nonuniform.get_tec_results(0)
    for field in ("J", "UE", "UC", "phiE", "phiC", "Vd"):
        values = np.asarray(nonuniform_results[field], dtype=float)
        assert values.shape == (N_NODES,)
        assert np.all(np.isfinite(values))
    _assert_fvm_joule_outputs(nonuniform, nonuniform_results)
    assert not np.allclose(nonuniform_results["UE"], uniform_results["UE"])


if __name__ == "__main__":
    test_input_shape_validation()
    test_fixed_current_is_explicitly_rejected()
    test_uniform_and_nonuniform_runtime_interfaces()
    print("ThermoCalc interface checks passed.")
