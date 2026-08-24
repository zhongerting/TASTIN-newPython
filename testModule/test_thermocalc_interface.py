import os
import sys
import multiprocessing as mp

import numpy as np


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel


N_NODES = 5


def _run_low_temperature_fixed_voltage_case(queue, disable_guard=False):
    os.environ.pop("THERMOCALC_ENABLE_LOOKUP", None)
    if disable_guard:
        os.environ["THERMOCALC_DISABLE_ZERO_EMISSION_GUARD"] = "1"
    else:
        os.environ.pop("THERMOCALC_DISABLE_ZERO_EMISSION_GUARD", None)
    model = ThermoCalcModel(n_elements=1, n_nodes=37)
    model.setup_circuit_mode("fixed_u", 0.8, I_guess=1.0)
    model.set_temperatures(
        np.full((1, 37), 800.0),
        np.full((1, 37), 600.0),
    )
    model.set_tcs(np.full((1, 37), 520.0))
    elapsed_ms = model.calculate(verbose=False)
    global_results = model.get_global_results()
    tec_results = model.get_tec_results(0)
    queue.put(
        {
            "elapsed_ms": float(elapsed_ms),
            "global": global_results,
            "J": np.asarray(tec_results["J"], dtype=float),
            "joulePowerE": np.asarray(tec_results["joulePowerE"], dtype=float),
            "joulePowerC": np.asarray(tec_results["joulePowerC"], dtype=float),
        }
    )

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
    global_results = model.get_global_results()
    assert global_results["converged"] is True
    assert abs(float(global_results["Uout"]) - 0.8) <= 1.0e-9
    assert int(global_results["iteration_count"]) < 50
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


def test_fixed_current_rejects_invalid_target():
    model = ThermoCalcModel(n_elements=1, n_nodes=N_NODES)
    _expect_value_error(
        lambda: model.setup_circuit_mode("fixed_i", -1.0),
        "fixed_i accepted a negative target current.",
    )


def test_low_temperature_fixed_voltage_auto_skips_zero_emission_case():
    queue = mp.Queue()
    proc = mp.Process(target=_run_low_temperature_fixed_voltage_case, args=(queue,))
    proc.start()
    proc.join(2.0)
    if proc.is_alive():
        proc.terminate()
        proc.join(5.0)
        raise AssertionError("Low-temperature fixed-voltage ThermoCalc calculation did not return.")
    assert proc.exitcode == 0
    result = queue.get_nowait()
    global_results = result["global"]
    assert global_results["zero_emission_skipped"] is True
    assert global_results["converged"] is True
    assert global_results["Iout"] == 0.0
    assert global_results["Uout"] == 0.8
    assert np.all(result["J"] == 0.0)
    assert np.all(result["joulePowerE"] == 0.0)
    assert np.all(result["joulePowerC"] == 0.0)


def test_low_temperature_fixed_voltage_cpp_iteration_returns_when_guard_disabled():
    queue = mp.Queue()
    proc = mp.Process(target=_run_low_temperature_fixed_voltage_case, args=(queue, True))
    proc.start()
    proc.join(3.0)
    if proc.is_alive():
        proc.terminate()
        proc.join(5.0)
        raise AssertionError("C++ low-temperature fixed-voltage ThermoCalc calculation did not return.")
    assert proc.exitcode == 0
    result = queue.get_nowait()
    global_results = result["global"]
    assert global_results["zero_emission_skipped"] is False
    assert np.isfinite(global_results["Iout"])
    assert np.isfinite(global_results["Uout"])
    assert np.all(np.isfinite(result["J"]))
    assert np.all(np.isfinite(result["joulePowerE"]))
    assert np.all(np.isfinite(result["joulePowerC"]))

def test_fixed_resistance_zero_output_clears_failed_node_state():
    from types import SimpleNamespace

    fields = (
        "J", "V", "UE", "UC", "IEsecSingle", "ICsecSingle",
        "phiE", "phiC", "Vd", "joulePowerE", "joulePowerC",
    )
    tec = SimpleNamespace(**{
        field: np.full(N_NODES, np.nan) for field in fields
    })
    tec.I = tec.U = 1.0
    for name in (
        "terminalPointUE1", "terminalPointUE2",
        "terminalPointUC1", "terminalPointUC2",
    ):
        setattr(tec, name, np.nan)
    circuit = SimpleNamespace(
        isFixedR=True, isFixedU=False, isParallelFixedU=False,
        isFixedI=False, isParallelFixedI=False, isParallelLoadCurve=False,
        Iout=0.0, Uout=0.0, converged=False, iterationCount=7,
        TECs=[tec], calc=lambda: None,
    )
    model = ThermoCalcModel(n_elements=1, n_nodes=N_NODES)
    model._circuit = circuit
    model._should_skip_zero_emission = lambda: False
    model.calculate(verbose=False)

    assert circuit.Iout == 0.0
    assert circuit.Uout == 0.0
    assert circuit.converged is False
    assert circuit.iterationCount == 7
    assert model._zero_emission_skipped is False
    for field in fields:
        assert np.all(np.asarray(getattr(tec, field)) == 0.0)


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
    test_fixed_current_rejects_invalid_target()
    test_low_temperature_fixed_voltage_auto_skips_zero_emission_case()
    test_low_temperature_fixed_voltage_cpp_iteration_returns_when_guard_disabled()
    test_fixed_resistance_zero_output_clears_failed_node_state()
    test_uniform_and_nonuniform_runtime_interfaces()
    print("ThermoCalc interface checks passed.")
