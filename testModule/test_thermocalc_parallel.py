import os
import sys
from typing import List

import numpy as np


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

BUILD_PYD_DIR = os.path.join(ROOT_DIR, "ThermoCalc", "build_cp312", "Release")
if os.path.exists(BUILD_PYD_DIR) and BUILD_PYD_DIR not in sys.path:
    sys.path.insert(0, BUILD_PYD_DIR)


from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel


N_NODES = 12
N_ELEMENTS = 2

RTOL_BALANCE = 0.05
ATOL_BALANCE = 1.0e-8
RTOL_SUM = 0.02


def _base_temperatures() -> tuple[np.ndarray, np.ndarray]:
    emitter = np.linspace(1450.0, 1700.0, N_NODES, dtype=float)[np.newaxis, :]
    collector = np.linspace(750.0, 900.0, N_NODES, dtype=float)[np.newaxis, :]
    return np.repeat(emitter, N_ELEMENTS, axis=0), np.repeat(collector, N_ELEMENTS, axis=0)


def _build_parallel_model(mode: str, target_value: float, i_guess: float = 150.0) -> ThermoCalcModel:
    model = ThermoCalcModel(n_elements=N_ELEMENTS, n_nodes=N_NODES)
    Te, Tc = _base_temperatures()
    model.set_temperatures(Te, Tc)
    model.set_tcs(610.0)

    try:
        model.setup_circuit_mode(mode, target_value, I_guess=i_guess)
    except ValueError as exc:
        raise AssertionError(
            f"Pending implementation: parallel mode '{mode}' is not exposed in production API yet. "
            f"Observed setup error: {exc}"
        ) from exc

    return model


def _collect_tec_results(model: ThermoCalcModel) -> tuple[dict, List[dict]]:
    global_results = model.get_global_results()
    if global_results is None:
        raise AssertionError("ThermoCalcModel did not return global results after calculate().")

    tec_results = []
    for i in range(N_ELEMENTS):
        item = model.get_tec_results(i)
        if item is None:
            raise AssertionError(f"Missing TEC result payload for element {i}.")
        tec_results.append(item)

    return global_results, tec_results


def _read_currents_and_voltages(tec_results: List[dict]) -> tuple[np.ndarray, np.ndarray]:
    i_branch = np.asarray([float(r.get("I", np.nan)) for r in tec_results], dtype=float)
    u_branch = np.asarray([float(r.get("U", np.nan)) for r in tec_results], dtype=float)
    return i_branch, u_branch


def test_parallel_fixed_u_balanced_branches():
    model = _build_parallel_model("parallel_fixed_u", target_value=0.5, i_guess=180.0)
    model.calculate(verbose=False)

    global_results, tec_results = _collect_tec_results(model)
    i_branch, u_branch = _read_currents_and_voltages(tec_results)

    assert np.all(np.isfinite(i_branch))
    assert np.all(np.isfinite(u_branch))
    assert np.all(i_branch > 0.0)
    assert np.allclose(np.mean(i_branch), i_branch, rtol=RTOL_BALANCE, atol=ATOL_BALANCE)
    assert np.allclose(u_branch, global_results["Uout"], rtol=RTOL_BALANCE, atol=ATOL_BALANCE)
    assert np.isclose(i_branch.sum(), global_results["Iout"], rtol=RTOL_SUM, atol=ATOL_BALANCE)


def test_parallel_fixed_i_respects_target_current():
    target_i = 100.0
    model = _build_parallel_model("parallel_fixed_i", target_value=target_i, i_guess=target_i)
    model.calculate(verbose=False)

    global_results, tec_results = _collect_tec_results(model)
    i_branch, u_branch = _read_currents_and_voltages(tec_results)

    assert np.all(np.isfinite(i_branch))
    assert np.all(np.isfinite(u_branch))
    assert np.all(i_branch > 0.0)
    assert np.isclose(global_results["Iout"], target_i, rtol=0.05, atol=ATOL_BALANCE)
    assert np.allclose(i_branch.sum(), target_i, rtol=RTOL_SUM, atol=ATOL_BALANCE)
    assert np.allclose(np.mean(i_branch), i_branch, rtol=RTOL_BALANCE, atol=ATOL_BALANCE)
    assert np.allclose(u_branch, global_results["Uout"], rtol=RTOL_BALANCE, atol=ATOL_BALANCE)


def test_parallel_load_curve_follows_u_i_relation():
    load_resistance = 0.002
    model = _build_parallel_model("parallel_load_curve", target_value=load_resistance, i_guess=150.0)
    i_axis = np.array([0.0, 400.0], dtype=float)
    model.set_load_curve(i_axis, load_resistance * i_axis)
    model.calculate(verbose=False)

    global_results, tec_results = _collect_tec_results(model)
    i_branch, u_branch = _read_currents_and_voltages(tec_results)
    i_total = float(global_results["Iout"])
    u_total = float(global_results["Uout"])

    assert np.all(np.isfinite(i_branch))
    assert np.all(np.isfinite(u_branch))
    assert np.isfinite(i_total)
    assert np.isfinite(u_total)

    predicted_u = float(load_resistance) * i_total
    assert np.isclose(u_total, predicted_u, rtol=0.05, atol=ATOL_BALANCE)
    assert np.allclose(np.mean(i_branch), i_branch, rtol=RTOL_BALANCE, atol=ATOL_BALANCE)


def test_parallel_runtime_updates_for_tcs_and_temperatures():
    model = _build_parallel_model("parallel_fixed_u", target_value=0.5, i_guess=170.0)
    model.calculate(verbose=False)
    baseline_global, baseline_tec_results = _collect_tec_results(model)

    model.set_tcs(np.full((N_ELEMENTS, N_NODES), 700.0, dtype=float))

    te_update, tc_update = _base_temperatures()
    te_update = te_update + 40.0
    tc_update = tc_update + 30.0
    model.set_temperatures(te_update, tc_update)

    model.calculate(verbose=False)
    updated_global, updated_tec_results = _collect_tec_results(model)
    updated_branch_temps = np.asarray(updated_tec_results[0]["TE"], dtype=float)

    assert np.all(np.isfinite(updated_branch_temps))
    assert not np.array_equal(
        np.asarray(baseline_tec_results[0]["TE"], dtype=float),
        updated_branch_temps,
    )
    assert np.isfinite(updated_global["Iout"])
    assert np.isfinite(baseline_global["Iout"])
    assert np.isfinite(updated_global["Uout"])
    if updated_global.get("Rload") is not None:
        assert np.isfinite(updated_global["Rload"])


if __name__ == "__main__":
    test_parallel_fixed_u_balanced_branches()
    test_parallel_fixed_i_respects_target_current()
    test_parallel_load_curve_follows_u_i_relation()
    test_parallel_runtime_updates_for_tcs_and_temperatures()
    print("ThermoCalc parallel circuit checks passed.")
