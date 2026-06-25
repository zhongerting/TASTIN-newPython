import os
import sys
import types
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
os.chdir(ROOT_DIR)

from Components.ReactorCore import ReactorCore, TecCircuitGroup
from testModule.run_v13_caseA_closed_loop import load_tec_load_curve
from testModule.test_core_assemble_v7_caseA import _case_a_electric_diagnostics


class DummyThermoCalc:
    def __init__(self):
        self.mode_calls = []
        self.load_curve_calls = []

    def setup_circuit_mode(self, mode_str, target_value, I_guess=150.0):
        self.mode_calls.append((mode_str, float(target_value), float(I_guess)))

    def set_load_curve(self, current_a, voltage_v):
        self.load_curve_calls.append((
            np.asarray(current_a, dtype=float),
            np.asarray(voltage_v, dtype=float),
        ))


def _core_with_dummy_thermo():
    core = ReactorCore.__new__(ReactorCore)
    core.thermo_calc = DummyThermoCalc()
    core.tec_topology = "series"
    core.tec_circuit_mode = "fixed_u"
    return core


def test_series_fixed_u_keeps_legacy_mode():
    core = _core_with_dummy_thermo()
    core.setup_tec_circuit("fixed_u", 27.2, I_guess=123.0)
    assert core.tec_topology == "series"
    assert core.tec_circuit_mode == "fixed_u"
    assert core.thermo_calc.mode_calls == [("fixed_u", 27.2, 123.0)]


def test_parallel_fixed_u_maps_to_cpp_parallel_mode():
    core = _core_with_dummy_thermo()
    core.setup_tec_circuit("fixed_u", 0.5, I_guess=12.0, topology="parallel")
    assert core.tec_topology == "parallel"
    assert core.tec_circuit_mode == "fixed_u"
    assert core.thermo_calc.mode_calls == [("parallel_fixed_u", 0.5, 12.0)]


def test_parallel_fixed_i_maps_to_cpp_parallel_mode():
    core = _core_with_dummy_thermo()
    core.setup_tec_circuit("fixed_i", 200.0, I_guess=200.0, topology="parallel")
    assert core.thermo_calc.mode_calls == [("parallel_fixed_i", 200.0, 200.0)]


def test_parallel_fixed_r_uses_load_curve_solver():
    core = _core_with_dummy_thermo()
    core.setup_tec_circuit("fixed_r", 0.002, I_guess=100.0, topology="parallel")
    assert core.tec_circuit_mode == "fixed_r"
    assert core.thermo_calc.mode_calls == [("parallel_load_curve", 0.002, 100.0)]


def test_parallel_load_curve_is_forwarded_before_mode_setup():
    core = _core_with_dummy_thermo()
    current = np.array([0.0, 100.0, 200.0])
    voltage = np.array([0.0, 10.0, 30.0])
    core.setup_tec_circuit(
        "load_curve",
        0.0,
        I_guess=100.0,
        topology="parallel",
        load_curve=(current, voltage),
    )
    assert core.tec_circuit_mode == "load_curve"
    assert core.thermo_calc.mode_calls == [("parallel_load_curve", 0.0, 100.0)]
    np.testing.assert_allclose(core.thermo_calc.load_curve_calls[0][0], current)
    np.testing.assert_allclose(core.thermo_calc.load_curve_calls[0][1], voltage)


def test_series_fixed_i_still_reaches_wrapper_rejection_mode():
    mode, topology, logical = ReactorCore._resolve_tec_mode("fixed_i", "series")
    assert (mode, topology, logical) == ("fixed_i", "series", "fixed_i")


def test_invalid_topology_raises():
    try:
        ReactorCore._resolve_tec_mode("fixed_u", "mixed")
    except ValueError as exc:
        assert "topology" in str(exc)
    else:
        raise AssertionError("invalid topology did not raise")


def test_csv_load_curve_reader(tmp_path):
    path = tmp_path / "load_curve.csv"
    path.write_text("current_a,voltage_v\n0,0\n10,2.5\n20,6\n", encoding="utf-8")
    current, voltage = load_tec_load_curve(str(path))
    np.testing.assert_allclose(current, [0.0, 10.0, 20.0])
    np.testing.assert_allclose(voltage, [0.0, 2.5, 6.0])


def test_reserved_parallel_group_uses_ring3_open_only():
    core = ReactorCore.__new__(ReactorCore)
    core.enable_tec_coupled = True
    core.tfes = {"Ring3_TEC": object(), "Ring3_Open": object()}
    core.tfe_multipliers = {"Ring3_TEC": 15, "Ring3_Open": 3}
    core.tec_multipliers = {"Ring3_TEC": 15, "Ring3_Open": 0}
    core.tec_circuit_groups = {}
    core._last_thermo_update_time = 10.0
    created = {}

    def fake_create(self, name, multipliers, topology="series", circuit_mode="fixed_u"):
        group = TecCircuitGroup(
            name=name,
            multipliers={key: int(value) for key, value in multipliers.items() if int(value) > 0},
            thermo_calc=DummyThermoCalc(),
            total_virtual_elements=sum(int(value) for value in multipliers.values() if int(value) > 0),
            topology=topology,
            circuit_mode=circuit_mode,
            last_update_time=self._last_thermo_update_time,
        )
        created[name] = group
        return group

    core._create_tec_circuit_group = types.MethodType(fake_create, core)
    core.setup_reserved_parallel_tec_circuit(mode_str="fixed_u", target_value=0.8, I_guess=6000.0)

    group = core.tec_circuit_groups["reserved_parallel"]
    assert group.multipliers == {"Ring3_Open": 3}
    assert group.total_virtual_elements == 3
    assert group.topology == "parallel"
    assert group.circuit_mode == "fixed_u"
    assert group.thermo_calc.mode_calls == [("parallel_fixed_u", 0.8, 6000.0)]
    assert core.reserved_parallel_tec_enabled is True
    assert created["reserved_parallel"] is group


def test_reserved_parallel_rejects_overlap_with_main_group():
    core = ReactorCore.__new__(ReactorCore)
    core.tfes = {"Ring3_TEC": object()}
    core.tfe_multipliers = {"Ring3_TEC": 15}
    core.tec_multipliers = {"Ring3_TEC": 15}
    try:
        core._validate_aux_tec_multipliers({"Ring3_TEC": 1})
    except ValueError as exc:
        assert "exceeding thermal multiplier" in str(exc)
    else:
        raise AssertionError("overlapping reserved TEC multipliers did not raise")


def test_connected_tec_tfe_names_includes_reserved_parallel_group():
    core = ReactorCore.__new__(ReactorCore)
    core.enable_tec_coupled = True
    core.tec_circuit_groups = {
        "main": TecCircuitGroup(
            name="main",
            multipliers={"Ring3_TEC": 15},
            thermo_calc=DummyThermoCalc(),
            total_virtual_elements=15,
        ),
        "reserved_parallel": TecCircuitGroup(
            name="reserved_parallel",
            multipliers={"Ring3_Open": 3},
            thermo_calc=DummyThermoCalc(),
            total_virtual_elements=3,
            topology="parallel",
        ),
    }
    assert core._connected_tec_tfe_names() == {"Ring3_TEC", "Ring3_Open"}


def test_case_a_electric_diagnostics_sums_main_and_reserved_power():
    class DummyCore:
        def get_tec_circuit_global_results(self):
            return {
                "main": {"Uout": 27.2, "Iout": 200.0, "Rload": 0.136},
                "reserved_parallel": {"Uout": 0.8, "Iout": 6000.0, "Rload": 0.0},
            }

    diag = _case_a_electric_diagnostics(DummyCore())
    assert diag["tec_main_electric_power_w"] == 5440.0
    assert diag["tec_reserved_parallel_electric_power_w"] == 4800.0
    assert diag["tec_total_electric_power_w"] == 10240.0
    assert diag["tec_total_voltage_v"] == 27.2
    assert diag["tec_total_current_a"] == 200.0


if __name__ == "__main__":
    test_series_fixed_u_keeps_legacy_mode()
    test_parallel_fixed_u_maps_to_cpp_parallel_mode()
    test_parallel_fixed_i_maps_to_cpp_parallel_mode()
    test_parallel_fixed_r_uses_load_curve_solver()
    test_parallel_load_curve_is_forwarded_before_mode_setup()
    test_series_fixed_i_still_reaches_wrapper_rejection_mode()
    test_invalid_topology_raises()
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        test_csv_load_curve_reader(Path(tmp))
    test_reserved_parallel_group_uses_ring3_open_only()
    test_reserved_parallel_rejects_overlap_with_main_group()
    test_connected_tec_tfe_names_includes_reserved_parallel_group()
    test_case_a_electric_diagnostics_sums_main_and_reserved_power()
    print("ReactorCore TEC topology checks passed.")
