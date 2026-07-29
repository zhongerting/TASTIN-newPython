import unittest
from types import SimpleNamespace

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_low_electric_power_fixed_I.run_v14_low_power_fixed_i import (
    _collect_heat_rejection_diagnostics,
    _refresh_ring_hp_diagnostic_state,
)


class _Boundary:
    def __init__(self, temperature, flux=None):
        self.T_surface = np.asarray(temperature, dtype=float)
        self.current_flux = np.asarray(
            np.zeros_like(self.T_surface) if flux is None else flux, dtype=float)

    def get_coupling_surface_snapshot(self):
        return self.T_surface.copy(), None


class _Wick:
    @staticmethod
    def conductivity_axial(temperature):
        return 0.1 * np.asarray(temperature, dtype=float)


class _HPUnit:
    fin_thickness = 0.2
    fin_height = 2.0
    n_fin_height = 2
    fin_width_array = np.array([1.0])
    last_fin_temperature = np.array([[400.0, 600.0]])

    def __init__(self):
        hp = SimpleNamespace()
        hp.shape_nodes = (2, 3)
        hp.T = np.array([100.0, 200.0, 300.0, 300.0, 400.0, 500.0])
        hp.mesh = SimpleNamespace(
            geom_data=SimpleNamespace(volumes=np.array([1., 1., 1., 3., 3., 3.])))
        hp.n_wick = 1
        hp.wick_mat = _Wick()
        hp._slice_eva = slice(0, 1)
        hp._slice_aba = slice(1, 2)
        hp._slice_con = slice(2, 3)
        hp.boundaries = {
            "outer_eva": _Boundary([300.0], [5.0]),
            "outer_con": _Boundary([500.0]),
        }
        self.hp = hp

    def get_heat_exchange_breakdown(self):
        return {"gross_rejection": np.array([11.0]), "net_rejection": np.array([7.0])}

    def get_heat_rejection_distribution(self):
        return np.array([2.0]), np.array([11.0])

    def get_external_heat_absorption_distribution(self, current_time):
        return np.array([1.0]), np.array([3.0]), np.array([4.0])


class _Ring:
    def __init__(self, multiplier=2.0, wall_flux=(2.0, -1.0)):
        self.hp = _HPUnit()
        self.solid_header = SimpleNamespace(
            T=np.array([700.0, 900.0]),
            mesh=SimpleNamespace(geom_data=SimpleNamespace(volumes=np.array([1.0, 3.0]))),
            boundaries={"right": _Boundary([700.0, 900.0], wall_flux)},
        )
        self.multiplier = multiplier
        self.pre_step_calls = []

    def _iter_present_hp_units_with_multiplier(self):
        yield 0, self.hp, self.multiplier

    def pre_step(self, dt, current_time):
        self.pre_step_calls.append((dt, current_time))


class HeatRejectionDiagnosticsTests(unittest.TestCase):
    def test_uses_physical_volume_weights_and_signed_wall_flux(self):
        ring = _Ring()
        result = _collect_heat_rejection_diagnostics(
            {"ring_hps": [ring]}, current_time_s=25.0,
            external_heat_period_s=100.0, external_heat_time_origin_s=10.0,
        )

        self.assertAlmostEqual(result["hp_evaporator_temperature_mean_K"], 250.0)
        self.assertAlmostEqual(result["hp_adiabatic_temperature_mean_K"], 350.0)
        self.assertAlmostEqual(result["hp_condenser_temperature_mean_K"], 450.0)
        self.assertAlmostEqual(result["hp_wick_axial_conductivity_mean_W_mK"], 20.0)
        self.assertAlmostEqual(result["hp_evaporator_minus_condenser_mean_K"], -200.0)
        self.assertAlmostEqual(result["hp_evaporator_heat_input_W"], 10.0)
        self.assertAlmostEqual(result["radiator_gross_heat_rejection_W"], 26.0)
        self.assertAlmostEqual(result["radiator_external_heat_absorption_W"], 8.0)
        self.assertAlmostEqual(result["radiator_net_heat_rejection_W"], 18.0)
        self.assertAlmostEqual(result["collector_ring_wall_temperature_mean_K"], 850.0)
        self.assertAlmostEqual(result["collector_ring_wall_outward_rejection_W"], -1.0)
        self.assertAlmostEqual(result["radiator_fin_temperature_mean_K"], 500.0)
        self.assertAlmostEqual(result["external_heat_phase_s"], 15.0)
        self.assertAlmostEqual(result["external_heat_phase_fraction"], 0.15)

    def test_refreshes_each_ring_once_at_checkpoint_time(self):
        rings = [_Ring(), _Ring()]
        _refresh_ring_hp_diagnostic_state({"ring_hps": rings}, 19265.0)
        self.assertEqual(rings[0].pre_step_calls, [(0.0, 19265.0)])
        self.assertEqual(rings[1].pre_step_calls, [(0.0, 19265.0)])


if __name__ == "__main__":
    unittest.main()
