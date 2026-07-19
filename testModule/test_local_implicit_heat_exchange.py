import unittest

import numpy as np

from Solvers.Couplers import FluidSolidCouple


class _DummyMaterial:
    def heat_capacity(self, T, P):
        return np.full(np.asarray(T, dtype=float).shape, 900.0, dtype=float)

    def conductivity(self, T, P):
        return np.full(np.asarray(T, dtype=float).shape, 75.0, dtype=float)

    def viscosity(self, T, P):
        return np.full(np.asarray(T, dtype=float).shape, 4.0e-4, dtype=float)

    def prandtl_number(self, T, P):
        return np.full(np.asarray(T, dtype=float).shape, 0.004, dtype=float)


class _DummyFluxBC:
    def __init__(self, n_nodes: int):
        self.q_flux = np.zeros(n_nodes, dtype=float)

    def update_params(self, q_flux):
        self.q_flux = np.asarray(q_flux, dtype=float)


class _DummyResistanceBC:
    def __init__(self, n_nodes: int):
        self.T_ext = np.zeros(n_nodes, dtype=float)
        self.R_ext = np.zeros(n_nodes, dtype=float)


class _DummyBoundaryRegion:
    def __init__(self, n_nodes: int):
        self.shape = (n_nodes,)
        self.T_surface = np.full(n_nodes, 900.0, dtype=float)
        self.current_flux = np.zeros(n_nodes, dtype=float)

    def add_resistance_condition(self, T_ext, R_ext):
        bc = _DummyResistanceBC(self.shape[0])
        bc.T_ext[:] = np.asarray(T_ext, dtype=float)
        bc.R_ext[:] = np.asarray(R_ext, dtype=float)
        return bc

    def add_flux_condition(self, q_flux):
        bc = _DummyFluxBC(self.shape[0])
        bc.update_params(q_flux)
        return bc

    def compute_net_flux_for_solver(self):
        return self.current_flux


class _DummyFluid:
    def __init__(self, n_nodes: int):
        self.n_nodes = n_nodes
        self.temperature_vector = np.array([730.0, 740.0, 760.0], dtype=float)
        self.pressure_vector = np.full(n_nodes, 101325.0, dtype=float)
        self.density_vector = np.full(n_nodes, 900.0, dtype=float)
        self.velocity_vector = np.full(n_nodes, 1.0, dtype=float)
        self.d_h = np.full(n_nodes, 0.1, dtype=float)
        self.node_length = np.array([1.0, 1.0, 1.0], dtype=float)
        self.area = np.full(n_nodes, 0.05, dtype=float)
        self.material = _DummyMaterial()
        self.last_explicit = None
        self.last_implicit = None

    def add_coupling_source_distribution(self, explicit_arr, implicit_arr):
        self.last_explicit = np.asarray(explicit_arr, dtype=float)
        self.last_implicit = np.asarray(implicit_arr, dtype=float)


class LocalImplicitHeatExchangeTests(unittest.TestCase):
    def test_coupler_defaults_to_local_implicit(self):
        coupler = FluidSolidCouple(
            name="default-coupler",
            fluid=_DummyFluid(n_nodes=3),
            solid_boundary_region=_DummyBoundaryRegion(n_nodes=3),
            heated_perimeter=2.0,
            correlation_func=lambda Re, Pr, ratio: np.full_like(np.asarray(Re), 25.0),
            solid_node_capacitance=np.full(3, 5.0, dtype=float),
        )
        self.assertEqual(coupler.coupling_time_scheme, "local_implicit")

    def test_wall_hotter_reduces_delta_without_sign_reversal(self):
        q_to_fluid, delta_new, c_eff = FluidSolidCouple.compute_local_implicit_exchange(
            delta_old=np.array([100.0]),
            lambda_vals=np.array([10.0]),
            C_solid=np.array([20.0]),
            C_fluid=np.array([5.0]),
            dt=100.0,
        )

        self.assertGreater(q_to_fluid[0], 0.0)
        self.assertGreater(delta_new[0], 0.0)
        self.assertLess(delta_new[0], 100.0)
        self.assertAlmostEqual(c_eff[0], 4.0)

    def test_fluid_hotter_is_symmetric(self):
        q_to_fluid, delta_new, _ = FluidSolidCouple.compute_local_implicit_exchange(
            delta_old=np.array([-80.0]),
            lambda_vals=np.array([7.0]),
            C_solid=np.array([30.0]),
            C_fluid=np.array([10.0]),
            dt=50.0,
        )

        self.assertLess(q_to_fluid[0], 0.0)
        self.assertLess(delta_new[0], 0.0)
        self.assertLess(abs(delta_new[0]), 80.0)

    def test_small_dt_approaches_explicit_heat_rate(self):
        delta_old = np.array([12.0])
        lambda_vals = np.array([3.5])
        q_to_fluid, _, _ = FluidSolidCouple.compute_local_implicit_exchange(
            delta_old=delta_old,
            lambda_vals=lambda_vals,
            C_solid=np.array([1000.0]),
            C_fluid=np.array([900.0]),
            dt=1.0e-6,
        )

        self.assertAlmostEqual(q_to_fluid[0], (lambda_vals * delta_old)[0], places=4)

    def test_rejects_non_positive_dt(self):
        with self.assertRaises(ValueError):
            FluidSolidCouple.compute_local_implicit_exchange(
                delta_old=np.array([1.0]),
                lambda_vals=np.array([1.0]),
                C_solid=np.array([1.0]),
                C_fluid=np.array([1.0]),
                dt=0.0,
            )

    def test_local_implicit_diagnostics_reflect_energy_closure(self):
        fluid = _DummyFluid(n_nodes=3)
        boundary = _DummyBoundaryRegion(n_nodes=3)

        def corr(Re, Pr, ratio):
            return np.full_like(np.asarray(Re, dtype=float), 25.0)

        coupler = FluidSolidCouple(
            name="fake-coupler",
            fluid=fluid,
            solid_boundary_region=boundary,
            heated_perimeter=2.0,
            correlation_func=corr,
            solid_node_capacitance=np.full(3, 5.0, dtype=float),
            coupling_time_scheme="local_implicit",
        )

        boundary.T_surface[:] = np.array([950.0, 970.0, 990.0], dtype=float)
        fluid.temperature_vector[:] = np.array([900.0, 910.0, 920.0], dtype=float)
        coupler.execute(dt=0.2)

        diagnostics = coupler.get_coupling_diagnostics()
        self.assertIsNotNone(diagnostics)
        self.assertIn("local_implicit_q_to_fluid_sum_w", diagnostics)
        self.assertIn("local_implicit_fluid_source_sum_w", diagnostics)
        self.assertIn("local_implicit_solid_boundary_sum_w", diagnostics)
        self.assertIn("local_implicit_energy_mismatch_w", diagnostics)
        self.assertIn("coupling_tau_min_s", diagnostics)
        self.assertIn("coupling_dt_limit_s", diagnostics)
        self.assertIn("dt_over_coupling_tau_max", diagnostics)
        self.assertNotEqual(diagnostics["local_implicit_q_to_fluid_sum_w"], 0.0)
        self.assertAlmostEqual(
            diagnostics["local_implicit_fluid_source_sum_w"],
            float(np.sum(fluid.last_explicit)),
        )
        self.assertAlmostEqual(
            diagnostics["local_implicit_energy_mismatch_w"],
            diagnostics["local_implicit_fluid_source_sum_w"]
            + diagnostics["local_implicit_solid_boundary_sum_w"],
            places=12,
        )
        self.assertAlmostEqual(diagnostics["local_implicit_q_to_fluid_sum_w"], diagnostics["local_implicit_fluid_source_sum_w"], places=12)
        self.assertAlmostEqual(diagnostics["local_implicit_energy_mismatch_w"], 0.0, delta=1.0e-9)
        self.assertGreater(diagnostics["coupling_tau_min_s"], 0.0)
        self.assertAlmostEqual(diagnostics["coupling_dt_limit_s"], diagnostics["coupling_tau_min_s"])
        self.assertAlmostEqual(diagnostics["dt_over_coupling_tau_max"], diagnostics["local_implicit_dt_over_tau_max"])

    def test_local_implicit_stable_dt_does_not_limit_backward_euler_exchange(self):
        fluid = _DummyFluid(n_nodes=3)
        boundary = _DummyBoundaryRegion(n_nodes=3)

        def corr(Re, Pr, ratio):
            return np.full_like(np.asarray(Re, dtype=float), 25.0)

        solid_capacitance = np.array([5.0, 10.0, 20.0], dtype=float)
        coupler = FluidSolidCouple(
            name="fake-coupler",
            fluid=fluid,
            solid_boundary_region=boundary,
            heated_perimeter=2.0,
            correlation_func=corr,
            solid_node_capacitance=solid_capacitance,
            coupling_time_scheme="local_implicit",
        )

        coupler.execute(dt=0.2)

        self.assertAlmostEqual(
            coupler.get_max_stable_dt(safety_factor=0.8, max_limit=10.0),
            10.0,
        )

    def test_local_implicit_stable_dt_keeps_max_limit_before_first_execute(self):
        fluid = _DummyFluid(n_nodes=3)
        boundary = _DummyBoundaryRegion(n_nodes=3)

        def corr(Re, Pr, ratio):
            return np.full_like(np.asarray(Re, dtype=float), 25.0)

        coupler = FluidSolidCouple(
            name="fake-coupler",
            fluid=fluid,
            solid_boundary_region=boundary,
            heated_perimeter=2.0,
            correlation_func=corr,
            solid_node_capacitance=np.full(3, 5.0, dtype=float),
            coupling_time_scheme="local_implicit",
        )

        self.assertAlmostEqual(
            coupler.get_max_stable_dt(safety_factor=0.8, max_limit=10.0),
            10.0,
        )

if __name__ == "__main__":
    unittest.main()
