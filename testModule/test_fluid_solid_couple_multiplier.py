"""Focused regression tests for effective heat-transfer failure scaling."""

import unittest

import numpy as np

from Solvers.Couplers import FluidSolidCouple


class _FakeResistanceBC:
    def __init__(self, size):
        self.R_ext = np.full(size, 1.0e-3, dtype=float)
        self.T_ext = np.full(size, 300.0, dtype=float)
        self.q_flux = np.zeros(size, dtype=float)

    def update_params(self, T_ext=None, R_ext=None, q_flux=None):
        if T_ext is not None:
            self.T_ext[:] = np.asarray(T_ext, dtype=float)
        if R_ext is not None:
            self.R_ext[:] = np.asarray(R_ext, dtype=float)
        if q_flux is not None:
            self.q_flux[:] = np.asarray(q_flux, dtype=float)


class _FakeBoundary:
    shape = (1,)
    area = np.array([1.0], dtype=float)

    def __init__(self, wall_temperature=500.0):
        self.T_surface = np.array([wall_temperature], dtype=float)
        self.resistance_bc = None
        self.flux_bc = None

    def compute_net_flux_for_solver(self):
        return None

    def add_resistance_condition(self, T_ext, R_ext, R_add=0.0):
        del R_add
        self.resistance_bc = _FakeResistanceBC(1)
        self.resistance_bc.update_params(T_ext=T_ext, R_ext=R_ext)
        return self.resistance_bc

    def add_flux_condition(self, q_flux):
        self.flux_bc = _FakeResistanceBC(1)
        self.flux_bc.update_params(q_flux=q_flux)
        return self.flux_bc


class _FakeMaterial:
    @staticmethod
    def viscosity(T, P):
        return np.ones_like(np.asarray(T, dtype=float))

    @staticmethod
    def conductivity(T, P):
        return np.full_like(np.asarray(T, dtype=float), 2.0)

    @staticmethod
    def heat_capacity(T, P):
        return np.full_like(np.asarray(T, dtype=float), 5.0)

    @staticmethod
    def prandtl_number(T, P):
        return np.ones_like(np.asarray(T, dtype=float))


class _FakeFluid:
    n_nodes = 1
    node_length = 4.0
    d_h = 2.0
    material = _FakeMaterial()

    def __init__(self):
        self.temperature_vector = np.array([300.0], dtype=float)
        self.pressure_vector = np.array([1.0e5], dtype=float)
        self.density_vector = np.array([1.0], dtype=float)
        self.velocity_vector = np.array([2.0], dtype=float)
        self.sources = []

    def add_coupling_source_distribution(self, explicit, implicit):
        self.sources.append((np.asarray(explicit, dtype=float).copy(),
                             np.asarray(implicit, dtype=float).copy()))


def _make_coupler(multiplier):
    fluid = _FakeFluid()
    boundary = _FakeBoundary()
    coupler = FluidSolidCouple(
        name="test_coupler",
        fluid=fluid,
        solid_boundary_region=boundary,
        heated_perimeter=3.0,
        correlation_func=lambda Re, Pr, ratio: np.ones_like(np.asarray(Re, dtype=float)),
        coupling_time_scheme="current",
        coupling_multiplier=multiplier,
    )
    return coupler, fluid, boundary


class FluidSolidCoupleMultiplierTests(unittest.TestCase):
    def test_half_transfer_scales_lambda_and_source(self):
        coupler, fluid, boundary = _make_coupler(0.5)
        coupler.execute()

        # h=Nu*k/d_h=1, A=perimeter*node_length=12, lambda=6.
        self.assertAlmostEqual(float(boundary.resistance_bc.R_ext[0]), 1.0 / 6.0)
        explicit, implicit = fluid.sources[-1]
        self.assertAlmostEqual(float(explicit[0]), 6.0 * 500.0)
        self.assertAlmostEqual(float(implicit[0]), 6.0)
        self.assertEqual(coupler.get_coupling_diagnostics()["coupling_multiplier"], 0.5)

    def test_zero_transfer_clears_source_and_boundary_exchange(self):
        coupler, fluid, boundary = _make_coupler(0.0)
        coupler.execute()

        explicit, implicit = fluid.sources[-1]
        np.testing.assert_array_equal(explicit, np.zeros(1))
        np.testing.assert_array_equal(implicit, np.zeros(1))
        self.assertEqual(float(boundary.resistance_bc.R_ext[0]), 1.0e30)
        self.assertTrue(coupler.get_coupling_diagnostics()["coupling_disabled"])


if __name__ == "__main__":
    unittest.main()
