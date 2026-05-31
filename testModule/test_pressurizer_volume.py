import unittest

import numpy as np

from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidVolume,
    PressurizerVolume,
    PumpJunction,
)
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork


class ConstantLiquid:
    def __init__(self, rho=1000.0, mu=1.0e-3, cp=1000.0):
        self.rho = rho
        self.mu = mu
        self.cp = cp

    def enthalpy(self, T, P):
        return self.cp * T

    def temperature_from_enthalpy(self, h, P):
        return np.asarray(h, dtype=float) / self.cp

    def heat_capacity(self, T, P):
        return np.asarray(T, dtype=float) * 0.0 + self.cp

    def density(self, T, P):
        return np.asarray(T, dtype=float) * 0.0 + self.rho

    def viscosity(self, T, P):
        return np.asarray(T, dtype=float) * 0.0 + self.mu


class PressurizerVolumeTests(unittest.TestCase):
    def setUp(self):
        self.material = ConstantLiquid()

    def _pressurizer(self, pressure=160000.0):
        return PressurizerVolume(
            name="pressurizer",
            volume=1.0,
            length=1.0,
            flow_area=0.01,
            hydraulic_diam=0.1,
            initial_P=pressure,
            initial_T=300.0,
            material=self.material,
        )

    def test_pressurizer_defaults_to_passive_pressure_reference(self):
        vol = self._pressurizer(160000.0)

        self.assertTrue(vol.is_pressure_reference)
        self.assertFalse(vol.is_pressure_boundary)
        self.assertAlmostEqual(vol.P, 160000.0)
        self.assertAlmostEqual(vol.target_P, 160000.0)

        vol.set_pressure(165000.0)
        self.assertAlmostEqual(vol.target_P, 165000.0)
        self.assertAlmostEqual(vol.compute_target_pressure(), 165000.0)

    def test_pressure_table_interpolates_clamps_and_clears_to_fixed_target(self):
        vol = self._pressurizer(160000.0)
        vol.set_pressure(165000.0)
        vol.set_pressure_table(
            times=[0.0, 5.0, 10.0],
            pressures=[160000.0, 170000.0, 165000.0],
        )

        self.assertAlmostEqual(vol.compute_target_pressure(0.0), 160000.0)
        self.assertAlmostEqual(vol.compute_target_pressure(2.5), 165000.0)
        self.assertAlmostEqual(vol.compute_target_pressure(5.0), 170000.0)
        self.assertAlmostEqual(vol.compute_target_pressure(20.0), 165000.0)

        vol.clear_pressure_table()
        self.assertAlmostEqual(vol.compute_target_pressure(), 165000.0)
        self.assertAlmostEqual(vol.target_P, 165000.0)


class PressurizerHydraulicNetworkTests(unittest.TestCase):
    def setUp(self):
        self.material = ConstantLiquid()
        self.area = 0.01
        self.hydraulic_diam = 0.1
        self.dt = 0.1

    def _volume(self, name, pressure=160000.0):
        return IncompressibleFluidVolume(
            name=name,
            volume=1.0,
            length=1.0,
            flow_area=self.area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=pressure,
            initial_T=300.0,
            material=self.material,
        )

    def _pressurizer(self, pressure=160000.0):
        return PressurizerVolume(
            name="pressurizer",
            volume=1.0,
            length=1.0,
            flow_area=self.area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=pressure,
            initial_T=300.0,
            material=self.material,
        )

    def _build_closed_loop(self, reference_pressure=160000.0, pump_delta=3000.0):
        ref = self._pressurizer(reference_pressure)
        node_a = self._volume("node_a", reference_pressure)
        node_b = self._volume("node_b", reference_pressure)

        pump = PumpJunction(
            "pump",
            ref,
            node_a,
            flow_area=self.area,
            k_loss=10.0,
            custom_length=1.0,
            delta_p=pump_delta,
        )
        pipe_ab = FlowJunction(
            "pipe_ab",
            node_a,
            node_b,
            flow_area=self.area,
            k_loss=10.0,
            custom_length=1.0,
        )
        pipe_br = FlowJunction(
            "pipe_br",
            node_b,
            ref,
            flow_area=self.area,
            k_loss=10.0,
            custom_length=1.0,
        )

        network = HydraulicNetwork([ref, node_a, node_b], [pump, pipe_ab, pipe_br], gravity_vector=0.0)
        return network, ref, pump

    def _junction_pressure_drops(self, network):
        return network.P_vec[network.idx_from_vec] - network.P_vec[network.idx_to_vec]

    def _equivalent_source(self, network, junction_index=0):
        network._update_fluid_properties()
        network._calc_momentum_coeffs(self.dt)
        self.assertNotEqual(network.A_coeffs[junction_index], 0.0)
        return network.B_coeffs[junction_index] / network.A_coeffs[junction_index]

    def test_reference_is_not_fixed_boundary_and_pressure_shift_preserves_deltas(self):
        network, ref, _ = self._build_closed_loop()

        self.assertEqual(network.pressure_reference_idx, 0)
        self.assertIs(network.pressure_reference_volume, ref)
        self.assertNotIn(0, network.fixed_pressure_indices)
        self.assertFalse(network.fixed_pressure_mask[0])

        network.P_vec[:] = np.array([150000.0, 151250.0, 149400.0])
        drops_before = self._junction_pressure_drops(network).copy()

        offset = network._apply_pressure_reference_shift()
        drops_after = self._junction_pressure_drops(network)

        self.assertAlmostEqual(offset, 10000.0)
        self.assertAlmostEqual(network.P_vec[0], 160000.0)
        np.testing.assert_allclose(drops_after, drops_before, atol=1e-9)

    def test_closed_loop_solve_anchors_reference_and_keeps_pump_source(self):
        network, _, _ = self._build_closed_loop()
        initial_mass_proxy = float(np.sum(network.rho_vec * network.V_vec))

        self.assertAlmostEqual(self._equivalent_source(network), 3000.0)

        A, B = network._assemble_pressure_system(self.dt, include_thermal_expansion=False)
        network._solve_linear_system(A, B)
        network._update_flow_rates()

        self.assertAlmostEqual(network.P_vec[0], 160000.0, places=6)
        self.assertGreater(network.W_vec[0], 0.0)

        for _ in range(5):
            network.step_hydraulic(self.dt)

        final_mass_proxy = float(np.sum(network.rho_vec * network.V_vec))
        self.assertAlmostEqual(final_mass_proxy, initial_mass_proxy)

    def test_invalid_reference_configurations_raise(self):
        ref_a = self._pressurizer()
        ref_b = self._pressurizer()
        FlowJunction("link", ref_a, ref_b, flow_area=self.area, k_loss=0.0, custom_length=1.0)

        with self.assertRaises(ValueError):
            HydraulicNetwork([ref_a, ref_b], list(ref_a.outlet_junctions), gravity_vector=0.0)

        ref = self._pressurizer()
        ref.is_pressure_boundary = True
        plain = self._volume("plain")
        junction = FlowJunction("bad_link", ref, plain, flow_area=self.area, k_loss=0.0, custom_length=1.0)

        with self.assertRaises(ValueError):
            HydraulicNetwork([ref, plain], [junction], gravity_vector=0.0)

    def test_set_time_reaches_pressurizer_and_pump(self):
        network, ref, pump = self._build_closed_loop()
        ref.set_pressure_table([0.0, 5.0], [160000.0, 170000.0])
        pump.set_pressure_table([0.0, 5.0], [3000.0, 4000.0])

        network.set_time(2.5)

        self.assertAlmostEqual(network.current_time, 2.5)
        self.assertAlmostEqual(ref.current_time, 2.5)
        self.assertAlmostEqual(ref.target_P, 165000.0)
        self.assertAlmostEqual(pump.current_time, 2.5)
        self.assertAlmostEqual(pump.compute_pump_head(), 3500.0)


if __name__ == "__main__":
    unittest.main()
