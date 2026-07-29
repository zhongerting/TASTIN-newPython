import unittest

import numpy as np

from Solvers.Hydrodynamics.Components import (
    FlowJunction,
    IncompressibleFluidVolume,
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
        return h / self.cp

    def heat_capacity(self, T, P):
        return np.asarray(T, dtype=float) * 0.0 + self.cp

    def density(self, T, P):
        return np.asarray(T, dtype=float) * 0.0 + self.rho

    def viscosity(self, T, P):
        return np.asarray(T, dtype=float) * 0.0 + self.mu


class PumpJunctionHydraulicNetworkTests(unittest.TestCase):
    def setUp(self):
        self.material = ConstantLiquid()
        self.area = 1.0
        self.hydraulic_diam = 1.0
        self.initial_pressure = 100000.0
        self.dt = 0.1

    def _volume(self, name):
        vol = IncompressibleFluidVolume(
            name=name,
            volume=1.0,
            length=1.0,
            flow_area=self.area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=self.initial_pressure,
            initial_T=300.0,
            material=self.material,
        )
        vol.is_pressure_boundary = True
        vol.target_P = self.initial_pressure
        return vol

    def _build_network(self, junction_factory):
        upstream = self._volume("upstream")
        downstream = self._volume("downstream")
        junction = junction_factory(upstream, downstream)
        network = HydraulicNetwork([upstream, downstream], [junction], gravity_vector=0.0)
        return network, junction

    def _equivalent_source(self, network, junction_index=0):
        network._calc_momentum_coeffs(self.dt)
        self.assertNotEqual(network.A_coeffs[junction_index], 0.0)
        return network.B_coeffs[junction_index] / network.A_coeffs[junction_index]

    def test_fixed_pump_pressure_enters_momentum_source(self):
        pump_net, _ = self._build_network(
            lambda upstream, downstream: PumpJunction(
                "pump",
                upstream,
                downstream,
                flow_area=self.area,
                k_loss=0.0,
                custom_length=1.0,
                delta_p=5000.0,
            )
        )
        plain_net, _ = self._build_network(
            lambda upstream, downstream: FlowJunction(
                "plain",
                upstream,
                downstream,
                flow_area=self.area,
                k_loss=0.0,
                custom_length=1.0,
            )
        )

        self.assertAlmostEqual(self._equivalent_source(pump_net), 5000.0)
        self.assertAlmostEqual(self._equivalent_source(plain_net), 0.0)

        pump_net._update_flow_rates()
        plain_net._update_flow_rates()
        self.assertGreater(pump_net.W_vec[0], 0.0)
        self.assertAlmostEqual(plain_net.W_vec[0], 0.0)

    def test_pressure_table_interpolates_clamps_and_feeds_network_time(self):
        network, pump = self._build_network(
            lambda upstream, downstream: PumpJunction(
                "pump",
                upstream,
                downstream,
                flow_area=self.area,
                k_loss=0.0,
                custom_length=1.0,
                delta_p=750.0,
            )
        )

        pump.set_pressure_table(
            times=[0.0, 5.0, 10.0],
            pressures=[1000.0, 3000.0, 2000.0],
        )

        self.assertAlmostEqual(pump.compute_pump_head(0.0), 1000.0)
        self.assertAlmostEqual(pump.compute_pump_head(2.5), 2000.0)
        self.assertAlmostEqual(pump.compute_pump_head(5.0), 3000.0)
        self.assertAlmostEqual(pump.compute_pump_head(20.0), 2000.0)

        network.set_time(2.5)
        self.assertAlmostEqual(network.current_time, 2.5)
        self.assertAlmostEqual(pump.current_time, 2.5)
        self.assertAlmostEqual(self._equivalent_source(network), 2000.0)

        pump.clear_pressure_table()
        self.assertAlmostEqual(self._equivalent_source(network), 750.0)

    def test_prescribed_flow_relaxation_does_not_overshoot_target(self):
        network, junction = self._build_network(
            lambda upstream, downstream: FlowJunction(
                'controlled',
                upstream,
                downstream,
                flow_area=self.area,
                k_loss=0.0,
                custom_length=1.0,
            )
        )
        junction.target_W = 2.4
        junction.W = 2.46
        network = HydraulicNetwork(
            network.volumes_obj,
            [junction],
            gravity_vector=0.0,
        )
        network._calc_momentum_coeffs(0.05)

        self.assertAlmostEqual(network.B_coeffs[0], 2.4)


if __name__ == "__main__":
    unittest.main()
