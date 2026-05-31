import unittest
import os
import sys

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

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


class PressurizerPumpedClosedLoopTests(unittest.TestCase):
    def setUp(self):
        self.material = ConstantLiquid()
        self.area = 0.01
        self.hydraulic_diam = 0.1
        self.dt = 0.05
        self.t_end = 100.0

    def _volume(self, name, pressure):
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

    def _build_loop(self, reference_pressure=160000.0):
        ref = PressurizerVolume(
            name="pressurizer",
            volume=1.0,
            length=1.0,
            flow_area=self.area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=reference_pressure,
            initial_T=300.0,
            material=self.material,
        )
        nodes = [self._volume(f"node_{i}", reference_pressure) for i in range(3)]
        volumes = [ref] + nodes

        junctions = [
            PumpJunction(
                "pump",
                ref,
                nodes[0],
                flow_area=self.area,
                k_loss=40.0,
                custom_length=1.0,
                delta_p=3000.0,
            ),
            FlowJunction(
                "pipe_01",
                nodes[0],
                nodes[1],
                flow_area=self.area,
                k_loss=40.0,
                custom_length=1.0,
            ),
            FlowJunction(
                "pipe_12",
                nodes[1],
                nodes[2],
                flow_area=self.area,
                k_loss=40.0,
                custom_length=1.0,
            ),
            FlowJunction(
                "pipe_2r",
                nodes[2],
                ref,
                flow_area=self.area,
                k_loss=40.0,
                custom_length=1.0,
            ),
        ]

        network = HydraulicNetwork(volumes, junctions, gravity_vector=0.0)
        return network, ref, junctions[0]

    def _run_loop(self, network):
        n_steps = int(round(self.t_end / self.dt))
        flow_history = []
        ref_pressure_history = []
        mass_history = []

        for step in range(n_steps + 1):
            t = step * self.dt
            network.set_time(t)
            network.step_hydraulic(self.dt)
            flow_history.append(float(network.W_vec[0]))
            ref_pressure_history.append(float(network.P_vec[network.pressure_reference_idx]))
            mass_history.append(float(np.sum(network.rho_vec * network.V_vec)))

        return (
            np.asarray(flow_history),
            np.asarray(ref_pressure_history),
            np.asarray(mass_history),
        )

    def test_fixed_pressurizer_holds_absolute_pressure_and_flow_converges(self):
        network, _, _ = self._build_loop(reference_pressure=160000.0)
        initial_mass = float(np.sum(network.rho_vec * network.V_vec))

        flow, ref_pressure, mass = self._run_loop(network)
        tail = flow[-200:]

        self.assertTrue(np.allclose(ref_pressure, 160000.0, atol=1e-6))
        self.assertGreater(float(np.mean(tail)), 0.0)
        self.assertLess(float(np.ptp(tail)), max(1.0e-5, 1.0e-3 * abs(float(np.mean(tail)))))
        self.assertAlmostEqual(float(mass[-1]), initial_mass)
        self.assertTrue(np.allclose(mass, initial_mass))

    def test_reference_pressure_shift_does_not_change_stable_flow(self):
        low_network, _, _ = self._build_loop(reference_pressure=160000.0)
        high_network, _, _ = self._build_loop(reference_pressure=170000.0)

        low_flow, _, _ = self._run_loop(low_network)
        high_flow, _, _ = self._run_loop(high_network)

        self.assertAlmostEqual(float(low_flow[-1]), float(high_flow[-1]), places=9)

    def test_tabulated_pressurizer_pressure_tracks_table_without_changing_pump_drive(self):
        table_network, ref, pump = self._build_loop(reference_pressure=160000.0)
        ref.set_pressure_table(
            times=[0.0, 50.0, 100.0],
            pressures=[160000.0, 170000.0, 155000.0],
        )
        pump.set_delta_p(3000.0)

        fixed_network, _, _ = self._build_loop(reference_pressure=160000.0)

        table_flow, table_ref_pressure, _ = self._run_loop(table_network)
        fixed_flow, _, _ = self._run_loop(fixed_network)

        self.assertAlmostEqual(table_ref_pressure[0], 160000.0, places=6)
        self.assertAlmostEqual(table_ref_pressure[len(table_ref_pressure) // 2], 170000.0, places=6)
        self.assertAlmostEqual(table_ref_pressure[-1], 155000.0, places=6)
        self.assertGreater(float(table_flow[-1]), 0.0)
        self.assertAlmostEqual(float(table_flow[-1]), float(fixed_flow[-1]), places=8)


if __name__ == "__main__":
    unittest.main()
