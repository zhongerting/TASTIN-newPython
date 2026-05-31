import unittest

from Solvers.Hydrodynamics.Components import (
    FixedPressurePumpVolume,
    FlowJunction,
    IncompressibleFluidChannel,
    IncompressibleFluidVolume,
)


class ConstantLiquid:
    def __init__(self, rho=1000.0, mu=1.0e-3, cp=1000.0):
        self.rho = rho
        self.mu = mu
        self.cp = cp

    def enthalpy(self, T, P):
        return self.cp * T

    def temperature_from_enthalpy(self, h, P):
        return h / self.cp

    def density(self, T, P):
        return self.rho

    def viscosity(self, T, P):
        return self.mu


class FixedPressurePumpVolumeTests(unittest.TestCase):
    def setUp(self):
        self.material = ConstantLiquid()
        self.area = 1.0
        self.hydraulic_diam = 1.0
        self.mass_flow = 0.2

    def _plain_volume(self, name):
        return IncompressibleFluidVolume(
            name=name,
            volume=1.0,
            length=1.0,
            flow_area=self.area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=100000.0,
            initial_T=300.0,
            material=self.material,
        )

    def _pump_volume(self, name, delta_p):
        return FixedPressurePumpVolume(
            name=name,
            volume=1.0,
            length=1.0,
            flow_area=self.area,
            hydraulic_diam=self.hydraulic_diam,
            initial_P=100000.0,
            initial_T=300.0,
            material=self.material,
            delta_p=delta_p,
        )

    def _channel_with_middle_volume(self, middle_volume):
        channel = IncompressibleFluidChannel.__new__(IncompressibleFluidChannel)
        channel.name = "pump_test_channel"
        channel.material = self.material
        channel.volumes = [
            self._plain_volume("upstream"),
            middle_volume,
            self._plain_volume("downstream"),
        ]
        channel.internal_junctions = [
            FlowJunction(
                "j0",
                channel.volumes[0],
                channel.volumes[1],
                flow_area=self.area,
                k_loss=0.0,
                custom_length=1.0,
            ),
            FlowJunction(
                "j1",
                channel.volumes[1],
                channel.volumes[2],
                flow_area=self.area,
                k_loss=0.0,
                custom_length=1.0,
            ),
        ]
        for junction in channel.internal_junctions:
            junction.W = self.mass_flow
        return channel

    @staticmethod
    def _pressures(channel):
        return [vol.P for vol in channel.volumes]

    def test_downstream_distribution_applies_fixed_pressure_rise(self):
        baseline = self._channel_with_middle_volume(self._plain_volume("middle"))
        pumped = self._channel_with_middle_volume(self._pump_volume("pump", 5000.0))

        baseline.update_pressure_distribution_downstream(100000.0)
        pumped.update_pressure_distribution_downstream(100000.0)

        self.assertAlmostEqual(pumped.volumes[1].P - baseline.volumes[1].P, 5000.0)
        self.assertAlmostEqual(pumped.volumes[2].P - baseline.volumes[2].P, 5000.0)

    def test_upstream_distribution_reverses_fixed_pressure_rise(self):
        channel = self._channel_with_middle_volume(self._pump_volume("pump", 5000.0))
        channel.update_pressure_distribution_downstream(100000.0)
        downstream_pressures = self._pressures(channel)

        channel.update_pressure_distribution_upstream(downstream_pressures[-1])

        for actual, expected in zip(self._pressures(channel), downstream_pressures):
            self.assertAlmostEqual(actual, expected)

    def test_zero_delta_p_matches_plain_volume_distribution(self):
        baseline = self._channel_with_middle_volume(self._plain_volume("middle"))
        zero_pump = self._channel_with_middle_volume(self._pump_volume("pump", 0.0))

        baseline.update_pressure_distribution_downstream(100000.0)
        zero_pump.update_pressure_distribution_downstream(100000.0)

        for actual, expected in zip(self._pressures(zero_pump), self._pressures(baseline)):
            self.assertAlmostEqual(actual, expected)

        baseline.update_pressure_distribution_upstream(95000.0)
        zero_pump.update_pressure_distribution_upstream(95000.0)

        for actual, expected in zip(self._pressures(zero_pump), self._pressures(baseline)):
            self.assertAlmostEqual(actual, expected)

    def test_set_delta_p_updates_next_pressure_distribution(self):
        channel = self._channel_with_middle_volume(self._pump_volume("pump", 1000.0))
        channel.update_pressure_distribution_downstream(100000.0)
        initial_outlet_pressure = channel.volumes[2].P

        channel.volumes[1].set_delta_p(3500.0)
        channel.update_pressure_distribution_downstream(100000.0)

        self.assertAlmostEqual(channel.volumes[2].P - initial_outlet_pressure, 2500.0)

    def test_pump_volume_derivative_remains_scalar_enthalpy_derivative(self):
        pump = self._pump_volume("pump", 5000.0)

        derivative = pump.get_volume_derivatives(self.material)

        self.assertIsInstance(derivative, float)
        self.assertNotIsInstance(derivative, tuple)


if __name__ == "__main__":
    unittest.main()
