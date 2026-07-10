import os
import sys

import numpy as np


ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)


from Components.RadiatorThermalShield import RadiatorThermalShield


class FakeRadiatorUnit:
    def __init__(self, temperatures, background=3.0):
        self._temperatures = np.asarray(temperatures, dtype=float)
        self.background_updates = []
        self.radiation_background_temperature = np.full_like(self._temperatures, float(background))

    def get_radiation_surface_temperature(self):
        return np.array(self._temperatures, copy=True)

    def set_radiation_background_temperature(self, value):
        arr = np.asarray(value, dtype=float)
        if arr.ndim == 0:
            arr = np.full_like(self._temperatures, float(arr))
        self.radiation_background_temperature = np.array(arr, copy=True)
        self.background_updates.append(np.array(arr, copy=True))


def test_active_shield_raises_effective_background_temperature():
    unit = FakeRadiatorUnit(np.linspace(760.0, 820.0, 8))
    shield = RadiatorThermalShield(
        name="shield",
        radiator_units=[unit],
        active_until_s=100.0,
        background_temperature_k=3.0,
        shield_view_factor=0.85,
        solar_heat_flux_w_m2=0.0,
    )

    shield.pre_step(dt=0.1, current_time=10.0)

    assert shield.last_active is True
    assert len(unit.background_updates) == 1
    updated = unit.background_updates[-1]
    assert updated.shape == (8,)
    assert np.all(updated > 3.0)
    assert np.all(updated < unit.get_radiation_surface_temperature())
    diagnostics = shield.get_diagnostics()
    assert diagnostics["radiation_shield_active"] is True
    assert diagnostics["radiation_shield_effective_background_mean_k"] > 3.0


def test_active_shield_closes_quasi_steady_energy_balance():
    unit = FakeRadiatorUnit(np.array([780.0, 800.0, 820.0]))
    unit.tube_bare_area = np.array([1.0, 1.0, 1.0])
    unit.fin_radiating_area = np.zeros(3)
    shield = RadiatorThermalShield(
        name="shield",
        radiator_units=[unit],
        active_until_s=100.0,
        background_temperature_k=3.0,
        shield_view_factor=0.85,
        inner_emissivity=0.8,
        outer_emissivity=0.75,
        conductivity_w_m_k=1.0,
        thickness_m=0.002,
        solar_heat_flux_w_m2=120.0,
    )

    shield.pre_step(dt=0.1, current_time=10.0)

    updated = unit.background_updates[-1]
    diagnostics = shield.get_diagnostics()
    assert np.all(updated > 3.0)
    assert np.all(updated < unit.get_radiation_surface_temperature())
    assert diagnostics["radiation_shield_outer_temperature_mean_k"] < diagnostics["radiation_shield_inner_temperature_mean_k"]
    assert diagnostics["radiation_shield_q_solar_w"] > 0.0
    assert diagnostics["radiation_shield_solver_failures"] == 0
    assert abs(diagnostics["radiation_shield_energy_residual_w"]) < 1.0e-6
    assert diagnostics["radiation_shield_energy_residual_rel"] < 1.0e-8


def test_inactive_shield_restores_deep_space_background():
    unit = FakeRadiatorUnit(np.linspace(760.0, 820.0, 8), background=450.0)
    shield = RadiatorThermalShield(
        name="shield",
        radiator_units=[unit],
        active_until_s=5.0,
        background_temperature_k=3.0,
        shield_view_factor=0.85,
    )

    shield.pre_step(dt=0.1, current_time=10.0)

    assert shield.last_active is False
    np.testing.assert_allclose(unit.background_updates[-1], np.full(8, 3.0))
    diagnostics = shield.get_diagnostics()
    assert diagnostics["radiation_shield_active"] is False
    assert diagnostics["radiation_shield_effective_background_mean_k"] == 3.0


def test_fortran_shield2_groups_78_tubes_into_12_circumferential_sectors():
    units = []
    fin_t4_by_tube = []
    for tube_index in range(78):
        surface_t = 720.0 + 0.25 * tube_index
        unit = FakeRadiatorUnit(np.full(8, surface_t))
        fin_t = 680.0 + float(tube_index)
        unit.last_fin_temperature = np.array(
            [
                [fin_t, fin_t + 1.0],
                [fin_t + 2.0, fin_t + 3.0],
            ],
            dtype=float,
        )
        units.append(unit)
        fin_t4_by_tube.append(float(np.mean(unit.last_fin_temperature ** 4)))

    shield = RadiatorThermalShield(
        name="shield",
        radiator_units=units,
        active_until_s=100.0,
        model="fortran_shield2",
        background_temperature_k=3.0,
        outer_emissivity=0.1,
        conductivity_w_m_k=0.0008,
        thickness_m=0.01,
    )

    shield.pre_step(dt=0.1, current_time=10.0)

    diagnostics = shield.get_diagnostics()
    assert diagnostics["radiation_shield_model"] == "fortran_shield2"
    assert diagnostics["radiation_shield_converged"] is True
    assert diagnostics["radiation_shield_iteration_count"] > 0
    assert len(diagnostics["radiation_shield_tc4_rad_12"]) == 12
    assert len(diagnostics["radiation_shield_qrr_12"]) == 12
    assert len(diagnostics["radiation_shield_qrr_weight_12"]) == 12
    assert len(diagnostics["radiation_shield_qrr_bg4_12"]) == 12
    assert len(diagnostics["radiation_shield_qrr_bgT_12"]) == 12
    assert len(diagnostics["radiation_shield_inner_t4_8"]) == 8
    assert len(diagnostics["radiation_shield_inner_temperature_8_k"]) == 8
    assert len(diagnostics["radiation_shield_outer_temperature_8_k"]) == 8
    np.testing.assert_allclose(
        diagnostics["radiation_shield_qrr_bg4_12"],
        np.asarray(diagnostics["radiation_shield_qrr_12"]) / np.asarray(diagnostics["radiation_shield_qrr_weight_12"]),
    )

    expected_sector_t4 = []
    for sector in range(12):
        values = [
            fin_t4_by_tube[i]
            for i in range(78)
            if int(np.floor(i * 12 / 78)) == sector
        ]
        expected_sector_t4.append(float(np.mean(values)))
    np.testing.assert_allclose(
        diagnostics["radiation_shield_tc4_rad_12"],
        expected_sector_t4,
        rtol=1.0e-12,
        atol=0.0,
    )

    for unit in units:
        assert len(unit.background_updates) == 1
        assert unit.background_updates[-1].shape == (8,)
        assert np.all(np.isfinite(unit.background_updates[-1]))
        assert np.all(unit.background_updates[-1] > 0.0)


def test_fortran_shield2_normalizes_raw_qrr_before_background_writeback():
    temperature = 760.0
    units = []
    for _ in range(78):
        unit = FakeRadiatorUnit(np.full(8, temperature))
        unit.last_fin_temperature = np.full((2, 2), temperature)
        units.append(unit)

    shield = RadiatorThermalShield(
        name="shield",
        radiator_units=units,
        active_until_s=100.0,
        model="fortran_shield2",
        background_temperature_k=3.0,
        outer_emissivity=0.1,
        conductivity_w_m_k=0.0008,
        thickness_m=0.01,
    )

    shield.pre_step(dt=0.1, current_time=10.0)

    diagnostics = shield.get_diagnostics()
    raw_qrr = np.asarray(diagnostics["radiation_shield_qrr_12"], dtype=float)
    weights = np.asarray(diagnostics["radiation_shield_qrr_weight_12"], dtype=float)
    bg4 = np.asarray(diagnostics["radiation_shield_qrr_bg4_12"], dtype=float)
    bgT = np.asarray(diagnostics["radiation_shield_qrr_bgT_12"], dtype=float)

    assert diagnostics["radiation_shield_converged"] is True
    assert np.all(weights > 1.0)
    np.testing.assert_allclose(bg4, raw_qrr / weights, rtol=1.0e-12)
    np.testing.assert_allclose(bgT, bg4 ** 0.25, rtol=1.0e-12)
    assert np.all(bgT <= temperature * (1.0 + 1.0e-12))

    sectors = np.floor(np.arange(78, dtype=float) * 12.0 / 78.0).astype(int)
    for index, unit in enumerate(units):
        expected = np.full(8, bgT[sectors[index]])
        np.testing.assert_allclose(unit.background_updates[-1], expected)


if __name__ == "__main__":
    test_active_shield_raises_effective_background_temperature()
    test_active_shield_closes_quasi_steady_energy_balance()
    test_inactive_shield_restores_deep_space_background()
    test_fortran_shield2_groups_78_tubes_into_12_circumferential_sectors()
    test_fortran_shield2_normalizes_raw_qrr_before_background_writeback()
    print("Radiator thermal shield checks passed.")
