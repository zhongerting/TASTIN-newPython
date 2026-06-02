import numpy as np
import os
import sys

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Components.tec_electric import (
    distribute_axial_power_by_volume,
    electric_field_from_node_potential,
    joule_power_from_electric_field,
)


def _centers(faces):
    faces = np.asarray(faces, dtype=float)
    return 0.5 * (faces[:-1] + faces[1:])


def test_uniform_linear_matches_legacy_total_joule():
    n = 5
    faces = np.linspace(0.0, 1.0, n + 1)
    centers = _centers(faces)
    slope = 3.0
    potential = 1.2 + slope * centers
    rho = np.full(n, 2.0)
    volumes = np.ones((2, n)).reshape(-1) * 0.1

    field = electric_field_from_node_potential(potential, y_faces=faces)
    watts, _ = joule_power_from_electric_field(field, rho, volumes, (2, n))

    legacy_du = np.gradient(potential)
    legacy_l_node = (faces[-1] - faces[0]) / n
    legacy_q_vol = legacy_du**2 / (rho * legacy_l_node**2)
    legacy_watts = np.broadcast_to(legacy_q_vol[np.newaxis, :], (2, n)).reshape(-1) * volumes

    assert np.allclose(field, slope)
    assert np.allclose(np.sum(watts), np.sum(legacy_watts))


def test_nonuniform_linear_gives_constant_field():
    faces = np.array([0.0, 0.05, 0.20, 0.55, 1.0])
    centers = _centers(faces)
    slope = -4.0
    potential = 0.7 + slope * centers
    rho = np.full(centers.size, 5.0)
    volumes = np.array([[0.01, 0.03, 0.07, 0.09], [0.02, 0.04, 0.08, 0.10]]).reshape(-1)

    field = electric_field_from_node_potential(potential, y_faces=faces)
    watts, q_vol = joule_power_from_electric_field(field, rho, volumes, (2, centers.size))

    expected_q_vol = np.full(centers.size, slope**2 / rho[0])
    assert np.allclose(field, slope)
    assert np.allclose(q_vol, expected_q_vol)
    assert np.allclose(np.sum(watts), np.sum(volumes * np.tile(expected_q_vol, 2)))


def test_nonuniform_quadratic_differs_from_uniform_index_gradient():
    faces = np.array([0.0, 0.05, 0.20, 0.55, 1.0])
    centers = _centers(faces)
    potential = centers**2

    field = electric_field_from_node_potential(potential, y_faces=faces)
    legacy_like = np.gradient(potential) / ((faces[-1] - faces[0]) / centers.size)

    assert np.allclose(field, 2.0 * centers)
    assert not np.allclose(field, legacy_like)


def test_plasma_flux_power_uses_emitter_area_basis():
    q_e = np.array([-10.0, -20.0])
    q_c = np.array([8.0, 12.0])
    side_e = np.array([0.4, 0.6])
    side_c = np.array([0.5, 0.7])

    q_emitter = q_e * side_e
    q_collector = q_c * side_e
    q_collector_if_collector_area = q_c * side_c

    assert np.allclose(q_emitter, np.array([-4.0, -12.0]))
    assert np.allclose(q_collector, np.array([3.2, 7.2]))
    assert np.allclose(q_collector_if_collector_area - q_collector, q_c * (side_c - side_e))


def test_distribute_axial_power_by_volume_conserves_each_column():
    axial_power = np.array([3.0, 7.0, 11.0])
    volumes = np.array([[1.0, 1.0, 3.0], [2.0, 3.0, 1.0]]).reshape(-1)

    watts = distribute_axial_power_by_volume(axial_power, volumes, (2, 3)).reshape(2, 3)

    assert np.allclose(np.sum(watts, axis=0), axial_power)
    assert np.isclose(np.sum(watts), np.sum(axial_power))
    assert np.allclose(watts[:, 0], np.array([1.0, 2.0]))
    assert np.allclose(watts[:, 1], np.array([1.75, 5.25]))
    assert np.allclose(watts[:, 2], np.array([8.25, 2.75]))


def _assert_value_error(callback):
    try:
        callback()
    except ValueError:
        return
    raise AssertionError("Expected ValueError")


def test_distribute_axial_power_by_volume_rejects_invalid_inputs():
    _assert_value_error(
        lambda: distribute_axial_power_by_volume(np.array([1.0]), np.ones(4), (2, 2))
    )
    _assert_value_error(
        lambda: distribute_axial_power_by_volume(np.ones(2), np.ones(3), (2, 2))
    )
    _assert_value_error(
        lambda: distribute_axial_power_by_volume(np.array([1.0, np.nan]), np.ones(4), (2, 2))
    )
    _assert_value_error(
        lambda: distribute_axial_power_by_volume(np.ones(2), np.array([1.0, 0.0, 1.0, 0.0]), (2, 2))
    )
    _assert_value_error(
        lambda: distribute_axial_power_by_volume(np.ones(2), np.array([1.0, 1.0, -1.0, 1.0]), (2, 2))
    )


if __name__ == "__main__":
    test_uniform_linear_matches_legacy_total_joule()
    test_nonuniform_linear_gives_constant_field()
    test_nonuniform_quadratic_differs_from_uniform_index_gradient()
    test_plasma_flux_power_uses_emitter_area_basis()
    test_distribute_axial_power_by_volume_conserves_each_column()
    test_distribute_axial_power_by_volume_rejects_invalid_inputs()
