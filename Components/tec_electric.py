"""Shared TEC electric-field and Joule-heating helpers."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def node_centers_from_lengths(node_lengths: np.ndarray) -> np.ndarray:
    """Return cell-center coordinates from axial cell lengths."""
    lengths = np.asarray(node_lengths, dtype=float)
    if lengths.ndim != 1:
        raise ValueError(f"node_lengths must be 1D, got shape {lengths.shape}")
    if np.any(lengths <= 0.0):
        raise ValueError("node_lengths must be positive")

    faces = np.concatenate(([0.0], np.cumsum(lengths)))
    return 0.5 * (faces[:-1] + faces[1:])


def node_centers_from_faces(y_faces: np.ndarray) -> np.ndarray:
    """Return cell-center coordinates from axial face coordinates."""
    faces = np.asarray(y_faces, dtype=float)
    if faces.ndim != 1:
        raise ValueError(f"y_faces must be 1D, got shape {faces.shape}")
    if faces.size < 2:
        raise ValueError("y_faces must contain at least two entries")
    if np.any(np.diff(faces) <= 0.0):
        raise ValueError("y_faces must be strictly increasing")
    return 0.5 * (faces[:-1] + faces[1:])


def electric_field_from_node_potential(
    potential: np.ndarray,
    *,
    y_faces: Optional[np.ndarray] = None,
    node_lengths: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Compute dU/dy at axial cell centers.

    ThermoCalc returns UE/UC at axial nodes. For non-uniform axial grids the
    derivative must be taken against the node-center coordinates, not against
    implicit unit-spaced indices.
    """
    values = np.asarray(potential, dtype=float)
    if values.ndim != 1:
        raise ValueError(f"potential must be 1D, got shape {values.shape}")

    if y_faces is not None:
        centers = node_centers_from_faces(y_faces)
    elif node_lengths is not None:
        centers = node_centers_from_lengths(node_lengths)
    else:
        raise ValueError("either y_faces or node_lengths must be provided")

    if centers.shape != values.shape:
        raise ValueError(
            f"potential shape {values.shape} does not match axial grid shape {centers.shape}"
        )
    if values.size == 1:
        return np.zeros_like(values)

    edge_order = 2 if values.size >= 3 else 1
    return np.gradient(values, centers, edge_order=edge_order)


def joule_power_from_electric_field(
    electric_field: np.ndarray,
    resistivity: np.ndarray,
    volumes_flat: np.ndarray,
    shape_nodes: Tuple[int, int],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Map 1D axial electric field to flat 2D per-cell Joule power.

    Returns (q_watts_flat, q_vol_1d). q_vol_1d is useful for diagnostics.
    """
    field = np.asarray(electric_field, dtype=float)
    rho = np.asarray(resistivity, dtype=float)
    if field.ndim != 1 or rho.ndim != 1:
        raise ValueError(f"field and resistivity must be 1D, got {field.shape}/{rho.shape}")
    if field.shape != rho.shape:
        raise ValueError(f"field shape {field.shape} does not match resistivity shape {rho.shape}")

    nx, ny = shape_nodes
    if field.shape != (ny,):
        raise ValueError(f"field shape {field.shape} does not match axial node count {(ny,)}")

    safe_rho = np.maximum(rho, 1.0e-12)
    q_vol_1d = field**2 / safe_rho
    q_vol_2d = np.broadcast_to(q_vol_1d[np.newaxis, :], (nx, ny))
    q_watts_flat = q_vol_2d.reshape(-1) * np.asarray(volumes_flat, dtype=float)
    return q_watts_flat, q_vol_1d


def distribute_axial_power_by_volume(
    axial_power_w: np.ndarray,
    volumes_flat: np.ndarray,
    shape_nodes: Tuple[int, int],
) -> np.ndarray:
    """Distribute per-axial-node power to a 2D thermal mesh without changing totals."""
    axial_power = np.asarray(axial_power_w, dtype=float)
    volumes = np.asarray(volumes_flat, dtype=float)
    nx, ny = shape_nodes

    if axial_power.shape != (ny,):
        raise ValueError(f"axial_power_w shape {axial_power.shape} does not match axial node count {(ny,)}")
    if volumes.shape != (nx * ny,):
        raise ValueError(f"volumes_flat shape {volumes.shape} does not match flattened mesh size {(nx * ny,)}")
    if not np.all(np.isfinite(axial_power)):
        raise ValueError("axial_power_w must contain only finite values")
    if not np.all(np.isfinite(volumes)):
        raise ValueError("volumes_flat must contain only finite values")
    if np.any(volumes < 0.0):
        raise ValueError("volumes_flat must not contain negative values")

    volumes_2d = volumes.reshape(nx, ny)
    column_volumes = np.sum(volumes_2d, axis=0)
    if np.any(column_volumes <= 0.0):
        raise ValueError("each axial column must have positive total volume")

    return (volumes_2d / column_volumes[np.newaxis, :] * axial_power[np.newaxis, :]).reshape(-1)


def joule_power_from_node_potential(
    potential: np.ndarray,
    resistivity: np.ndarray,
    volumes_flat: np.ndarray,
    shape_nodes: Tuple[int, int],
    *,
    y_faces: Optional[np.ndarray] = None,
    node_lengths: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute flat per-cell Joule power from axial node potentials."""
    field = electric_field_from_node_potential(
        potential,
        y_faces=y_faces,
        node_lengths=node_lengths,
    )
    watts, q_vol = joule_power_from_electric_field(
        field,
        resistivity,
        volumes_flat,
        shape_nodes,
    )
    return watts, q_vol, field
