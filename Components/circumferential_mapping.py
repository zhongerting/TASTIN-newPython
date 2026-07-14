import operator

import numpy as np


_PERIOD_DEG = 360.0


def _positive_segment_count(name: str, value) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer")
    try:
        count = operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if count <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return count


def _wrapped_interval_parts(start_deg: float, width_deg: float):
    start = float(start_deg) % _PERIOD_DEG
    end = start + float(width_deg)
    if end <= _PERIOD_DEG:
        return ((start, end),)
    return ((start, _PERIOD_DEG), (0.0, end - _PERIOD_DEG))


def build_uniform_circumferential_mapping(
    n_source: int,
    n_target: int,
    source_offset_deg: float = 0.0,
    target_offset_deg: float = 0.0,
) -> np.ndarray:
    """Return a column-conservative periodic source-to-target allocation matrix."""
    n_source = _positive_segment_count("n_source", n_source)
    n_target = _positive_segment_count("n_target", n_target)
    source_offset_deg = float(source_offset_deg)
    target_offset_deg = float(target_offset_deg)
    if not np.isfinite(source_offset_deg) or not np.isfinite(target_offset_deg):
        raise ValueError("circumferential offsets must be finite")

    source_width = _PERIOD_DEG / n_source
    target_width = _PERIOD_DEG / n_target
    source_parts = [
        _wrapped_interval_parts(source_offset_deg + i * source_width, source_width)
        for i in range(n_source)
    ]
    target_parts = [
        _wrapped_interval_parts(target_offset_deg + j * target_width, target_width)
        for j in range(n_target)
    ]

    mapping = np.zeros((n_target, n_source), dtype=float)
    for j, target in enumerate(target_parts):
        for i, source in enumerate(source_parts):
            overlap = sum(
                max(0.0, min(source_end, target_end) - max(source_start, target_start))
                for source_start, source_end in source
                for target_start, target_end in target
            )
            mapping[j, i] = overlap / source_width

    mapping[np.abs(mapping) < 1.0e-15] = 0.0
    column_sums = mapping.sum(axis=0)
    if np.any(column_sums <= 0.0):
        raise RuntimeError("circumferential mapping left a source segment unallocated")
    mapping /= column_sums
    return mapping


def map_circumferential_intensive(
    source_values,
    mapping,
    source_weights=None,
) -> np.ndarray:
    """Map intensive values using physical source weights and angular allocation."""
    values = np.asarray(source_values, dtype=float)
    allocation = np.asarray(mapping, dtype=float)
    if allocation.ndim != 2 or not allocation.shape[0] or not allocation.shape[1]:
        raise ValueError("mapping must be a non-empty 2D array")
    if not np.all(np.isfinite(allocation)) or np.any(allocation < 0.0):
        raise ValueError("mapping must contain finite non-negative values")
    if values.ndim == 0 or values.shape[-1] != allocation.shape[1]:
        raise ValueError("source_values final axis must match mapping source segments")
    if not np.all(np.isfinite(values)):
        raise ValueError("source_values must be finite")

    if source_weights is None:
        weights = np.ones(allocation.shape[1], dtype=float)
    else:
        weights = np.asarray(source_weights, dtype=float)
        if weights.shape != (allocation.shape[1],):
            raise ValueError("source_weights must match mapping source segments")
        if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
            raise ValueError("source_weights must be finite and non-negative")

    denominator = allocation @ weights
    if np.any(denominator <= 0.0):
        raise ValueError("every target segment must receive positive physical weight")
    numerator = np.matmul(values * weights, allocation.T)
    return numerator / denominator
