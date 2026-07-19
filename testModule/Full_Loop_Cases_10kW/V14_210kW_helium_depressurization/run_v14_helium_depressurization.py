"""Run the V14 210 kW all-TFE helium depressurization accident."""

from __future__ import annotations

import math
from typing import Any, Dict


REPRESENTATIVE_NAMES = ('Center', 'Ring1', 'Ring2', 'Ring3', 'Ring4')
EXPECTED_MULTIPLIERS = (1, 6, 9, 18, 24)
HELIUM_GAP_KEY = 'collector_iclad_gap'
HELIUM_H_INITIAL_W_M2K = 5678.0
HELIUM_H_FINAL_W_M2K = 0.0
HELIUM_GAP_WIDTH_M = 5.0e-5


def collect_helium_gaps(build: Dict[str, Any]) -> Dict[str, tuple[Any, int]]:
    tfes = build['tfes']
    multipliers = build['ring_multipliers']
    if tuple(tfes) != REPRESENTATIVE_NAMES:
        raise ValueError(f'unexpected TFE names/order: {tuple(tfes)}')
    actual_multipliers = tuple(
        int(multipliers[name]) for name in REPRESENTATIVE_NAMES
    )
    if actual_multipliers != EXPECTED_MULTIPLIERS:
        raise ValueError(f'unexpected TFE multipliers: {actual_multipliers}')

    result: Dict[str, tuple[Any, int]] = {}
    for name, multiplier in zip(REPRESENTATIVE_NAMES, EXPECTED_MULTIPLIERS):
        gap = tfes[name].couplers.get(HELIUM_GAP_KEY)
        if gap is None:
            raise ValueError(f'{name} missing {HELIUM_GAP_KEY}')
        if not math.isclose(
                float(gap.gap), HELIUM_GAP_WIDTH_M,
                rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f'{name} helium gap width is {float(gap.gap)} m')
        h_eq = float(gap.k_gas) / float(gap.gap)
        if not math.isclose(
                h_eq, HELIUM_H_INITIAL_W_M2K,
                rel_tol=0.0, abs_tol=1.0e-9):
            raise ValueError(
                f'{name} initial helium h_eq is {h_eq} W/(m2*K)'
            )
        result[name] = (gap, multiplier)
    return result


def set_helium_h_eq(
        gaps: Dict[str, tuple[Any, int]], h_eq_w_m2k: float) -> None:
    h_eq = float(h_eq_w_m2k)
    if not math.isfinite(h_eq) or h_eq < 0.0:
        raise ValueError('helium h_eq must be finite and non-negative')
    for gap, _ in gaps.values():
        gap.k_gas = h_eq * float(gap.gap)


def read_source_accident_state(source_config: Dict[str, Any]) -> bool:
    value = source_config.get('helium_accident_active')
    if not isinstance(value, bool):
        raise ValueError(
            'run_config.json must contain boolean helium_accident_active'
        )
    return value
