"""Run the V14 210 kW all-TFE helium depressurization accident."""

from __future__ import annotations

import math
from typing import Any, Dict

import numpy as np


REPRESENTATIVE_NAMES = ('Center', 'Ring1', 'Ring2', 'Ring3', 'Ring4')
EXPECTED_MULTIPLIERS = (1, 6, 9, 18, 24)
HELIUM_GAP_KEY = 'collector_iclad_gap'
HELIUM_H_INITIAL_W_M2K = 5678.0
HELIUM_H_FINAL_W_M2K = 0.0
HELIUM_GAP_WIDTH_M = 5.0e-5
TEMPERATURE_LIMITS_K = {
    'channel_wall': 1058.0,
    'pellet': 2700.0,
    'collector': 1023.0,
    'moderator': 930.0,
    'reflector': 1000.0,
}


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


def _temperature_peak(
        solid: Any,
        *,
        component: str,
        representative: str,
        limits_k: Dict[str, float]) -> Dict[str, Any]:
    values = np.asarray(solid.T, dtype=float).ravel()
    if values.size == 0:
        raise ValueError(f'{representative} {component} has no temperature nodes')
    axial = np.asarray(
        solid.mesh.geom_data.node_centers_y, dtype=float
    ).ravel()
    if axial.size != values.size:
        raise ValueError(
            f'{representative} {component} temperature/mesh size mismatch'
        )
    nonfinite = np.flatnonzero(~np.isfinite(values))
    if nonfinite.size:
        index = int(nonfinite[0])
        return {
            'component': 'nonfinite_temperature',
            'source_component': component,
            'representative': representative,
            'actual_k': float(values[index]),
            'limit_k': float(limits_k[component]),
            'axial_position_m': float(axial[index]),
        }
    index = int(np.argmax(values))
    return {
        'component': component,
        'representative': representative,
        'actual_k': float(values[index]),
        'limit_k': float(limits_k[component]),
        'axial_position_m': float(axial[index]),
    }


def collect_temperature_peaks(
        core: Any,
        limits_k: Dict[str, float] = TEMPERATURE_LIMITS_K,
) -> list[Dict[str, Any]]:
    peaks = []
    for name in REPRESENTATIVE_NAMES:
        solids = core.tfes[name].solids
        peaks.append(_temperature_peak(
            solids['inner_clad'],
            component='channel_wall',
            representative=f'{name}:inner_clad',
            limits_k=limits_k,
        ))
        peaks.append(_temperature_peak(
            solids['outer_clad'],
            component='channel_wall',
            representative=f'{name}:outer_clad',
            limits_k=limits_k,
        ))
        for component in ('pellet', 'collector'):
            peaks.append(_temperature_peak(
                solids[component],
                component=component,
                representative=name,
                limits_k=limits_k,
            ))
        if 'moderator' in solids:
            peaks.append(_temperature_peak(
                solids['moderator'],
                component='moderator',
                representative=name,
                limits_k=limits_k,
            ))
    for index, solid in enumerate(core.mod_rings):
        peaks.append(_temperature_peak(
            solid,
            component='moderator',
            representative=f'global_mod_ring_{index}',
            limits_k=limits_k,
        ))
    if core.reflector is None:
        raise ValueError('V14 core missing global reflector')
    peaks.append(_temperature_peak(
        core.reflector,
        component='reflector',
        representative='global_reflector',
        limits_k=limits_k,
    ))
    return peaks


def find_limit_trip(peaks: list[Dict[str, Any]]) -> Dict[str, Any] | None:
    for peak in peaks:
        if peak['component'] == 'nonfinite_temperature':
            return peak
    violations = [
        peak for peak in peaks
        if float(peak['actual_k']) > float(peak['limit_k'])
    ]
    if not violations:
        return None
    return max(
        violations,
        key=lambda item: float(item['actual_k']) / float(item['limit_k']),
    )


def collect_helium_metrics(
        build: Dict[str, Any],
        gaps: Dict[str, tuple[Any, int]],
        *,
        accident_time_s: float,
        active: bool,
) -> Dict[str, Any]:
    row: Dict[str, Any] = {
        'accident_elapsed_s': (
            float(build['system'].global_time) - float(accident_time_s)
        ),
        'helium_accident_active': bool(active),
        'helium_h_eq_W_m2K': (
            HELIUM_H_FINAL_W_M2K if active else HELIUM_H_INITIAL_W_M2K
        ),
        'helium_conduction_fraction': 0.0 if active else 1.0,
    }
    total_scaled = 0.0
    resistance_values = []
    for name in REPRESENTATIVE_NAMES:
        gap, multiplier = gaps[name]
        tfe = build['tfes'][name]
        collector = np.asarray(tfe.solids['collector'].T, dtype=float)
        inner_clad = np.asarray(tfe.solids['inner_clad'].T, dtype=float)
        q_out = -float(np.sum(np.asarray(gap.bc1.current_flux, dtype=float)))
        row[f'{name}_collector_mean_T_K'] = float(np.mean(collector))
        row[f'{name}_collector_max_T_K'] = float(np.max(collector))
        row[f'{name}_inner_clad_mean_T_K'] = float(np.mean(inner_clad))
        row[f'{name}_inner_clad_max_T_K'] = float(np.max(inner_clad))
        row[f'{name}_helium_gap_heat_out_W'] = q_out
        total_scaled += int(multiplier) * q_out
        resistance_values.extend(
            np.asarray(gap.R_gap_total, dtype=float).ravel()
        )
    row['helium_gap_heat_out_scaled_W'] = total_scaled
    row['helium_gap_R_total_min_K_W'] = float(np.min(resistance_values))
    row['helium_gap_R_total_max_K_W'] = float(np.max(resistance_values))
    return row
