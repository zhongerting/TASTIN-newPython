"""Run the V14 210 kW all-TFE helium depressurization accident."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    build_debug_case,
)
from testModule.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    REACTIVITY_HISTORY_FIELDS,
    ReactivityControlRunConfig,
    _next_step_dt,
    _validate_runtime,
    _write_json,
    collect_reactivity_metrics,
    load_baseline_debug_config,
    prepare_reactivity_control,
)


CASE_NAME = 'V14_10kW_210kW_all_tfe_helium_depressurization'
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
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / 'runs' / 'default'
ACCIDENT_FIELDS = [
    'accident_elapsed_s', 'helium_accident_active',
    'helium_h_eq_W_m2K', 'helium_conduction_fraction',
    'helium_gap_heat_out_scaled_W',
    'helium_gap_R_total_min_K_W', 'helium_gap_R_total_max_K_W',
    'channel_wall_max_T_K', 'pellet_max_T_K', 'collector_max_T_K',
    'moderator_max_T_K', 'reflector_max_T_K',
]
for _name in REPRESENTATIVE_NAMES:
    ACCIDENT_FIELDS.extend([
        f'{_name}_collector_mean_T_K', f'{_name}_collector_max_T_K',
        f'{_name}_inner_clad_mean_T_K', f'{_name}_inner_clad_max_T_K',
        f'{_name}_helium_gap_heat_out_W',
    ])
ACCIDENT_HISTORY_FIELDS = REACTIVITY_HISTORY_FIELDS + ACCIDENT_FIELDS


@dataclass(frozen=True)
class HeliumAccidentRunConfig:
    restart_in: Path
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = 100.0
    dt_s: float = 0.05
    tec_update_interval_s: float = 0.05
    record_interval_s: float = 0.1
    checkpoint_interval_s: float = 10.0
    min_fluid_temperature_stop_k: Optional[float] = 500.0
    max_power_factor: float = 2.0
    wall_limit_k: float = 1058.0
    pellet_limit_k: float = 2700.0
    collector_limit_k: float = 1023.0
    moderator_limit_k: float = 930.0
    reflector_limit_k: float = 1000.0


def set_tec_update_interval(core: Any, interval_s: float) -> float:
    interval = float(interval_s)
    if not math.isfinite(interval) or interval <= 0.0:
        raise ValueError('TEC update interval must be finite and positive')
    if not hasattr(core, 'thermo_update_interval'):
        raise ValueError('core does not expose thermo_update_interval')
    tolerance = min(0.5 * interval, max(1.0e-12, interval * 1.0e-9))
    scheduler_threshold = interval - tolerance
    core.thermo_update_interval = scheduler_threshold
    return scheduler_threshold


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


def restore_or_trigger_accident(
        gaps: Dict[str, tuple[Any, int]],
        source_config: Dict[str, Any],
        current_time_s: float,
) -> Dict[str, Any]:
    active = read_source_accident_state(source_config)
    if active:
        if 'helium_accident_time_absolute_s' not in source_config:
            raise ValueError('active helium accident missing absolute event time')
        event_time = float(source_config['helium_accident_time_absolute_s'])
        triggered_now = False
    else:
        event_time = float(current_time_s)
        triggered_now = True
    if not math.isfinite(event_time):
        raise ValueError('helium accident time must be finite')
    set_helium_h_eq(gaps, HELIUM_H_FINAL_W_M2K)
    return {
        'helium_accident_active': True,
        'triggered_now': triggered_now,
        'accident_time_absolute_s': event_time,
        'h_before_W_m2K': HELIUM_H_INITIAL_W_M2K,
        'h_after_W_m2K': HELIUM_H_FINAL_W_M2K,
        'affected_representatives': list(REPRESENTATIVE_NAMES),
        'physical_tfe_count': sum(EXPECTED_MULTIPLIERS),
    }


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
        collector_surface, _ = gap.bound1.get_coupling_surface_snapshot()
        inner_clad_surface, _ = gap.bound2.get_coupling_surface_snapshot()
        q_out = float(np.sum(
            (
                np.asarray(collector_surface, dtype=float)
                - np.asarray(inner_clad_surface, dtype=float)
            ) / np.asarray(gap.R_gap_total, dtype=float)
        ))
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


def _limits_from_config(config: HeliumAccidentRunConfig) -> Dict[str, float]:
    return {
        'channel_wall': float(config.wall_limit_k),
        'pellet': float(config.pellet_limit_k),
        'collector': float(config.collector_limit_k),
        'moderator': float(config.moderator_limit_k),
        'reflector': float(config.reflector_limit_k),
    }


def _refresh_gap_diagnostics(
        build: Dict[str, Any], gaps: Dict[str, tuple[Any, int]]) -> None:
    for gap, _ in gaps.values():
        gap.sync()
    build['system']._refresh_solid_boundary_cache(
        update_flux=True,
        current_time=float(build['system'].global_time),
    )


def collect_all_metrics(
        build: Dict[str, Any],
        gaps: Dict[str, tuple[Any, int]],
        *,
        handoff_type: str,
        initial_power_w: float,
        accident_time_s: float,
        active: bool,
        dt_s: float,
        limits_k: Dict[str, float],
) -> tuple[Dict[str, Any], list[Dict[str, Any]]]:
    row = collect_reactivity_metrics(
        build,
        handoff_type=handoff_type,
        initial_power_w=initial_power_w,
        dt_s=dt_s,
    )
    row.update(collect_helium_metrics(
        build,
        gaps,
        accident_time_s=accident_time_s,
        active=active,
    ))
    peaks = collect_temperature_peaks(build['core'], limits_k)
    for component in (
            'channel_wall', 'pellet', 'collector', 'moderator', 'reflector'):
        values = [
            float(peak['actual_k']) for peak in peaks
            if peak.get('source_component', peak['component']) == component
        ]
        row[f'{component}_max_T_K'] = float(np.max(values))
    return row, peaks


def _append_history(path: Path, row: Dict[str, Any]) -> None:
    write_header = not path.exists()
    with path.open('a', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=ACCIDENT_HISTORY_FIELDS,
            extrasaction='ignore',
        )
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _print_progress(row: Dict[str, Any]) -> None:
    print(
        ('[t={time_s:.3f}s accident={accident_elapsed_s:.3f}s] '
         'P={core_total_power_W:.3f}W '
         'rho_T={effective_temperature_feedback:.6e} '
         'Twall={channel_wall_max_T_K:.3f}K '
         'Tpellet={pellet_max_T_K:.3f}K '
         'Tcollector={collector_max_T_K:.3f}K '
         'Tmod={moderator_max_T_K:.3f}K '
         'Tref={reflector_max_T_K:.3f}K').format(**row),
        flush=True,
    )


def _nonfinite_metric_trip(row: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    for field, raw_value in row.items():
        if (
                raw_value is None
                or isinstance(raw_value, (str, bool, np.bool_))
                or not np.isscalar(raw_value)):
            continue
        try:
            value = float(raw_value)
        except (TypeError, ValueError, OverflowError):
            continue
        if not math.isfinite(value):
            return {
                'component': 'nonfinite_metric',
                'field': field,
                'actual': value,
            }
    return None


def evaluate_state_trip(
        row: Dict[str, Any],
        *,
        peaks: Sequence[Dict[str, Any]],
        config: HeliumAccidentRunConfig,
        baseline_power_w: float,
        require_tec_convergence: bool,
) -> tuple[str, Optional[Dict[str, Any]]]:
    trip = _nonfinite_metric_trip(row)
    if trip is not None:
        return 'nonfinite_metric', trip
    if not bool(row.get('fluid_converged', False)):
        return 'hydraulic_nonconvergence', {
            'component': 'fluid_solver',
            'actual': row.get('fluid_converged'),
            'required': True,
        }
    if (
            require_tec_convergence
            and not bool(row.get('tec_main_converged', False))):
        return 'tec_nonconvergence', {
            'component': 'tec_main',
            'actual': row.get('tec_main_converged'),
            'required': True,
        }
    trip = find_limit_trip(peaks)
    if trip is not None:
        return 'temperature_limit', trip
    power_limit_w = baseline_power_w * float(config.max_power_factor)
    if float(row['core_total_power_W']) > power_limit_w:
        return 'maximum_power_factor', {
            'component': 'core_total_power',
            'actual_w': float(row['core_total_power_W']),
            'limit_w': power_limit_w,
        }
    if (
            config.min_fluid_temperature_stop_k is not None
            and float(config.min_fluid_temperature_stop_k) > 0.0
            and float(row['min_fluid_T_K'])
            < float(config.min_fluid_temperature_stop_k)):
        return 'low_fluid_temperature', {
            'component': 'minimum_fluid_temperature',
            'actual_k': float(row['min_fluid_T_K']),
            'limit_k': float(config.min_fluid_temperature_stop_k),
        }
    return 'completed', None


def _validate_accident_config(config: HeliumAccidentRunConfig) -> None:
    runtime = ReactivityControlRunConfig(
        restart_in=config.restart_in,
        output_dir=config.output_dir,
        duration_s=config.duration_s,
        dt_s=config.dt_s,
        record_interval_s=config.record_interval_s,
        checkpoint_interval_s=config.checkpoint_interval_s,
        min_fluid_temperature_stop_k=config.min_fluid_temperature_stop_k,
        max_power_factor=config.max_power_factor,
    )
    _validate_runtime(runtime)
    if (
            not math.isfinite(float(config.tec_update_interval_s))
            or float(config.tec_update_interval_s) <= 0.0):
        raise ValueError('TEC update interval must be finite and positive')
    for name, value in _limits_from_config(config).items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f'{name} temperature limit must be positive')


def run_helium_accident(config: HeliumAccidentRunConfig) -> Dict[str, Any]:
    _validate_accident_config(config)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / 'history.csv'
    if history_path.exists():
        raise FileExistsError(f'output history already exists: {history_path}')

    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in),
        output_dir=out_dir,
        duration_s=float(config.duration_s),
        dt_s=float(config.dt_s),
        record_interval_s=float(config.record_interval_s),
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=config.min_fluid_temperature_stop_k,
        max_power_factor=float(config.max_power_factor),
    )
    debug, source_config = load_baseline_debug_config(runtime)
    if source_config.get('external_heat_enabled') is not False:
        raise ValueError('helium accident baseline must have external heat disabled')
    build = build_debug_case(debug, apply_fixed_power=False)
    system = build['system']
    core = build['core']
    tec_scheduler_threshold_s = set_tec_update_interval(
        core, config.tec_update_interval_s)
    handoff_type = prepare_reactivity_control(
        core,
        source_point_kinetics_enabled=bool(source_config['point_kinetics_enabled']),
        expected_power_w=float(debug.power_w),
    )
    gaps = collect_helium_gaps(build)
    limits_k = _limits_from_config(config)
    baseline_power_w = float(debug.power_w)
    start_time = float(system.global_time)
    end_time = start_time + float(config.duration_s)
    source_active = read_source_accident_state(source_config)
    _refresh_gap_diagnostics(build, gaps)

    latest: Dict[str, Any]
    if not source_active:
        latest, initial_peaks = collect_all_metrics(
            build,
            gaps,
            handoff_type=handoff_type,
            initial_power_w=baseline_power_w,
            accident_time_s=start_time,
            active=False,
            dt_s=0.0,
            limits_k=limits_k,
        )
        initial_reason, initial_trip = evaluate_state_trip(
            latest,
            peaks=initial_peaks,
            config=config,
            baseline_power_w=baseline_power_w,
            require_tec_convergence=False,
        )
        if initial_reason != 'completed':
            _write_json(out_dir / 'limit_trip.json', {
                **initial_trip,
                'phase': 'initial_preflight',
                'stop_reason': initial_reason,
                'time_s': start_time,
            })
            raise RuntimeError(f'initial state violates limit: {initial_trip}')
        _append_history(history_path, latest)
        _print_progress(latest)

    event = restore_or_trigger_accident(
        gaps,
        source_config=source_config,
        current_time_s=start_time,
    )
    _refresh_gap_diagnostics(build, gaps)
    _write_json(out_dir / 'accident_event.json', event)

    run_config = dict(source_config)
    run_config.update({
        'case': CASE_NAME,
        'duration_s': float(config.duration_s),
        'stage_durations_s': [float(config.duration_s)],
        'dt_s': float(config.dt_s),
        'tec_update_interval_s': float(config.tec_update_interval_s),
        'tec_scheduler_threshold_s': float(tec_scheduler_threshold_s),
        'record_interval_s': float(config.record_interval_s),
        'checkpoint_interval_s': float(config.checkpoint_interval_s),
        'min_fluid_temperature_stop_k': config.min_fluid_temperature_stop_k,
        'max_power_factor': float(config.max_power_factor),
        'restart_in': str(config.restart_in),
        'source_run_config': str(Path(config.restart_in).parent / 'run_config.json'),
        'handoff_type': handoff_type,
        'point_kinetics_enabled': True,
        'reactivity_control_mode': 'temperature_feedback_only',
        'external_reactivity': 0.0,
        'control_drum_enabled': False,
        'external_heat_enabled': False,
        'helium_accident_active': True,
        'helium_accident_model': 'instantaneous_total_loss_of_gas_conduction',
        'helium_accident_time_absolute_s': event['accident_time_absolute_s'],
        'helium_h_initial_w_m2k': HELIUM_H_INITIAL_W_M2K,
        'helium_h_final_w_m2k': HELIUM_H_FINAL_W_M2K,
        'temperature_limits_k': limits_k,
    })
    _write_json(out_dir / 'run_config.json', run_config)

    last_record_time = start_time
    last_checkpoint_time = start_time
    stop_reason = 'completed'
    trip_payload = None

    if source_active:
        latest, restart_peaks = collect_all_metrics(
            build,
            gaps,
            handoff_type=handoff_type,
            initial_power_w=baseline_power_w,
            accident_time_s=float(event['accident_time_absolute_s']),
            active=True,
            dt_s=0.0,
            limits_k=limits_k,
        )
        _append_history(history_path, latest)
        _print_progress(latest)
        stop_reason, trip_payload = evaluate_state_trip(
            latest,
            peaks=restart_peaks,
            config=config,
            baseline_power_w=baseline_power_w,
            require_tec_convergence=False,
        )
        if stop_reason != 'completed':
            emergency_path = out_dir / 'emergency_restart.npz'
            system.save_global_state(str(emergency_path))
            _write_json(out_dir / 'limit_trip.json', {
                **trip_payload,
                'phase': 'restart_preflight',
                'stop_reason': stop_reason,
                'time_s': float(system.global_time),
                'accident_elapsed_s': (
                    float(system.global_time)
                    - float(event['accident_time_absolute_s'])
                ),
            })

    while stop_reason == 'completed':
        dt = _next_step_dt(
            current_time=float(system.global_time),
            end_time=end_time,
            requested_dt=float(config.dt_s),
        )
        if dt is None:
            break
        system.step(
            dt,
            inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=False,
            fluid_max_iter=int(debug.fluid_max_iter),
            reactivity_control=0.0,
        )
        latest, peaks = collect_all_metrics(
            build,
            gaps,
            handoff_type=handoff_type,
            initial_power_w=baseline_power_w,
            accident_time_s=float(event['accident_time_absolute_s']),
            active=True,
            dt_s=dt,
            limits_k=limits_k,
        )

        stop_reason, trip_payload = evaluate_state_trip(
            latest,
            peaks=peaks,
            config=config,
            baseline_power_w=baseline_power_w,
            require_tec_convergence=True,
        )

        should_record = (
            float(system.global_time) - last_record_time
            >= float(config.record_interval_s) - 1.0e-9
            or _next_step_dt(
                current_time=float(system.global_time),
                end_time=end_time,
                requested_dt=float(config.dt_s),
            ) is None
            or stop_reason != 'completed'
        )
        if should_record:
            _append_history(history_path, latest)
            _print_progress(latest)
            last_record_time = float(system.global_time)

        if stop_reason != 'completed':
            emergency_path = out_dir / 'emergency_restart.npz'
            system.save_global_state(str(emergency_path))
            _write_json(out_dir / 'limit_trip.json', {
                **trip_payload,
                'stop_reason': stop_reason,
                'time_s': float(system.global_time),
                'accident_elapsed_s': (
                    float(system.global_time)
                    - float(event['accident_time_absolute_s'])
                ),
            })
            break

        if (
                float(config.checkpoint_interval_s) > 0.0
                and float(system.global_time) - last_checkpoint_time
                >= float(config.checkpoint_interval_s) - 1.0e-9):
            checkpoint_path = out_dir / f'checkpoint_t{system.global_time:.3f}s.npz'
            system.save_global_state(str(checkpoint_path))
            last_checkpoint_time = float(system.global_time)

    restart_path = out_dir / 'stage_01_restart.npz'
    system.save_global_state(str(restart_path))
    result = {
        'case': CASE_NAME,
        'output_dir': str(out_dir),
        'history_path': str(history_path),
        'restart_path': str(restart_path),
        'source_restart_path': str(config.restart_in),
        'handoff_type': handoff_type,
        'start_time_s': start_time,
        'end_time_s': float(system.global_time),
        'accident_time_absolute_s': float(event['accident_time_absolute_s']),
        'stop_reason': stop_reason,
        'latest_metrics': latest,
    }
    summary_path = out_dir / 'run_summary.json'
    _write_json(summary_path, result)
    _write_json(out_dir / 'latest_state.json', {
        'case': CASE_NAME,
        'latest_restart_path': str(restart_path),
        'latest_summary_path': str(summary_path),
        'latest_metrics': latest,
        'stop_reason': stop_reason,
    })
    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restart-in', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--duration', type=float, default=100.0)
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--tec-update-interval', type=float, default=0.05)
    parser.add_argument('--record-interval', type=float, default=0.1)
    parser.add_argument('--checkpoint-interval', type=float, default=10.0)
    parser.add_argument('--min-fluid-temperature-stop', type=float, default=500.0)
    parser.add_argument('--max-power-factor', type=float, default=2.0)
    parser.add_argument('--wall-limit-k', type=float, default=1058.0)
    parser.add_argument('--pellet-limit-k', type=float, default=2700.0)
    parser.add_argument('--collector-limit-k', type=float, default=1023.0)
    parser.add_argument('--moderator-limit-k', type=float, default=930.0)
    parser.add_argument('--reflector-limit-k', type=float, default=1000.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    min_fluid = (
        None
        if float(args.min_fluid_temperature_stop) <= 0.0
        else float(args.min_fluid_temperature_stop)
    )
    result = run_helium_accident(HeliumAccidentRunConfig(
        restart_in=args.restart_in,
        output_dir=args.output_dir,
        duration_s=float(args.duration),
        dt_s=float(args.dt),
        tec_update_interval_s=float(args.tec_update_interval),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        min_fluid_temperature_stop_k=min_fluid,
        max_power_factor=float(args.max_power_factor),
        wall_limit_k=float(args.wall_limit_k),
        pellet_limit_k=float(args.pellet_limit_k),
        collector_limit_k=float(args.collector_limit_k),
        moderator_limit_k=float(args.moderator_limit_k),
        reflector_limit_k=float(args.reflector_limit_k),
    ))
    print(json.dumps(
        result['latest_metrics'],
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
    ))
    print(f"Stop reason: {result['stop_reason']}")
    print(f"Saved outputs to: {result['output_dir']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
