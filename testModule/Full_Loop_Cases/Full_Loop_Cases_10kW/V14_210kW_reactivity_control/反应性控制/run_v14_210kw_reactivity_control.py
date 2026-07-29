'''Run the V14 210 kW baseline with zero external reactivity.'''

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, Optional, Sequence, Tuple

REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    DebugRunConfig,
    HISTORY_FIELDS,
    build_debug_case,
    collect_metrics,
)

CASE_NAME = 'V14_10kW_210kW_temperature_feedback_only'
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / 'runs' / 'default'
REACTIVITY_FIELDS = [
    'handoff_type', 'point_kinetics_enabled', 'core_total_power_W',
    'fission_power_W', 'decay_power_W', 'power_relative_change',
    'external_reactivity', 'control_drum_enabled', 'control_drum_reactivity',
    'feedback_fuel', 'feedback_electrode', 'feedback_moderator',
    'feedback_reflector', 'feedback_total_absolute',
    'feedback_reference_total', 'effective_temperature_feedback',
    'total_reactivity',
]
REACTIVITY_HISTORY_FIELDS = HISTORY_FIELDS + REACTIVITY_FIELDS


@dataclass(frozen=True)
class ReactivityControlRunConfig:
    restart_in: Path
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = 10.0
    dt_s: float = 0.05
    record_interval_s: float = 1.0
    checkpoint_interval_s: float = 10.0
    min_fluid_temperature_stop_k: Optional[float] = 500.0
    max_power_factor: float = 2.0


_REQUIRED_SOURCE_KEYS = (
    'fluid_max_iter',
    'hp_up_view_factor',
    'initial_temperature_k',
    'inner_iter',
    'lower_hp_down_view_factor',
    'point_kinetics_enabled',
    'power_w',
    'radiator_emissivity',
    'space_temperature_k',
    'target_flow_kg_s',
    'tec_current_guess_a',
    'tec_electrical_calculation_enabled',
    'tec_lookup_enabled',
    'tec_lookup_regions',
    'tec_voltage_v',
    'upper_hp_down_view_factor',
    'wire_resistance_scale',
)


def load_baseline_debug_config(
        runtime: ReactivityControlRunConfig) -> Tuple[DebugRunConfig, Dict[str, Any]]:
    restart = Path(runtime.restart_in)
    if not restart.is_file():
        raise FileNotFoundError(f'restart not found: {restart}')
    config_path = restart.parent / 'run_config.json'
    if not config_path.is_file():
        raise FileNotFoundError(f'run_config.json not found beside restart: {config_path}')

    source = json.loads(config_path.read_text(encoding='utf-8'))
    missing = [key for key in _REQUIRED_SOURCE_KEYS if key not in source]
    if missing:
        raise ValueError(f'baseline run_config.json missing keys: {missing}')
    power_w = float(source['power_w'])
    if not math.isfinite(power_w) or power_w <= 0.0:
        raise ValueError('baseline power_w must be finite and positive')

    debug = DebugRunConfig(
        output_dir=Path(runtime.output_dir),
        stage_durations_s=(float(runtime.duration_s),),
        dt_s=float(runtime.dt_s),
        record_interval_s=float(runtime.record_interval_s),
        checkpoint_interval_s=float(runtime.checkpoint_interval_s),
        min_fluid_temperature_stop_k=runtime.min_fluid_temperature_stop_k,
        tec_electrical_enabled=bool(source['tec_electrical_calculation_enabled']),
        tec_voltage_v=float(source['tec_voltage_v']),
        tec_current_guess_a=float(source['tec_current_guess_a']),
        tec_lookup_enabled=bool(source['tec_lookup_enabled']),
        tec_lookup_db=source.get('tec_lookup_db'),
        tec_lookup_regions=tuple(source['tec_lookup_regions']),
        wire_resistance_scale=float(source['wire_resistance_scale']),
        radiator_emissivity=float(source['radiator_emissivity']),
        hp_up_view_factor=float(source['hp_up_view_factor']),
        upper_hp_down_view_factor=float(source['upper_hp_down_view_factor']),
        lower_hp_down_view_factor=float(source['lower_hp_down_view_factor']),
        inner_iter=int(source['inner_iter']),
        fluid_max_iter=int(source['fluid_max_iter']),
        power_w=power_w,
        target_flow_kg_s=float(source['target_flow_kg_s']),
        initial_temperature_k=float(source['initial_temperature_k']),
        space_temperature_k=float(source['space_temperature_k']),
        restart_in=restart,
    )
    return debug, source


def prepare_reactivity_control(
        core: Any,
        *,
        source_point_kinetics_enabled: bool,
        expected_power_w: float) -> str:
    has_point_reactor = bool(core.has_point_reactor)
    if has_point_reactor != bool(source_point_kinetics_enabled):
        raise ValueError('restart point-kinetics state does not match run_config.json')
    if bool(core.control_drum_reactivity_model.enabled):
        raise ValueError('control drum must be disabled for this runner')

    if has_point_reactor:
        return 'reactivity_continuation'

    loaded_power = float(core.last_total_core_power)
    if not math.isfinite(loaded_power) or not math.isclose(
            loaded_power, float(expected_power_w), rel_tol=1.0e-6, abs_tol=1.0e-6):
        raise ValueError(
            f'restart core power {loaded_power} W does not match run_config power '
            f'{float(expected_power_w)} W'
        )
    # 点堆初态设置
    core.initialize_point_reactor(total_power_initial=loaded_power)
    effective_feedback = float(core.get_effective_reactivity_feedback())
    if not math.isfinite(effective_feedback) or abs(effective_feedback) > 1.0e-12:
        raise RuntimeError(
            f'point-kinetics handoff did not zero feedback: {effective_feedback}'
        )
    return 'fixed_power_handoff'


def collect_reactivity_diagnostics(
        core: Any,
        *,
        handoff_type: str,
        initial_power_w: float) -> Dict[str, Any]:
    point = core.point_reactor
    feedback = core.compute_reactivity_feedback()
    reference_total = float(core.feedback_reference_result.total)
    effective_feedback = float(feedback.total) - reference_total
    drum_reactivity = float(core.get_control_drum_reactivity())
    total_power = float(point.total_power)
    external_reactivity = 0.0
    return {
        'handoff_type': str(handoff_type),
        'point_kinetics_enabled': True,
        'core_total_power_W': total_power,
        'fission_power_W': float(point.fission_power),
        'decay_power_W': float(point.decay_power),
        'power_relative_change': (
            total_power - float(initial_power_w)
        ) / float(initial_power_w),
        'external_reactivity': external_reactivity,
        'control_drum_enabled': bool(core.control_drum_reactivity_model.enabled),
        'control_drum_reactivity': drum_reactivity,
        'feedback_fuel': float(feedback.fuel),
        'feedback_electrode': float(feedback.electrode),
        'feedback_moderator': float(feedback.moderator),
        'feedback_reflector': float(feedback.reflector),
        'feedback_total_absolute': float(feedback.total),
        'feedback_reference_total': reference_total,
        'effective_temperature_feedback': effective_feedback,
        'total_reactivity': external_reactivity + drum_reactivity + effective_feedback,
    }


def collect_reactivity_metrics(
        build: Dict[str, Any],
        *,
        handoff_type: str,
        initial_power_w: float,
        dt_s: float) -> Dict[str, Any]:
    row = collect_metrics(build, stage_index=1, dt_s=dt_s)
    row.update(collect_reactivity_diagnostics(
        build['core'],
        handoff_type=handoff_type,
        initial_power_w=initial_power_w,
    ))
    return row


def _validate_runtime(config: ReactivityControlRunConfig) -> None:
    for name in ('duration_s', 'dt_s', 'record_interval_s', 'max_power_factor'):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f'{name} must be finite and positive')
    if not math.isfinite(float(config.checkpoint_interval_s)):
        raise ValueError('checkpoint_interval_s must be finite')


def _next_step_dt(
        *,
        current_time: float,
        end_time: float,
        requested_dt: float) -> Optional[float]:
    remaining = float(end_time) - float(current_time)
    if remaining <= 1.0e-9:
        return None
    return min(float(requested_dt), remaining)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False),
        encoding='utf-8',
    )


def _append_history(path: Path, row: Dict[str, Any]) -> None:
    write_header = not path.exists()
    with path.open('a', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=REACTIVITY_HISTORY_FIELDS,
            extrasaction='ignore',
        )
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _print_progress(row: Dict[str, Any]) -> None:
    print(
        ('[t={time_s:.3f}s] P={core_total_power_W:.3f} W '
         'dP/P0={power_relative_change:.6e} '
         'rho_T={effective_temperature_feedback:.6e} '
         'rho_ext={external_reactivity:.1f} '
         'rho_drum={control_drum_reactivity:.1f} '
         'Tfluid=[{min_fluid_T_K:.3f}, {max_fluid_T_K:.3f}] K').format(**row),
        flush=True,
    )


def run_reactivity_control(config: ReactivityControlRunConfig) -> Dict[str, Any]:
    _validate_runtime(config)
    debug, source_config = load_baseline_debug_config(config)
    # 不再每一步强制给堆芯施加固定的 210 kW，而是让点堆动力学决定功率。
    build = build_debug_case(debug, apply_fixed_power=False)
    system = build['system']
    core = build['core']
    handoff_type = prepare_reactivity_control(
        core,
        source_point_kinetics_enabled=bool(source_config['point_kinetics_enabled']),
        expected_power_w=float(debug.power_w),
    )

    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    history_path = out_dir / 'history.csv'
    if history_path.exists():
        raise FileExistsError(f'output history already exists: {history_path}')

    run_config = dict(source_config)
    run_config.update({
        'case': CASE_NAME,
        'duration_s': float(config.duration_s),
        'stage_durations_s': [float(config.duration_s)],
        'dt_s': float(config.dt_s),
        'record_interval_s': float(config.record_interval_s),
        'checkpoint_interval_s': float(config.checkpoint_interval_s),
        'min_fluid_temperature_stop_k': config.min_fluid_temperature_stop_k,
        'max_power_factor': float(config.max_power_factor),
        'restart_in': str(config.restart_in),
        'source_run_config': str(
            Path(config.restart_in).parent / 'run_config.json'
        ),
        'handoff_type': handoff_type,
        'point_kinetics_enabled': True,
        'reactivity_control_mode': 'temperature_feedback_only',
        'external_reactivity': 0.0,
        'control_drum_enabled': False,
    })
    _write_json(out_dir / 'run_config.json', run_config)

    baseline_power_w = float(debug.power_w)
    start_time = float(system.global_time)
    end_time = start_time + float(config.duration_s)
    last_record_time = start_time
    last_checkpoint_time = start_time
    stop_reason = 'completed'

    latest = collect_reactivity_metrics(
        build,
        handoff_type=handoff_type,
        initial_power_w=baseline_power_w,
        dt_s=0.0,
    )
    _append_history(history_path, latest)
    _print_progress(latest)

    while True:
        dt = _next_step_dt(
            current_time=float(system.global_time),
            end_time=end_time,
            requested_dt=float(config.dt_s),
        )
        if dt is None:
            break
        # 点堆反应性推进功率
        system.step(
            dt,
            inner_iter=int(debug.inner_iter),
            fail_on_fluid_nonconvergence=False,
            fluid_max_iter=int(debug.fluid_max_iter),
            reactivity_control=0.0,
        )
        latest = collect_reactivity_metrics(
            build,
            handoff_type=handoff_type,
            initial_power_w=baseline_power_w,
            dt_s=dt,
        )

        power = float(latest['core_total_power_W'])
        if not math.isfinite(power):
            stop_reason = 'nonfinite_core_power'
        elif power > baseline_power_w * float(config.max_power_factor):
            stop_reason = 'maximum_power_factor'
        else:
            threshold = config.min_fluid_temperature_stop_k
            if threshold is not None and float(threshold) > 0.0:
                if float(latest['min_fluid_T_K']) < float(threshold):
                    stop_reason = 'low_fluid_temperature'

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
            break

        if (
                float(config.checkpoint_interval_s) > 0.0
                and float(system.global_time) - last_checkpoint_time
                >= float(config.checkpoint_interval_s) - 1.0e-9):
            checkpoint_path = out_dir / (
                f'checkpoint_t{system.global_time:.3f}s.npz'
            )
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
    })
    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restart-in', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--duration', type=float, default=10.0)
    parser.add_argument('--dt', type=float, default=0.05)
    parser.add_argument('--record-interval', type=float, default=1.0)
    parser.add_argument('--checkpoint-interval', type=float, default=10.0)
    parser.add_argument(
        '--min-fluid-temperature-stop',
        type=float,
        default=500.0,
        help='Emergency stop threshold in K; <=0 disables it.',
    )
    parser.add_argument('--max-power-factor', type=float, default=2.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    threshold = (
        None
        if float(args.min_fluid_temperature_stop) <= 0.0
        else float(args.min_fluid_temperature_stop)
    )
    result = run_reactivity_control(ReactivityControlRunConfig(
        restart_in=args.restart_in,
        output_dir=args.output_dir,
        duration_s=float(args.duration),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        min_fluid_temperature_stop_k=threshold,
        max_power_factor=float(args.max_power_factor),
    ))
    print(json.dumps(
        result['latest_metrics'],
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
    ))
    print('Stop reason: {}'.format(result['stop_reason']))
    print('Saved outputs to: {}'.format(result['output_dir']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
