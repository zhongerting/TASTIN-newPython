"""Restartable V15 operating-point calibration and orbital-heat evaluation."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from testModule.Full_Loop_Cases import (  # noqa: E402
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    ReservedParallelTecConfig,
    V15PipeFinRadiatorConfig,
    build_v15_case_a_system,
)


from Solvers.Couplers import FluidSolidCouple  # noqa: E402

BASE_WIRE_RESISTANCE_OHM = np.array([0.001552, 0.001024, 0.000336, 0.000608])
LOOKUP_DB = REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr"
LOOKUP_REGIONS = ("core", "startup", "high_power", "accident")
DEFAULT_ROOT = Path(__file__).resolve().parent / "runs"


@dataclass(frozen=True)
class RunConfig:
    phase: str
    output_dir: Path
    restart_in: Optional[Path] = None
    duration_s: float = 12000.0
    dt_s: float = 0.5
    record_interval_s: float = 10.0
    checkpoint_interval_s: float = 200.0
    initial_temperature_k: float = 727.0
    space_temperature_k: float = 4.0
    core_power_w: float = 115000.0
    target_flow_kg_s: float = 1.3
    pump_total_head_pa: float = 7900.0
    surface_emissivity: float = 0.8
    wire_resistance_scale: float = 1.0
    main_voltage_v: float = 27.2
    main_current_guess_a: float = 185.0
    parallel_voltage_v: float = 0.35
    parallel_current_guess_a: float = 500.0
    thermo_update_interval_s: float = 0.8
    steady_window_s: float = 1000.0
    steady_minimum_time_s: float = 2000.0
    steady_temperature_span_k: float = 0.5
    steady_rejection_relative_span: float = 0.01

    @property
    def tec_enabled(self) -> bool:
        return self.phase in ("coupled", "orbit")

    @property
    def external_heat_enabled(self) -> bool:
        return self.phase == "orbit"

    @property
    def fixed_flow(self) -> bool:
        return self.phase == "thermal-flow"


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def _finite(value: Any, default: float = float("nan")) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _set_uniform_state(build: Dict[str, Any], temperature_k: float) -> None:
    system = build["system"]
    net = system.fluid_solver
    for volume in net.volumes_obj:
        volume.T = float(temperature_k)
        volume.h = float(volume.material.enthalpy(volume.T, volume.P))
        volume.update_properties(volume.material)
    net._initialize_state_from_objects()
    net._update_fluid_properties()
    net._sync_vectors_to_objects()
    for solid in system.solid_components.values():
        solid.T[...] = float(temperature_k)
        if hasattr(solid, "dTdt"):
            solid.dTdt[...] = 0.0
        solid.current_time = float(system.global_time)
        solid._update_properties()
        solid._update_boundaries_state(current_time=float(system.global_time))
        solid.set_ode_method("implicit_euler")
    for unit in build["radiator_units"]:
        unit.last_fin_temperature[...] = float(temperature_k)
        unit.last_fin_effective_temperature_distribution[...] = float(temperature_k)


def _apply_core_power(build: Dict[str, Any], power_w: float) -> None:
    build["core"].update_neutronic_power(
        p_total=float(power_w), p_fiss=float(power_w), p_decay=0.0, alpha=1.0
    )


def _set_fluid_solid_scheme(system: Any, scheme: str) -> None:
    for coupler in system.couplers:
        if isinstance(coupler, FluidSolidCouple) and coupler.solid_node_capacitance is not None:
            coupler.set_coupling_time_scheme(scheme)


def _apply_wire_resistance(core: Any, scale: float) -> list[float]:
    values = BASE_WIRE_RESISTANCE_OHM * float(scale)
    for group in core.iter_tec_circuit_groups():
        model = group.thermo_calc
        model._input_data.resistanceWire = np.tile(values, (model.N_elem, 1))
        model.build()
    return values.tolist()


def _configure_tec(build: Dict[str, Any], config: RunConfig) -> list[float]:
    core = build["core"]
    core.enable_tec_coupled = True
    core.tec_lookup_enabled = True
    core.tec_lookup_db = str(LOOKUP_DB)
    core.tec_lookup_regions = LOOKUP_REGIONS
    core._build_thermo_calc()
    core.setup_tec_circuit(
        "fixed_u", config.main_voltage_v, I_guess=config.main_current_guess_a, topology="series"
    )
    core.setup_reserved_parallel_tec_circuit(
        "fixed_u",
        config.parallel_voltage_v,
        I_guess=config.parallel_current_guess_a,
        multipliers={"Ring3_Open": 3},
    )
    core.thermo_update_interval = float(config.thermo_update_interval_s)
    core.post_step(0.0, float(build["system"].global_time))
    values = _apply_wire_resistance(core, config.wire_resistance_scale)
    core.setup_tec_circuit(
        "fixed_u", config.main_voltage_v, I_guess=config.main_current_guess_a, topology="series"
    )
    reserved = core.tec_circuit_groups["reserved_parallel"]
    reserved.thermo_calc.setup_circuit_mode(
        "parallel_fixed_u",
        config.parallel_voltage_v,
        config.parallel_current_guess_a,
    )
    core.post_step(0.0, float(build["system"].global_time))
    core.set_thermo_update_time(float(build["system"].global_time) - core.thermo_update_interval)
    return values


def build_case(config: RunConfig) -> Dict[str, Any]:
    core_config = FullLoopCoreConfig(
        inlet_temperature_k=config.initial_temperature_k,
        main_tec_enabled=False,
        reserved_parallel_tec=ReservedParallelTecConfig(enabled=False),
    )
    build = build_v15_case_a_system(
        core_config=core_config,
        flow_config=FullLoopFlowConfig(total_flow_kg_s=config.target_flow_kg_s),
        pump_config=FullLoopPumpConfig(
            pump_total_head_pa=max(1.0, config.pump_total_head_pa),
            pump_flow_control=config.fixed_flow,
            target_flow_kg_s=config.target_flow_kg_s if config.fixed_flow else None,
        ),
        radiator_config=V15PipeFinRadiatorConfig(
            tube_emissivity=config.surface_emissivity,
            fin_emissivity=config.surface_emissivity,
            t_space_k=config.space_temperature_k,
            external_heat_enabled=config.external_heat_enabled,
            external_heat_scale_factor=1.0,
            solid_ode_method="implicit_euler",
            fluid_solid_coupling_scheme="local_implicit",
        ),
    )
    system = build["system"]
    _set_uniform_state(build, config.initial_temperature_k)
    system.initialize_system(dt_init=0.01, tol=1.0e-4, max_iter=1000)
    if config.restart_in is not None:
        _set_fluid_solid_scheme(system, "current")
        system.load_global_state(str(config.restart_in))
    _set_fluid_solid_scheme(system, "local_implicit")
    for solid in system.solid_components.values():
        solid.set_ode_method("implicit_euler")
    _apply_core_power(build, config.core_power_w)
    build["wire_resistance_ohm"] = (
        _configure_tec(build, config)
        if config.tec_enabled
        else (BASE_WIRE_RESISTANCE_OHM * config.wire_resistance_scale).tolist()
    )
    return build


def _pump_metrics(build: Dict[str, Any]) -> Dict[str, float]:
    pump_a = build["pump_a"]
    pump_b = build["pump_b"]
    head_a = _finite(pump_a.to_vol.P) - _finite(pump_a.from_vol.P)
    head_b = _finite(pump_b.to_vol.P) - _finite(pump_b.from_vol.P)
    return {
        "pump_a_flow_kg_s": _finite(pump_a.W),
        "pump_b_flow_kg_s": _finite(pump_b.W),
        "pump_a_required_head_Pa": head_a,
        "pump_b_required_head_Pa": head_b,
        "pump_required_head_total_Pa": head_a + head_b,
    }


def _tec_metrics(build: Dict[str, Any]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    groups = build["core"].get_tec_circuit_global_results() if build["core"].enable_tec_coupled else {}
    for output_name, group_name in (("main", "main"), ("parallel", "reserved_parallel")):
        values = groups.get(group_name) or {}
        current = _finite(values.get("Iout"), 0.0)
        voltage = _finite(values.get("Uout"), 0.0)
        result.update({
            f"tec_{output_name}_current_A": current,
            f"tec_{output_name}_voltage_V": voltage,
            f"tec_{output_name}_power_W": current * voltage,
            f"tec_{output_name}_converged": bool(values.get("converged", not build["core"].enable_tec_coupled)),
            f"tec_{output_name}_iterations": int(values.get("iteration_count", 0)),
        })
    return result


def collect_metrics(build: Dict[str, Any], config: RunConfig, dt_s: float) -> Dict[str, Any]:
    system = build["system"]
    net = system.fluid_solver
    solids = np.concatenate([
        np.asarray(solid.T, dtype=float).reshape(-1) for solid in system.solid_components.values()
    ])
    radiator_q = 0.0
    external_q = 0.0
    for unit in build["radiator_units"]:
        radiator_q += float(np.sum(unit.get_heat_exchange_breakdown()["gross_rejection"]))
        if hasattr(unit, "get_external_heat_absorption_distribution"):
            external_q += float(np.sum(
                unit.get_external_heat_absorption_distribution(float(system.global_time))[2]
            ))
    diag = system.last_step_diagnostics or {}
    flags = diag.get("fluid_converged_by_iteration", [])
    return {
        "time_s": float(system.global_time),
        "dt_s": float(dt_s),
        "phase": config.phase,
        "surface_emissivity": float(config.surface_emissivity),
        "effective_emissivity": 0.585 * float(config.surface_emissivity),
        "wire_resistance_scale": float(config.wire_resistance_scale),
        "core_power_W": float(build["core"].last_total_core_power),
        "core_inlet_T_K": float(build["core_inlet_connector"].T),
        "core_outlet_T_K": float(build["core_outlet_connector"].T),
        "core_delta_T_K": float(build["core_outlet_connector"].T - build["core_inlet_connector"].T),
        "radiator_gross_rejection_W": radiator_q,
        "radiator_external_absorption_W": external_q,
        "min_fluid_T_K": float(np.min(net.T_vec)),
        "max_fluid_T_K": float(np.max(net.T_vec)),
        "min_solid_T_K": float(np.min(solids)),
        "max_solid_T_K": float(np.max(solids)),
        "mean_solid_T_K": float(np.mean(solids)),
        "fluid_converged": bool(all(flags)) if flags else True,
        **_pump_metrics(build),
        **_tec_metrics(build),
    }


class HistoryWriter:
    def __init__(self, output_dir: Path, build: Dict[str, Any], samples_per_chunk: int = 60):
        self.path = output_dir / "state_history"
        self.path.mkdir(parents=True, exist_ok=True)
        self.build = build
        self.samples_per_chunk = max(1, int(samples_per_chunk))
        self.samples: list[Dict[str, Any]] = []
        self.index = len(list(self.path.glob("state_*.npz")))

    def append(self, metrics: Dict[str, Any]) -> None:
        units = self.build["radiator_units"]
        self.samples.append({
            "metrics": metrics,
            "fin_temperature_K": np.stack([unit.last_fin_temperature for unit in units]),
            "fin_radiation_W": np.stack([unit.last_fin_radiation_distribution for unit in units]),
            "fin_absorption_W": np.stack([unit.last_fin_absorption_distribution for unit in units]),
            "fin_net_from_root_W": np.stack([unit.last_fin_net_from_root_distribution for unit in units]),
            "radiation_background_K": np.stack([unit.radiation_background_temperature for unit in units]),
            "tube_wall_temperature_K": np.stack([unit.get_temperature_distribution() for unit in units]),
            "tube_fluid_temperature_K": np.stack([
                np.asarray([volume.T for volume in channel.volumes], dtype=float)
                for channel in self.build["radiator_tube_channels"]
            ]),
            "tube_flow_kg_s": np.asarray([
                channel.internal_junctions[0].W for channel in self.build["radiator_tube_channels"]
            ]),
        })
        if len(self.samples) >= self.samples_per_chunk:
            self.flush()

    def flush(self) -> None:
        if not self.samples:
            return
        payload: Dict[str, Any] = {
            "metrics_json": np.asarray([json.dumps(row["metrics"], sort_keys=True) for row in self.samples]),
            "radiator_unit_names": np.asarray([unit.name for unit in self.build["radiator_units"]]),
        }
        for key in self.samples[0]:
            if key != "metrics":
                payload[key] = np.stack([row[key] for row in self.samples])
        start = int(round(self.samples[0]["metrics"]["time_s"]))
        end = int(round(self.samples[-1]["metrics"]["time_s"]))
        np.savez_compressed(self.path / f"state_{self.index:04d}_t{start:06d}_to_t{end:06d}.npz", **payload)
        self.index += 1
        self.samples.clear()


def _steady(rows: Iterable[Dict[str, Any]], config: RunConfig, phase_start: float) -> bool:
    rows = list(rows)
    if not rows or rows[-1]["time_s"] - phase_start < config.steady_minimum_time_s:
        return False
    if rows[-1]["time_s"] - rows[0]["time_s"] < config.steady_window_s:
        return False
    for key in ("core_inlet_T_K", "core_outlet_T_K", "mean_solid_T_K"):
        values = np.asarray([row[key] for row in rows], dtype=float)
        if np.ptp(values) > config.steady_temperature_span_k:
            return False
    q = np.asarray([row["radiator_gross_rejection_W"] for row in rows], dtype=float)
    return float(np.ptp(q) / max(abs(float(np.mean(q))), 1.0)) <= config.steady_rejection_relative_span


def run_phase(config: RunConfig) -> Dict[str, Any]:
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    serial_config = asdict(config)
    serial_config["output_dir"] = str(config.output_dir)
    serial_config["restart_in"] = str(config.restart_in) if config.restart_in else None
    serial_config.update({
        "solid_ode_method": "implicit_euler",
        "fluid_solid_coupling_scheme": "local_implicit",
        "lookup_db": str(LOOKUP_DB),
        "wire_resistance_ohm": (BASE_WIRE_RESISTANCE_OHM * config.wire_resistance_scale).tolist(),
    })
    _write_json(output_dir / "run_config.json", serial_config)
    build = build_case(config)
    system = build["system"]
    phase_start = float(system.global_time)
    end_time = phase_start + float(config.duration_s)
    history_path = output_dir / "history.csv"
    writer = HistoryWriter(output_dir, build)
    window: deque[Dict[str, Any]] = deque()
    last_record = -float("inf")
    last_checkpoint = phase_start
    latest = collect_metrics(build, config, 0.0)
    fields = list(latest)
    with history_path.open("w", newline="", encoding="utf-8") as handle:
        csv.DictWriter(handle, fieldnames=fields).writeheader()
    stop_reason = "duration"
    while system.global_time < end_time - 1.0e-12:
        dt = min(config.dt_s, end_time - float(system.global_time))
        _apply_core_power(build, config.core_power_w)
        system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=False, fluid_max_iter=100)
        _apply_core_power(build, config.core_power_w)
        if system.global_time - last_record >= config.record_interval_s - 1.0e-12 or system.global_time >= end_time - 1.0e-12:
            latest = collect_metrics(build, config, dt)
            if not all(math.isfinite(float(value)) for key, value in latest.items() if key not in ("phase", "fluid_converged", "tec_main_converged", "tec_parallel_converged")):
                raise FloatingPointError(f"Non-finite V15 state: {latest}")
            with history_path.open("a", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=fields).writerow(latest)
            writer.append(latest)
            print(json.dumps(latest, sort_keys=True), flush=True)
            last_record = float(system.global_time)
            window.append(latest)
            while window and latest["time_s"] - window[0]["time_s"] > config.steady_window_s:
                window.popleft()
            if config.phase != "orbit" and _steady(window, config, phase_start):
                stop_reason = "near_steady"
                break
        if config.checkpoint_interval_s > 0.0 and system.global_time - last_checkpoint >= config.checkpoint_interval_s - 1.0e-12:
            checkpoint = output_dir / f"checkpoint_t{int(round(system.global_time)):07d}s.npz"
            system.save_global_state(str(checkpoint))
            last_checkpoint = float(system.global_time)
            _write_json(output_dir / "latest_checkpoint.json", {"path": str(checkpoint), "metrics": latest})
    writer.flush()
    restart = output_dir / "final_restart.npz"
    system.save_global_state(str(restart))
    summary = {
        "phase": config.phase,
        "stop_reason": stop_reason,
        "start_time_s": phase_start,
        "end_time_s": float(system.global_time),
        "restart_path": str(restart),
        "latest_metrics": latest,
        "wire_resistance_ohm": build["wire_resistance_ohm"],
        "lookup_enabled": bool(config.tec_enabled and all(
            group.thermo_calc.lookup_enabled for group in build["core"].iter_tec_circuit_groups()
        )),
    }
    _write_json(output_dir / "summary.json", summary)
    return summary


def _score(row: Dict[str, Any]) -> float:
    return (
        abs(row["core_inlet_T_K"] - 727.0) / 8.0
        + abs(row["core_outlet_T_K"] - 823.0) / 8.0
        + abs(row["tec_main_power_W"] - 5000.0) / 400.0
    )


def run_auto(root: Path, max_iterations: int, coupled_duration_s: float) -> Dict[str, Any]:
    root.mkdir(parents=True, exist_ok=True)
    thermal = run_phase(RunConfig(phase="thermal-flow", output_dir=root / "01_thermal_fixed_flow"))
    required_head = float(thermal["latest_metrics"]["pump_required_head_total_Pa"])
    if not math.isfinite(required_head) or required_head <= 100.0:
        raise RuntimeError(f"Invalid measured pump head: {required_head}")
    fixed = run_phase(RunConfig(
        phase="fixed-head",
        output_dir=root / "02_fixed_head_confirmation",
        restart_in=Path(thermal["restart_path"]),
        duration_s=3000.0,
        pump_total_head_pa=required_head,
    ))
    restart = Path(fixed["restart_path"])
    emissivity = 0.8
    wire_scale = 1.0
    candidates = []
    for index in range(1, max_iterations + 1):
        candidate = run_phase(RunConfig(
            phase="coupled",
            output_dir=root / f"03_coupled_iter_{index:02d}_eps{emissivity:.4f}_wire{wire_scale:.4f}",
            restart_in=restart,
            duration_s=coupled_duration_s,
            dt_s=0.05,
            pump_total_head_pa=required_head,
            surface_emissivity=emissivity,
            wire_resistance_scale=wire_scale,
        ))
        row = candidate["latest_metrics"]
        candidate.update({"surface_emissivity": emissivity, "wire_resistance_scale": wire_scale, "score": _score(row)})
        candidates.append(candidate)
        restart = Path(candidate["restart_path"])
        temperature_ok = abs(row["core_inlet_T_K"] - 727.0) <= 8.0 and abs(row["core_outlet_T_K"] - 823.0) <= 8.0
        power_ok = abs(row["tec_main_power_W"] - 5000.0) <= 400.0
        parallel_ok = abs(row["tec_parallel_voltage_V"] - 0.35) <= 0.01 and row["tec_parallel_converged"]
        if temperature_ok and power_ok and parallel_ok and row["tec_main_converged"]:
            break
        emissivity = float(np.clip(emissivity + np.clip((row["core_inlet_T_K"] - 727.0) / 120.0, -0.05, 0.05), 0.2, 1.0))
        if row["tec_main_power_W"] > 0.0:
            wire_scale = float(np.clip(wire_scale * np.clip(row["tec_main_power_W"] / 5000.0, 0.85, 1.15), 0.1, 3.0))
    best = min(candidates, key=lambda item: item["score"])
    orbit = run_phase(RunConfig(
        phase="orbit",
        output_dir=root / "04_orbit_N78_two_periods",
        restart_in=Path(best["restart_path"]),
        duration_s=13104.0,
        dt_s=0.1,
        pump_total_head_pa=required_head,
        surface_emissivity=float(best["surface_emissivity"]),
        wire_resistance_scale=float(best["wire_resistance_scale"]),
    ))
    result = {
        "measured_total_pump_head_Pa": required_head,
        "pump_head_each_Pa": 0.5 * required_head,
        "candidates": candidates,
        "best": best,
        "orbit": orbit,
    }
    _write_json(root / "auto_summary.json", result)
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("phase", choices=("thermal-flow", "fixed-head", "coupled", "orbit", "auto"))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_ROOT / "manual")
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--duration", type=float, default=12000.0)
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--checkpoint-interval", type=float, default=200.0)
    parser.add_argument("--pump-head", type=float, default=7900.0)
    parser.add_argument("--emissivity", type=float, default=0.8)
    parser.add_argument("--wire-scale", type=float, default=1.0)
    parser.add_argument("--max-iterations", type=int, default=6)
    parser.add_argument("--coupled-duration", type=float, default=5000.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.phase == "auto":
        result = run_auto(Path(args.output_dir), int(args.max_iterations), float(args.coupled_duration))
    else:
        result = run_phase(RunConfig(
            phase=args.phase,
            output_dir=Path(args.output_dir),
            restart_in=args.restart,
            duration_s=float(args.duration),
            dt_s=float(args.dt),
            record_interval_s=float(args.record_interval),
            checkpoint_interval_s=float(args.checkpoint_interval),
            pump_total_head_pa=float(args.pump_head),
            surface_emissivity=float(args.emissivity),
            wire_resistance_scale=float(args.wire_scale),
        ))
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
