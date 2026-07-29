"""Restartable V14 operating-point calibration runner."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Solvers.Couplers import FluidSolidCouple  # noqa: E402
from testModule.Full_Loop_Cases import (  # noqa: E402
    FullLoopCoreConfig,
    FullLoopFlowConfig,
    FullLoopPumpConfig,
    V14HeatPipeRadiatorConfig,
    build_v14_case_a_system,
)
from testModule.Full_Loop_Cases.common_config import ReservedParallelTecConfig  # noqa: E402


INITIAL_TEMPERATURE_K = 727.0
SPACE_TEMPERATURE_K = 4.0
CORE_POWER_W = 115000.0
TARGET_FLOW_KG_S = 1.3
MAIN_VOLTAGE_V = 27.2
PARALLEL_VOLTAGE_V = 0.35
WIRE_RESISTANCE_BASE_OHM = np.array(
    [0.001552, 0.001024, 0.000336, 0.000608], dtype=float
)
LOOKUP_DB = REPO_ROOT / "ThermoCalc" / "emission_runtime_db_v2" / "pcs_0p02_5torr"
ORBIT_PERIOD_S = 6552.0


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )


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


def _configure_numerics(build: Dict[str, Any]) -> int:
    system = build["system"]
    for solid in system.solid_components.values():
        solid.set_ode_method("implicit_euler")
        solid.implicit_fallback_to_solve_ivp = False
    count = 0
    for coupler in system.couplers:
        if isinstance(coupler, FluidSolidCouple) and coupler.solid_node_capacitance is not None:
            coupler.set_coupling_time_scheme("local_implicit")
            count += 1
    return count


def _suppress_area_change_acceleration(build: Dict[str, Any]) -> None:
    """Remove the unsigned W**2 area-change source for this V14 runner only."""
    net = build["system"].fluid_solver
    original = net._calc_momentum_coeffs

    def calculate_without_acceleration(dt: float, W_old_frozen=None):
        original(dt, W_old_frozen=W_old_frozen)
        normal = ~net.is_inlet_junction_mask
        if not np.any(normal):
            return
        flow = net.W_vec if W_old_frozen is None else np.asarray(W_old_frozen, dtype=float)
        idx_from = net.idx_from_vec[normal]
        idx_to = net.idx_to_vec[normal]
        multiplier_from = net.M_from_vec[normal]
        multiplier_to = net.M_to_vec[normal]
        area_from = np.maximum(net.A_from_node_vec[normal], 1.0e-10)
        area_to = np.maximum(net.A_to_node_vec[normal], 1.0e-10)
        rho_from = np.maximum(net.rho_vec[idx_from], 1.0e-10)
        rho_to = np.maximum(net.rho_vec[idx_to], 1.0e-10)
        acceleration_dp = 0.5 * flow[normal] ** 2 * (
            multiplier_from ** 2 / (rho_from * area_from ** 2)
            - multiplier_to ** 2 / (rho_to * area_to ** 2)
        )
        net.B_coeffs[normal] -= net.A_coeffs[normal] * acceleration_dp

    net._calc_momentum_coeffs = calculate_without_acceleration


def _apply_core_power(build: Dict[str, Any]) -> None:
    build["core"].update_neutronic_power(
        p_total=CORE_POWER_W, p_fiss=CORE_POWER_W, p_decay=0.0, alpha=1.0
    )


def build_case(
    *,
    emissivity: float,
    pump_head_pa: float,
    flow_control: bool,
    external_heat: bool,
) -> Dict[str, Any]:
    if not 0.0 <= float(emissivity) <= 1.0:
        raise ValueError("emissivity must be within [0, 1]")
    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(
            inlet_temperature_k=INITIAL_TEMPERATURE_K,
            main_tec_enabled=False,
            reserved_parallel_tec=ReservedParallelTecConfig(enabled=False),
        ),
        flow_config=FullLoopFlowConfig(total_flow_kg_s=TARGET_FLOW_KG_S),
        pump_config=FullLoopPumpConfig(
            pump_total_head_pa=float(pump_head_pa),
            pump_flow_control=bool(flow_control),
            target_flow_kg_s=TARGET_FLOW_KG_S if flow_control else None,
        ),
        radiator_config=V14HeatPipeRadiatorConfig(
            t_space_k=SPACE_TEMPERATURE_K,
            hp_initial_temp_k=INITIAL_TEMPERATURE_K,
            hp_emissivity=float(emissivity),
            fin_emissivity=float(emissivity),
            ring_emissivity=0.2,
            external_heat_enabled=bool(external_heat),
            external_heat_scale_factor=1.0,
        ),
    )
    build["core"].point_reactor = None
    build["core"].enable_tec_coupled = False
    build["operating_emissivity"] = float(emissivity)
    build["flow_control"] = bool(flow_control)
    for solid in build["system"].solid_components.values():
        solid.set_ode_method("implicit_euler")
    _apply_core_power(build)
    return build


def _set_pump_head(build: Dict[str, Any], total_head_pa: float) -> None:
    single = 0.5 * float(total_head_pa)
    build["pump_a"].set_delta_p(single)
    build["pump_b"].set_delta_p(single)
    build["pump_total_head_pa"] = float(total_head_pa)
    build["pump_single_head_pa"] = single


def _solve_flow(build: Dict[str, Any], total_head_pa: float) -> float:
    _set_pump_head(build, total_head_pa)
    net = build["system"].fluid_solver
    net.initialize_hydraulics(dt=0.1, tol=1.0e-8, max_iter=1500, omega=0.5)
    return 0.5 * (float(build["pump_a"].W) + float(build["pump_b"].W))


def calibrate_pump_head(build: Dict[str, Any], initial_head_pa: float) -> Dict[str, float]:
    low = max(10.0, 0.5 * float(initial_head_pa))
    high = max(100.0, 1.5 * float(initial_head_pa))
    flow_low = _solve_flow(build, low)
    flow_high = _solve_flow(build, high)
    for _ in range(30):
        if flow_low <= TARGET_FLOW_KG_S:
            break
        low *= 0.5
        flow_low = _solve_flow(build, low)
    for _ in range(30):
        if flow_high >= TARGET_FLOW_KG_S:
            break
        high *= 1.5
        flow_high = _solve_flow(build, high)
    if not flow_low <= TARGET_FLOW_KG_S <= flow_high:
        raise RuntimeError(
            f"Could not bracket 1.3 kg/s: ({low}, {flow_low}) to ({high}, {flow_high})"
        )
    for _ in range(36):
        mid = 0.5 * (low + high)
        if _solve_flow(build, mid) < TARGET_FLOW_KG_S:
            low = mid
        else:
            high = mid
    head = 0.5 * (low + high)
    flow = _solve_flow(build, head)
    return {
        "pump_total_head_pa": float(head),
        "pump_single_head_pa": 0.5 * float(head),
        "flow_kg_s": float(flow),
        "flow_error_kg_s": float(flow - TARGET_FLOW_KG_S),
    }


def configure_tec(
    build: Dict[str, Any], wire_scale: float, update_interval_s: float
) -> Dict[str, Any]:
    if float(wire_scale) < 0.0:
        raise ValueError("wire_scale must be non-negative")
    if not LOOKUP_DB.exists():
        raise FileNotFoundError(f"TEC lookup database not found: {LOOKUP_DB}")
    core = build["core"]
    core.tec_lookup_enabled = True
    core.tec_lookup_db = str(LOOKUP_DB)
    core.tec_lookup_regions = ("core", "startup", "high_power", "accident")
    core.enable_tec_coupled = True
    core._build_thermo_calc()
    core.setup_tec_circuit("fixed_u", MAIN_VOLTAGE_V, I_guess=185.0, topology="series")
    core.setup_reserved_parallel_tec_circuit(
        mode_str="fixed_u",
        target_value=PARALLEL_VOLTAGE_V,
        I_guess=500.0,
        multipliers={"Ring3_Open": 3},
    )
    wire = WIRE_RESISTANCE_BASE_OHM * float(wire_scale)
    for group in core.iter_tec_circuit_groups():
        thermo = group.thermo_calc
        thermo._input_data.resistanceWire = np.tile(wire, (thermo.N_elem, 1))
        thermo.build()
    core.thermo_update_interval = float(update_interval_s)
    core.post_step(0.0, float(build["system"].global_time))
    for group in core.iter_tec_circuit_groups():
        group.thermo_calc.calculate(verbose=False)
    return {
        "wire_scale": float(wire_scale),
        "wire_resistance_ohm": wire.tolist(),
        "lookup_db": str(LOOKUP_DB),
        "update_interval_s": float(update_interval_s),
    }


def _tec_metrics(build: Dict[str, Any], name: str) -> Dict[str, Any]:
    result = build["core"].get_tec_circuit_global_results().get(name) or {}
    voltage = float(result.get("Uout", 0.0))
    current = float(result.get("Iout", 0.0))
    return {
        f"tec_{name}_voltage_V": voltage,
        f"tec_{name}_current_A": current,
        f"tec_{name}_power_W": voltage * current,
        f"tec_{name}_converged": bool(result.get("converged", False)),
        f"tec_{name}_iterations": int(result.get("iteration_count", 0)),
        f"tec_{name}_zero_emission_skipped": bool(
            result.get("zero_emission_skipped", False)
        ),
    }


def _v14_radiator_metrics(build: Dict[str, Any]) -> Dict[str, float]:
    symmetric = float(build["radiator_config"].symmetric_ring_multiplier)
    hp_rejection = symmetric * sum(
        ring.get_total_heat_rejection_scaled() for ring in build["ring_hps"]
    )
    wall_single = sum(
        max(0.0, -float(np.sum(solid.boundaries["right"].current_flux)))
        for solid in build["ring_solids"]
    )
    external = 0.0
    if build.get("external_heat_enabled", False):
        external = symmetric * sum(
            ring.get_total_external_heat_absorption_scaled(build["system"].global_time)
            for ring in build["ring_hps"]
        )
    return {
        "radiator_heatpipe_rejection_W": float(hp_rejection),
        "radiator_ring_wall_rejection_W": float(symmetric * wall_single),
        "radiator_total_rejection_W": float(hp_rejection + symmetric * wall_single),
        "radiator_external_heat_W": float(external),
    }


def collect_metrics(build: Dict[str, Any], dt_s: float) -> Dict[str, Any]:
    system = build["system"]
    net = system.fluid_solver
    solids = np.concatenate(
        [np.asarray(solid.T, dtype=float).reshape(-1) for solid in system.solid_components.values()]
    )
    radiator_out = build["radiator_outlet_header"].volumes[-1]
    pump_out = build["pump_outlet_node"]
    row: Dict[str, Any] = {
        "time_s": float(system.global_time),
        "dt_s": float(dt_s),
        "core_power_W": float(build["core"].last_total_core_power),
        "core_inlet_T_K": float(build["core_inlet_connector"].T),
        "core_outlet_T_K": float(build["core_outlet_connector"].T),
        "core_delta_T_K": float(
            build["core_outlet_connector"].T - build["core_inlet_connector"].T
        ),
        "pump_flow_kg_s": float(build["pump_a"].W),
        "pump_total_head_setting_Pa": float(build["pump_a"].delta_p + build["pump_b"].delta_p),
        "pump_section_pressure_rise_Pa": float(pump_out.P - radiator_out.P),
        "loop_pressure_range_Pa": float(np.max(net.P_vec) - np.min(net.P_vec)),
        "min_fluid_T_K": float(np.min(net.T_vec)),
        "max_fluid_T_K": float(np.max(net.T_vec)),
        "min_solid_T_K": float(np.min(solids)),
        "max_solid_T_K": float(np.max(solids)),
        "mean_solid_T_K": float(np.mean(solids)),
        "radiator_emissivity": float(build["operating_emissivity"]),
        "external_heat_enabled": bool(build.get("external_heat_enabled", False)),
        **_v14_radiator_metrics(build),
    }
    row.update(_tec_metrics(build, "main"))
    row.update(_tec_metrics(build, "reserved_parallel"))
    return row


def _steady_diagnostic(history: list[Dict[str, Any]], window_s: float) -> Dict[str, Any]:
    if len(history) < 3:
        return {"steady": False, "window_s": 0.0}
    end = float(history[-1]["time_s"])
    rows = [row for row in history if float(row["time_s"]) >= end - float(window_s)]
    if len(rows) < 3:
        return {"steady": False, "window_s": 0.0}
    times = np.asarray([row["time_s"] for row in rows], dtype=float)
    result: Dict[str, Any] = {"window_s": float(times[-1] - times[0])}
    steady = result["window_s"] >= 0.8 * float(window_s)
    for key in ("core_inlet_T_K", "core_outlet_T_K", "mean_solid_T_K"):
        values = np.asarray([row[key] for row in rows], dtype=float)
        slope = float(np.polyfit(times - times[0], values, 1)[0])
        span = float(np.max(values) - np.min(values))
        result[f"{key}_slope_K_s"] = slope
        result[f"{key}_span_K"] = span
        steady = steady and abs(slope) <= 5.0e-4 and span <= 0.75
    result["steady"] = bool(steady)
    return result


def run_stage(args: argparse.Namespace) -> Dict[str, Any]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    flow_control = args.stage == "thermal-flow"
    external_heat = args.stage == "external"
    tec_enabled = args.stage in ("coupled", "external")
    build = build_case(
        emissivity=float(args.emissivity),
        pump_head_pa=float(args.pump_head),
        flow_control=flow_control,
        external_heat=external_heat,
    )
    system = build["system"]
    if args.restart is None:
        if args.stage != "thermal-flow":
            raise ValueError(f"--restart is required for stage {args.stage}")
        _set_uniform_state(build, INITIAL_TEMPERATURE_K)
        system.initialize_system(dt_init=0.01, tol=1.0e-4, max_iter=2000)
        _configure_numerics(build)
    else:
        for coupler in system.couplers:
            if isinstance(coupler, FluidSolidCouple):
                coupler.set_coupling_time_scheme("current")
        system.load_global_state(str(args.restart))
        _configure_numerics(build)
    _apply_core_power(build)

    if args.suppress_area_change_acceleration:
        _suppress_area_change_acceleration(build)

    pump_calibration = None
    measured_flow_control_head = None
    if flow_control:
        measured_flow_control_head = float(
            build["pump_outlet_node"].P - build["radiator_outlet_header"].volumes[-1].P
        )
    else:
        _set_pump_head(build, float(args.pump_head))
        pump_calibration = {
            "pump_total_head_pa": float(args.pump_head),
            "pump_single_head_pa": 0.5 * float(args.pump_head),
            "flow_kg_s_at_stage_start": float(build["pump_a"].W),
        }
    hydraulic_reconditioning = None
    if args.recondition_hydraulics:
        net = system.fluid_solver
        converged = net.initialize_hydraulics(
            dt=0.02, tol=1.0e-6, max_iter=2000, omega=0.5
        )
        if not converged:
            raise RuntimeError("V14 hydraulic restart reconditioning did not converge")
        hydraulic_reconditioning = {
            "converged": True,
            "pump_a_flow_kg_s": float(build["pump_a"].W),
            "pump_b_flow_kg_s": float(build["pump_b"].W),
            "pressure_range_Pa": float(np.ptp(net.P_vec)),
        }

    wire = None
    if tec_enabled:
        wire = configure_tec(
            build, float(args.wire_scale), float(args.tec_update_interval)
        )

    history_path = output_dir / "history.csv"
    latest_path = output_dir / "latest_state.json"
    restart_path = output_dir / "latest_restart.npz"
    start_time = float(system.global_time)
    end_time = start_time + float(args.duration)
    next_record = start_time
    next_restart = start_time + float(args.restart_interval)
    fields: Optional[list[str]] = None
    history: list[Dict[str, Any]] = []
    stop_reason = "duration_limit"

    while system.global_time < end_time - 1.0e-10:
        _apply_core_power(build)
        dt = system.compute_adaptive_dt(
            min_dt=1.0e-4,
            max_dt=float(args.max_dt),
            safety_factor=0.8,
            respect_fluid_cfl=bool(args.respect_fluid_cfl),
        )
        dt = min(float(dt), float(args.max_dt), end_time - float(system.global_time))
        if next_record > system.global_time:
            dt = min(dt, next_record - float(system.global_time))
        while True:
            try:
                system.step(dt, inner_iter=1, fail_on_fluid_nonconvergence=True, fluid_max_iter=300)
                break
            except RuntimeError as exc:
                if "Fluid solver NOT converged" not in str(exc) or dt <= 1.0e-4:
                    raise
                dt *= 0.5
        _apply_core_power(build)

        if system.global_time >= next_record - 1.0e-10:
            row = collect_metrics(build, dt)
            numeric = [float(value) for value in row.values() if not isinstance(value, (bool, str))]
            if not all(math.isfinite(value) for value in numeric):
                raise FloatingPointError(f"Non-finite state: {row}")
            history.append(row)
            fields = list(row) if fields is None else fields
            write_header = not history_path.exists()
            with history_path.open("a", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                if write_header:
                    writer.writeheader()
                writer.writerow(row)
            _write_json(latest_path, row)
            print(json.dumps(row, sort_keys=True), flush=True)
            next_record += float(args.record_interval)

            if flow_control:
                measured_flow_control_head = float(row["pump_section_pressure_rise_Pa"])
            elapsed = float(system.global_time) - start_time
            steady = _steady_diagnostic(history, float(args.steady_window))
            if (
                args.stop_when_steady
                and not external_heat
                and elapsed >= float(args.minimum_duration)
                and steady["steady"]
            ):
                stop_reason = "near_steady"
                break

        if system.global_time >= next_restart - 1.0e-10:
            checkpoint = output_dir / f"restart_t{int(round(system.global_time)):06d}s.npz"
            system.save_global_state(str(checkpoint))
            system.save_global_state(str(restart_path))
            next_restart += float(args.restart_interval)

    if not history or history[-1]["time_s"] != float(system.global_time):
        history.append(collect_metrics(build, 0.0))
    system.save_global_state(str(restart_path))
    steady = _steady_diagnostic(history, float(args.steady_window))
    summary = {
        "case": "V14_operating_steady",
        "stage": args.stage,
        "restart_in": None if args.restart is None else str(args.restart),
        "restart_out": str(restart_path),
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "stop_reason": stop_reason,
        "initial_temperature_K": INITIAL_TEMPERATURE_K,
        "core_power_W": CORE_POWER_W,
        "target_flow_kg_s": TARGET_FLOW_KG_S,
        "radiator_emissivity": float(args.emissivity),
        "ring_wall_emissivity": 0.2,
        "external_heat_enabled": external_heat,
        "external_heat_period_s": ORBIT_PERIOD_S if external_heat else None,
        "pump_calibration": pump_calibration,
        "hydraulic_reconditioning": hydraulic_reconditioning,
        "flow_control_measured_total_head_pa": measured_flow_control_head,
        "tec_enabled": tec_enabled,
        "tec_lookup_enabled": tec_enabled,
        "tec_main_target_voltage_V": MAIN_VOLTAGE_V if tec_enabled else None,
        "tec_parallel_target_voltage_V": PARALLEL_VOLTAGE_V if tec_enabled else None,
        "wire": wire,
        "numerics": {
            "solid_ode_method": "implicit_euler",
            "fluid_solid_coupling": "local_implicit",
            "area_change_acceleration_suppressed": bool(
                args.suppress_area_change_acceleration
            ),
            "respect_fluid_cfl": bool(args.respect_fluid_cfl),
            "max_dt_s": float(args.max_dt),
        },
        "steady_diagnostic": steady,
        "latest": history[-1],
    }
    _write_json(output_dir / "summary.json", summary)
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("thermal-flow", "thermal-head", "coupled", "external"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--duration", type=float, default=20000.0)
    parser.add_argument("--max-dt", type=float, default=0.5)
    parser.add_argument("--record-interval", type=float, default=10.0)
    parser.add_argument("--restart-interval", type=float, default=200.0)
    parser.add_argument("--minimum-duration", type=float, default=2000.0)
    parser.add_argument("--steady-window", type=float, default=1000.0)
    parser.add_argument("--stop-when-steady", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--pump-head", type=float, default=7900.0)
    parser.add_argument("--emissivity", type=float, default=0.75)
    parser.add_argument("--wire-scale", type=float, default=1.0)
    parser.add_argument("--tec-update-interval", type=float, default=0.8)
    parser.add_argument(
        "--respect-fluid-cfl",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Apply the hydraulic-network fluid CFL limit when choosing the time step.",
    )
    parser.add_argument(
        "--suppress-area-change-acceleration",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Suppress the unsigned W**2 area-change source in this V14 runner.",
    )
    parser.add_argument(
        "--recondition-hydraulics",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Re-equilibrate restart pressure and flow before coupled advancement.",
    )
    return parser


if __name__ == "__main__":
    result = run_stage(_parser().parse_args())
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False), flush=True)
