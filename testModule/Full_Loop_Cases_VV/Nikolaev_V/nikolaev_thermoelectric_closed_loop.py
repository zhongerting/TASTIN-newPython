from __future__ import annotations

import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

CASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CASE_DIR.parents[2]
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nikolaev_source_data import Table2OperatingPoint
from nikolaev_thermocalc_model import (
    NikolaevThermalNetworkConfig,
    NikolaevTfeGeometry,
    build_thermocalc_case,
)
from nikolaev_thermocalc_runner import apply_case_to_thermocalc

BOLTZMANN_EV_PER_K = 8.617e-5


@dataclass(frozen=True)
class ClosedLoopConfig:
    max_iterations: int = 30
    relaxation: float = 0.35
    temperature_tolerance_k: float = 0.25
    collector_to_boundary_resistance_k_per_w: float = 0.010
    i_guess_a: float = 300.0
    min_temperature_k: float = 500.0
    max_temperature_k: float = 2600.0


@dataclass(frozen=True)
class ThermoElectricFeedback:
    electron_emitter_power_w: np.ndarray
    electron_collector_power_w: np.ndarray
    electron_emitter_flux_w_m2: np.ndarray
    electron_collector_flux_w_m2: np.ndarray
    joule_power_emitter_w: np.ndarray
    joule_power_collector_w: np.ndarray


@dataclass(frozen=True)
class ClosedLoopIteration:
    iteration: int
    current_a: float
    uout_v: float
    max_temperature_change_k: float
    emitter_temperature_mean_k: float
    collector_temperature_mean_k: float
    electron_cooling_power_w: float
    collector_electron_heating_power_w: float
    joule_power_emitter_w: float
    joule_power_collector_w: float
    converged: bool


@dataclass(frozen=True)
class ClosedLoopResult:
    voltage_v: float
    thermal_power_kw: float
    current_exp_a: float
    current_calc_a: float
    current_error_a: float
    electric_power_exp_w: float
    electric_power_calc_w: float
    electric_power_error_w: float
    efficiency_exp_percent: float
    efficiency_calc_percent: float
    efficiency_error_percent: float
    emitter_temperature_exp_k: float
    initial_emitter_temperature_mean_k: float
    emitter_temperature_mean_k: float
    emitter_temperature_min_k: float
    emitter_temperature_max_k: float
    collector_temperature_mean_k: float
    collector_temperature_min_k: float
    collector_temperature_max_k: float
    electron_cooling_power_w: float
    collector_electron_heating_power_w: float
    joule_power_emitter_w: float
    joule_power_collector_w: float
    outer_iterations: int
    max_temperature_change_k: float
    closed_loop_converged: bool
    thermocalc_converged: bool
    finite: bool
    iteration_history: list[ClosedLoopIteration] = field(default_factory=list)


def _default_model_factory(n_elements: int, n_nodes: int):
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    return ThermoCalcModel(n_elements=n_elements, n_nodes=n_nodes)


def _as_node_array(name: str, value, n_nodes: int, default: float = 0.0) -> np.ndarray:
    if value is None:
        arr = np.full(n_nodes, default, dtype=float)
    else:
        arr = np.asarray(value, dtype=float)
    if arr.shape != (n_nodes,):
        raise ValueError(f"{name} must have shape {(n_nodes,)}, got {arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def extract_thermoelectric_feedback(tec_result: dict, geometry: NikolaevTfeGeometry) -> ThermoElectricFeedback:
    n = geometry.n_nodes
    j_a_cm2 = _as_node_array("J", tec_result.get("J"), n)
    current_density_a_m2 = j_a_cm2 * 1.0e4
    phi_e = _as_node_array("phiE", tec_result.get("phiE"), n)
    te = _as_node_array("TE", tec_result.get("TE"), n)
    ue = _as_node_array("UE", tec_result.get("UE"), n)
    uc = _as_node_array("UC", tec_result.get("UC"), n)
    joule_e = _as_node_array("joulePowerE", tec_result.get("joulePowerE"), n)
    joule_c = _as_node_array("joulePowerC", tec_result.get("joulePowerC"), n)

    electron_emitter_flux = -current_density_a_m2 * (phi_e + 2.0 * BOLTZMANN_EV_PER_K * te)
    electron_collector_flux = current_density_a_m2 * (phi_e + 2.0 * BOLTZMANN_EV_PER_K * te - (ue - uc))
    emitter_area = np.full(n, geometry.emitter_side_area_total_m2 / n, dtype=float)

    return ThermoElectricFeedback(
        electron_emitter_power_w=electron_emitter_flux * emitter_area,
        electron_collector_power_w=electron_collector_flux * emitter_area,
        electron_emitter_flux_w_m2=electron_emitter_flux,
        electron_collector_flux_w_m2=electron_collector_flux,
        joule_power_emitter_w=joule_e.copy(),
        joule_power_collector_w=joule_c.copy(),
    )


def update_temperature_fields(
    *,
    heat_source_w: np.ndarray,
    feedback: ThermoElectricFeedback,
    thermal_config: NikolaevThermalNetworkConfig,
    closed_loop_config: ClosedLoopConfig,
    old_emitter_k: np.ndarray,
    old_collector_k: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    n = heat_source_w.size
    r_ec_node = float(thermal_config.emitter_to_collector_resistance_k_per_w) * n
    r_cb_node = float(closed_loop_config.collector_to_boundary_resistance_k_per_w) * n

    emitter_net_to_collector_w = heat_source_w + feedback.electron_emitter_power_w + feedback.joule_power_emitter_w
    collector_net_to_boundary_w = (
        emitter_net_to_collector_w
        + feedback.electron_collector_power_w
        + feedback.joule_power_collector_w
    )

    target_collector = thermal_config.collector_boundary_temperature_k + collector_net_to_boundary_w * r_cb_node
    target_emitter = target_collector + emitter_net_to_collector_w * r_ec_node

    target_emitter = np.clip(target_emitter, closed_loop_config.min_temperature_k, closed_loop_config.max_temperature_k)
    target_collector = np.clip(target_collector, closed_loop_config.min_temperature_k, closed_loop_config.max_temperature_k)

    alpha = float(closed_loop_config.relaxation)
    if not 0.0 < alpha <= 1.0:
        raise ValueError("relaxation must be in (0, 1].")
    new_emitter = old_emitter_k + alpha * (target_emitter - old_emitter_k)
    new_collector = old_collector_k + alpha * (target_collector - old_collector_k)
    max_change = float(max(np.max(np.abs(new_emitter - old_emitter_k)), np.max(np.abs(new_collector - old_collector_k))))
    return new_emitter, new_collector, max_change


def solve_closed_loop_point(
    point: Table2OperatingPoint,
    thermal_config: NikolaevThermalNetworkConfig = NikolaevThermalNetworkConfig(),
    closed_loop_config: ClosedLoopConfig = ClosedLoopConfig(),
    model_factory: Callable[[int, int], object] = _default_model_factory,
) -> ClosedLoopResult:
    if closed_loop_config.max_iterations <= 0:
        raise ValueError("max_iterations must be positive.")
    case_model = build_thermocalc_case(point, thermal_config)
    geometry = case_model.geometry
    thermo_model = model_factory(1, geometry.n_nodes)
    apply_case_to_thermocalc(thermo_model, case_model, point.voltage_v, closed_loop_config.i_guess_a)

    emitter = np.asarray(case_model.arrays.temitter_k[0], dtype=float).copy()
    collector = np.asarray(case_model.arrays.tcollector_k[0], dtype=float).copy()
    initial_emitter_mean = float(np.mean(emitter))
    history: list[ClosedLoopIteration] = []
    last_feedback = ThermoElectricFeedback(
        electron_emitter_power_w=np.zeros(geometry.n_nodes),
        electron_collector_power_w=np.zeros(geometry.n_nodes),
        electron_emitter_flux_w_m2=np.zeros(geometry.n_nodes),
        electron_collector_flux_w_m2=np.zeros(geometry.n_nodes),
        joule_power_emitter_w=np.zeros(geometry.n_nodes),
        joule_power_collector_w=np.zeros(geometry.n_nodes),
    )
    global_results = {}
    max_change = math.inf
    closed_converged = False

    for iteration in range(1, closed_loop_config.max_iterations + 1):
        thermo_model.set_temperatures(emitter.reshape(1, geometry.n_nodes), collector.reshape(1, geometry.n_nodes))
        thermo_model.set_tcs(np.full((1, geometry.n_nodes), thermal_config.cesium_reservoir_temperature_k, dtype=float))
        thermo_model.calculate(verbose=False)
        global_results = thermo_model.get_global_results() or {}
        tec_result = thermo_model.get_tec_results(0) or {}
        last_feedback = extract_thermoelectric_feedback(tec_result, geometry)
        new_emitter, new_collector, max_change = update_temperature_fields(
            heat_source_w=case_model.heat_source_w,
            feedback=last_feedback,
            thermal_config=thermal_config,
            closed_loop_config=closed_loop_config,
            old_emitter_k=emitter,
            old_collector_k=collector,
        )
        emitter = new_emitter
        collector = new_collector
        current = float(global_results.get("Iout", math.nan))
        uout = float(global_results.get("Uout", point.voltage_v))
        thermo_converged = bool(global_results.get("converged", False))
        closed_converged = bool(max_change <= closed_loop_config.temperature_tolerance_k and thermo_converged)
        history.append(
            ClosedLoopIteration(
                iteration=iteration,
                current_a=current,
                uout_v=uout,
                max_temperature_change_k=max_change,
                emitter_temperature_mean_k=float(np.mean(emitter)),
                collector_temperature_mean_k=float(np.mean(collector)),
                electron_cooling_power_w=float(np.sum(last_feedback.electron_emitter_power_w)),
                collector_electron_heating_power_w=float(np.sum(last_feedback.electron_collector_power_w)),
                joule_power_emitter_w=float(np.sum(last_feedback.joule_power_emitter_w)),
                joule_power_collector_w=float(np.sum(last_feedback.joule_power_collector_w)),
                converged=closed_converged,
            )
        )
        if closed_converged:
            break

    current = float(global_results.get("Iout", math.nan))
    uout = float(global_results.get("Uout", point.voltage_v))
    power = current * uout
    finite = bool(np.isfinite([current, uout, power, max_change]).all())
    efficiency = 100.0 * power / (point.thermal_power_kw * 1000.0) if finite else math.nan
    return ClosedLoopResult(
        voltage_v=point.voltage_v,
        thermal_power_kw=point.thermal_power_kw,
        current_exp_a=point.current_a,
        current_calc_a=current,
        current_error_a=current - point.current_a if finite else math.nan,
        electric_power_exp_w=point.electric_power_w,
        electric_power_calc_w=power if finite else math.nan,
        electric_power_error_w=power - point.electric_power_w if finite else math.nan,
        efficiency_exp_percent=point.efficiency_percent,
        efficiency_calc_percent=efficiency,
        efficiency_error_percent=efficiency - point.efficiency_percent if finite else math.nan,
        emitter_temperature_exp_k=point.emitter_temperature_k,
        initial_emitter_temperature_mean_k=initial_emitter_mean,
        emitter_temperature_mean_k=float(np.mean(emitter)),
        emitter_temperature_min_k=float(np.min(emitter)),
        emitter_temperature_max_k=float(np.max(emitter)),
        collector_temperature_mean_k=float(np.mean(collector)),
        collector_temperature_min_k=float(np.min(collector)),
        collector_temperature_max_k=float(np.max(collector)),
        electron_cooling_power_w=float(np.sum(last_feedback.electron_emitter_power_w)),
        collector_electron_heating_power_w=float(np.sum(last_feedback.electron_collector_power_w)),
        joule_power_emitter_w=float(np.sum(last_feedback.joule_power_emitter_w)),
        joule_power_collector_w=float(np.sum(last_feedback.joule_power_collector_w)),
        outer_iterations=len(history),
        max_temperature_change_k=float(max_change),
        closed_loop_converged=closed_converged,
        thermocalc_converged=bool(global_results.get("converged", False)),
        finite=finite,
        iteration_history=history,
    )


def summarize_closed_loop_results(results: Sequence[ClosedLoopResult]) -> dict:
    current_errors = np.asarray([r.current_error_a for r in results], dtype=float)
    power_errors = np.asarray([r.electric_power_error_w for r in results], dtype=float)
    te_errors = np.asarray([r.emitter_temperature_mean_k - r.emitter_temperature_exp_k for r in results], dtype=float)
    return {
        "case_count": len(results),
        "finite_all": all(r.finite for r in results),
        "thermocalc_converged_all": all(r.thermocalc_converged for r in results),
        "closed_loop_converged_all": all(r.closed_loop_converged for r in results),
        "current_mae_a": float(np.nanmean(np.abs(current_errors))),
        "current_max_abs_a": float(np.nanmax(np.abs(current_errors))),
        "electric_power_mae_w": float(np.nanmean(np.abs(power_errors))),
        "emitter_temperature_mae_k": float(np.nanmean(np.abs(te_errors))),
        "max_outer_iterations": int(max((r.outer_iterations for r in results), default=0)),
    }


def result_to_dict(result: ClosedLoopResult) -> dict:
    data = asdict(result)
    data["iteration_history"] = [asdict(item) for item in result.iteration_history]
    return data
