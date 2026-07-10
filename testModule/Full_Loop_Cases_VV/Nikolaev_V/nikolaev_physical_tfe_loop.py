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
from nikolaev_thermocalc_model import NikolaevThermalNetworkConfig, NikolaevTfeGeometry, centered_heater_power_profile, build_thermocalc_case
from nikolaev_thermocalc_runner import apply_case_to_thermocalc
from nikolaev_thermoelectric_closed_loop import ThermoElectricFeedback, extract_thermoelectric_feedback
from Materials.Solids.MoNb import MoNb
from Materials.Solids.Molybdenum import Molybdenum


@dataclass(frozen=True)
class PhysicalLoopConfig:
    n_nodes: int = 50
    coolant_inlet_temperature_k: float = 770.0
    coolant_mass_flow_kg_s: float = 0.040
    coolant_heat_capacity_j_kg_k: float = 1000.0
    collector_convective_h_w_m2_k: float = 8500.0
    emitter_to_collector_resistance_k_per_w: float = 0.340
    cesium_reservoir_temperature_k: float = 560.0
    wire_resistance_ohm: float = 0.0005
    axial_shape_amplitude: float = 0.18
    axial_conduction_enabled: bool = False
    axial_conduction_smoothing: float = 0.0
    axial_conduction_passes: int = 1
    max_iterations: int = 40
    relaxation: float = 0.25
    temperature_tolerance_k: float = 0.75
    i_guess_a: float = 300.0
    min_temperature_k: float = 500.0
    max_temperature_k: float = 2600.0


@dataclass(frozen=True)
class ThermalHydraulicState:
    emitter_temperature_k: np.ndarray
    collector_temperature_k: np.ndarray
    coolant_temperature_faces_k: np.ndarray
    coolant_bulk_temperature_k: np.ndarray
    heat_to_coolant_w: np.ndarray
    emitter_to_collector_heat_w: np.ndarray
    coolant_inlet_temperature_k: float
    coolant_outlet_temperature_k: float
    coolant_heat_gain_w: float
    thermal_balance_residual_w: float


@dataclass(frozen=True)
class PhysicalLoopIteration:
    iteration: int
    current_a: float
    uout_v: float
    max_temperature_change_k: float
    emitter_temperature_mean_k: float
    collector_temperature_mean_k: float
    coolant_outlet_temperature_k: float
    coolant_heat_gain_w: float
    electron_cooling_power_w: float
    collector_electron_heating_power_w: float
    joule_power_emitter_w: float
    joule_power_collector_w: float
    converged: bool


@dataclass(frozen=True)
class PhysicalLoopResult:
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
    emitter_temperature_mean_k: float
    emitter_temperature_min_k: float
    emitter_temperature_max_k: float
    collector_temperature_mean_k: float
    collector_temperature_min_k: float
    collector_temperature_max_k: float
    coolant_inlet_temperature_k: float
    coolant_outlet_temperature_k: float
    coolant_heat_gain_w: float
    heat_to_coolant_w: float
    thermal_balance_residual_w: float
    electron_cooling_power_w: float
    collector_electron_heating_power_w: float
    joule_power_emitter_w: float
    joule_power_collector_w: float
    outer_iterations: int
    max_temperature_change_k: float
    physical_loop_converged: bool
    thermocalc_converged: bool
    finite: bool
    iteration_history: list[PhysicalLoopIteration] = field(default_factory=list)


def _default_model_factory(n_elements: int, n_nodes: int):
    from ThermoCalc.ThermoCalcWrapper import ThermoCalcModel

    return ThermoCalcModel(n_elements=n_elements, n_nodes=n_nodes)


def _validate_config(config: PhysicalLoopConfig) -> None:
    if config.n_nodes <= 0:
        raise ValueError("n_nodes must be positive")
    if config.coolant_mass_flow_kg_s <= 0.0:
        raise ValueError("coolant_mass_flow_kg_s must be positive")
    if config.coolant_heat_capacity_j_kg_k <= 0.0:
        raise ValueError("coolant_heat_capacity_j_kg_k must be positive")
    if config.collector_convective_h_w_m2_k <= 0.0:
        raise ValueError("collector_convective_h_w_m2_k must be positive")
    if not 0.0 < config.relaxation <= 1.0:
        raise ValueError("relaxation must be in (0, 1]")
    if not 0.0 <= config.axial_conduction_smoothing <= 0.5:
        raise ValueError("axial_conduction_smoothing must be in [0, 0.5]")
    if config.axial_conduction_passes <= 0:
        raise ValueError("axial_conduction_passes must be positive")

def zero_feedback(n_nodes: int) -> ThermoElectricFeedback:
    z = np.zeros(n_nodes, dtype=float)
    return ThermoElectricFeedback(
        electron_emitter_power_w=z.copy(),
        electron_collector_power_w=z.copy(),
        electron_emitter_flux_w_m2=z.copy(),
        electron_collector_flux_w_m2=z.copy(),
        joule_power_emitter_w=z.copy(),
        joule_power_collector_w=z.copy(),
    )



def apply_axial_conduction_smoothing(temperature: np.ndarray, strength: float, passes: int = 1) -> np.ndarray:
    """Apply explicit 1D axial diffusion with zero end heat flux while preserving the mean.

    Kept only for compatibility with historical runs. New validation runs should use
    ``axial_conduction_enabled=True`` so axial redistribution is solved from material
    conductivity and geometry instead of this dimensionless smoothing operator.
    """
    values = np.asarray(temperature, dtype=float)
    if values.ndim != 1:
        raise ValueError(f"temperature must be 1D, got {values.shape}")
    if not 0.0 <= strength <= 0.5:
        raise ValueError("strength must be in [0, 0.5]")
    if passes <= 0:
        raise ValueError("passes must be positive")
    result = values.copy()
    if result.size < 2 or strength == 0.0:
        return result
    original_mean = float(np.mean(result))
    for _ in range(int(passes)):
        old = result.copy()
        result[0] = old[0] + strength * (old[1] - old[0])
        result[-1] = old[-1] + strength * (old[-2] - old[-1])
        if result.size > 2:
            result[1:-1] = old[1:-1] + strength * (old[:-2] - 2.0 * old[1:-1] + old[2:])
    result += original_mean - float(np.mean(result))
    return result


def axial_face_conductance_w_per_k(temperature: np.ndarray, cross_area_m2: float, dz_m: float, material) -> np.ndarray:
    """Return finite-volume axial face conductance ``k(T_face) * A / dz``."""
    values = np.asarray(temperature, dtype=float)
    if values.ndim != 1:
        raise ValueError(f"temperature must be 1D, got {values.shape}")
    if values.size < 2:
        return np.zeros(0, dtype=float)
    if cross_area_m2 <= 0.0:
        raise ValueError("cross_area_m2 must be positive")
    if dz_m <= 0.0:
        raise ValueError("dz_m must be positive")
    face_temperature = 0.5 * (values[:-1] + values[1:])
    conductivity = np.asarray(material.conductivity(face_temperature), dtype=float)
    return conductivity * float(cross_area_m2) / float(dz_m)


def _add_axial_conductance(matrix: np.ndarray, row: int, var: int, conductance_faces: np.ndarray) -> None:
    n = conductance_faces.size + 1
    if var > 0:
        g_left = float(conductance_faces[var - 1])
        matrix[row, row] += g_left
        matrix[row, row - 1] -= g_left
    if var < n - 1:
        g_right = float(conductance_faces[var])
        matrix[row, row] += g_right
        matrix[row, row + 1] -= g_right


def _solve_solid_temperatures_with_axial_conduction(
    emitter_source_w: np.ndarray,
    collector_source_w: np.ndarray,
    coolant_bulk_k: np.ndarray,
    config: PhysicalLoopConfig,
    geometry: NikolaevTfeGeometry,
    initial_emitter_k: np.ndarray,
    initial_collector_k: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n = config.n_nodes
    dz = geometry.active_length_m / n
    radial_resistance = config.emitter_to_collector_resistance_k_per_w * n
    radial_g = 1.0 / radial_resistance
    collector_area_node = geometry.collector_side_area_total_m2 / n
    conv_g = config.collector_convective_h_w_m2_k * collector_area_node
    emitter_material = MoNb()
    collector_material = Molybdenum()
    emitter = np.asarray(initial_emitter_k, dtype=float).copy()
    collector = np.asarray(initial_collector_k, dtype=float).copy()

    for _ in range(30):
        emitter_axial_g = axial_face_conductance_w_per_k(emitter, geometry.emitter_cross_area_m2, dz, emitter_material)
        collector_axial_g = axial_face_conductance_w_per_k(collector, geometry.collector_cross_area_m2, dz, collector_material)
        matrix = np.zeros((2 * n, 2 * n), dtype=float)
        rhs = np.zeros(2 * n, dtype=float)

        for idx in range(n):
            row_e = idx
            matrix[row_e, idx] += radial_g
            matrix[row_e, n + idx] -= radial_g
            _add_axial_conductance(matrix, row_e, idx, emitter_axial_g)
            rhs[row_e] = emitter_source_w[idx]

            row_c = n + idx
            matrix[row_c, idx] -= radial_g
            matrix[row_c, n + idx] += radial_g + conv_g
            if idx > 0:
                g_left = float(collector_axial_g[idx - 1])
                matrix[row_c, n + idx] += g_left
                matrix[row_c, n + idx - 1] -= g_left
            if idx < n - 1:
                g_right = float(collector_axial_g[idx])
                matrix[row_c, n + idx] += g_right
                matrix[row_c, n + idx + 1] -= g_right
            rhs[row_c] = collector_source_w[idx] + conv_g * coolant_bulk_k[idx]

        solved = np.linalg.solve(matrix, rhs)
        new_emitter = solved[:n]
        new_collector = solved[n:]
        max_change = max(float(np.max(np.abs(new_emitter - emitter))), float(np.max(np.abs(new_collector - collector))))
        emitter = new_emitter
        collector = new_collector
        if max_change < 1.0e-8:
            break
    return emitter, collector


def _solve_thermal_hydraulic_update_with_axial_conduction(
    emitter_source_w: np.ndarray,
    collector_source_w: np.ndarray,
    config: PhysicalLoopConfig,
    geometry: NikolaevTfeGeometry,
) -> ThermalHydraulicState:
    n = config.n_nodes
    mcp = config.coolant_mass_flow_kg_s * config.coolant_heat_capacity_j_kg_k
    collector_area_node = geometry.collector_side_area_total_m2 / n
    conv_g = config.collector_convective_h_w_m2_k * collector_area_node
    total_source = emitter_source_w + collector_source_w

    coolant_faces = np.empty(n + 1, dtype=float)
    coolant_faces[0] = config.coolant_inlet_temperature_k
    for idx, q_node in enumerate(total_source):
        coolant_faces[idx + 1] = coolant_faces[idx] + q_node / mcp
    coolant_bulk = 0.5 * (coolant_faces[:-1] + coolant_faces[1:])

    local_state = solve_thermal_hydraulic_update(
        emitter_source_w,
        ThermoElectricFeedback(
            electron_emitter_power_w=np.zeros(n, dtype=float),
            electron_collector_power_w=collector_source_w.copy(),
            electron_emitter_flux_w_m2=np.zeros(n, dtype=float),
            electron_collector_flux_w_m2=np.zeros(n, dtype=float),
            joule_power_emitter_w=np.zeros(n, dtype=float),
            joule_power_collector_w=np.zeros(n, dtype=float),
        ),
        PhysicalLoopConfig(
            n_nodes=n,
            coolant_inlet_temperature_k=config.coolant_inlet_temperature_k,
            coolant_mass_flow_kg_s=config.coolant_mass_flow_kg_s,
            coolant_heat_capacity_j_kg_k=config.coolant_heat_capacity_j_kg_k,
            collector_convective_h_w_m2_k=config.collector_convective_h_w_m2_k,
            emitter_to_collector_resistance_k_per_w=config.emitter_to_collector_resistance_k_per_w,
            cesium_reservoir_temperature_k=config.cesium_reservoir_temperature_k,
            wire_resistance_ohm=config.wire_resistance_ohm,
            axial_shape_amplitude=config.axial_shape_amplitude,
            axial_conduction_enabled=False,
            max_iterations=config.max_iterations,
            relaxation=config.relaxation,
            temperature_tolerance_k=config.temperature_tolerance_k,
            i_guess_a=config.i_guess_a,
            min_temperature_k=config.min_temperature_k,
            max_temperature_k=config.max_temperature_k,
        ),
    )
    emitter = local_state.emitter_temperature_k.copy()
    collector = local_state.collector_temperature_k.copy()

    for _ in range(80):
        new_emitter, new_collector = _solve_solid_temperatures_with_axial_conduction(
            emitter_source_w,
            collector_source_w,
            coolant_bulk,
            config,
            geometry,
            emitter,
            collector,
        )
        heat_to_coolant = conv_g * (new_collector - coolant_bulk)
        new_faces = np.empty(n + 1, dtype=float)
        new_faces[0] = config.coolant_inlet_temperature_k
        for idx, q_node in enumerate(heat_to_coolant):
            new_faces[idx + 1] = new_faces[idx] + q_node / mcp
        new_bulk = 0.5 * (new_faces[:-1] + new_faces[1:])
        max_change = max(
            float(np.max(np.abs(new_emitter - emitter))),
            float(np.max(np.abs(new_collector - collector))),
            float(np.max(np.abs(new_faces - coolant_faces))),
        )
        emitter = new_emitter
        collector = new_collector
        coolant_faces = 0.5 * coolant_faces + 0.5 * new_faces
        coolant_bulk = 0.5 * (coolant_faces[:-1] + coolant_faces[1:])
        if max_change < 1.0e-7:
            coolant_faces = new_faces
            coolant_bulk = new_bulk
            break

    heat_to_coolant = conv_g * (collector - coolant_bulk)
    coolant_faces[0] = config.coolant_inlet_temperature_k
    for idx, q_node in enumerate(heat_to_coolant):
        coolant_faces[idx + 1] = coolant_faces[idx] + q_node / mcp
    coolant_bulk = 0.5 * (coolant_faces[:-1] + coolant_faces[1:])
    coolant_heat_gain = mcp * (coolant_faces[-1] - coolant_faces[0])
    radial_g = 1.0 / (config.emitter_to_collector_resistance_k_per_w * n)
    emitter_to_collector_heat = radial_g * (emitter - collector)
    residual = float(np.sum(total_source) - coolant_heat_gain)
    return ThermalHydraulicState(
        emitter_temperature_k=emitter,
        collector_temperature_k=collector,
        coolant_temperature_faces_k=coolant_faces,
        coolant_bulk_temperature_k=coolant_bulk,
        heat_to_coolant_w=heat_to_coolant,
        emitter_to_collector_heat_w=emitter_to_collector_heat,
        coolant_inlet_temperature_k=float(coolant_faces[0]),
        coolant_outlet_temperature_k=float(coolant_faces[-1]),
        coolant_heat_gain_w=float(coolant_heat_gain),
        thermal_balance_residual_w=residual,
    )


def solve_thermal_hydraulic_update(
    heat_source_w: np.ndarray,
    feedback: ThermoElectricFeedback,
    config: PhysicalLoopConfig = PhysicalLoopConfig(),
) -> ThermalHydraulicState:
    _validate_config(config)
    geometry = NikolaevTfeGeometry(n_nodes=config.n_nodes)
    heat_source = np.asarray(heat_source_w, dtype=float)
    if heat_source.shape != (config.n_nodes,):
        raise ValueError(f"heat_source_w must have shape {(config.n_nodes,)}, got {heat_source.shape}")

    emitter_source = heat_source + feedback.electron_emitter_power_w + feedback.joule_power_emitter_w
    collector_source = feedback.electron_collector_power_w + feedback.joule_power_collector_w
    if config.axial_conduction_enabled:
        return _solve_thermal_hydraulic_update_with_axial_conduction(emitter_source, collector_source, config, geometry)

    emitter_to_collector_heat = emitter_source
    heat_to_coolant = emitter_to_collector_heat + collector_source

    mcp = config.coolant_mass_flow_kg_s * config.coolant_heat_capacity_j_kg_k
    coolant_faces = np.empty(config.n_nodes + 1, dtype=float)
    coolant_faces[0] = config.coolant_inlet_temperature_k
    for idx, q_node in enumerate(heat_to_coolant):
        coolant_faces[idx + 1] = coolant_faces[idx] + q_node / mcp
    coolant_bulk = 0.5 * (coolant_faces[:-1] + coolant_faces[1:])

    collector_area_node = geometry.collector_side_area_total_m2 / config.n_nodes
    r_collector_to_coolant_node = 1.0 / (config.collector_convective_h_w_m2_k * collector_area_node)
    r_emitter_to_collector_node = config.emitter_to_collector_resistance_k_per_w * config.n_nodes

    collector = coolant_bulk + heat_to_coolant * r_collector_to_coolant_node
    emitter = collector + emitter_to_collector_heat * r_emitter_to_collector_node
    if config.axial_conduction_smoothing > 0.0:
        collector = apply_axial_conduction_smoothing(collector, config.axial_conduction_smoothing, config.axial_conduction_passes)
        emitter = apply_axial_conduction_smoothing(emitter, config.axial_conduction_smoothing, config.axial_conduction_passes)
    collector = np.clip(collector, config.min_temperature_k, config.max_temperature_k)
    emitter = np.clip(emitter, config.min_temperature_k, config.max_temperature_k)
    coolant_heat_gain = mcp * (coolant_faces[-1] - coolant_faces[0])
    residual = float(np.sum(heat_to_coolant) - coolant_heat_gain)
    return ThermalHydraulicState(
        emitter_temperature_k=emitter,
        collector_temperature_k=collector,
        coolant_temperature_faces_k=coolant_faces,
        coolant_bulk_temperature_k=coolant_bulk,
        heat_to_coolant_w=heat_to_coolant,
        emitter_to_collector_heat_w=emitter_to_collector_heat,
        coolant_inlet_temperature_k=float(coolant_faces[0]),
        coolant_outlet_temperature_k=float(coolant_faces[-1]),
        coolant_heat_gain_w=float(coolant_heat_gain),
        thermal_balance_residual_w=residual,
    )

def _make_thermal_config(config: PhysicalLoopConfig) -> NikolaevThermalNetworkConfig:
    return NikolaevThermalNetworkConfig(
        n_nodes=config.n_nodes,
        collector_boundary_temperature_k=config.coolant_inlet_temperature_k,
        emitter_to_collector_resistance_k_per_w=config.emitter_to_collector_resistance_k_per_w,
        axial_shape_amplitude=config.axial_shape_amplitude,
        cesium_reservoir_temperature_k=config.cesium_reservoir_temperature_k,
        wire_resistance_ohm=config.wire_resistance_ohm,
    )


def solve_physical_tfe_point(
    point: Table2OperatingPoint,
    config: PhysicalLoopConfig = PhysicalLoopConfig(),
    model_factory: Callable[[int, int], object] = _default_model_factory,
) -> PhysicalLoopResult:
    _validate_config(config)
    thermal_config = _make_thermal_config(config)
    case_model = build_thermocalc_case(point, thermal_config)
    geometry = case_model.geometry
    heat_source = centered_heater_power_profile(point.thermal_power_kw * 1000.0, geometry)

    initial_state = solve_thermal_hydraulic_update(heat_source, zero_feedback(config.n_nodes), config)
    emitter = initial_state.emitter_temperature_k.copy()
    collector = initial_state.collector_temperature_k.copy()
    state = initial_state
    feedback = zero_feedback(config.n_nodes)
    thermo_model = model_factory(1, geometry.n_nodes)
    apply_case_to_thermocalc(thermo_model, case_model, point.voltage_v, config.i_guess_a)

    history: list[PhysicalLoopIteration] = []
    global_results = {}
    max_change = math.inf
    physical_converged = False

    for iteration in range(1, config.max_iterations + 1):
        thermo_model.set_temperatures(emitter.reshape(1, config.n_nodes), collector.reshape(1, config.n_nodes))
        thermo_model.set_tcs(np.full((1, config.n_nodes), config.cesium_reservoir_temperature_k, dtype=float))
        thermo_model.calculate(verbose=False)
        global_results = thermo_model.get_global_results() or {}
        tec_result = thermo_model.get_tec_results(0) or {}
        feedback = extract_thermoelectric_feedback(tec_result, geometry)
        target_state = solve_thermal_hydraulic_update(heat_source, feedback, config)
        new_emitter = emitter + config.relaxation * (target_state.emitter_temperature_k - emitter)
        new_collector = collector + config.relaxation * (target_state.collector_temperature_k - collector)
        max_change = float(max(np.max(np.abs(new_emitter - emitter)), np.max(np.abs(new_collector - collector))))
        emitter = new_emitter
        collector = new_collector
        state = ThermalHydraulicState(
            emitter_temperature_k=emitter.copy(),
            collector_temperature_k=collector.copy(),
            coolant_temperature_faces_k=target_state.coolant_temperature_faces_k.copy(),
            coolant_bulk_temperature_k=target_state.coolant_bulk_temperature_k.copy(),
            heat_to_coolant_w=target_state.heat_to_coolant_w.copy(),
            emitter_to_collector_heat_w=target_state.emitter_to_collector_heat_w.copy(),
            coolant_inlet_temperature_k=target_state.coolant_inlet_temperature_k,
            coolant_outlet_temperature_k=target_state.coolant_outlet_temperature_k,
            coolant_heat_gain_w=target_state.coolant_heat_gain_w,
            thermal_balance_residual_w=target_state.thermal_balance_residual_w,
        )
        current = float(global_results.get("Iout", math.nan))
        uout = float(global_results.get("Uout", point.voltage_v))
        thermo_converged = bool(global_results.get("converged", False))
        physical_converged = bool(max_change <= config.temperature_tolerance_k and thermo_converged)
        history.append(
            PhysicalLoopIteration(
                iteration=iteration,
                current_a=current,
                uout_v=uout,
                max_temperature_change_k=max_change,
                emitter_temperature_mean_k=float(np.mean(emitter)),
                collector_temperature_mean_k=float(np.mean(collector)),
                coolant_outlet_temperature_k=state.coolant_outlet_temperature_k,
                coolant_heat_gain_w=state.coolant_heat_gain_w,
                electron_cooling_power_w=float(np.sum(feedback.electron_emitter_power_w)),
                collector_electron_heating_power_w=float(np.sum(feedback.electron_collector_power_w)),
                joule_power_emitter_w=float(np.sum(feedback.joule_power_emitter_w)),
                joule_power_collector_w=float(np.sum(feedback.joule_power_collector_w)),
                converged=physical_converged,
            )
        )
        if physical_converged:
            break

    current = float(global_results.get("Iout", math.nan))
    uout = float(global_results.get("Uout", point.voltage_v))
    electric_power = current * uout
    finite = bool(np.isfinite([current, uout, electric_power, max_change, state.coolant_outlet_temperature_k]).all())
    efficiency = 100.0 * electric_power / (point.thermal_power_kw * 1000.0) if finite else math.nan
    return PhysicalLoopResult(
        voltage_v=point.voltage_v,
        thermal_power_kw=point.thermal_power_kw,
        current_exp_a=point.current_a,
        current_calc_a=current,
        current_error_a=current - point.current_a if finite else math.nan,
        electric_power_exp_w=point.electric_power_w,
        electric_power_calc_w=electric_power if finite else math.nan,
        electric_power_error_w=electric_power - point.electric_power_w if finite else math.nan,
        efficiency_exp_percent=point.efficiency_percent,
        efficiency_calc_percent=efficiency,
        efficiency_error_percent=efficiency - point.efficiency_percent if finite else math.nan,
        emitter_temperature_exp_k=point.emitter_temperature_k,
        emitter_temperature_mean_k=float(np.mean(emitter)),
        emitter_temperature_min_k=float(np.min(emitter)),
        emitter_temperature_max_k=float(np.max(emitter)),
        collector_temperature_mean_k=float(np.mean(collector)),
        collector_temperature_min_k=float(np.min(collector)),
        collector_temperature_max_k=float(np.max(collector)),
        coolant_inlet_temperature_k=state.coolant_inlet_temperature_k,
        coolant_outlet_temperature_k=state.coolant_outlet_temperature_k,
        coolant_heat_gain_w=state.coolant_heat_gain_w,
        heat_to_coolant_w=float(np.sum(state.heat_to_coolant_w)),
        thermal_balance_residual_w=state.thermal_balance_residual_w,
        electron_cooling_power_w=float(np.sum(feedback.electron_emitter_power_w)),
        collector_electron_heating_power_w=float(np.sum(feedback.electron_collector_power_w)),
        joule_power_emitter_w=float(np.sum(feedback.joule_power_emitter_w)),
        joule_power_collector_w=float(np.sum(feedback.joule_power_collector_w)),
        outer_iterations=len(history),
        max_temperature_change_k=float(max_change),
        physical_loop_converged=physical_converged,
        thermocalc_converged=bool(global_results.get("converged", False)),
        finite=finite,
        iteration_history=history,
    )


def summarize_physical_results(results: Sequence[PhysicalLoopResult]) -> dict:
    current_errors = np.asarray([r.current_error_a for r in results], dtype=float)
    power_errors = np.asarray([r.electric_power_error_w for r in results], dtype=float)
    te_errors = np.asarray([r.emitter_temperature_mean_k - r.emitter_temperature_exp_k for r in results], dtype=float)
    residuals = np.asarray([r.thermal_balance_residual_w for r in results], dtype=float)
    return {
        "case_count": len(results),
        "finite_all": all(r.finite for r in results),
        "thermocalc_converged_all": all(r.thermocalc_converged for r in results),
        "physical_loop_converged_all": all(r.physical_loop_converged for r in results),
        "current_mae_a": float(np.nanmean(np.abs(current_errors))),
        "current_max_abs_a": float(np.nanmax(np.abs(current_errors))),
        "electric_power_mae_w": float(np.nanmean(np.abs(power_errors))),
        "emitter_temperature_mae_k": float(np.nanmean(np.abs(te_errors))),
        "thermal_balance_max_abs_w": float(np.nanmax(np.abs(residuals))),
        "max_outer_iterations": int(max((r.outer_iterations for r in results), default=0)),
    }


def result_to_dict(result: PhysicalLoopResult) -> dict:
    data = asdict(result)
    data["iteration_history"] = [asdict(item) for item in result.iteration_history]
    return data
