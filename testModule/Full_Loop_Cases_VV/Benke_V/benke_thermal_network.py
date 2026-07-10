from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

BENKE_ACTIVE_ZONE_FRACTION = 0.88
BENKE_THERMOCOUPLE_POSITIONS_MM: tuple[float | None, ...] = (
    -205.0,
    -163.0,
    -108.0,
    -55.0,
    -55.0,
    0.0,
    0.0,
    55.0,
    None,
    108.0,
    163.0,
    205.0,
)
BENKE_AVERAGE_THERMOCOUPLE_INDICES: tuple[int, ...] = (2, 3, 4, 5, 6, 7, 8, 10, 11)
BENKE_THERMOCOUPLE_DEPTH_FROM_SLEEVE_INNER_M = 0.0018


@dataclass(frozen=True)
class BenkeThermalCase:
    name: str
    tisa_power_w: float
    regulated_he_pressure_torr: float = 10.0


@dataclass(frozen=True)
class BenkeThermalNetworkConfig:
    active_length_m: float = 0.410
    tisa_heated_length_m: float = 0.300
    n_nodes: int = 60
    water_inlet_temperature_k: float = 310.0
    water_mass_flow_kg_s: float = 0.03
    water_cp_j_kg_k: float = 4180.0
    water_h_w_m2_k: float = 800.0
    regulated_he_effective_k_w_m_k: float = 0.08
    regulated_he_gap_m: float = 0.0005
    unregulated_he_effective_k_w_m_k: float = 0.276
    unregulated_he_gap_m: float = 0.00005
    cs_gap_effective_k_w_m_k: float = 0.12
    collector_k_w_m_k: float = 115.0
    alumina_k_w_m_k: float = 25.0
    sleeve_k_w_m_k: float = 16.0
    collector_inner_radius_m: float = 0.0206 / 2.0
    collector_wall_thickness_m: float = 0.0014
    emitter_outer_radius_m: float = 0.0196 / 2.0
    alumina_thickness_m: float = 0.00015
    sleeve_outer_radius_m: float = 0.0299 / 2.0
    extra_resistance_k_per_w: float = 0.0
    coolant_heat_fraction: float = 1.0


@dataclass(frozen=True)
class BenkeThermalNetworkResult:
    case: BenkeThermalCase
    config: BenkeThermalNetworkConfig
    active_zone_power_w: float
    heat_source_w: np.ndarray
    water_bulk_temperature_k: np.ndarray
    water_bulk_outlet_k: float
    collector_inner_temperature_k: np.ndarray
    sleeve_outer_temperature_k: np.ndarray
    sleeve_thermocouple_temperature_k: np.ndarray
    sleeve_thermocouple_z_mm: tuple[float | None, ...]
    sleeve_thermocouple_radius_m: float
    sleeve_thermocouple_included_in_benke_average: tuple[bool, ...]
    r_collector_to_water_k_per_w: np.ndarray
    r_sleeve_to_water_k_per_w: np.ndarray
    energy_balance_error_w: float


BENKE_TYPICAL_CASE = BenkeThermalCase(
    name="benke_typical_insulation_case_3003w_10torr",
    tisa_power_w=3003.0 / BENKE_ACTIVE_ZONE_FRACTION,
    regulated_he_pressure_torr=10.0,
)


def active_zone_power_w(tisa_power_w: float) -> float:
    if tisa_power_w < 0.0:
        raise ValueError("TISA power cannot be negative.")
    return float(tisa_power_w) * BENKE_ACTIVE_ZONE_FRACTION


def cylindrical_resistance_k_per_w(r_inner_m: float, r_outer_m: float, k_w_m_k: float, length_m: float) -> float:
    if r_inner_m <= 0.0 or r_outer_m <= 0.0:
        raise ValueError("Radii must be positive.")
    if r_outer_m < r_inner_m:
        raise ValueError("Outer radius must be no smaller than inner radius.")
    if k_w_m_k <= 0.0:
        raise ValueError("Thermal conductivity must be positive.")
    if length_m <= 0.0:
        raise ValueError("Length must be positive.")
    if r_outer_m == r_inner_m:
        return 0.0
    return math.log(r_outer_m / r_inner_m) / (2.0 * math.pi * k_w_m_k * length_m)


def convection_resistance_k_per_w(radius_m: float, h_w_m2_k: float, length_m: float) -> float:
    if radius_m <= 0.0:
        raise ValueError("Radius must be positive.")
    if h_w_m2_k <= 0.0:
        raise ValueError("Heat-transfer coefficient must be positive.")
    if length_m <= 0.0:
        raise ValueError("Length must be positive.")
    return 1.0 / (h_w_m2_k * 2.0 * math.pi * radius_m * length_m)


def build_centered_tisa_heat_source_w(
    active_zone_power_w: float,
    n_nodes: int,
    active_length_m: float,
    heated_length_m: float,
) -> np.ndarray:
    if active_zone_power_w < 0.0:
        raise ValueError("Active-zone power cannot be negative.")
    if n_nodes <= 0:
        raise ValueError("n_nodes must be positive.")
    if active_length_m <= 0.0:
        raise ValueError("active_length_m must be positive.")
    if heated_length_m <= 0.0 or heated_length_m > active_length_m:
        raise ValueError("heated_length_m must be positive and no longer than active_length_m.")

    centers_m = (np.arange(n_nodes, dtype=float) + 0.5) * active_length_m / n_nodes
    margin_m = 0.5 * (active_length_m - heated_length_m)
    shape = np.zeros(n_nodes, dtype=float)
    heated = (centers_m >= margin_m) & (centers_m <= margin_m + heated_length_m)
    shape[heated] = 1.0
    if margin_m > 0.0:
        left = centers_m < margin_m
        right = centers_m > margin_m + heated_length_m
        shape[left] = centers_m[left] / margin_m
        shape[right] = (active_length_m - centers_m[right]) / margin_m
    shape = np.clip(shape, 0.0, None)
    total_shape = float(shape.sum())
    if total_shape <= 0.0:
        raise ValueError("TISA heat source shape has zero integral.")
    return float(active_zone_power_w) * shape / total_shape


def _validate_config(config: BenkeThermalNetworkConfig) -> None:
    positive = (
        "active_length_m",
        "tisa_heated_length_m",
        "water_mass_flow_kg_s",
        "water_cp_j_kg_k",
        "water_h_w_m2_k",
        "regulated_he_effective_k_w_m_k",
        "unregulated_he_effective_k_w_m_k",
        "cs_gap_effective_k_w_m_k",
        "collector_k_w_m_k",
        "alumina_k_w_m_k",
        "sleeve_k_w_m_k",
    )
    for name in positive:
        if getattr(config, name) <= 0.0:
            raise ValueError(f"{name} must be positive.")
    if config.n_nodes <= 0:
        raise ValueError("n_nodes must be positive.")
    if config.tisa_heated_length_m > config.active_length_m:
        raise ValueError("TISA heated length cannot exceed active length.")
    if config.regulated_he_gap_m < 0.0 or config.unregulated_he_gap_m < 0.0:
        raise ValueError("Helium gap widths cannot be negative.")
    if config.extra_resistance_k_per_w < 0.0:
        raise ValueError("extra_resistance_k_per_w cannot be negative.")
    if not (0.0 < config.coolant_heat_fraction <= 1.0):
        raise ValueError("coolant_heat_fraction must be in (0, 1].")


def _centered_axial_coordinates_m(n_nodes: int, active_length_m: float) -> np.ndarray:
    return (np.arange(n_nodes, dtype=float) + 0.5) * active_length_m / n_nodes - 0.5 * active_length_m


def _sample_benke_thermocouple_values(values: np.ndarray, active_length_m: float) -> np.ndarray:
    z_m = _centered_axial_coordinates_m(values.size, active_length_m)
    sampled = []
    for position_mm in BENKE_THERMOCOUPLE_POSITIONS_MM:
        if position_mm is None:
            sampled.append(float("nan"))
        else:
            sampled.append(float(np.interp(position_mm / 1000.0, z_m, values, left=values[0], right=values[-1])))
    return np.asarray(sampled, dtype=float)


def benke_average_inclusion_flags() -> tuple[bool, ...]:
    included = set(BENKE_AVERAGE_THERMOCOUPLE_INDICES)
    return tuple((idx + 1) in included for idx in range(len(BENKE_THERMOCOUPLE_POSITIONS_MM)))


def solve_benke_thermal_network(
    case: BenkeThermalCase = BENKE_TYPICAL_CASE,
    config: BenkeThermalNetworkConfig | None = None,
) -> BenkeThermalNetworkResult:
    config = BenkeThermalNetworkConfig() if config is None else config
    _validate_config(config)
    p_az_w = active_zone_power_w(case.tisa_power_w)
    coolant_heat_w = p_az_w * config.coolant_heat_fraction
    heat_source_w = build_centered_tisa_heat_source_w(
        coolant_heat_w,
        n_nodes=config.n_nodes,
        active_length_m=config.active_length_m,
        heated_length_m=config.tisa_heated_length_m,
    )
    dz_m = config.active_length_m / config.n_nodes

    collector_inner_r = config.collector_inner_radius_m
    collector_outer_r = collector_inner_r + config.collector_wall_thickness_m
    alumina_outer_r = collector_outer_r + config.alumina_thickness_m
    unregulated_outer_r = alumina_outer_r + config.unregulated_he_gap_m
    sleeve_outer_r = max(config.sleeve_outer_radius_m, unregulated_outer_r)
    regulated_outer_r = sleeve_outer_r + config.regulated_he_gap_m
    sleeve_thermocouple_r = unregulated_outer_r + BENKE_THERMOCOUPLE_DEPTH_FROM_SLEEVE_INNER_M
    if sleeve_thermocouple_r >= sleeve_outer_r:
        raise ValueError("Benke sleeve thermocouple radius must be inside the collector sleeve.")

    r_collector = cylindrical_resistance_k_per_w(
        collector_inner_r,
        collector_outer_r,
        config.collector_k_w_m_k,
        dz_m,
    )
    r_alumina = cylindrical_resistance_k_per_w(
        collector_outer_r,
        alumina_outer_r,
        config.alumina_k_w_m_k,
        dz_m,
    )
    r_he_unreg = cylindrical_resistance_k_per_w(
        alumina_outer_r,
        unregulated_outer_r,
        config.unregulated_he_effective_k_w_m_k,
        dz_m,
    )
    r_sleeve = cylindrical_resistance_k_per_w(
        unregulated_outer_r,
        sleeve_outer_r,
        config.sleeve_k_w_m_k,
        dz_m,
    )
    r_sleeve_tc_to_outer = cylindrical_resistance_k_per_w(
        sleeve_thermocouple_r,
        sleeve_outer_r,
        config.sleeve_k_w_m_k,
        dz_m,
    )
    r_he_reg = cylindrical_resistance_k_per_w(
        sleeve_outer_r,
        regulated_outer_r,
        config.regulated_he_effective_k_w_m_k,
        dz_m,
    )
    r_water = convection_resistance_k_per_w(regulated_outer_r, config.water_h_w_m2_k, dz_m)
    r_extra_node = config.extra_resistance_k_per_w * config.n_nodes

    r_sleeve_to_water = r_he_reg + r_water + r_extra_node
    r_sleeve_tc_to_water = r_sleeve_tc_to_outer + r_sleeve_to_water
    r_collector_to_water = r_collector + r_alumina + r_he_unreg + r_sleeve + r_sleeve_to_water

    water_bulk = np.empty(config.n_nodes, dtype=float)
    sleeve_outer = np.empty(config.n_nodes, dtype=float)
    sleeve_thermocouple_profile = np.empty(config.n_nodes, dtype=float)
    collector_inner = np.empty(config.n_nodes, dtype=float)
    water_inlet = float(config.water_inlet_temperature_k)
    mcp = config.water_mass_flow_kg_s * config.water_cp_j_kg_k

    for idx, q_node_w in enumerate(heat_source_w):
        delta_t_water = float(q_node_w) / mcp
        water_mean = water_inlet + 0.5 * delta_t_water
        water_bulk[idx] = water_mean
        sleeve_outer[idx] = water_mean + float(q_node_w) * r_sleeve_to_water
        sleeve_thermocouple_profile[idx] = water_mean + float(q_node_w) * r_sleeve_tc_to_water
        collector_inner[idx] = water_mean + float(q_node_w) * r_collector_to_water
        water_inlet += delta_t_water

    water_outlet = water_inlet
    energy_error = float(heat_source_w.sum()) - mcp * (water_outlet - config.water_inlet_temperature_k)
    return BenkeThermalNetworkResult(
        case=case,
        config=config,
        active_zone_power_w=p_az_w,
        heat_source_w=heat_source_w,
        water_bulk_temperature_k=water_bulk,
        water_bulk_outlet_k=water_outlet,
        collector_inner_temperature_k=collector_inner,
        sleeve_outer_temperature_k=sleeve_outer,
        sleeve_thermocouple_temperature_k=_sample_benke_thermocouple_values(sleeve_thermocouple_profile, config.active_length_m),
        sleeve_thermocouple_z_mm=BENKE_THERMOCOUPLE_POSITIONS_MM,
        sleeve_thermocouple_radius_m=sleeve_thermocouple_r,
        sleeve_thermocouple_included_in_benke_average=benke_average_inclusion_flags(),
        r_collector_to_water_k_per_w=np.full(config.n_nodes, r_collector_to_water, dtype=float),
        r_sleeve_to_water_k_per_w=np.full(config.n_nodes, r_sleeve_to_water, dtype=float),
        energy_balance_error_w=float(energy_error),
    )
