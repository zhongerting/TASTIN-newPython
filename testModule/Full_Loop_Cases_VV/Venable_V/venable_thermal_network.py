from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class VenableThermalNetworkConfig:
    cooling_water_inlet_temperature_k: float = 310.0
    cooling_water_mass_flow_kg_s: float = 0.03
    cooling_water_cp_j_kg_k: float = 4180.0
    water_heat_transfer_coefficient_w_m2_k: float = 800.0
    regulated_he_gap_effective_k_w_m_k: float = 0.08
    regulated_he_gap_m: float = 0.0005
    unregulated_he_gap_effective_k_w_m_k: float = 0.276
    unregulated_he_gap_m: float = 0.00005
    tisa_heated_length_m: float = 0.300
    heat_pickup_fraction: float = 1.0
    collector_extra_resistance_k_per_w: float = 0.0
    cs_gap_effective_k_w_m_k: float = 0.12
    collector_k_w_m_k: float = 115.0
    alumina_k_w_m_k: float = 25.0
    sleeve_k_w_m_k: float = 16.0
    alumina_thickness_m: float = 0.00015
    sleeve_outer_radius_m: float = 0.0299 / 2.0
    emitter_outer_radius_m: float = 0.0196 / 2.0
    collector_inner_radius_m: float = 0.0206 / 2.0
    collector_wall_thickness_m: float = 0.0014


@dataclass(frozen=True)
class VenableThermalNetworkResult:
    config: VenableThermalNetworkConfig
    heat_source_w: np.ndarray
    water_bulk_temperature_k: np.ndarray
    water_bulk_outlet_k: float
    collector_temperature_k: np.ndarray
    emitter_temperature_k: np.ndarray
    r_cs_gap_k_per_w: np.ndarray
    r_collector_to_water_k_per_w: np.ndarray
    energy_balance_error_w: float


def cylindrical_conduction_resistance_k_per_w(
    inner_radius_m: float,
    outer_radius_m: float,
    conductivity_w_m_k: float,
    length_m: float,
) -> float:
    if inner_radius_m <= 0.0 or outer_radius_m <= 0.0:
        raise ValueError("Radii must be positive.")
    if outer_radius_m < inner_radius_m:
        raise ValueError("Outer radius must be no smaller than inner radius.")
    if conductivity_w_m_k <= 0.0:
        raise ValueError("Thermal conductivity must be positive.")
    if length_m <= 0.0:
        raise ValueError("Length must be positive.")
    if outer_radius_m == inner_radius_m:
        return 0.0
    return math.log(outer_radius_m / inner_radius_m) / (2.0 * math.pi * conductivity_w_m_k * length_m)


def convection_resistance_k_per_w(radius_m: float, h_w_m2_k: float, length_m: float) -> float:
    if radius_m <= 0.0:
        raise ValueError("Radius must be positive.")
    if h_w_m2_k <= 0.0:
        raise ValueError("Heat-transfer coefficient must be positive.")
    if length_m <= 0.0:
        raise ValueError("Length must be positive.")
    return 1.0 / (h_w_m2_k * 2.0 * math.pi * radius_m * length_m)


def build_tisa_heat_source_w(
    q_total_w: float,
    n_nodes: int,
    active_length_m: float,
    heated_length_m: float,
) -> np.ndarray:
    if q_total_w < 0.0:
        raise ValueError("Total heat source cannot be negative.")
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
    return float(q_total_w) * shape / total_shape


def _validate_config(config: VenableThermalNetworkConfig) -> None:
    positive_fields = (
        "cooling_water_mass_flow_kg_s",
        "cooling_water_cp_j_kg_k",
        "water_heat_transfer_coefficient_w_m2_k",
        "regulated_he_gap_effective_k_w_m_k",
        "unregulated_he_gap_effective_k_w_m_k",
        "cs_gap_effective_k_w_m_k",
        "collector_k_w_m_k",
        "alumina_k_w_m_k",
        "sleeve_k_w_m_k",
    )
    for field in positive_fields:
        if getattr(config, field) <= 0.0:
            raise ValueError(f"{field} must be positive.")
    if not (0.0 < config.heat_pickup_fraction <= 1.0):
        raise ValueError("heat_pickup_fraction must be in (0, 1].")
    if config.regulated_he_gap_m < 0.0 or config.unregulated_he_gap_m < 0.0:
        raise ValueError("Helium gap widths cannot be negative.")
    if config.collector_extra_resistance_k_per_w < 0.0:
        raise ValueError("collector_extra_resistance_k_per_w cannot be negative.")


def solve_thermal_network(case, geometry, config: VenableThermalNetworkConfig | None = None) -> VenableThermalNetworkResult:
    config = VenableThermalNetworkConfig() if config is None else config
    _validate_config(config)
    if geometry.n_nodes <= 0:
        raise ValueError("n_nodes must be positive.")
    if geometry.active_length_m <= 0.0:
        raise ValueError("active_length_m must be positive.")

    q_to_water_total_w = float(case.q_az_w) * config.heat_pickup_fraction
    heat_source_w = build_tisa_heat_source_w(
        q_to_water_total_w,
        geometry.n_nodes,
        geometry.active_length_m,
        config.tisa_heated_length_m,
    )
    dz_m = geometry.active_length_m / geometry.n_nodes
    n_nodes = geometry.n_nodes

    collector_inner_radius_m = config.collector_inner_radius_m
    collector_outer_radius_m = collector_inner_radius_m + config.collector_wall_thickness_m
    alumina_outer_radius_m = collector_outer_radius_m + config.alumina_thickness_m
    unregulated_he_outer_radius_m = alumina_outer_radius_m + config.unregulated_he_gap_m
    sleeve_outer_radius_m = max(config.sleeve_outer_radius_m, unregulated_he_outer_radius_m)
    regulated_he_outer_radius_m = sleeve_outer_radius_m + config.regulated_he_gap_m

    r_cs = cylindrical_conduction_resistance_k_per_w(
        config.emitter_outer_radius_m,
        collector_inner_radius_m,
        config.cs_gap_effective_k_w_m_k,
        dz_m,
    )
    r_collector = cylindrical_conduction_resistance_k_per_w(
        collector_inner_radius_m,
        collector_outer_radius_m,
        config.collector_k_w_m_k,
        dz_m,
    )
    r_alumina = cylindrical_conduction_resistance_k_per_w(
        collector_outer_radius_m,
        alumina_outer_radius_m,
        config.alumina_k_w_m_k,
        dz_m,
    )
    r_he_unreg = cylindrical_conduction_resistance_k_per_w(
        alumina_outer_radius_m,
        unregulated_he_outer_radius_m,
        config.unregulated_he_gap_effective_k_w_m_k,
        dz_m,
    )
    r_sleeve = cylindrical_conduction_resistance_k_per_w(
        unregulated_he_outer_radius_m,
        sleeve_outer_radius_m,
        config.sleeve_k_w_m_k,
        dz_m,
    )
    r_he_reg = cylindrical_conduction_resistance_k_per_w(
        sleeve_outer_radius_m,
        regulated_he_outer_radius_m,
        config.regulated_he_gap_effective_k_w_m_k,
        dz_m,
    )
    r_water = convection_resistance_k_per_w(
        regulated_he_outer_radius_m,
        config.water_heat_transfer_coefficient_w_m2_k,
        dz_m,
    )
    r_extra_node = config.collector_extra_resistance_k_per_w * n_nodes
    r_collector_to_water = r_collector + r_alumina + r_he_unreg + r_sleeve + r_he_reg + r_water + r_extra_node

    water_bulk_temperature_k = np.empty(n_nodes, dtype=float)
    collector_temperature_k = np.empty(n_nodes, dtype=float)
    emitter_temperature_k = np.empty(n_nodes, dtype=float)
    water_inlet_k = float(config.cooling_water_inlet_temperature_k)
    mcp_w_k = config.cooling_water_mass_flow_kg_s * config.cooling_water_cp_j_kg_k
    for idx, q_node_w in enumerate(heat_source_w):
        water_delta_node_k = float(q_node_w) / mcp_w_k
        water_mean_k = water_inlet_k + 0.5 * water_delta_node_k
        water_bulk_temperature_k[idx] = water_mean_k
        collector_temperature_k[idx] = water_mean_k + float(q_node_w) * r_collector_to_water
        emitter_temperature_k[idx] = collector_temperature_k[idx] + float(q_node_w) * r_cs
        water_inlet_k += water_delta_node_k

    water_bulk_outlet_k = water_inlet_k
    energy_balance_error_w = q_to_water_total_w - mcp_w_k * (
        water_bulk_outlet_k - config.cooling_water_inlet_temperature_k
    )
    return VenableThermalNetworkResult(
        config=config,
        heat_source_w=heat_source_w,
        water_bulk_temperature_k=water_bulk_temperature_k,
        water_bulk_outlet_k=water_bulk_outlet_k,
        collector_temperature_k=collector_temperature_k,
        emitter_temperature_k=emitter_temperature_k,
        r_cs_gap_k_per_w=np.full(n_nodes, r_cs, dtype=float),
        r_collector_to_water_k_per_w=np.full(n_nodes, r_collector_to_water, dtype=float),
        energy_balance_error_w=float(energy_balance_error_w),
    )
