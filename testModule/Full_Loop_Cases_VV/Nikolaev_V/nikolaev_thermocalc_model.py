from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from nikolaev_source_data import Table2OperatingPoint


@dataclass(frozen=True)
class NikolaevTfeGeometry:
    active_length_m: float = 0.500
    heater_length_m: float = 0.350
    n_nodes: int = 50
    emitter_outer_diameter_m: float = 0.0173
    emitter_wall_thickness_m: float = 0.0023
    interelectrode_gap_m: float = 0.0005
    collector_wall_thickness_m: float = 0.00185

    @property
    def emitter_outer_radius_m(self) -> float:
        return 0.5 * self.emitter_outer_diameter_m

    @property
    def emitter_inner_radius_m(self) -> float:
        return self.emitter_outer_radius_m - self.emitter_wall_thickness_m

    @property
    def collector_inner_radius_m(self) -> float:
        return self.emitter_outer_radius_m + self.interelectrode_gap_m

    @property
    def collector_outer_radius_m(self) -> float:
        return self.collector_inner_radius_m + self.collector_wall_thickness_m

    @property
    def emitter_cross_area_m2(self) -> float:
        return math.pi * (self.emitter_outer_radius_m**2 - self.emitter_inner_radius_m**2)

    @property
    def collector_cross_area_m2(self) -> float:
        return math.pi * (self.collector_outer_radius_m**2 - self.collector_inner_radius_m**2)

    @property
    def emitter_side_area_total_m2(self) -> float:
        return math.pi * self.emitter_outer_diameter_m * self.active_length_m

    @property
    def collector_side_area_total_m2(self) -> float:
        return 2.0 * math.pi * self.collector_inner_radius_m * self.active_length_m


@dataclass(frozen=True)
class NikolaevThermalNetworkConfig:
    n_nodes: int = 50
    collector_boundary_temperature_k: float = 870.0
    emitter_to_collector_resistance_k_per_w: float = 0.248
    axial_shape_amplitude: float = 0.18
    cesium_reservoir_temperature_k: float = 610.0
    wire_resistance_ohm: float = 0.0


@dataclass(frozen=True)
class NikolaevThermoCalcArrays:
    temitter_k: np.ndarray
    tcollector_k: np.ndarray
    tcs_k: np.ndarray
    dl_emitter_m: np.ndarray
    dl_collector_m: np.ndarray
    side_area_emitter_m2: np.ndarray
    side_area_collector_m2: np.ndarray
    emitter_cross_area_m2: np.ndarray
    collector_cross_area_m2: np.ndarray


@dataclass(frozen=True)
class NikolaevThermoCalcCase:
    source_point: Table2OperatingPoint
    geometry: NikolaevTfeGeometry
    thermal_config: NikolaevThermalNetworkConfig
    arrays: NikolaevThermoCalcArrays
    heat_source_w: np.ndarray


def centered_heater_power_profile(total_power_w: float, geometry: NikolaevTfeGeometry) -> np.ndarray:
    if total_power_w <= 0.0:
        raise ValueError("total_power_w must be positive.")
    if geometry.n_nodes <= 0:
        raise ValueError("n_nodes must be positive.")
    if geometry.heater_length_m <= 0.0 or geometry.heater_length_m > geometry.active_length_m:
        raise ValueError("heater length must be positive and no longer than active length.")
    centers = (np.arange(geometry.n_nodes, dtype=float) + 0.5) * geometry.active_length_m / geometry.n_nodes
    margin = 0.5 * (geometry.active_length_m - geometry.heater_length_m)
    shape = np.zeros(geometry.n_nodes, dtype=float)
    heated = (centers >= margin) & (centers <= margin + geometry.heater_length_m)
    shape[heated] = 1.0
    if margin > 0.0:
        left = centers < margin
        right = centers > margin + geometry.heater_length_m
        shape[left] = centers[left] / margin
        shape[right] = (geometry.active_length_m - centers[right]) / margin
    shape = np.clip(shape, 0.0, None)
    if float(shape.sum()) <= 0.0:
        raise ValueError("heater power profile has zero integral.")
    return float(total_power_w) * shape / float(shape.sum())


def build_thermocalc_case(
    point: Table2OperatingPoint,
    thermal_config: NikolaevThermalNetworkConfig = NikolaevThermalNetworkConfig(),
) -> NikolaevThermoCalcCase:
    geometry = NikolaevTfeGeometry(n_nodes=thermal_config.n_nodes)
    if geometry.emitter_inner_radius_m <= 0.0:
        raise ValueError("emitter wall is too thick for the specified outer diameter.")
    heat_source_w = centered_heater_power_profile(point.thermal_power_kw * 1000.0, geometry)
    node_resistance_k_per_w = thermal_config.emitter_to_collector_resistance_k_per_w * geometry.n_nodes
    local_delta_k = heat_source_w * node_resistance_k_per_w
    mean_delta_k = float(np.mean(local_delta_k))
    shaped_delta_k = mean_delta_k + (local_delta_k - mean_delta_k) * thermal_config.axial_shape_amplitude
    temitter = thermal_config.collector_boundary_temperature_k + shaped_delta_k
    tcollector = np.full(geometry.n_nodes, thermal_config.collector_boundary_temperature_k, dtype=float)
    shape = (1, geometry.n_nodes)
    arrays = NikolaevThermoCalcArrays(
        temitter_k=temitter.reshape(shape),
        tcollector_k=tcollector.reshape(shape),
        tcs_k=np.full(shape, thermal_config.cesium_reservoir_temperature_k, dtype=float),
        dl_emitter_m=np.full(shape, geometry.active_length_m / geometry.n_nodes, dtype=float),
        dl_collector_m=np.full(shape, geometry.active_length_m / geometry.n_nodes, dtype=float),
        side_area_emitter_m2=np.full(shape, geometry.emitter_side_area_total_m2 / geometry.n_nodes, dtype=float),
        side_area_collector_m2=np.full(shape, geometry.collector_side_area_total_m2 / geometry.n_nodes, dtype=float),
        emitter_cross_area_m2=np.full(shape, geometry.emitter_cross_area_m2, dtype=float),
        collector_cross_area_m2=np.full(shape, geometry.collector_cross_area_m2, dtype=float),
    )
    return NikolaevThermoCalcCase(
        source_point=point,
        geometry=geometry,
        thermal_config=thermal_config,
        arrays=arrays,
        heat_source_w=heat_source_w,
    )
