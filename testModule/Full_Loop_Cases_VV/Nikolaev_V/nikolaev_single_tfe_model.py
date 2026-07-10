from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Sequence

from nikolaev_source_data import CapillaryLimitPoint, FuelTemperaturePoint, Table2OperatingPoint


@dataclass(frozen=True)
class NikolaevModelConfig:
    """Compact TOPAZ-II-like single TFE mock-up model.

    The paper gives calculated integral operating points rather than enough
    boundary conditions to rebuild the original Russian electro-thermal code.
    This model keeps the explicit table values traceable and exposes the
    missing closure parameters as adjustable coefficients.
    """

    nominal_output_power_w: float = 300.0
    nominal_efficiency_percent: float = 7.3
    reference_voltage_v: float = 0.8
    thermal_power_reference_kw: float = 4.1
    low_voltage_thermal_power_slope_kw_per_v: float = 1.0
    emitter_temp_reference_k: float = 1890.0
    emitter_temp_linear_k_per_v: float = 150.0
    emitter_temp_quadratic_k_per_v2: float = 500.0
    collector_temperature_k: float = 870.0
    active_core_length_m: float = 0.40
    effective_height_m: float = 0.70
    heater_length_m: float = 0.350
    interelectrode_gap_m: float = 0.0005


@dataclass(frozen=True)
class NikolaevOperatingResult:
    voltage_v: float
    current_a: float
    electric_power_w: float
    thermal_power_kw: float
    emitter_temperature_k: float
    collector_temperature_k: float
    efficiency_percent: float


def calculate_operating_point(voltage_v: float, config: NikolaevModelConfig = NikolaevModelConfig()) -> NikolaevOperatingResult:
    if voltage_v <= 0.0:
        raise ValueError("voltage_v must be positive.")
    dv = float(voltage_v) - config.reference_voltage_v
    current_a = config.nominal_output_power_w / float(voltage_v)
    thermal_power_kw = config.thermal_power_reference_kw + max(0.0, -dv) * config.low_voltage_thermal_power_slope_kw_per_v
    emitter_temperature_k = (
        config.emitter_temp_reference_k
        + config.emitter_temp_linear_k_per_v * dv
        + config.emitter_temp_quadratic_k_per_v2 * dv * dv
    )
    electric_power_w = float(voltage_v) * current_a
    efficiency_percent = 100.0 * electric_power_w / (thermal_power_kw * 1000.0)
    return NikolaevOperatingResult(
        voltage_v=float(voltage_v),
        current_a=current_a,
        electric_power_w=electric_power_w,
        thermal_power_kw=thermal_power_kw,
        emitter_temperature_k=emitter_temperature_k,
        collector_temperature_k=config.collector_temperature_k,
        efficiency_percent=efficiency_percent,
    )


def _linear_interpolate(x: float, x0: float, x1: float, y0: float, y1: float) -> float:
    if x1 == x0:
        return float(y0)
    weight = (float(x) - float(x0)) / (float(x1) - float(x0))
    return float(y0) + weight * (float(y1) - float(y0))


def bilinear_table_value(points: Sequence, free_volume_percent: float, radial_factor: float, value_attr: str) -> float:
    volumes = sorted({float(point.free_volume_percent) for point in points})
    radials = sorted({float(point.radial_factor) for point in points})
    if free_volume_percent < volumes[0] or free_volume_percent > volumes[-1]:
        raise ValueError("free_volume_percent is outside the Nikolaev table range.")
    if radial_factor < radials[0] or radial_factor > radials[-1]:
        raise ValueError("radial_factor is outside the Nikolaev table range.")

    def bracket(values: list[float], x: float) -> tuple[float, float]:
        for idx, value in enumerate(values):
            if math.isclose(x, value):
                return value, value
            if x < value:
                return values[idx - 1], value
        return values[-1], values[-1]

    v0, v1 = bracket(volumes, float(free_volume_percent))
    r0, r1 = bracket(radials, float(radial_factor))
    lookup = {(float(point.free_volume_percent), float(point.radial_factor)): float(getattr(point, value_attr)) for point in points}
    y00 = lookup[(v0, r0)]
    y10 = lookup[(v1, r0)]
    y01 = lookup[(v0, r1)]
    y11 = lookup[(v1, r1)]
    y_r0 = _linear_interpolate(free_volume_percent, v0, v1, y00, y10)
    y_r1 = _linear_interpolate(free_volume_percent, v0, v1, y01, y11)
    return _linear_interpolate(radial_factor, r0, r1, y_r0, y_r1)


def fuel_max_temperature_k(
    free_volume_percent: float,
    radial_factor: float,
    points: Iterable[FuelTemperaturePoint],
) -> float:
    return bilinear_table_value(tuple(points), free_volume_percent, radial_factor, "max_fuel_temperature_k")


def capillary_limit_diameter_mm(
    free_volume_percent: float,
    radial_factor: float,
    points: Iterable[CapillaryLimitPoint],
) -> float:
    return bilinear_table_value(tuple(points), free_volume_percent, radial_factor, "max_capillary_diameter_mm")


def compare_operating_point(target: Table2OperatingPoint, calculated: NikolaevOperatingResult) -> dict:
    return {
        "voltage_v": target.voltage_v,
        "current_exp_a": target.current_a,
        "current_calc_a": calculated.current_a,
        "current_error_a": calculated.current_a - target.current_a,
        "electric_power_exp_w": target.electric_power_w,
        "electric_power_calc_w": calculated.electric_power_w,
        "electric_power_error_w": calculated.electric_power_w - target.electric_power_w,
        "thermal_power_exp_kw": target.thermal_power_kw,
        "thermal_power_calc_kw": calculated.thermal_power_kw,
        "thermal_power_error_kw": calculated.thermal_power_kw - target.thermal_power_kw,
        "emitter_temperature_exp_k": target.emitter_temperature_k,
        "emitter_temperature_calc_k": calculated.emitter_temperature_k,
        "emitter_temperature_error_k": calculated.emitter_temperature_k - target.emitter_temperature_k,
        "collector_temperature_calc_k": calculated.collector_temperature_k,
        "efficiency_exp_percent": target.efficiency_percent,
        "efficiency_calc_percent": calculated.efficiency_percent,
        "efficiency_error_percent": calculated.efficiency_percent - target.efficiency_percent,
    }
