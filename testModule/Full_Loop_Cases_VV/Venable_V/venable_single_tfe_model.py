from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np

from venable_table71_data import VenableTable71Case
from venable_thermal_network import (
    VenableThermalNetworkConfig,
    VenableThermalNetworkResult,
    solve_thermal_network,
)


TCS_MODE_PLACEHOLDER = "placeholder_table"
TCS_MODE_PRESSURE_FORMULA = "pressure_formula"
VALID_TCS_MODES = (TCS_MODE_PLACEHOLDER, TCS_MODE_PRESSURE_FORMULA)
THERMAL_MODEL_EMPIRICAL = "empirical"
THERMAL_MODEL_THERMAL_NETWORK_V1 = "thermal_network_v1"
VALID_THERMAL_MODEL_MODES = (THERMAL_MODEL_EMPIRICAL, THERMAL_MODEL_THERMAL_NETWORK_V1)
AXIAL_PROFILE_COSINE = "cosine"
AXIAL_PROFILE_TISA_300MM = "tisa_300mm"
VALID_AXIAL_PROFILE_MODES = (AXIAL_PROFILE_COSINE, AXIAL_PROFILE_TISA_300MM)
COLLECTOR_BOUNDARY_LINEAR = "linear"
COLLECTOR_BOUNDARY_BENKE_WATER_JACKET = "benke_water_jacket"
VALID_COLLECTOR_BOUNDARY_MODES = (COLLECTOR_BOUNDARY_LINEAR, COLLECTOR_BOUNDARY_BENKE_WATER_JACKET)

SOURCE_ACTIVE_LENGTH_M = 0.375
SOURCE_EMITTER_OUTER_DIAMETER_M = 0.0196
SOURCE_EMITTER_WALL_THICKNESS_M = 0.00115
SOURCE_COLLECTOR_INNER_DIAMETER_M = 0.0206
SOURCE_COLLECTOR_WALL_THICKNESS_M = 0.0014
SOURCE_CS_GAP_M = 0.0005
SOURCE_SLEEVE_OUTER_DIAMETER_M = 0.0299


def _annulus_area_m2(inner_radius_m: float, outer_radius_m: float) -> float:
    return math.pi * (outer_radius_m**2 - inner_radius_m**2)


def _source_emitter_cross_area_m2() -> float:
    outer_radius_m = SOURCE_EMITTER_OUTER_DIAMETER_M / 2.0
    inner_radius_m = outer_radius_m - SOURCE_EMITTER_WALL_THICKNESS_M
    return _annulus_area_m2(inner_radius_m, outer_radius_m)


def _source_collector_cross_area_m2() -> float:
    inner_radius_m = SOURCE_COLLECTOR_INNER_DIAMETER_M / 2.0
    outer_radius_m = inner_radius_m + SOURCE_COLLECTOR_WALL_THICKNESS_M
    return _annulus_area_m2(inner_radius_m, outer_radius_m)


@dataclass(frozen=True)
class VenableTfeGeometry:
    active_length_m: float = SOURCE_ACTIVE_LENGTH_M
    n_nodes: int = 50
    gap_m: float = SOURCE_CS_GAP_M
    emitter_cross_area_m2: float = _source_emitter_cross_area_m2()
    collector_cross_area_m2: float = _source_collector_cross_area_m2()
    emitter_side_area_total_m2: float = math.pi * SOURCE_EMITTER_OUTER_DIAMETER_M * SOURCE_ACTIVE_LENGTH_M
    collector_side_area_total_m2: float = math.pi * SOURCE_COLLECTOR_INNER_DIAMETER_M * SOURCE_ACTIVE_LENGTH_M


@dataclass(frozen=True)
class VenableThermalClosure:
    name: str = "placeholder_qaz_scaled_temperature_fields"
    thermal_model_mode: str = THERMAL_MODEL_EMPIRICAL
    emitter_mean_min_k: float = 1450.0
    emitter_mean_max_k: float = 1900.0
    collector_mean_min_k: float = 760.0
    collector_mean_max_k: float = 870.0
    axial_shape_amplitude: float = 0.04
    axial_profile_mode: str = AXIAL_PROFILE_COSINE
    tisa_heated_length_m: float = 0.300
    emitter_quadratic_peak_k: float = 0.0
    collector_quadratic_peak_k: float = 0.0
    collector_boundary_mode: str = COLLECTOR_BOUNDARY_LINEAR
    cooling_water_inlet_temperature_k: float = 310.0
    cooling_water_mass_flow_kg_s: float = 0.03
    cooling_water_cp_j_kg_k: float = 4180.0
    coolant_heat_pickup_fraction: float = 0.4
    water_heat_transfer_coefficient_w_m2_k: float = 800.0
    regulated_he_gap_effective_k_w_m_k: float = 0.08
    regulated_he_gap_m: float = 0.0005
    unregulated_he_gap_effective_k_w_m_k: float = 0.276
    unregulated_he_gap_m: float = 0.00005
    thermal_network_heat_pickup_fraction: float = 1.0
    cs_gap_effective_k_w_m_k: float = 0.12
    collector_extra_resistance_k_per_w: float = 0.0


@dataclass(frozen=True)
class ThermoCalcInputArrays:
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
class VenableCaseModel:
    case: VenableTable71Case
    geometry: VenableTfeGeometry
    thermal_closure: VenableThermalClosure
    arrays: ThermoCalcInputArrays
    tcs_mode: str = TCS_MODE_PLACEHOLDER
    runs_thermocalc_calculation: bool = False
    thermal_network_result: VenableThermalNetworkResult | None = None


DEFAULT_GEOMETRY = VenableTfeGeometry()
DEFAULT_THERMAL_CLOSURE = VenableThermalClosure()


_TCS_BY_PCS_TORR: Dict[float, float] = {
    0.4: 560.0,
    0.5: 570.0,
    0.8: 590.0,
    1.0: 600.0,
}


def cesium_pressure_from_tcs(tcs_k: float) -> float:
    """Return Cs pressure in torr from Cs reservoir temperature in K.

    Source: ThermoCalc production lookup workflow, same relation documented in
    ``ThermoCalc/EMISSION_SCAN_GUIDE.md`` and used by
    ``ThermoCalc/tools/scan_emission_map.py``.
    """
    if tcs_k <= 0.0:
        raise ValueError("Cs reservoir temperature must be positive.")
    return 2.45e8 / math.sqrt(float(tcs_k)) * math.exp(-8910.0 / float(tcs_k))


def tcs_from_cesium_pressure(pcs_torr: float) -> float:
    """Invert the production Cs pressure relation to obtain Tcs in K."""
    if pcs_torr <= 0.0:
        raise ValueError("Cs pressure must be positive.")
    target = float(pcs_torr)
    lo = 300.0
    hi = 1200.0
    while cesium_pressure_from_tcs(lo) > target:
        lo *= 0.8
    while cesium_pressure_from_tcs(hi) < target:
        hi *= 1.2
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if cesium_pressure_from_tcs(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cesium_reservoir_temperature_k(pcs_torr: float, tcs_mode: str = TCS_MODE_PLACEHOLDER) -> float:
    """Return Cs reservoir temperature for the selected input closure."""
    if tcs_mode == TCS_MODE_PRESSURE_FORMULA:
        return tcs_from_cesium_pressure(float(pcs_torr))
    if tcs_mode != TCS_MODE_PLACEHOLDER:
        raise ValueError(f"Unknown Tcs mode: {tcs_mode!r}; expected one of {VALID_TCS_MODES}")
    try:
        return _TCS_BY_PCS_TORR[float(pcs_torr)]
    except KeyError as exc:
        raise ValueError(f"No placeholder Tcs value for Pcs={pcs_torr:g} torr") from exc


def _qaz_fraction(case: VenableTable71Case) -> float:
    q_min = 892.0
    q_max = 3162.0
    return (case.q_az_w - q_min) / (q_max - q_min)


def _temperature_profile(
    mean_k: float,
    n_nodes: int,
    amplitude: float,
    profile_mode: str = AXIAL_PROFILE_COSINE,
    active_length_m: float = SOURCE_ACTIVE_LENGTH_M,
    heated_length_m: float = 0.300,
) -> np.ndarray:
    centers = (np.arange(n_nodes, dtype=float) + 0.5) / n_nodes
    if profile_mode == AXIAL_PROFILE_COSINE:
        shape = 1.0 + amplitude * np.cos(np.pi * (centers - 0.5))
        return (mean_k * shape).reshape(1, n_nodes)
    if profile_mode != AXIAL_PROFILE_TISA_300MM:
        raise ValueError(f"Unknown axial profile mode: {profile_mode!r}; expected one of {VALID_AXIAL_PROFILE_MODES}")
    if heated_length_m <= 0.0 or heated_length_m > active_length_m:
        raise ValueError("TISA heated length must be positive and no longer than active length.")

    x_m = centers * active_length_m
    margin_m = 0.5 * (active_length_m - heated_length_m)
    raw = np.zeros(n_nodes, dtype=float)
    heated = (x_m >= margin_m) & (x_m <= margin_m + heated_length_m)
    raw[heated] = 1.0
    if margin_m > 0.0:
        left = x_m < margin_m
        right = x_m > margin_m + heated_length_m
        raw[left] = x_m[left] / margin_m
        raw[right] = (active_length_m - x_m[right]) / margin_m
    raw = np.clip(raw, 0.0, 1.0)
    raw_mean = float(np.mean(raw))
    if raw_mean <= 0.0:
        raise ValueError("TISA axial profile has zero mean.")
    normalized = raw / raw_mean
    shape = 1.0 + amplitude * (normalized - 1.0)
    return (mean_k * shape).reshape(1, n_nodes)



def _collector_mean_from_benke_water_jacket(
    case: VenableTable71Case,
    geometry: VenableTfeGeometry,
    thermal_closure: VenableThermalClosure,
) -> float:
    if thermal_closure.cooling_water_mass_flow_kg_s <= 0.0:
        raise ValueError("Cooling water mass flow must be positive.")
    if thermal_closure.cooling_water_cp_j_kg_k <= 0.0:
        raise ValueError("Cooling water heat capacity must be positive.")
    if thermal_closure.water_heat_transfer_coefficient_w_m2_k <= 0.0:
        raise ValueError("Water heat-transfer coefficient must be positive.")
    if thermal_closure.regulated_he_gap_effective_k_w_m_k <= 0.0:
        raise ValueError("Regulated He gap effective conductivity must be positive.")
    if thermal_closure.regulated_he_gap_m < 0.0:
        raise ValueError("Regulated He gap width cannot be negative.")
    if not (0.0 < thermal_closure.coolant_heat_pickup_fraction <= 1.0):
        raise ValueError("Coolant heat pickup fraction must be in (0, 1].")

    q_to_water_w = case.q_az_w * thermal_closure.coolant_heat_pickup_fraction
    water_delta_k = q_to_water_w / (
        thermal_closure.cooling_water_mass_flow_kg_s * thermal_closure.cooling_water_cp_j_kg_k
    )
    water_bulk_mean_k = thermal_closure.cooling_water_inlet_temperature_k + 0.5 * water_delta_k

    sleeve_outer_radius_m = SOURCE_SLEEVE_OUTER_DIAMETER_M / 2.0
    gap_outer_radius_m = sleeve_outer_radius_m + thermal_closure.regulated_he_gap_m
    heat_area_m2 = 2.0 * math.pi * gap_outer_radius_m * geometry.active_length_m
    r_water_k_per_w = 1.0 / (thermal_closure.water_heat_transfer_coefficient_w_m2_k * heat_area_m2)
    if thermal_closure.regulated_he_gap_m == 0.0:
        r_he_k_per_w = 0.0
    else:
        r_he_k_per_w = math.log(gap_outer_radius_m / sleeve_outer_radius_m) / (
            2.0
            * math.pi
            * thermal_closure.regulated_he_gap_effective_k_w_m_k
            * geometry.active_length_m
        )
    return water_bulk_mean_k + q_to_water_w * (
        r_water_k_per_w + r_he_k_per_w + thermal_closure.collector_extra_resistance_k_per_w
    )
def _thermal_network_config_from_closure(thermal_closure: VenableThermalClosure) -> VenableThermalNetworkConfig:
    return VenableThermalNetworkConfig(
        cooling_water_inlet_temperature_k=thermal_closure.cooling_water_inlet_temperature_k,
        cooling_water_mass_flow_kg_s=thermal_closure.cooling_water_mass_flow_kg_s,
        cooling_water_cp_j_kg_k=thermal_closure.cooling_water_cp_j_kg_k,
        water_heat_transfer_coefficient_w_m2_k=thermal_closure.water_heat_transfer_coefficient_w_m2_k,
        regulated_he_gap_effective_k_w_m_k=thermal_closure.regulated_he_gap_effective_k_w_m_k,
        regulated_he_gap_m=thermal_closure.regulated_he_gap_m,
        unregulated_he_gap_effective_k_w_m_k=thermal_closure.unregulated_he_gap_effective_k_w_m_k,
        unregulated_he_gap_m=thermal_closure.unregulated_he_gap_m,
        tisa_heated_length_m=thermal_closure.tisa_heated_length_m,
        heat_pickup_fraction=thermal_closure.thermal_network_heat_pickup_fraction,
        collector_extra_resistance_k_per_w=thermal_closure.collector_extra_resistance_k_per_w,
        cs_gap_effective_k_w_m_k=thermal_closure.cs_gap_effective_k_w_m_k,
    )


def _empirical_temperature_arrays(
    case: VenableTable71Case,
    geometry: VenableTfeGeometry,
    thermal_closure: VenableThermalClosure,
) -> tuple[np.ndarray, np.ndarray]:
    q_fraction = _qaz_fraction(case)
    quadratic_shape = 4.0 * q_fraction * (1.0 - q_fraction)
    emitter_mean_k = (
        thermal_closure.emitter_mean_min_k
        + q_fraction
        * (thermal_closure.emitter_mean_max_k - thermal_closure.emitter_mean_min_k)
        + thermal_closure.emitter_quadratic_peak_k * quadratic_shape
    )
    if thermal_closure.collector_boundary_mode == COLLECTOR_BOUNDARY_LINEAR:
        collector_mean_k = (
            thermal_closure.collector_mean_min_k
            + q_fraction
            * (thermal_closure.collector_mean_max_k - thermal_closure.collector_mean_min_k)
            + thermal_closure.collector_quadratic_peak_k * quadratic_shape
        )
    elif thermal_closure.collector_boundary_mode == COLLECTOR_BOUNDARY_BENKE_WATER_JACKET:
        collector_mean_k = _collector_mean_from_benke_water_jacket(case, geometry, thermal_closure)
    else:
        raise ValueError(
            f"Unknown collector boundary mode: {thermal_closure.collector_boundary_mode!r}; "
            f"expected one of {VALID_COLLECTOR_BOUNDARY_MODES}"
        )
    return (
        _temperature_profile(
            emitter_mean_k,
            geometry.n_nodes,
            thermal_closure.axial_shape_amplitude,
            thermal_closure.axial_profile_mode,
            geometry.active_length_m,
            thermal_closure.tisa_heated_length_m,
        ),
        _temperature_profile(
            collector_mean_k,
            geometry.n_nodes,
            thermal_closure.axial_shape_amplitude * 0.5,
            thermal_closure.axial_profile_mode,
            geometry.active_length_m,
            thermal_closure.tisa_heated_length_m,
        ),
    )


def build_case_model(
    case: VenableTable71Case,
    geometry: VenableTfeGeometry = DEFAULT_GEOMETRY,
    thermal_closure: VenableThermalClosure = DEFAULT_THERMAL_CLOSURE,
    tcs_mode: str = TCS_MODE_PLACEHOLDER,
) -> VenableCaseModel:
    if geometry.n_nodes <= 0:
        raise ValueError("n_nodes must be positive")
    if geometry.active_length_m <= 0.0:
        raise ValueError("active_length_m must be positive")
    if thermal_closure.thermal_model_mode not in VALID_THERMAL_MODEL_MODES:
        raise ValueError(
            f"Unknown thermal model mode: {thermal_closure.thermal_model_mode!r}; "
            f"expected one of {VALID_THERMAL_MODEL_MODES}"
        )

    thermal_network_result = None
    if thermal_closure.thermal_model_mode == THERMAL_MODEL_EMPIRICAL:
        temitter_k, tcollector_k = _empirical_temperature_arrays(case, geometry, thermal_closure)
    else:
        thermal_network_result = solve_thermal_network(
            case,
            geometry,
            _thermal_network_config_from_closure(thermal_closure),
        )
        temitter_k = thermal_network_result.emitter_temperature_k.reshape(1, geometry.n_nodes)
        tcollector_k = thermal_network_result.collector_temperature_k.reshape(1, geometry.n_nodes)

    shape = (1, geometry.n_nodes)
    arrays = ThermoCalcInputArrays(
        temitter_k=temitter_k,
        tcollector_k=tcollector_k,
        tcs_k=np.full(shape, cesium_reservoir_temperature_k(case.pcs_torr, tcs_mode=tcs_mode), dtype=float),
        dl_emitter_m=np.full(shape, geometry.active_length_m / geometry.n_nodes, dtype=float),
        dl_collector_m=np.full(shape, geometry.active_length_m / geometry.n_nodes, dtype=float),
        side_area_emitter_m2=np.full(
            shape, geometry.emitter_side_area_total_m2 / geometry.n_nodes, dtype=float
        ),
        side_area_collector_m2=np.full(
            shape, geometry.collector_side_area_total_m2 / geometry.n_nodes, dtype=float
        ),
        emitter_cross_area_m2=np.full(shape, geometry.emitter_cross_area_m2, dtype=float),
        collector_cross_area_m2=np.full(shape, geometry.collector_cross_area_m2, dtype=float),
    )
    return VenableCaseModel(
        case=case,
        geometry=geometry,
        thermal_closure=thermal_closure,
        arrays=arrays,
        tcs_mode=tcs_mode,
        thermal_network_result=thermal_network_result,
    )

def build_all_case_models(
    geometry: VenableTfeGeometry = DEFAULT_GEOMETRY,
    thermal_closure: VenableThermalClosure = DEFAULT_THERMAL_CLOSURE,
    tcs_mode: str = TCS_MODE_PLACEHOLDER,
) -> List[VenableCaseModel]:
    from venable_table71_data import iter_table71_cases

    return [
        build_case_model(
            case, geometry=geometry, thermal_closure=thermal_closure, tcs_mode=tcs_mode
        )
        for case in iter_table71_cases()
    ]


def model_summary_row(model: VenableCaseModel) -> Dict[str, object]:
    arrays = model.arrays
    return {
        "case_id": model.case.case_id,
        "q_az_w": model.case.q_az_w,
        "p_out_exp_w": model.case.p_out_exp_w,
        "eta_exp_percent": model.case.eta_exp_percent,
        "pcs_torr": model.case.pcs_torr,
        "tcs_k": float(arrays.tcs_k[0, 0]),
        "tcs_mode": model.tcs_mode,
        "n_nodes": model.geometry.n_nodes,
        "active_length_m": model.geometry.active_length_m,
        "gap_m": model.geometry.gap_m,
        "temitter_min_k": float(np.min(arrays.temitter_k)),
        "temitter_max_k": float(np.max(arrays.temitter_k)),
        "tcollector_min_k": float(np.min(arrays.tcollector_k)),
        "tcollector_max_k": float(np.max(arrays.tcollector_k)),
        "runs_thermocalc_calculation": model.runs_thermocalc_calculation,
    }


def model_config_summary(models: List[VenableCaseModel]) -> Dict[str, object]:
    if not models:
        raise ValueError("At least one model is required")
    first = models[0]
    return {
        "case_family": "Venable 1995 Table 7-1 single TFE maximum output power setup",
        "table_data_source": "TOPAZII_VV_public_experimental_data.md, Venable Table 7-1",
        "no_validation_calculation": True,
        "qaz_note": "Table 7-1 active-zone power is used directly; no additional 0.88 factor.",
        "geometry": {
            **asdict(first.geometry),
            "source": "TOPAZII_VV_public_experimental_data.md Benke/Venable single-cell dimensions: active length 375 mm, emitter OD 19.6 mm, emitter thickness 1.15 mm, collector ID 20.6 mm, collector thickness 1.4 mm, Cs gap 0.5 mm",
        },
        "thermal_closure": {
            **asdict(first.thermal_closure),
            "status": "empirical baseline or thermal_network_v1 case-local steady heat-transfer closure",
            "tcs_mode": first.tcs_mode,
            "tcs_by_pcs_torr_k": _TCS_BY_PCS_TORR,
            "pressure_formula": "Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)",
        },
        "thermocalc_input_shape": [
            1,
            first.geometry.n_nodes,
        ],
        "case_count": len(models),
    }
