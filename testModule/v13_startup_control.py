from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional

import numpy as np

from Components.ReactorCore import ControlDrumReactivityModel
from Components.ExternalHeatSources import BaseExternalHeatSource, ExternalHeatFluxBC
from Components.ExternalHeatSources.embedded_flux_tables import W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY


@dataclass
class V13StartupControlConfig:
    """Simplified TITAM-style TOPAZ-II cold-start schedule for V13-start."""

    source_power_w: float = 1.0
    fixed_power_w: Optional[float] = None
    initial_temperature_k: float = 373.0
    safety_drum_duration_s: float = 8.0
    safety_drum_worth_dollars: float = 2.0
    control_drum_speed_deg_s: float = 1.4
    critical_control_angle_deg: float = 125.0
    supercritical_control_angle_deg: float = 154.0
    pullback_control_angle_deg: float = 145.0
    final_control_angle_deg: float = 88.0
    initial_supercritical_duration_s: float = 20.0
    pullback_duration_s: float = 10.0
    low_power_w: float = 5_000.0
    fast_ramp_target_w: float = 35_000.0
    steady_power_w: float = 110_000.0
    fast_power_ramp_w_s: float = 600.0
    slow_power_ramp_w_s: float = 80.0
    shield_jettison_temperature_k: float = 400.0
    tfe_start_after_critical_s: float = 1500.0
    tfe_start_emitter_temperature_k: float = 1050.0
    tec_electrical_start_after_cesium_s: float = 0.0
    tec_electrical_start_cs_fraction: float = 0.0
    tec_electrical_start_emitter_temperature_k: float = 0.0
    helium_gap_h_eq_w_m2_k: float = 1200.0
    cesium_gap_h_eq_w_m2_k: float = 29.0
    cs_transition_tau_s: float = 120.0
    beta_total: float = ControlDrumReactivityModel.default_beta_total

    def __post_init__(self) -> None:
        self.safety_drum_duration_s = max(1.0e-12, float(self.safety_drum_duration_s))
        self.control_drum_speed_deg_s = max(1.0e-12, float(self.control_drum_speed_deg_s))
        self.initial_supercritical_duration_s = max(0.0, float(self.initial_supercritical_duration_s))
        self.pullback_duration_s = max(0.0, float(self.pullback_duration_s))
        self.fast_power_ramp_w_s = max(1.0e-12, float(self.fast_power_ramp_w_s))
        self.slow_power_ramp_w_s = max(1.0e-12, float(self.slow_power_ramp_w_s))
        self.cs_transition_tau_s = max(1.0e-12, float(self.cs_transition_tau_s))
        self.tec_electrical_start_after_cesium_s = max(0.0, float(self.tec_electrical_start_after_cesium_s))
        self.tec_electrical_start_cs_fraction = float(np.clip(self.tec_electrical_start_cs_fraction, 0.0, 1.0))
        self.critical_time_s = (
            self.safety_drum_duration_s
            + self.critical_control_angle_deg / self.control_drum_speed_deg_s
        )
        self.supercritical_end_s = self.critical_time_s + self.initial_supercritical_duration_s
        self.low_power_hold_start_s = self.supercritical_end_s + self.pullback_duration_s
        self.fast_ramp_end_s = (
            self.low_power_hold_start_s
            + (self.fast_ramp_target_w - self.low_power_w) / self.fast_power_ramp_w_s
        )
        self.slow_ramp_end_s = (
            self.fast_ramp_end_s
            + (self.steady_power_w - self.fast_ramp_target_w) / self.slow_power_ramp_w_s
        )


@dataclass
class V13StartupCommand:
    absolute_time_s: float
    time_after_critical_s: float
    phase: str
    thermal_power_w: float
    fission_power_w: float
    decay_power_w: float
    safety_drum_angle_deg: float
    control_drum_angle_deg: float
    safety_reactivity_dollars: float
    control_reactivity_dollars: float
    total_startup_reactivity_dollars: float
    total_startup_reactivity: float
    radiation_shield_active: bool
    shield_jettisoned: bool
    cesium_conditioning_started: bool
    tec_enabled: bool
    cs_fraction: float
    tec_gap_h_eq_w_m2_k: float

    def diagnostics(self) -> Dict[str, Any]:
        return {
            "startup_phase": self.phase,
            "startup_time_after_critical_s": float(self.time_after_critical_s),
            "startup_thermal_power_w": float(self.thermal_power_w),
            "startup_fission_power_w": float(self.fission_power_w),
            "startup_decay_power_w": float(self.decay_power_w),
            "startup_safety_drum_angle_deg": float(self.safety_drum_angle_deg),
            "startup_control_drum_angle_deg": float(self.control_drum_angle_deg),
            "startup_safety_reactivity_dollars": float(self.safety_reactivity_dollars),
            "startup_control_reactivity_dollars": float(self.control_reactivity_dollars),
            "startup_total_reactivity_dollars": float(self.total_startup_reactivity_dollars),
            "startup_total_reactivity": float(self.total_startup_reactivity),
            "startup_shield_jettisoned": bool(self.shield_jettisoned),
            "startup_cesium_conditioning_started": bool(self.cesium_conditioning_started),
            "startup_tec_enabled": bool(self.tec_enabled),
            "startup_cs_fraction": float(self.cs_fraction),
            "startup_tec_gap_h_eq_w_m2_k": float(self.tec_gap_h_eq_w_m2_k),
        }


class V13StartupController:
    """Stateful startup sequencer used by the V13-start runner."""

    def __init__(self, config: Optional[V13StartupControlConfig] = None):
        self.config = config or V13StartupControlConfig()
        self.shield_jettisoned = False
        self.cesium_conditioning_started = False
        self.cesium_start_time_s: Optional[float] = None
        self.tec_electrical_started = False
        self.tec_electrical_start_time_s: Optional[float] = None

    def seed_cesium_conditioning(self, absolute_time_s: float, cs_fraction: float) -> None:
        fraction = float(np.clip(cs_fraction, 0.0, 1.0))
        if fraction <= 0.0:
            return
        t = float(absolute_time_s)
        self.cesium_conditioning_started = True
        if fraction >= 1.0:
            elapsed = 50.0 * self.config.cs_transition_tau_s
        else:
            elapsed = -self.config.cs_transition_tau_s * float(np.log(max(1.0 - fraction, 1.0e-12)))
        self.cesium_start_time_s = t - elapsed

    def _phase(self, t: float) -> str:
        c = self.config
        if t < c.safety_drum_duration_s:
            return "SAFETY_DRUM_WITHDRAWAL"
        if t < c.critical_time_s:
            return "CONTROL_DRUM_APPROACH"
        if t < c.supercritical_end_s:
            return "INITIAL_SUPERCRITICAL_RAMP"
        if t < c.low_power_hold_start_s:
            return "REACTIVITY_PULLBACK"
        if abs(t - c.low_power_hold_start_s) <= 1.0e-10:
            return "LOW_POWER_HOLD"
        if t < c.fast_ramp_end_s:
            return "FAST_POWER_RAMP"
        if t < c.slow_ramp_end_s:
            return "SLOW_POWER_RAMP"
        return "CRITICAL_POWER_HOLD"

    def _safety_angle(self, t: float) -> float:
        c = self.config
        return float(np.clip(180.0 * t / c.safety_drum_duration_s, 0.0, 180.0))

    def _control_angle(self, t: float) -> float:
        c = self.config
        if t <= c.safety_drum_duration_s:
            return 0.0
        if t < c.critical_time_s:
            return float(np.clip((t - c.safety_drum_duration_s) * c.control_drum_speed_deg_s, 0.0, c.critical_control_angle_deg))
        if t < c.supercritical_end_s and c.initial_supercritical_duration_s > 0.0:
            frac = (t - c.critical_time_s) / c.initial_supercritical_duration_s
            return float(c.critical_control_angle_deg + frac * (c.supercritical_control_angle_deg - c.critical_control_angle_deg))
        if t < c.low_power_hold_start_s and c.pullback_duration_s > 0.0:
            frac = (t - c.supercritical_end_s) / c.pullback_duration_s
            return float(c.supercritical_control_angle_deg + frac * (c.pullback_control_angle_deg - c.supercritical_control_angle_deg))
        if t < c.slow_ramp_end_s:
            denom = max(c.slow_ramp_end_s - c.low_power_hold_start_s, 1.0e-12)
            frac = (t - c.low_power_hold_start_s) / denom
            return float(c.pullback_control_angle_deg + frac * (c.final_control_angle_deg - c.pullback_control_angle_deg))
        return float(c.final_control_angle_deg)

    def _thermal_power(self, t: float) -> float:
        c = self.config
        if c.fixed_power_w is not None:
            return float(c.fixed_power_w)
        if t < c.critical_time_s:
            return float(c.source_power_w)
        if t < c.low_power_hold_start_s:
            denom = max(c.low_power_hold_start_s - c.critical_time_s, 1.0e-12)
            frac = (t - c.critical_time_s) / denom
            return float(c.source_power_w + frac * (c.low_power_w - c.source_power_w))
        if t < c.fast_ramp_end_s:
            return float(c.low_power_w + (t - c.low_power_hold_start_s) * c.fast_power_ramp_w_s)
        if t < c.slow_ramp_end_s:
            return float(c.fast_ramp_target_w + (t - c.fast_ramp_end_s) * c.slow_power_ramp_w_s)
        return float(c.steady_power_w)

    def _safety_reactivity_dollars(self, safety_angle_deg: float) -> float:
        c = self.config
        withdrawn_fraction = float(np.clip(safety_angle_deg / 180.0, 0.0, 1.0))
        return -float(c.safety_drum_worth_dollars) * (1.0 - withdrawn_fraction)

    def _control_reactivity_dollars(self, control_angle_deg: float) -> float:
        model = ControlDrumReactivityModel(
            enabled=True,
            theta_deg=control_angle_deg,
            reference_theta_deg=0.0,
            cold_reference_keff=0.952,
        )
        return float(model.total_reactivity_dollars(self.config.beta_total))

    def _update_latches(self, t: float, core_inlet_temperature_k: float, emitter_temperature_k: float) -> None:
        c = self.config
        if not self.shield_jettisoned and core_inlet_temperature_k >= c.shield_jettison_temperature_k:
            self.shield_jettisoned = True
        if (
            not self.cesium_conditioning_started
            and t - c.critical_time_s >= c.tfe_start_after_critical_s - 1.0e-9
        ):
            self.cesium_conditioning_started = True
            self.cesium_start_time_s = float(t)
        if self.cesium_conditioning_started and not self.tec_electrical_started:
            elapsed = 0.0 if self.cesium_start_time_s is None else float(t) - self.cesium_start_time_s
            if (
                elapsed >= c.tec_electrical_start_after_cesium_s
                and self._cs_fraction(t) >= c.tec_electrical_start_cs_fraction
                and emitter_temperature_k >= c.tec_electrical_start_emitter_temperature_k
            ):
                self.tec_electrical_started = True
                self.tec_electrical_start_time_s = float(t)

    def _cs_fraction(self, t: float) -> float:
        if not self.cesium_conditioning_started or self.cesium_start_time_s is None:
            return 0.0
        return 1.0

    def evaluate(
            self,
            absolute_time_s: float,
            *,
            core_inlet_temperature_k: float,
            emitter_temperature_k: float) -> V13StartupCommand:
        t = max(0.0, float(absolute_time_s))
        c = self.config
        self._update_latches(
            t,
            float(core_inlet_temperature_k),
            float(emitter_temperature_k),
        )
        safety_angle = self._safety_angle(t)
        control_angle = self._control_angle(t)
        safety_rho_dollars = self._safety_reactivity_dollars(safety_angle)
        control_rho_dollars = self._control_reactivity_dollars(control_angle)
        total_rho_dollars = safety_rho_dollars + control_rho_dollars
        cs_fraction = self._cs_fraction(t)
        h_eq = (
            (1.0 - cs_fraction) * c.helium_gap_h_eq_w_m2_k
            + cs_fraction * c.cesium_gap_h_eq_w_m2_k
        )
        thermal_power = self._thermal_power(t)
        return V13StartupCommand(
            absolute_time_s=t,
            time_after_critical_s=t - c.critical_time_s,
            phase=self._phase(t),
            thermal_power_w=thermal_power,
            fission_power_w=thermal_power,
            decay_power_w=0.0,
            safety_drum_angle_deg=safety_angle,
            control_drum_angle_deg=control_angle,
            safety_reactivity_dollars=safety_rho_dollars,
            control_reactivity_dollars=control_rho_dollars,
            total_startup_reactivity_dollars=total_rho_dollars,
            total_startup_reactivity=total_rho_dollars * c.beta_total,
            radiation_shield_active=not self.shield_jettisoned,
            shield_jettisoned=self.shield_jettisoned,
            cesium_conditioning_started=self.cesium_conditioning_started,
            tec_enabled=self.tec_electrical_started,
            cs_fraction=cs_fraction,
            tec_gap_h_eq_w_m2_k=float(h_eq),
        )


def _iter_tfes(tfes_or_core: Any) -> Iterable[Any]:
    if hasattr(tfes_or_core, "tfes"):
        tfes_or_core = getattr(tfes_or_core, "tfes")
    if isinstance(tfes_or_core, dict):
        return tfes_or_core.values()
    return tfes_or_core


def apply_tec_gap_h_eq(tfes_or_core: Any, h_eq_w_m2_k: float) -> int:
    """Update simplified emitter-collector gap gas conductance on all TFE units."""
    count = 0
    for tfe in _iter_tfes(tfes_or_core):
        couplers = getattr(tfe, "couplers", {})
        coupler = couplers.get("tec_couple") if isinstance(couplers, dict) else None
        if coupler is None:
            continue
        gap_width = max(float(getattr(coupler, "gap", 0.0)), 0.0)
        setattr(coupler, "k_gas", float(h_eq_w_m2_k) * gap_width)
        count += 1
    return count
def reset_fluid_temperatures(system_or_network: Any, temperature_k: float) -> int:
    """Reset all hydraulic control-volume temperatures and enthalpies to a uniform value."""
    network = getattr(system_or_network, "fluid_solver", system_or_network)
    volumes = list(getattr(network, "volumes_obj", []))
    temperature = float(temperature_k)
    count = 0
    for vol in volumes:
        material = getattr(vol, "material", None)
        if material is None:
            continue
        vol.T = temperature
        vol.h = float(material.enthalpy(vol.T, vol.P))
        if hasattr(vol, "update_properties"):
            vol.update_properties(material)
        for attr in ("Q_wall", "Q_vol", "implicit_coeff", "source_explicit", "source_implicit"):
            if hasattr(vol, attr):
                setattr(vol, attr, 0.0)
        count += 1
    if hasattr(network, "_initialize_state_from_objects"):
        network._initialize_state_from_objects()
    else:
        for idx, vol in enumerate(volumes):
            if hasattr(network, "T_vec"):
                network.T_vec[idx] = vol.T
            if hasattr(network, "h_vec"):
                network.h_vec[idx] = vol.h
            if hasattr(network, "P_vec"):
                network.P_vec[idx] = vol.P
    if hasattr(network, "_update_fluid_properties"):
        network._update_fluid_properties()
    return count

def reset_solid_temperatures(solids_or_system: Any, temperature_k: float, *, current_time_s: float = 0.0) -> int:
    """Reset registered solid heat-conduction objects to a uniform cold-start temperature."""
    if hasattr(solids_or_system, "solid_components"):
        solids = getattr(solids_or_system, "solid_components").values()
    elif isinstance(solids_or_system, dict):
        solids = solids_or_system.values()
    else:
        solids = solids_or_system

    temperature = float(temperature_k)
    count = 0
    for solid in solids:
        if solid is None or not hasattr(solid, "T"):
            continue
        solid.T[...] = temperature
        if hasattr(solid, "dTdt"):
            solid.dTdt[...] = 0.0
        if hasattr(solid, "current_time"):
            solid.current_time = float(current_time_s)
        if hasattr(solid, "last_trial_temperature_min"):
            solid.last_trial_temperature_min = temperature
        if hasattr(solid, "last_trial_temperature_max"):
            solid.last_trial_temperature_max = temperature
        if hasattr(solid, "last_trial_temperature_time"):
            solid.last_trial_temperature_time = float(current_time_s)
        if hasattr(solid, "last_step_success"):
            solid.last_step_success = True
        if hasattr(solid, "last_step_failure_message"):
            solid.last_step_failure_message = ""
        if hasattr(solid, "_update_properties"):
            solid._update_properties()
        for boundary in getattr(solid, "boundaries", {}).values():
            if hasattr(boundary, "T_surface"):
                boundary.T_surface[...] = temperature
            if hasattr(boundary, "T_adj_node"):
                boundary.T_adj_node[...] = temperature
            if hasattr(boundary, "current_flux"):
                boundary.current_flux[...] = 0.0
        count += 1
    return count

class MatrixColumnHeatSource(BaseExternalHeatSource):
    """Broadcast one embedded matrix column as heat-flux density [W/m2]."""

    def __init__(
            self,
            shape,
            matrix_key: str,
            column_index: int,
            scale_factor: float = 1.0,
            offset: float = 0.0,
            periodic: bool = True,
            matrix_library=W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY):
        super().__init__(tuple(shape))
        self.matrix_key = str(matrix_key)
        self.column_index = int(column_index)
        self.scale_factor = float(scale_factor)
        self.offset = float(offset)
        self.periodic = bool(periodic)
        self.matrix_library = matrix_library
        self.matrix = self.matrix_library.get_matrix(self.matrix_key)
        if self.column_index < 0 or self.column_index >= self.matrix.values.shape[1]:
            raise IndexError(
                f"column_index {self.column_index} out of range for {self.matrix_key} "
                f"with {self.matrix.values.shape[1]} columns."
            )

    @staticmethod
    def _wrap_time(time: float, sample_time: np.ndarray) -> float:
        start = float(sample_time[0])
        end = float(sample_time[-1])
        step = float(sample_time[1] - sample_time[0]) if sample_time.size > 1 else 0.0
        period = (end - start) + step
        if period <= 0.0:
            return float(time)
        return ((float(time) - start) % period) + start

    def _sample_column(self, time: float) -> float:
        sample_time = float(time)
        if self.periodic and self.matrix.periodic:
            sample_time = self._wrap_time(sample_time, self.matrix.time)
        value = float(np.interp(sample_time, self.matrix.time, self.matrix.values[:, self.column_index]))
        return value * self.scale_factor + self.offset

    def get_heat_flux(self, time: float) -> np.ndarray:
        return self._broadcast_flux(self._sample_column(time))

    def update_params(self, **kwargs):
        if "scale_factor" in kwargs and kwargs["scale_factor"] is not None:
            self.scale_factor = float(kwargs["scale_factor"])
        if "offset" in kwargs and kwargs["offset"] is not None:
            self.offset = float(kwargs["offset"])
        if "periodic" in kwargs and kwargs["periodic"] is not None:
            self.periodic = bool(kwargs["periodic"])

class MatrixColumnLineHeatSource(MatrixColumnHeatSource):
    """Convert one embedded matrix column from line load [W/m] to equivalent [W/m2]."""

    def __init__(
            self,
            shape,
            matrix_key: str,
            column_index: int,
            node_lengths_m,
            area_array_m2,
            scale_factor: float = 1.0,
            offset: float = 0.0,
            periodic: bool = True,
            matrix_library=W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY):
        super().__init__(
            shape=shape,
            matrix_key=matrix_key,
            column_index=column_index,
            scale_factor=scale_factor,
            offset=offset,
            periodic=periodic,
            matrix_library=matrix_library,
        )
        self.node_lengths_m = np.asarray(node_lengths_m, dtype=float)
        if self.node_lengths_m.shape != self.shape:
            self.node_lengths_m = np.broadcast_to(self.node_lengths_m, self.shape).astype(float, copy=True)
        self.area_array_m2 = np.asarray(area_array_m2, dtype=float)
        if self.area_array_m2.shape != self.shape:
            self.area_array_m2 = np.broadcast_to(self.area_array_m2, self.shape).astype(float, copy=True)
        self.area_array_m2 = np.nan_to_num(self.area_array_m2, nan=0.0, posinf=0.0, neginf=0.0)
        self.node_lengths_m = np.nan_to_num(self.node_lengths_m, nan=0.0, posinf=0.0, neginf=0.0)

    def get_heat_flux(self, time: float) -> np.ndarray:
        q_line_w_m = self._sample_column(time)
        q_node_w = q_line_w_m * self.node_lengths_m
        q_density = np.zeros(self.shape, dtype=float)
        np.divide(q_node_w, self.area_array_m2, out=q_density, where=self.area_array_m2 > 0.0)
        return q_density

def matrix_values_at_time(matrix_key: str, time_s: float, *, scale_factor: float = 1.0, offset: float = 0.0, periodic: bool = True) -> np.ndarray:
    matrix = W0_8P12_ORBITAL_HEAT_MATRIX_LIBRARY.get_matrix(str(matrix_key))
    sample_time = float(time_s)
    if periodic and matrix.periodic:
        sample_time = MatrixColumnHeatSource._wrap_time(sample_time, matrix.time)
    values = np.empty(matrix.values.shape[1], dtype=float)
    for idx in range(matrix.values.shape[1]):
        values[idx] = np.interp(sample_time, matrix.time, matrix.values[:, idx])
    values *= float(scale_factor)
    if offset != 0.0:
        values += float(offset)
    return values


def shield_qsss_from_matrix(matrix_key: str, time_s: float, *, scale_factor: float = 1.0, offset: float = 0.0, periodic: bool = True) -> np.ndarray:
    values = matrix_values_at_time(
        matrix_key,
        time_s,
        scale_factor=scale_factor,
        offset=offset,
        periodic=periodic,
    )
    if values.size != 6:
        raise ValueError(f"Shield qsss matrix must have 6 columns, got {values.size}.")
    qsss = np.zeros(8, dtype=float)
    qsss[:6] = values
    return qsss

def radiator_external_heat_power_w(radiator_units: Iterable[Any]) -> float:
    """Return the currently applied radiator-tube external heat power [W]."""
    total = 0.0
    for unit in radiator_units:
        wall = getattr(unit, "wall", None)
        if wall is None or "right" not in getattr(wall, "boundaries", {}):
            continue
        boundary = wall.boundaries["right"]
        for condition in getattr(boundary, "conditions", []):
            source = getattr(condition, "heat_source", None)
            if isinstance(source, (MatrixColumnHeatSource, MatrixColumnLineHeatSource)) and hasattr(condition, "q_flux"):
                total += float(np.sum(np.asarray(condition.q_flux, dtype=float)))
    return total

def attach_radiator_tube_external_heat(
        radiator_units: Iterable[Any],
        *,
        matrix_key: str = "is58p5_w0_8p12_N78_sum",
        scale_factor: float = 1.0,
        offset: float = 0.0,
        periodic: bool = True,
        area_fraction: float = 1.0,
        input_units: str = "W/m") -> int:
    count = 0
    for index, unit in enumerate(radiator_units):
        wall = getattr(unit, "wall", None)
        if wall is None or "right" not in getattr(wall, "boundaries", {}):
            continue
        boundary = wall.boundaries["right"]
        area = np.asarray(boundary.area, dtype=float) * float(area_fraction)
        units = str(input_units).lower().replace(" ", "")
        if units in {"w/m", "wperm", "w_per_m"}:
            node_length = getattr(unit, "node_length", None)
            if node_length is None:
                tube_length = float(getattr(unit, "tube_length"))
                n_axial = int(getattr(unit, "n_axial"))
                node_length = tube_length / n_axial
            node_lengths = np.full(boundary.shape, float(node_length), dtype=float)
            source = MatrixColumnLineHeatSource(
                shape=boundary.shape,
                matrix_key=matrix_key,
                column_index=index,
                node_lengths_m=node_lengths,
                area_array_m2=area,
                scale_factor=scale_factor,
                offset=offset,
                periodic=periodic,
            )
        elif units in {"w/m2", "w/m^2", "wperm2", "w_per_m2"}:
            source = MatrixColumnHeatSource(
                shape=boundary.shape,
                matrix_key=matrix_key,
                column_index=index,
                scale_factor=scale_factor,
                offset=offset,
                periodic=periodic,
            )
        else:
            raise ValueError(f"Unsupported external heat input_units: {input_units!r}")
        boundary.conditions.append(ExternalHeatFluxBC(source, area))
        count += 1
    return count




