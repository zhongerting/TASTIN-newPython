"""Run V14 after instantaneous complete NaK loss without hydraulics."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from types import MethodType
from typing import Any, Dict, Optional, Sequence

import numpy as np

REPO_ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / "testModule").is_dir())
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Solvers.Couplers import FluidSolidCouple, GapCouple2D
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_debug.run_v14_210kw_debug import (
    _apply_fixed_core_power,
    _ring_rejection,
    _ring_wall_rejection,
    _tec_main_metrics,
    build_debug_case,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.run_v14_helium_depressurization import (
    refresh_tec_now,
    set_tec_update_interval,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.run_v14_210kw_reactivity_control import (
    ReactivityControlRunConfig,
    load_baseline_debug_config,
    prepare_reactivity_control,
)


CASE_DIR = Path(__file__).resolve().parent
CASE_NAME = "V14_10kW_210kW_fixed_power_LOCA_1"
DEFAULT_RESTART = (
    CASE_DIR.parent
    / "V14_210kW_fixed_power_external_heat_2orbits"
    / "runs"
    / "two_orbits_from13864_20260720"
    / "checkpoint_t016864s.npz"
)
DEFAULT_OUTPUT_DIR = CASE_DIR / "runs" / "default"
DEFAULT_EMISSIVITY = 0.8
DEFAULT_RECORD_INTERVAL_S = 0.2

SUMMARY_FIELDS = [
    "time_s", "accident_elapsed_s", "dt_s", "coolant_present",
    "hydraulic_solve_enabled", "core_power_W", "radiation_emissivity",
    "vacuum_gap_count", "removed_fluid_solid_coupler_count",
    "core_structure_min_T_K", "core_structure_max_T_K",
    "pipe_wall_min_T_K", "pipe_wall_max_T_K",
    "heat_pipe_min_T_K", "heat_pipe_max_T_K",
    "coolant_max_T_K", "collector_max_T_K", "emitter_max_T_K",
    "moderator_max_T_K", "reflector_max_T_K", "failure_reason",
    "point_kinetics_enabled", "fission_power_W", "decay_power_W",
    "external_reactivity", "external_reactivity_dollars",
    "effective_temperature_feedback", "total_reactivity",
    "tec_main_current_A", "tec_main_voltage_V",
    "tec_main_electric_power_W", "tec_main_converged",
    "tec_open_circuit_active", "tec_open_circuit_time_s",
    "tec_open_circuit_trigger_current_A",
    "feedback_fuel", "feedback_electrode", "feedback_moderator",
    "feedback_reflector", "feedback_total_absolute",
    "feedback_total_change_from_accident",
    "inner_outer_clad_radiation_W", "radiator_net_rejection_W",
    "solid_relative_energy_J", "solid_energy_rate_W",
    "post_accident_energy_residual_W", "snapshot_path",
]

COOLANT_HISTORY_FIELDS = [
    "time_s", "accident_elapsed_s", "coolant_present", "category",
    "entity_type", "name", "temperature_K", "pressure_Pa", "enthalpy_J_kg",
    "mass_flow_kg_s", "velocity_m_s", "reference_temperature_K",
    "reference_pressure_Pa", "reference_enthalpy_J_kg",
    "reference_mass_flow_kg_s", "reference_velocity_m_s",
]
SOLID_HISTORY_FIELDS = [
    "time_s", "accident_elapsed_s", "category", "solid_name", "solid_shape",
    "flat_node_index", "temperature_K",
]
ELECTRICAL_HISTORY_FIELDS = [
    "time_s", "accident_elapsed_s", "tec_main_current_A", "tec_main_voltage_V",
    "tec_main_electric_power_W", "tec_main_converged", "tfe_name",
    "tfe_multiplier", "axial_node_index", "current_density_A_m2",
    "emitter_potential_V", "collector_potential_V",
    "emitter_collector_voltage_drop_V", "electron_cooling_flux_W_m2",
    "electron_heating_flux_W_m2", "electron_cooling_power_W",
    "electron_heating_power_W", "emitter_joule_power_axial_W",
    "collector_joule_power_axial_W",
]
REACTIVITY_HISTORY_FIELDS = [
    "time_s", "accident_elapsed_s", "feedback_fuel", "feedback_electrode",
    "feedback_moderator", "feedback_reflector", "feedback_total_absolute",
    "feedback_total_change_from_accident",
    "point_kinetics_enabled", "core_power_W", "fission_power_W", "decay_power_W",
    "external_reactivity", "external_reactivity_dollars",
    "effective_temperature_feedback", "total_reactivity",
    "tec_open_circuit_active", "tec_open_circuit_time_s",
]


@dataclass(frozen=True)
class LocARunConfig:
    restart_in: Path = DEFAULT_RESTART
    output_dir: Path = DEFAULT_OUTPUT_DIR
    duration_s: float = 0.4
    dt_s: float = 0.05
    record_interval_s: float = DEFAULT_RECORD_INTERVAL_S
    checkpoint_interval_s: float = 0.0
    tec_update_interval_s: float = 0.05
    radiation_emissivity: float = DEFAULT_EMISSIVITY
    collector_failure_temperature_k: float = 1500.0
    emitter_failure_temperature_k: float = 3000.0
    coolant_failure_temperature_k: float = 1058.0
    moderator_failure_temperature_k: float = 930.0
    reflector_failure_temperature_k: float = 1000.0
    enable_reactivity_feedback: bool = False
    scram_time_s: Optional[float] = None
    scram_reactivity_dollars: float = -2.0
    staged_recording: bool = False
    tec_open_circuit_current_threshold_a: float = 0.01


def _read_source_config(restart: Path) -> Dict[str, Any]:
    path = Path(restart).parent / "run_config.json"
    if not path.is_file():
        raise FileNotFoundError(f"run_config.json not found beside restart: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _validate_config(config: LocARunConfig) -> None:
    if not Path(config.restart_in).is_file():
        raise FileNotFoundError(f"restart not found: {config.restart_in}")
    for name in ("duration_s", "dt_s", "record_interval_s", "tec_update_interval_s"):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not math.isfinite(float(config.checkpoint_interval_s)):
        raise ValueError("checkpoint_interval_s must be finite")
    if not 0.0 < float(config.radiation_emissivity) <= 1.0:
        raise ValueError("radiation_emissivity must be within (0, 1]")
    for name in (
        "collector_failure_temperature_k", "emitter_failure_temperature_k",
        "coolant_failure_temperature_k", "moderator_failure_temperature_k",
        "reflector_failure_temperature_k",
    ):
        value = float(getattr(config, name))
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if config.scram_time_s is not None:
        if not config.enable_reactivity_feedback:
            raise ValueError("scram requires reactivity feedback/point kinetics")
        if not math.isfinite(float(config.scram_time_s)) or float(config.scram_time_s) < 0.0:
            raise ValueError("scram_time_s must be finite and non-negative")
        if not math.isfinite(float(config.scram_reactivity_dollars)):
            raise ValueError("scram_reactivity_dollars must be finite")
    if (
        not math.isfinite(float(config.tec_open_circuit_current_threshold_a))
        or float(config.tec_open_circuit_current_threshold_a) < 0.0
    ):
        raise ValueError("tec_open_circuit_current_threshold_a must be finite and non-negative")


def _remove_boundary_condition(boundary: Any, condition: Any) -> None:
    if condition is None:
        return
    boundary.conditions[:] = [bc for bc in boundary.conditions if bc is not condition]
    boundary.clear_boundary_conditions()


def detach_fluid_solid_couplers(system: Any) -> list[FluidSolidCouple]:
    """Remove every liquid-solid path and its installed boundary conditions."""
    removed = [c for c in system.couplers if isinstance(c, FluidSolidCouple)]
    for coupler in removed:
        _remove_boundary_condition(coupler.solid_bound, coupler.solid_bc)
        _remove_boundary_condition(
            coupler.solid_bound,
            getattr(coupler, "_local_implicit_flux_bc", None),
        )
    system.couplers[:] = [c for c in system.couplers if not isinstance(c, FluidSolidCouple)]
    return removed


def install_core_vacuum_radiation(build: Dict[str, Any], emissivity: float) -> list[GapCouple2D]:
    gaps: list[GapCouple2D] = []
    for name, tfe in build["tfes"].items():
        for key in ("iclad_coolant_fsc", "oclad_coolant_fsc"):
            tfe.couplers.pop(key, None)
        width = float(tfe.geom.r_coolant_outer - tfe.geom.r_inner_clad_outer)
        gap = GapCouple2D(
            obj1=tfe.solids["inner_clad"],
            obj2=tfe.solids["outer_clad"],
            direction="right",
            gap_width=width,
            gas_conductivity=0.0,
            emissivity1=float(emissivity),
            emissivity2=float(emissivity),
        )
        gap.name = f"{name}_inner_outer_clad_vacuum_radiation"
        tfe.couplers["iclad_oclad_vacuum_gap"] = gap
        build["system"].couplers.append(gap)
        gaps.append(gap)
    return gaps


def _zero_fluid_sources_and_flow(system: Any) -> None:
    net = system.fluid_solver
    for volume in net.volumes_obj:
        for attr in ("Q_wall", "Q_vol", "implicit_coeff", "source_explicit", "source_implicit"):
            if hasattr(volume, attr):
                setattr(volume, attr, 0.0)
    net.W_vec.fill(0.0)
    for attr in ("W_old", "W_iterate"):
        if hasattr(net, attr):
            getattr(net, attr).fill(0.0)
    for junction in net.junctions_obj:
        junction.W = 0.0
        junction.vel = 0.0


def freeze_hydraulics(system: Any) -> None:
    """Keep the restored liquid arrays only as references; never solve them."""
    net = system.fluid_solver

    def frozen_step(this: Any, dt: float, *args: Any, **kwargs: Any) -> bool:
        _zero_fluid_sources_and_flow(system)
        return True

    net.step_Picard = MethodType(frozen_step, net)
    net.loca_hydraulics_frozen = True
    _zero_fluid_sources_and_flow(system)


def take_over_accident(
    build: Dict[str, Any], emissivity: float, *, fixed_power: bool = True,
) -> Dict[str, Any]:
    system = build["system"]
    core = build["core"]
    refresh_tec_now(core, float(system.global_time))
    reference = {
        "T": np.asarray(system.fluid_solver.T_vec, dtype=float).copy(),
        "P": np.asarray(system.fluid_solver.P_vec, dtype=float).copy(),
        "h": np.asarray(system.fluid_solver.h_vec, dtype=float).copy(),
        "W": np.asarray(system.fluid_solver.W_vec, dtype=float).copy(),
    }
    removed = detach_fluid_solid_couplers(system)
    gaps = install_core_vacuum_radiation(build, emissivity)
    freeze_hydraulics(system)
    system._refresh_solid_boundary_cache(update_flux=False, current_time=system.global_time)
    for gap in gaps:
        gap.sync()
    system._refresh_solid_boundary_cache(update_flux=True, current_time=system.global_time)
    if fixed_power:
        _apply_fixed_core_power(build, 210000.0)
    feedback = core.compute_reactivity_feedback()
    return {
        "reference_fluid": reference,
        "removed_couplers": removed,
        "vacuum_gaps": gaps,
        "feedback_reference_total": float(feedback.total),
        "solid_reference_T": {
            name: np.asarray(solid.T, dtype=float).copy()
            for name, solid in system.solid_components.items()
            if hasattr(solid, "T")
        },
        "tec_open_circuit_active": False,
        "tec_open_circuit_time_s": float("nan"),
        "tec_open_circuit_trigger_current_A": float("nan"),
    }


def classify_fluid_name(name: str) -> str:
    if name.startswith("Chan_") or name.startswith("CoreInletConnector") or name.startswith("CoreOutletConnector"):
        return "core"
    if name.startswith(("Upper_", "Lower_", "OutletMix_", "Manifold_")):
        return "collector_ring"
    return "ordinary_pipe"


def classify_solid_name(name: str) -> str:
    lower = name.lower()
    if "heatpipe" in lower or "_hp" in lower or "hpwithfin" in lower:
        return "heat_pipe"
    if (
        "ringwall" in lower or "pipe" in lower or "header" in lower
        or (name.startswith(("Upper_", "Lower_")) and name.endswith("_Solid"))
    ):
        return "pipe_wall"
    return "core_structure"


def staged_record_interval(elapsed_s: float) -> float:
    if elapsed_s < 20.0 - 1.0e-9:
        return 0.5
    if elapsed_s < 100.0 - 1.0e-9:
        return 2.0
    if elapsed_s < 400.0 - 1.0e-9:
        return 5.0
    if elapsed_s < 600.0 - 1.0e-9:
        return 10.0
    return 20.0


def external_reactivity_for_elapsed(core: Any, config: LocARunConfig, elapsed_s: float) -> tuple[float, float]:
    if config.scram_time_s is None or elapsed_s < float(config.scram_time_s) - 1.0e-9:
        return 0.0, 0.0
    dollars = float(config.scram_reactivity_dollars)
    return dollars * float(core.point_reactor.beta_total), dollars


def _neutronics_metrics(core: Any, external_reactivity: float, external_dollars: float) -> Dict[str, Any]:
    if not core.has_point_reactor:
        return {
            "point_kinetics_enabled": False,
            "fission_power_W": float("nan"), "decay_power_W": float("nan"),
            "external_reactivity": 0.0, "external_reactivity_dollars": 0.0,
            "effective_temperature_feedback": float("nan"),
            "total_reactivity": float("nan"),
        }
    point = core.point_reactor
    effective = float(core.get_effective_reactivity_feedback())
    return {
        "point_kinetics_enabled": True,
        "fission_power_W": float(point.fission_power),
        "decay_power_W": float(point.decay_power),
        "external_reactivity": float(external_reactivity),
        "external_reactivity_dollars": float(external_dollars),
        "effective_temperature_feedback": effective,
        "total_reactivity": float(external_reactivity) + effective,
    }


def maybe_transition_tec_to_open_circuit(
    build: Dict[str, Any], accident: Dict[str, Any], config: LocARunConfig, elapsed_s: float,
) -> bool:
    core = build["core"]
    if (
        config.scram_time_s is None
        or elapsed_s < float(config.scram_time_s) - 1.0e-9
        or accident["tec_open_circuit_active"]
        or not core.enable_tec_coupled
    ):
        return False
    current = float(_tec_main_metrics(build)["tec_main_current_A"])
    if not math.isfinite(current) or current > float(config.tec_open_circuit_current_threshold_a):
        return False
    for tfe in build["tfes"].values():
        tfe.clear_tec_sources()
    core.enable_tec_coupled = False
    accident["tec_open_circuit_active"] = True
    accident["tec_open_circuit_time_s"] = float(elapsed_s)
    accident["tec_open_circuit_trigger_current_A"] = current
    return True


def _solid_temperature_payload(system: Any) -> Dict[str, np.ndarray]:
    names, categories, shapes, offsets, values = [], [], [], [0], []
    for name, solid in system.solid_components.items():
        if not hasattr(solid, "T"):
            continue
        arr = np.asarray(solid.T, dtype=float)
        names.append(name)
        categories.append(classify_solid_name(name))
        shapes.append(json.dumps(arr.shape))
        values.append(arr.reshape(-1))
        offsets.append(offsets[-1] + arr.size)
    return {
        "solid_names": np.asarray(names),
        "solid_categories": np.asarray(categories),
        "solid_shapes": np.asarray(shapes),
        "solid_offsets": np.asarray(offsets, dtype=np.int64),
        "solid_temperature_K": np.concatenate(values) if values else np.array([], dtype=float),
    }


def _fluid_payload(system: Any, reference: Dict[str, np.ndarray], *, coolant_present: bool) -> Dict[str, np.ndarray]:
    net = system.fluid_solver
    volume_names = np.asarray([v.name for v in net.volumes_obj])
    volume_categories = np.asarray([classify_fluid_name(v.name) for v in net.volumes_obj])
    junction_names = np.asarray([j.name for j in net.junctions_obj])
    junction_categories = np.asarray([classify_fluid_name(j.name) for j in net.junctions_obj])
    reference_velocity = np.asarray([
        reference["W"][i] / max(float(j.area) * float(j.from_vol.rho), 1.0e-30)
        for i, j in enumerate(net.junctions_obj)
    ])
    n_vol = len(net.volumes_obj)
    if coolant_present:
        temperature = np.asarray(net.T_vec, dtype=float).copy()
        pressure = np.asarray(net.P_vec, dtype=float).copy()
        enthalpy = np.asarray(net.h_vec, dtype=float).copy()
        mass_flow = np.asarray(net.W_vec, dtype=float).copy()
        velocity = np.asarray([
            mass_flow[i] / max(float(j.area) * float(j.from_vol.rho), 1.0e-30)
            for i, j in enumerate(net.junctions_obj)
        ])
    else:
        temperature = np.full(n_vol, np.nan)
        pressure = np.full(n_vol, np.nan)
        enthalpy = np.full(n_vol, np.nan)
        mass_flow = np.zeros(len(net.junctions_obj))
        velocity = np.zeros(len(net.junctions_obj))
    return {
        "coolant_present": np.asarray([int(coolant_present)], dtype=np.int8),
        "fluid_volume_names": volume_names,
        "fluid_volume_categories": volume_categories,
        "fluid_temperature_K": temperature,
        "fluid_pressure_Pa": pressure,
        "fluid_enthalpy_J_kg": enthalpy,
        "fluid_reference_temperature_K": reference["T"],
        "fluid_reference_pressure_Pa": reference["P"],
        "fluid_reference_enthalpy_J_kg": reference["h"],
        "fluid_junction_names": junction_names,
        "fluid_junction_categories": junction_categories,
        "fluid_mass_flow_kg_s": mass_flow,
        "fluid_velocity_m_s": velocity,
        "fluid_reference_mass_flow_kg_s": reference["W"],
        "fluid_reference_velocity_m_s": reference_velocity,
    }


def _tec_payload(build: Dict[str, Any]) -> Dict[str, np.ndarray]:
    names = list(build["tfes"])
    tfes = [build["tfes"][name] for name in names]
    multipliers = build["ring_multipliers"]

    def stack(domain: str, attr: str) -> np.ndarray:
        return np.stack([
            np.asarray(getattr(getattr(tfe, domain), attr), dtype=float)
            for tfe in tfes
        ])

    side_area = np.stack([
        np.asarray(tfe.solids["emitter"].boundaries["right"].area, dtype=float)
        for tfe in tfes
    ])
    q_e_flux = stack("plasma_data", "electron_cooling_flux")
    q_c_flux = stack("plasma_data", "electron_heating_flux")
    return {
        "tfe_names": np.asarray(names),
        "tfe_multipliers": np.asarray([
            multipliers[name] if isinstance(multipliers, dict) else multipliers[index]
            for index, name in enumerate(names)
        ], dtype=int),
        "tec_current_density_A_m2": stack("electric_data", "current_density"),
        "tec_emitter_potential_V": stack("electric_data", "emitter_potential"),
        "tec_collector_potential_V": stack("electric_data", "collector_potential"),
        "tec_emitter_collector_voltage_drop_V": stack("electric_data", "emitter_collector_voltage_drop"),
        "tec_electron_cooling_flux_W_m2": q_e_flux,
        "tec_electron_heating_flux_W_m2": q_c_flux,
        "tec_electron_cooling_power_W": q_e_flux * side_area,
        "tec_electron_heating_power_W": q_c_flux * side_area,
        "tec_emitter_joule_power_axial_W": stack("electric_data", "emitter_joule_heat"),
        "tec_collector_joule_power_axial_W": stack("electric_data", "collector_joule_heat"),
    }


def _solid_min_max_by_category(system: Any) -> Dict[str, float]:
    grouped: Dict[str, list[np.ndarray]] = {key: [] for key in ("core_structure", "pipe_wall", "heat_pipe")}
    for name, solid in system.solid_components.items():
        if hasattr(solid, "T"):
            grouped[classify_solid_name(name)].append(np.asarray(solid.T, dtype=float).reshape(-1))
    result: Dict[str, float] = {}
    for category, arrays in grouped.items():
        if arrays:
            values = np.concatenate(arrays)
            result[f"{category}_min_T_K"] = float(np.min(values))
            result[f"{category}_max_T_K"] = float(np.max(values))
        else:
            result[f"{category}_min_T_K"] = float("nan")
            result[f"{category}_max_T_K"] = float("nan")
    return result


def _failure_temperature_metrics(system: Any, *, coolant_present: bool) -> Dict[str, float]:
    grouped: Dict[str, list[np.ndarray]] = {
        "collector": [], "emitter": [], "moderator": [], "reflector": [],
    }
    for name, solid in system.solid_components.items():
        if not hasattr(solid, "T"):
            continue
        lower = name.lower()
        for category in grouped:
            if lower.endswith(f"_{category}"):
                grouped[category].append(np.asarray(solid.T, dtype=float).reshape(-1))
                break
    result = {
        f"{category}_max_T_K": (
            float(np.max(np.concatenate(arrays))) if arrays else float("nan")
        )
        for category, arrays in grouped.items()
    }
    result["coolant_max_T_K"] = (
        float(np.max(np.asarray(system.fluid_solver.T_vec, dtype=float)))
        if coolant_present else float("nan")
    )
    return result


def evaluate_failure_reason(system: Any, config: LocARunConfig, *, coolant_present: bool) -> tuple[str, Dict[str, float]]:
    metrics = _failure_temperature_metrics(system, coolant_present=coolant_present)
    limits = (
        ("collector", config.collector_failure_temperature_k),
        ("emitter", config.emitter_failure_temperature_k),
        ("coolant", config.coolant_failure_temperature_k),
        ("moderator", config.moderator_failure_temperature_k),
        ("reflector", config.reflector_failure_temperature_k),
    )
    for name, limit in limits:
        value = metrics[f"{name}_max_T_K"]
        if math.isfinite(value) and value >= float(limit):
            return f"{name}_temperature_limit", metrics
    return "", metrics


def _relative_solid_energy(system: Any, reference: Dict[str, np.ndarray]) -> float:
    total = 0.0
    for name, solid in system.solid_components.items():
        if name not in reference or not hasattr(solid, "T"):
            continue
        temperature = np.asarray(solid.T, dtype=float)
        delta = temperature - reference[name]
        material = solid.material
        rho = np.asarray(material.density(temperature), dtype=float)
        cp = np.asarray(material.heat_capacity(temperature), dtype=float)
        volumes = np.asarray(solid.mesh.geom_data.volumes, dtype=float)
        total += float(np.sum(rho * cp * volumes * delta))
    return total


def _vacuum_gap_heat(gaps: list[GapCouple2D]) -> float:
    total = 0.0
    for gap in gaps:
        t1, _ = gap.bound1.get_coupling_surface_snapshot()
        t2, _ = gap.bound2.get_coupling_surface_snapshot()
        total += float(np.sum((np.asarray(t1) - np.asarray(t2)) / np.asarray(gap.R_gap_total)))
    return total


def _radiator_rejection(build: Dict[str, Any]) -> float:
    ring_hps = build.get("ring_hps", [])
    return float(_ring_rejection(ring_hps) + _ring_wall_rejection(build))


def collect_summary(
    build: Dict[str, Any], accident: Dict[str, Any], *, start_time: float,
    dt_s: float, snapshot_path: str, previous_energy: Optional[tuple[float, float]],
    failure_reason: str = "",
    external_reactivity: float = 0.0,
    external_reactivity_dollars: float = 0.0,
) -> tuple[Dict[str, Any], tuple[float, float]]:
    system, core = build["system"], build["core"]
    feedback = core.compute_reactivity_feedback()
    tec = _tec_main_metrics(build)
    energy = _relative_solid_energy(system, accident["solid_reference_T"])
    if previous_energy is None or float(system.global_time) <= previous_energy[0]:
        energy_rate = float("nan")
    else:
        energy_rate = (energy - previous_energy[1]) / (float(system.global_time) - previous_energy[0])
    rejection = _radiator_rejection(build)
    residual = (
        float(core.last_total_core_power)
        - float(tec["tec_main_electric_power_W"])
        - rejection
        - energy_rate
        if math.isfinite(energy_rate) else float("nan")
    )
    row = {
        "time_s": float(system.global_time),
        "accident_elapsed_s": float(system.global_time) - float(start_time),
        "dt_s": float(dt_s),
        "coolant_present": False,
        "hydraulic_solve_enabled": False,
        "core_power_W": float(core.last_total_core_power),
        "radiation_emissivity": float(accident["vacuum_gaps"][0].eps1),
        "vacuum_gap_count": len(accident["vacuum_gaps"]),
        "removed_fluid_solid_coupler_count": len(accident["removed_couplers"]),
        **_solid_min_max_by_category(system),
        **_failure_temperature_metrics(system, coolant_present=False),
        "failure_reason": failure_reason,
        "tec_open_circuit_active": bool(accident["tec_open_circuit_active"]),
        "tec_open_circuit_time_s": float(accident["tec_open_circuit_time_s"]),
        "tec_open_circuit_trigger_current_A": float(
            accident["tec_open_circuit_trigger_current_A"]
        ),
        **_neutronics_metrics(core, external_reactivity, external_reactivity_dollars),
        **tec,
        "feedback_fuel": float(feedback.fuel),
        "feedback_electrode": float(feedback.electrode),
        "feedback_moderator": float(feedback.moderator),
        "feedback_reflector": float(feedback.reflector),
        "feedback_total_absolute": float(feedback.total),
        "feedback_total_change_from_accident": float(feedback.total) - accident["feedback_reference_total"],
        "inner_outer_clad_radiation_W": _vacuum_gap_heat(accident["vacuum_gaps"]),
        "radiator_net_rejection_W": rejection,
        "solid_relative_energy_J": energy,
        "solid_energy_rate_W": energy_rate,
        "post_accident_energy_residual_W": residual,
        "snapshot_path": snapshot_path,
    }
    return row, (float(system.global_time), energy)


def build_snapshot_payload(
    build: Dict[str, Any], accident: Dict[str, Any], *,
    start_time: float, coolant_present: bool,
    hydraulic_solve_enabled: bool = False,
    external_reactivity: float = 0.0, external_reactivity_dollars: float = 0.0,
) -> Dict[str, Any]:
    system = build["system"]
    feedback = build["core"].compute_reactivity_feedback()
    tec = _tec_main_metrics(build)
    return {
        "time_s": np.asarray([system.global_time]),
        "accident_elapsed_s": np.asarray([system.global_time - start_time]),
        "hydraulic_solve_enabled": np.asarray([
            int(hydraulic_solve_enabled)
        ], dtype=np.int8),
        "feedback_fuel": np.asarray([feedback.fuel]),
        "feedback_electrode": np.asarray([feedback.electrode]),
        "feedback_moderator": np.asarray([feedback.moderator]),
        "feedback_reflector": np.asarray([feedback.reflector]),
        "feedback_total_absolute": np.asarray([feedback.total]),
        "feedback_total_change_from_accident": np.asarray([
            feedback.total - accident["feedback_reference_total"]
        ]),
        "tec_open_circuit_active": np.asarray([
            int(accident.get("tec_open_circuit_active", False))
        ], dtype=np.int8),
        "tec_open_circuit_time_s": np.asarray([
            accident.get("tec_open_circuit_time_s", float("nan"))
        ]),
        "core_power_W": np.asarray([build["core"].last_total_core_power]),
        **{
            key: np.asarray([value])
            for key, value in _neutronics_metrics(
                build["core"], external_reactivity, external_reactivity_dollars,
            ).items()
        },
        **{key: np.asarray([value]) for key, value in tec.items()},
        **_fluid_payload(system, accident["reference_fluid"], coolant_present=coolant_present),
        **_solid_temperature_payload(system),
        **_tec_payload(build),
    }


def write_snapshot(
    path: Path, build: Dict[str, Any], accident: Dict[str, Any], *,
    start_time: float, coolant_present: bool,
    hydraulic_solve_enabled: bool = False,
    external_reactivity: float = 0.0, external_reactivity_dollars: float = 0.0,
) -> Dict[str, Any]:
    payload = build_snapshot_payload(
        build, accident, start_time=start_time, coolant_present=coolant_present,
        hydraulic_solve_enabled=hydraulic_solve_enabled,
        external_reactivity=external_reactivity,
        external_reactivity_dollars=external_reactivity_dollars,
    )
    np.savez_compressed(path, **payload)
    return payload


def _append_summary(path: Path, row: Dict[str, Any]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def _append_rows(path: Path, fieldnames: list[str], rows: list[Dict[str, Any]]) -> None:
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def _append_coolant_history(out_dir: Path, payload: Dict[str, Any]) -> None:
    common = {
        "time_s": float(payload["time_s"][0]),
        "accident_elapsed_s": float(payload["accident_elapsed_s"][0]),
        "coolant_present": bool(payload["coolant_present"][0]),
    }
    rows = []
    for index, name in enumerate(payload["fluid_volume_names"]):
        rows.append({
            **common, "category": payload["fluid_volume_categories"][index],
            "entity_type": "volume", "name": name,
            "temperature_K": payload["fluid_temperature_K"][index],
            "pressure_Pa": payload["fluid_pressure_Pa"][index],
            "enthalpy_J_kg": payload["fluid_enthalpy_J_kg"][index],
            "mass_flow_kg_s": "", "velocity_m_s": "",
            "reference_temperature_K": payload["fluid_reference_temperature_K"][index],
            "reference_pressure_Pa": payload["fluid_reference_pressure_Pa"][index],
            "reference_enthalpy_J_kg": payload["fluid_reference_enthalpy_J_kg"][index],
            "reference_mass_flow_kg_s": "", "reference_velocity_m_s": "",
        })
    for index, name in enumerate(payload["fluid_junction_names"]):
        rows.append({
            **common, "category": payload["fluid_junction_categories"][index],
            "entity_type": "junction", "name": name,
            "temperature_K": "", "pressure_Pa": "", "enthalpy_J_kg": "",
            "mass_flow_kg_s": payload["fluid_mass_flow_kg_s"][index],
            "velocity_m_s": payload["fluid_velocity_m_s"][index],
            "reference_temperature_K": "", "reference_pressure_Pa": "",
            "reference_enthalpy_J_kg": "",
            "reference_mass_flow_kg_s": payload["fluid_reference_mass_flow_kg_s"][index],
            "reference_velocity_m_s": payload["fluid_reference_velocity_m_s"][index],
        })
    _append_rows(out_dir / "history_coolant.csv", COOLANT_HISTORY_FIELDS, rows)


def _append_solid_history(out_dir: Path, payload: Dict[str, Any]) -> None:
    rows = []
    for index, name in enumerate(payload["solid_names"]):
        begin, end = map(int, payload["solid_offsets"][index:index + 2])
        for node, temperature in enumerate(payload["solid_temperature_K"][begin:end]):
            rows.append({
                "time_s": payload["time_s"][0],
                "accident_elapsed_s": payload["accident_elapsed_s"][0],
                "category": payload["solid_categories"][index], "solid_name": name,
                "solid_shape": payload["solid_shapes"][index],
                "flat_node_index": node, "temperature_K": temperature,
            })
    _append_rows(out_dir / "history_solids.csv", SOLID_HISTORY_FIELDS, rows)


def _append_electrical_history(out_dir: Path, payload: Dict[str, Any]) -> None:
    rows = []
    axial_count = payload["tec_current_density_A_m2"].shape[1]
    field_map = (
        ("current_density_A_m2", "tec_current_density_A_m2"),
        ("emitter_potential_V", "tec_emitter_potential_V"),
        ("collector_potential_V", "tec_collector_potential_V"),
        ("emitter_collector_voltage_drop_V", "tec_emitter_collector_voltage_drop_V"),
        ("electron_cooling_flux_W_m2", "tec_electron_cooling_flux_W_m2"),
        ("electron_heating_flux_W_m2", "tec_electron_heating_flux_W_m2"),
        ("electron_cooling_power_W", "tec_electron_cooling_power_W"),
        ("electron_heating_power_W", "tec_electron_heating_power_W"),
        ("emitter_joule_power_axial_W", "tec_emitter_joule_power_axial_W"),
        ("collector_joule_power_axial_W", "tec_collector_joule_power_axial_W"),
    )
    for tfe_index, name in enumerate(payload["tfe_names"]):
        for axial_index in range(axial_count):
            rows.append({
                "time_s": payload["time_s"][0],
                "accident_elapsed_s": payload["accident_elapsed_s"][0],
                "tec_main_current_A": payload["tec_main_current_A"][0],
                "tec_main_voltage_V": payload["tec_main_voltage_V"][0],
                "tec_main_electric_power_W": payload["tec_main_electric_power_W"][0],
                "tec_main_converged": bool(payload["tec_main_converged"][0]),
                "tfe_name": name, "tfe_multiplier": payload["tfe_multipliers"][tfe_index],
                "axial_node_index": axial_index,
                **{
                    field: payload[key][tfe_index, axial_index]
                    for field, key in field_map
                },
            })
    _append_rows(out_dir / "history_electrical.csv", ELECTRICAL_HISTORY_FIELDS, rows)


def _append_reactivity_history(out_dir: Path, payload: Dict[str, Any]) -> None:
    row = {
        field: (
            payload["time_s"][0] if field == "time_s"
            else payload["accident_elapsed_s"][0] if field == "accident_elapsed_s"
            else payload[field][0]
        )
        for field in REACTIVITY_HISTORY_FIELDS
    }
    _append_rows(out_dir / "history_reactivity.csv", REACTIVITY_HISTORY_FIELDS, [row])


def append_postprocessing_histories(out_dir: Path, payload: Dict[str, Any]) -> None:
    _append_coolant_history(out_dir, payload)
    _append_solid_history(out_dir, payload)
    _append_electrical_history(out_dir, payload)
    _append_reactivity_history(out_dir, payload)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False), encoding="utf-8")


def run_loca(config: LocARunConfig) -> Dict[str, Any]:
    _validate_config(config)
    out_dir = Path(config.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "history.csv"
    if summary_path.exists():
        raise FileExistsError(f"output history already exists: {summary_path}")

    source = _read_source_config(Path(config.restart_in))
    runtime = ReactivityControlRunConfig(
        restart_in=Path(config.restart_in),
        output_dir=out_dir,
        duration_s=float(config.duration_s),
        dt_s=float(config.dt_s),
        record_interval_s=float(config.record_interval_s),
        checkpoint_interval_s=float(config.checkpoint_interval_s),
        min_fluid_temperature_stop_k=None,
        external_heat_enabled=bool(source.get("external_heat_enabled", False)),
        external_heat_period_s=float(source.get("external_heat_period_s", 5668.144369)),
        external_heat_time_origin_s=float(source.get("external_heat_time_origin_s", 0.0)),
    )
    debug, source = load_baseline_debug_config(runtime)
    build = build_debug_case(
        debug, apply_fixed_power=not bool(config.enable_reactivity_feedback),
    )
    system, core = build["system"], build["core"]
    handoff_type = "fixed_power"
    if config.enable_reactivity_feedback:
        handoff_type = prepare_reactivity_control(
            core,
            source_point_kinetics_enabled=bool(source["point_kinetics_enabled"]),
            expected_power_w=float(source["power_w"]),
        )
    set_tec_update_interval(core, float(config.tec_update_interval_s))
    start_time = float(system.global_time)

    pre_path = out_dir / "snapshot_pre_accident.npz"
    pre_feedback = core.compute_reactivity_feedback()
    placeholder = {
        "reference_fluid": {
            "T": np.asarray(system.fluid_solver.T_vec).copy(),
            "P": np.asarray(system.fluid_solver.P_vec).copy(),
            "h": np.asarray(system.fluid_solver.h_vec).copy(),
            "W": np.asarray(system.fluid_solver.W_vec).copy(),
        },
        "feedback_reference_total": float(pre_feedback.total),
    }
    write_snapshot(pre_path, build, placeholder, start_time=start_time, coolant_present=True)

    accident = take_over_accident(
        build, float(config.radiation_emissivity),
        fixed_power=not bool(config.enable_reactivity_feedback),
    )
    event = {
        "case": CASE_NAME,
        "accident_time_absolute_s": start_time,
        "accident_model": "instantaneous_complete_NaK_loss",
        "hydraulic_solve_enabled": False,
        "fluid_solid_couplers_removed": len(accident["removed_couplers"]),
        "core_vacuum_radiation_gaps_added": len(accident["vacuum_gaps"]),
        "radiation_emissivity": float(config.radiation_emissivity),
    }
    _write_json(out_dir / "accident_event.json", event)
    run_config = dict(source)
    run_config.update({
        **event,
        "restart_in": str(config.restart_in),
        "source_run_config": str(Path(config.restart_in).parent / "run_config.json"),
        "duration_s": float(config.duration_s),
        "dt_s": float(config.dt_s),
        "record_interval_s": float(config.record_interval_s),
        "checkpoint_interval_s": float(config.checkpoint_interval_s),
        "tec_update_interval_s": float(config.tec_update_interval_s),
        "power_w": 210000.0,
        "point_kinetics_enabled": bool(config.enable_reactivity_feedback),
        "reactivity_control_mode": (
            "temperature_feedback_plus_scram"
            if config.scram_time_s is not None else "temperature_feedback_only"
        ),
        "handoff_type": handoff_type,
        "scram_time_s": config.scram_time_s,
        "scram_reactivity_dollars": (
            float(config.scram_reactivity_dollars)
            if config.scram_time_s is not None else 0.0
        ),
        "scram_reactivity": (
            float(config.scram_reactivity_dollars) * float(core.point_reactor.beta_total)
            if config.scram_time_s is not None else 0.0
        ),
        "record_schedule_s": (
            [[0.0, 20.0, 0.5], [20.0, 100.0, 2.0], [100.0, 400.0, 5.0],
             [400.0, 600.0, 10.0], [600.0, float(config.duration_s), 20.0]]
            if config.staged_recording else [[0.0, float(config.duration_s), float(config.record_interval_s)]]
        ),
        "tec_open_circuit_after_scram": config.scram_time_s is not None,
        "tec_open_circuit_current_threshold_A": float(
            config.tec_open_circuit_current_threshold_a
        ),
        "coolant_state_after_accident": "absent",
        "ordinary_pipe_internal_boundary": "adiabatic",
        "ringhp_fluid_solid_heat_transfer": "disabled",
        "cross_event_energy_balance_required": False,
        "failure_temperature_limits_K": {
            "collector": float(config.collector_failure_temperature_k),
            "emitter": float(config.emitter_failure_temperature_k),
            "coolant": float(config.coolant_failure_temperature_k),
            "moderator": float(config.moderator_failure_temperature_k),
            "reflector": float(config.reflector_failure_temperature_k),
        },
    })
    _write_json(out_dir / "run_config.json", run_config)

    post_path = out_dir / "snapshot_tplus_00000.000s.npz"
    post_payload = write_snapshot(
        post_path, build, accident, start_time=start_time, coolant_present=False,
    )
    append_postprocessing_histories(out_dir, post_payload)
    previous_energy: Optional[tuple[float, float]] = None
    latest, previous_energy = collect_summary(
        build, accident, start_time=start_time, dt_s=0.0,
        snapshot_path=str(post_path), previous_energy=previous_energy,
    )
    _append_summary(summary_path, latest)

    end_time = start_time + float(config.duration_s)
    next_record = (
        staged_record_interval(0.0)
        if config.staged_recording else float(config.record_interval_s)
    )
    last_checkpoint = 0.0
    stop_reason = "completed"
    while float(system.global_time) < end_time - 1.0e-9:
        dt = min(float(config.dt_s), end_time - float(system.global_time))
        elapsed_before = float(system.global_time) - start_time
        external_reactivity, external_dollars = external_reactivity_for_elapsed(
            core, config, elapsed_before,
        )
        if not config.enable_reactivity_feedback:
            _apply_fixed_core_power(build, 210000.0)
        system.step(
            dt, inner_iter=1, fail_on_fluid_nonconvergence=True, fluid_max_iter=1,
            reactivity_control=external_reactivity,
        )
        if not config.enable_reactivity_feedback:
            _apply_fixed_core_power(build, 210000.0)
        elapsed = float(system.global_time) - start_time
        tec_opened = maybe_transition_tec_to_open_circuit(
            build, accident, config, elapsed,
        )
        if tec_opened:
            print(
                f"[LOCA +{elapsed:.3f}s] TEC open circuit: "
                f"I={accident['tec_open_circuit_trigger_current_A']:.6g} A",
                flush=True,
            )
        external_reactivity, external_dollars = external_reactivity_for_elapsed(
            core, config, elapsed,
        )
        failure_reason, _ = evaluate_failure_reason(system, config, coolant_present=False)
        if (
            elapsed + 1.0e-9 >= next_record
            or float(system.global_time) >= end_time - 1.0e-9
            or failure_reason
        ):
            snap = out_dir / f"snapshot_tplus_{elapsed:09.3f}s.npz"
            payload = write_snapshot(
                snap, build, accident, start_time=start_time, coolant_present=False,
                external_reactivity=external_reactivity,
                external_reactivity_dollars=external_dollars,
            )
            append_postprocessing_histories(out_dir, payload)
            latest, previous_energy = collect_summary(
                build, accident, start_time=start_time, dt_s=dt,
                snapshot_path=str(snap), previous_energy=previous_energy,
                failure_reason=failure_reason,
                external_reactivity=external_reactivity,
                external_reactivity_dollars=external_dollars,
            )
            _append_summary(summary_path, latest)
            print(
                f"[LOCA +{elapsed:.3f}s] Tmax={latest['core_structure_max_T_K']:.3f} K "
                f"I={latest['tec_main_current_A']:.3f} A "
                f"rho_change={latest['feedback_total_change_from_accident']:.6e}",
                flush=True,
            )
            while next_record <= elapsed + 1.0e-9:
                next_record += (
                    staged_record_interval(next_record)
                    if config.staged_recording else float(config.record_interval_s)
                )
        if failure_reason:
            stop_reason = failure_reason
            break
        if (
            float(config.checkpoint_interval_s) > 0.0
            and elapsed - last_checkpoint >= float(config.checkpoint_interval_s) - 1.0e-9
        ):
            system.save_global_state(str(out_dir / f"checkpoint_tplus_{elapsed:.3f}s.npz"))
            last_checkpoint = elapsed

    restart_path = out_dir / "stage_01_restart.npz"
    system.save_global_state(str(restart_path))
    result = {
        "case": CASE_NAME,
        "output_dir": str(out_dir),
        "source_restart_path": str(config.restart_in),
        "start_time_s": start_time,
        "end_time_s": float(system.global_time),
        "duration_s": float(system.global_time) - start_time,
        "history_path": str(summary_path),
        "postprocessing_history_paths": {
            name: str(out_dir / f"history_{name}.csv")
            for name in ("coolant", "solids", "electrical", "reactivity")
        },
        "restart_path": str(restart_path),
        "latest_metrics": latest,
        "stop_reason": stop_reason,
    }
    _write_json(out_dir / "run_summary.json", result)
    _write_json(out_dir / "latest_state.json", result)
    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-in", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--duration", type=float, default=0.4)
    parser.add_argument("--dt", type=float, default=0.05)
    parser.add_argument("--record-interval", type=float, default=DEFAULT_RECORD_INTERVAL_S)
    parser.add_argument("--checkpoint-interval", type=float, default=0.0)
    parser.add_argument("--tec-update-interval", type=float, default=0.05)
    parser.add_argument("--radiation-emissivity", type=float, default=DEFAULT_EMISSIVITY)
    parser.add_argument("--collector-failure-temperature", type=float, default=1500.0)
    parser.add_argument("--emitter-failure-temperature", type=float, default=3000.0)
    parser.add_argument("--coolant-failure-temperature", type=float, default=1058.0)
    parser.add_argument("--moderator-failure-temperature", type=float, default=930.0)
    parser.add_argument("--reflector-failure-temperature", type=float, default=1000.0)
    parser.add_argument("--enable-reactivity-feedback", action="store_true")
    parser.add_argument("--scram-time", type=float)
    parser.add_argument("--scram-reactivity-dollars", type=float, default=-2.0)
    parser.add_argument("--staged-recording", action="store_true")
    parser.add_argument("--tec-open-circuit-current-threshold", type=float, default=0.01)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    result = run_loca(LocARunConfig(
        restart_in=args.restart_in,
        output_dir=args.output_dir,
        duration_s=float(args.duration),
        dt_s=float(args.dt),
        record_interval_s=float(args.record_interval),
        checkpoint_interval_s=float(args.checkpoint_interval),
        tec_update_interval_s=float(args.tec_update_interval),
        radiation_emissivity=float(args.radiation_emissivity),
        collector_failure_temperature_k=float(args.collector_failure_temperature),
        emitter_failure_temperature_k=float(args.emitter_failure_temperature),
        coolant_failure_temperature_k=float(args.coolant_failure_temperature),
        moderator_failure_temperature_k=float(args.moderator_failure_temperature),
        reflector_failure_temperature_k=float(args.reflector_failure_temperature),
        enable_reactivity_feedback=bool(args.enable_reactivity_feedback),
        scram_time_s=args.scram_time,
        scram_reactivity_dollars=float(args.scram_reactivity_dollars),
        staged_recording=bool(args.staged_recording),
        tec_open_circuit_current_threshold_a=float(
            args.tec_open_circuit_current_threshold
        ),
    ))
    print(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
