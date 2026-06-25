import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from run_v7_caseA_multipliers_short import json_default
from Solvers.HeatConduction.HeatConduction import BaseHeatConduction
from test_core_assemble_v8_caseA import (
    _case_a_reset_design_flows_after_restart,
    build_v8_case_a_system,
)


TOTAL_POWER_W = 115000.0
DEFAULT_SOLID_ODE_METHOD = "LSODA"
DEFAULT_COOLANT_MATERIAL = "SodiumPotassium78"
WIRE_RESISTANCE_OHM = [
    0.00155199999999970,
    0.00102400000000000,
    0.000336000000000000,
    0.000608000000000000,
]


def parse_solid_ode_method(method: str) -> str:
    method = str(method)
    method_by_lower = {
        valid.lower(): valid
        for valid in BaseHeatConduction.VALID_ODE_METHODS
    }
    method = method_by_lower.get(method.lower(), method)
    if method not in BaseHeatConduction.VALID_ODE_METHODS:
        valid = ", ".join(sorted(BaseHeatConduction.VALID_ODE_METHODS))
        raise argparse.ArgumentTypeError(
            f"Unsupported solid ODE method '{method}'. Valid methods: {valid}."
        )
    return method


def parse_v8_multipliers(text: str, *, allow_zero: bool = False) -> List[int]:
    values = [int(part.strip()) for part in str(text).split(",") if part.strip()]
    if len(values) != 5:
        raise argparse.ArgumentTypeError("Use five comma-separated integers, e.g. 1,6,12,15,3.")
    if allow_zero:
        invalid = any(value < 0 for value in values)
        message = "Multipliers must be non-negative."
    else:
        invalid = any(value <= 0 for value in values)
        message = "Multipliers must be positive."
    if invalid:
        raise argparse.ArgumentTypeError(message)
    return values


def passive_tec_source_totals(build: Dict[str, Any]) -> Dict[str, float]:
    totals = {
        "electron_cooling_flux_w_m2_sum": 0.0,
        "electron_heating_flux_w_m2_sum": 0.0,
        "emitter_joule_heat_w": 0.0,
        "collector_joule_heat_w": 0.0,
        "coupler_emitter_source_w": 0.0,
        "coupler_collector_source_w": 0.0,
    }
    for name in build.get("passive_tfe_names", []):
        tfe = build["tfes"][name]
        totals["electron_cooling_flux_w_m2_sum"] += float(np.sum(np.abs(tfe.plasma_data.electron_cooling_flux)))
        totals["electron_heating_flux_w_m2_sum"] += float(np.sum(np.abs(tfe.plasma_data.electron_heating_flux)))
        totals["emitter_joule_heat_w"] += float(np.sum(np.abs(tfe.electric_data.emitter_joule_heat)))
        totals["collector_joule_heat_w"] += float(np.sum(np.abs(tfe.electric_data.collector_joule_heat)))
        totals["coupler_emitter_source_w"] += float(np.sum(np.abs(tfe.couplers["tec_couple"].Q_source_1)))
        totals["coupler_collector_source_w"] += float(np.sum(np.abs(tfe.couplers["tec_couple"].Q_source_2)))
    return totals


def apply_solid_ode_method(build: Dict[str, Any], method: str) -> Dict[str, str]:
    method = parse_solid_ode_method(method)
    system = build["system"]
    applied = {}
    for name, solid in system.solid_components.items():
        solid.set_ode_method(method)
        applied[name] = solid.ode_method
    return applied


def get_solid_ode_methods(build: Dict[str, Any]) -> Dict[str, str]:
    system = build["system"]
    return {
        name: str(getattr(solid, "ode_method", ""))
        for name, solid in system.solid_components.items()
    }


def apply_wire_resistance(core: Any, scale: float = 1.0) -> None:
    thermos = list(_iter_core_thermo_calcs(core))
    if not thermos:
        return
    wire_res = np.asarray(WIRE_RESISTANCE_OHM, dtype=float) * float(scale)
    for _, thermo_calc in thermos:
        n_elem = thermo_calc.N_elem
        thermo_calc._input_data.resistanceWire = np.tile(wire_res, (n_elem, 1))
        thermo_calc.build()
        thermo_calc.calculate(verbose=False)


def get_wire_resistance(core: Any) -> List[float]:
    if core.thermo_calc is None:
        return list(WIRE_RESISTANCE_OHM)
    resistance = np.asarray(core.thermo_calc._input_data.resistanceWire, dtype=float)
    if resistance.size == 0:
        return list(WIRE_RESISTANCE_OHM)
    return [float(value) for value in resistance.reshape((-1, 4))[0]]


def _iter_core_thermo_calcs(core: Any):
    if hasattr(core, "iter_tec_circuit_groups"):
        for group in core.iter_tec_circuit_groups():
            if group.thermo_calc is not None:
                yield group.name, group.thermo_calc
        return
    thermo_calc = getattr(core, "thermo_calc", None)
    if thermo_calc is not None:
        yield "main", thermo_calc


def load_tec_load_curve(path: Optional[str]) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Load a two-column external TEC load curve as U_load=f(I_total)."""
    if path is None:
        return None

    curve_path = Path(path)
    if not curve_path.exists():
        raise FileNotFoundError(f"TEC load curve file does not exist: {path}")

    if curve_path.suffix.lower() == ".npz":
        with np.load(curve_path, allow_pickle=False) as data:
            current_key = "current_a" if "current_a" in data else "current"
            voltage_key = "voltage_v" if "voltage_v" in data else "voltage"
            if current_key not in data or voltage_key not in data:
                raise ValueError("TEC load curve npz must contain current_a/voltage_v or current/voltage.")
            return np.asarray(data[current_key], dtype=float), np.asarray(data[voltage_key], dtype=float)

    current: List[float] = []
    voltage: List[float] = []
    with curve_path.open("r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(2048)
        f.seek(0)
        has_header = csv.Sniffer().has_header(sample) if sample.strip() else False
        if has_header:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError("TEC load curve csv has no header.")
            lower_names = {name.lower(): name for name in reader.fieldnames}
            current_name = lower_names.get("current_a") or lower_names.get("current")
            voltage_name = lower_names.get("voltage_v") or lower_names.get("voltage")
            if current_name is None or voltage_name is None:
                raise ValueError("TEC load curve csv header must contain current_a/voltage_v or current/voltage.")
            for row in reader:
                current.append(float(row[current_name]))
                voltage.append(float(row[voltage_name]))
        else:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 2:
                    continue
                current.append(float(row[0]))
                voltage.append(float(row[1]))

    return np.asarray(current, dtype=float), np.asarray(voltage, dtype=float)


def build_loaded_case(args: argparse.Namespace) -> Dict[str, Any]:
    solid_ode_method = parse_solid_ode_method(
        getattr(args, "solid_ode_method", DEFAULT_SOLID_ODE_METHOD)
    )
    build = build_v8_case_a_system(
        pipe_n_nodes=args.pipe_n_nodes,
        solid_heat_capacity_scale=1.0,
        solid_heat_capacity_scale_scope="global_outer",
        coolant_material=getattr(args, "coolant_material", DEFAULT_COOLANT_MATERIAL),
        ring_multipliers=args.ring_multipliers,
        tec_ring_multipliers=args.tec_ring_multipliers,
    )
    system = build["system"]
    core = build["core"]
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    apply_solid_ode_method(build, solid_ode_method)
    system.initialize_system()
    system.load_global_state(args.restart_in)
    apply_solid_ode_method(build, solid_ode_method)
    core.point_reactor = None
    core.enable_tec_coupled = True
    core.thermo_update_interval = 0.0
    core.setup_tec_circuit("fixed_u", args.target_voltage, I_guess=150.0)
    core.update_neutronic_power(
        p_total=TOTAL_POWER_W,
        p_fiss=TOTAL_POWER_W,
        p_decay=0.0,
        alpha=1.0,
    )
    _case_a_reset_design_flows_after_restart(build)
    core.post_step(0.0, float(system.global_time))
    if core.thermo_calc is not None:
        apply_wire_resistance(core)
    core.pre_step(0.0, float(system.global_time))
    build["solid_ode_method"] = solid_ode_method
    build["solid_ode_methods"] = get_solid_ode_methods(build)
    build["wire_resistance_ohm"] = get_wire_resistance(core)
    build["coolant_material"] = build.get("coolant_material", DEFAULT_COOLANT_MATERIAL)
    return build


__all__ = [
    "DEFAULT_COOLANT_MATERIAL",
    "DEFAULT_SOLID_ODE_METHOD",
    "TOTAL_POWER_W",
    "WIRE_RESISTANCE_OHM",
    "apply_wire_resistance",
    "apply_solid_ode_method",
    "build_loaded_case",
    "get_solid_ode_methods",
    "get_wire_resistance",
    "json_default",
    "load_tec_load_curve",
    "parse_v8_multipliers",
    "parse_solid_ode_method",
    "passive_tec_source_totals",
]
