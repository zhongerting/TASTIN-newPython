"""Evaluate a saved V14 temperature field with the series TEC circuit."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.V14_210kW_start.run_v14_startup_stages import (
    CESIUM_GAP_H_W_M2_K,
    _apply_wire_resistance_without_calculate,
    build_case,
)
from testModule.v13_startup_control import apply_tec_gap_h_eq


CASE_DIR = Path(__file__).resolve().parent
DEFAULT_BASE = CASE_DIR / "startup_5000s_final"
DEFAULT_RESTART = (
    DEFAULT_BASE
    / "stage_2_ramp_60kw_to_210kw_hold_2800s"
    / "final_restart.npz"
)
DEFAULT_OUTPUT = DEFAULT_BASE / "t2800_fixed_i_200a_probe"


def run_probe(
    restart: Path,
    output_dir: Path,
    current_a: float,
    expected_time_s: float = 2800.0,
    mode: str = "fixed_i",
    resistance_ohm: float = 0.003,
) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    if not math.isfinite(current_a) or current_a < 0.0:
        raise ValueError("current_a must be finite and non-negative")
    if mode not in ("fixed_i", "fixed_r"):
        raise ValueError("mode must be fixed_i or fixed_r")
    if not math.isfinite(resistance_ohm) or resistance_ohm <= 0.0:
        raise ValueError("resistance_ohm must be finite and positive")

    build = build_case()
    system = build["system"]
    core = build["core"]
    system.load_global_state(str(restart))
    if not math.isclose(
        float(system.global_time), float(expected_time_s), abs_tol=2.0e-6
    ):
        raise ValueError(
            f"Expected t={expected_time_s:g} s restart, got {system.global_time}"
        )

    updated = apply_tec_gap_h_eq(core, CESIUM_GAP_H_W_M2_K)
    if updated != len(build["tfes"]):
        raise RuntimeError(f"Updated {updated} TEC gaps; expected {len(build['tfes'])}")
    target = current_a if mode == "fixed_i" else resistance_ohm
    core.setup_tec_circuit(mode, target, I_guess=current_a, topology="series")
    _apply_wire_resistance_without_calculate(core)
    core.enable_tec_coupled = True

    group = core.tec_circuit_groups["main"]
    core._sync_tec_group_temperatures(group)
    thermo = group.thermo_calc
    emitter_k = np.asarray(thermo._T_emitter, dtype=float).copy()
    collector_k = np.asarray(thermo._T_collector, dtype=float).copy()
    elapsed_ms = float(thermo.calculate(verbose=False))
    global_result = thermo.get_global_results()

    names = []
    copy_indices = []
    for name, multiplier in group.multipliers.items():
        for copy_index in range(int(multiplier)):
            names.append(name)
            copy_indices.append(copy_index)
    if len(names) != thermo.N_elem:
        raise RuntimeError("TEC multiplier expansion does not match circuit elements")

    fields = {
        key: np.stack([
            np.asarray(thermo.get_tec_results(i)[key], dtype=float)
            for i in range(thermo.N_elem)
        ])
        for key in (
            "J", "V", "UE", "UC", "phiE", "phiC", "Vd",
            "joulePowerE", "joulePowerC",
        )
    }
    element_current = np.asarray([
        float(thermo.get_tec_results(i)["I"]) for i in range(thermo.N_elem)
    ])
    element_voltage = np.asarray([
        float(thermo.get_tec_results(i)["U"]) for i in range(thermo.N_elem)
    ])

    output_dir.mkdir(parents=True)
    np.savez_compressed(
        output_dir / "result.npz",
        emitter_temperature_k=emitter_k,
        collector_temperature_k=collector_k,
        representative_name=np.asarray(names),
        copy_index=np.asarray(copy_indices, dtype=int),
        element_current_a=element_current,
        element_voltage_v=element_voltage,
        **fields,
    )

    with (output_dir / "result.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "virtual_element", "representative_name", "copy_index", "axial_node",
            "emitter_temperature_K", "collector_temperature_K", "element_current_A",
            "element_voltage_V", "current_density_A_m2", "gap_voltage_V",
            "emitter_joule_power_W", "collector_joule_power_W",
        ))
        for i, name in enumerate(names):
            for node in range(thermo.n_node):
                writer.writerow((
                    i, name, copy_indices[i], node, emitter_k[i, node],
                    collector_k[i, node], element_current[i], element_voltage[i],
                    fields["J"][i, node], fields["V"][i, node],
                    fields["joulePowerE"][i, node],
                    fields["joulePowerC"][i, node],
                ))

    summary = {
        "restart": str(restart.resolve()),
        "time_s": float(system.global_time),
        "mode": mode,
        "requested_current_a": float(current_a),
        "resistance_ohm": float(resistance_ohm) if mode == "fixed_r" else None,
        "series_element_count": int(thermo.N_elem),
        "gap_h_eq_w_m2_k": CESIUM_GAP_H_W_M2_K,
        "calculation_elapsed_ms": elapsed_ms,
        "converged": bool(global_result["converged"]),
        "current_a": float(global_result["Iout"]),
        "voltage_v": float(global_result["Uout"]),
        "electric_power_w": float(global_result["Iout"] * global_result["Uout"]),
        "iteration_count": int(global_result["iteration_count"]),
        "zero_emission_skipped": bool(global_result["zero_emission_skipped"]),
        "zero_emission_reason": global_result["zero_emission_reason"],
        "emitter_temperature_min_k": float(np.min(emitter_k)),
        "emitter_temperature_max_k": float(np.max(emitter_k)),
        "collector_temperature_min_k": float(np.min(collector_k)),
        "collector_temperature_max_k": float(np.max(collector_k)),
        "finite_outputs": bool(all(
            np.all(np.isfinite(value)) for value in (
                element_current, element_voltage, *fields.values()
            )
        )),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart", type=Path, default=DEFAULT_RESTART)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--current-a", type=float, default=200.0)
    parser.add_argument("--expected-time-s", type=float, default=2800.0)
    parser.add_argument("--mode", choices=("fixed_i", "fixed_r"), default="fixed_i")
    parser.add_argument("--resistance-ohm", type=float, default=0.003)
    args = parser.parse_args()
    run_probe(
        args.restart,
        args.output_dir,
        args.current_a,
        expected_time_s=args.expected_time_s,
        mode=args.mode,
        resistance_ohm=args.resistance_ohm,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
