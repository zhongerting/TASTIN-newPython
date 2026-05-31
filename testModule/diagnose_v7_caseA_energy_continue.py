import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)
os.chdir(root_dir)

from run_v7_caseA_multipliers_short import collect_tec_stats, parse_multipliers
from run_v7_caseA_newpyd_long200000 import advance_to, build_loaded_case
from test_core_assemble_v7_caseA import _case_a_electric_diagnostics
from test_core_assemble_v7_caseA_faststeady import compute_faststeady_energy_audit


DEFAULT_RESTART = (
    "testModule/v7_caseA_newpyd_long20000_from_t23800/"
    "v7_caseA_newpyd_long20000_from_t23800_latest_restart.npz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--restart-in", default=DEFAULT_RESTART)
    parser.add_argument("--output-csv", default=(
        "testModule/v7_caseA_newpyd_long20000_from_t23800/"
        "energy_diagnostic_continue.csv"
    ))
    parser.add_argument("--duration", type=float, default=1000.0)
    parser.add_argument("--record-interval", type=float, default=200.0)
    parser.add_argument("--max-dt", type=float, default=0.8)
    parser.add_argument("--safety-factor", type=float, default=5000.0)
    parser.add_argument("--pipe-n-nodes", type=int, default=8)
    parser.add_argument("--target-voltage", type=float, default=27.2)
    parser.add_argument("--ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,18"))
    parser.add_argument("--tec-ring-multipliers", type=parse_multipliers, default=parse_multipliers("1,6,12,15"))
    return parser.parse_args()


def applied_tec_heat_totals(build: Dict[str, Any]) -> Dict[str, float]:
    totals = {
        "applied_emitter_electron_heat_removed_positive_w": 0.0,
        "applied_collector_electron_heat_source_signed_w": 0.0,
        "applied_electron_boundary_heat_difference_w": 0.0,
        "applied_emitter_joule_heat_w": 0.0,
        "applied_collector_joule_heat_w": 0.0,
        "applied_total_joule_heat_w": 0.0,
    }
    for name, multiplier in build["tec_ring_multipliers"].items():
        tfe = build["tfes"][name]
        multiplier = float(multiplier)
        area = np.asarray(tfe.solids["emitter"].boundaries["right"].area, dtype=float)
        q_e = np.asarray(tfe.plasma_data.electron_cooling_flux, dtype=float)
        q_c = np.asarray(tfe.plasma_data.electron_heating_flux, dtype=float)
        emitter_removed = -float(np.sum(q_e * area)) * multiplier
        collector_source = float(np.sum(q_c * area)) * multiplier
        emitter_joule = float(np.sum(tfe.electric_data.emitter_joule_heat)) * multiplier
        collector_joule = float(np.sum(tfe.electric_data.collector_joule_heat)) * multiplier

        totals["applied_emitter_electron_heat_removed_positive_w"] += emitter_removed
        totals["applied_collector_electron_heat_source_signed_w"] += collector_source
        totals["applied_electron_boundary_heat_difference_w"] += emitter_removed - collector_source
        totals["applied_emitter_joule_heat_w"] += emitter_joule
        totals["applied_collector_joule_heat_w"] += collector_joule
        totals["applied_total_joule_heat_w"] += emitter_joule + collector_joule
    return totals


def scaled_solid_storage_w(build: Dict[str, Any]) -> float:
    total = 0.0
    ring_multipliers = build["ring_multipliers"]
    for name, solid in build["system"].solid_components.items():
        if hasattr(solid, "_update_properties"):
            solid._update_properties()
        storage = float(np.sum(np.asarray(solid.thermal_capacitance) * np.asarray(solid.dTdt)))
        multiplier = 1.0
        for ring_name, ring_multiplier in ring_multipliers.items():
            if name.startswith(f"{ring_name}_"):
                multiplier = float(ring_multiplier)
                break
        total += storage * multiplier
    return total


def collect_diagnostic_row(build: Dict[str, Any], start_time: float) -> Dict[str, float]:
    system = build["system"]
    core = build["core"]
    energy = compute_faststeady_energy_audit({"build": build, "system": system, "core": core})
    electric = _case_a_electric_diagnostics(core)
    tec_totals = collect_tec_stats(build).get("totals", {})
    applied = applied_tec_heat_totals(build)

    q_core = float(energy["core_heat_power_w"])
    q_cool = float(energy["coolant_heat_pickup_w"])
    q_rad = float(energy["outer_wall_radiation_w"])
    p_elec = float(electric["tec_total_electric_power_w"])
    instant_boundary = float(tec_totals["electron_boundary_heat_difference_w"])
    instant_joule = float(tec_totals["total_joule_heat_w"])
    applied_boundary = float(applied["applied_electron_boundary_heat_difference_w"])
    applied_joule = float(applied["applied_total_joule_heat_w"])
    instant_tec_removed = instant_boundary - instant_joule
    applied_tec_removed = applied_boundary - applied_joule

    return {
        "absolute_time_s": float(system.global_time),
        "relative_time_s": float(system.global_time) - start_time,
        "core_heat_power_w": q_core,
        "coolant_heat_pickup_w": q_cool,
        "outer_wall_radiation_w": q_rad,
        "terminal_electric_power_w": p_elec,
        "instant_electron_boundary_diff_w": instant_boundary,
        "instant_total_joule_w": instant_joule,
        "instant_tec_heat_removed_minus_joule_w": instant_tec_removed,
        "applied_electron_boundary_diff_w": applied_boundary,
        "applied_total_joule_w": applied_joule,
        "applied_tec_heat_removed_minus_joule_w": applied_tec_removed,
        "terminal_minus_instant_tec_heat_w": p_elec - instant_tec_removed,
        "terminal_minus_applied_tec_heat_w": p_elec - applied_tec_removed,
        "residual_using_terminal_power_w": q_core - q_cool - q_rad - p_elec,
        "residual_using_instant_tec_heat_w": q_core - q_cool - q_rad - instant_tec_removed,
        "residual_using_applied_tec_heat_w": q_core - q_cool - q_rad - applied_tec_removed,
        "scaled_solid_storage_w": scaled_solid_storage_w(build),
    }


def main() -> None:
    args = parse_args()
    build = build_loaded_case(args)
    system = build["system"]
    start_time = float(system.global_time)
    target_time = start_time + float(args.duration)
    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    next_record_time = start_time + float(args.record_interval)
    while next_record_time <= target_time + 1.0e-10:
        advance_to(build, next_record_time, args)
        row = collect_diagnostic_row(build, start_time)
        rows.append(row)
        print(
            f"t={row['absolute_time_s']:.1f}s "
            f"Pterm={row['terminal_electric_power_w']:.3f}W "
            f"Qtec_app={row['applied_tec_heat_removed_minus_joule_w']:.3f}W "
            f"Rterm={row['residual_using_terminal_power_w']:.3f}W "
            f"Rapp={row['residual_using_applied_tec_heat_w']:.3f}W",
            flush=True,
        )
        next_record_time += float(args.record_interval)

    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"CSV: {output_csv}", flush=True)


if __name__ == "__main__":
    main()
