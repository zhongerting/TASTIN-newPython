import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict

import numpy as np


CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
os.chdir(ROOT_DIR)

from run_v11_caseA_closed_loop import build_case, parse_args as parse_v11_args  # noqa: E402
from ThermoCalc.ThermoCalcWrapper import load_emission_lookup_database, te_solver  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Probe V11 ThermoCalc lookup hit rate.")
    parser.add_argument("--restart-in", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--lookup-db", required=True)
    parser.add_argument("--lookup-regions", default="core,startup,high_power,accident")
    parser.add_argument("--pump-total-head-pa", type=float, default=6483.548313292204)
    parser.add_argument("--ring-emissivity", type=float, default=0.24)
    parser.add_argument("--outer-header-emissivity", type=float, default=0.30)
    return parser.parse_args()


def _summary(values: np.ndarray) -> Dict[str, float]:
    values = np.asarray(values, dtype=float).ravel()
    return {
        "min": float(np.nanmin(values)),
        "mean": float(np.nanmean(values)),
        "max": float(np.nanmax(values)),
    }


def main() -> None:
    args = parse_args()
    regions = tuple(part.strip() for part in args.lookup_regions.split(",") if part.strip())
    os.environ["THERMOCALC_ENABLE_LOOKUP"] = "1"
    os.environ["THERMOCALC_LOOKUP_DB"] = str(Path(args.lookup_db).resolve())
    os.environ["THERMOCALC_LOOKUP_REGIONS"] = ",".join(regions)

    loaded = load_emission_lookup_database(args.lookup_db, enable=True, force=True, regions=regions)
    dense_regions = int(te_solver.emission_lookup_dense_region_count())

    old_argv = sys.argv[:]
    sys.argv = [
        "run_v11_caseA_closed_loop.py",
        "--restart-in",
        args.restart_in,
        "--duration",
        "0",
        "--record-interval",
        "50",
        "--restart-interval",
        "50",
        "--max-dt",
        "0.05",
        "--ring-emissivity",
        str(args.ring_emissivity),
        "--outer-header-emissivity",
        str(args.outer_header_emissivity),
        "--pump-total-head-pa",
        str(args.pump_total_head_pa),
        "--enable-pump-head-control",
        "--pump-control-interval",
        "50",
    ]
    try:
        v11_args = parse_v11_args()
    finally:
        sys.argv = old_argv

    build = build_case(v11_args)
    thermo = build["core"].thermo_calc
    if thermo is None:
        raise RuntimeError("V11 build did not create ThermoCalc model.")

    input_data = thermo._input_data
    te = np.asarray(input_data.Temitter, dtype=float).ravel()
    tc = np.asarray(input_data.Tcollector, dtype=float).ravel()
    tcs = np.asarray(input_data.Tcs, dtype=float).ravel()
    d_gap = np.asarray(input_data.d_gap, dtype=float).ravel()
    if d_gap.size == 0:
        d_gap_value = 0.5
    else:
        d_gap_value = float(d_gap[0])

    # The local emission voltage is not directly stored as a flattened lookup
    # input before circuit solve. For coverage probing, current V11 fixed-U
    # states use local Vo values in the same operating band seen after build.
    tic = time.perf_counter()
    returned_calculate_ms = float(thermo.calculate(verbose=False))
    calculate_wall_s = time.perf_counter() - tic

    vo_values = []
    for tec in thermo._circuit.TECs:
        vo_values.extend(float(value) for value in tec.V)
    vo = np.asarray(vo_values, dtype=float)
    if vo.size != te.size:
        vo = np.full_like(te, float(v11_args.target_voltage) / max(1, thermo.N_elem))

    query = te_solver.lookup_emission_points(te, tc, vo, tcs, d_gap_value)
    found = np.asarray(query["found"], dtype=np.uint8)
    n_points = int(found.size)
    n_found = int(found.sum())
    result: Dict[str, Any] = {
        "lookup_db": str(Path(args.lookup_db).resolve()),
        "lookup_regions": list(regions),
        "loaded_count": int(loaded),
        "dense_region_count": dense_regions,
        "restart_in": args.restart_in,
        "absolute_time_s": float(build["system"].global_time),
        "n_points": n_points,
        "found": n_found,
        "miss": int(n_points - n_found),
        "hit_frac": float(n_found / n_points) if n_points else None,
        "calculate_wall_s": float(calculate_wall_s),
        "returned_calculate_ms": returned_calculate_ms,
        "TE_K": _summary(te),
        "TC_K": _summary(tc),
        "Tcs_K": _summary(tcs),
        "Vo_V": _summary(vo),
        "d_gap_mm": d_gap_value,
    }
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

