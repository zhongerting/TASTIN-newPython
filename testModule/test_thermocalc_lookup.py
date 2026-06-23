import os
import sys
import time
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
BUILD_PYD_DIR = ROOT_DIR / "ThermoCalc" / "build_cp312" / "Release"
if str(BUILD_PYD_DIR) not in sys.path:
    sys.path.insert(0, str(BUILD_PYD_DIR))

import te_solver  # noqa: E402


def _load_block(path: Path, *, name: str, priority: int):
    with np.load(path, allow_pickle=False) as d:
        if "lookup_safe_flag" in d.files:
            safe = np.asarray(d["lookup_safe_flag"], dtype=np.uint8)
        else:
            safe = np.asarray(d["done"] & d["finite_flag"] & d["converged"], dtype=np.uint8)
        te_solver.add_emission_lookup_block(
            name,
            priority,
            np.asarray(d["TE_axis"], dtype=np.float64),
            np.asarray(d["TC_axis"], dtype=np.float64),
            np.asarray(d["Vo_axis"], dtype=np.float64),
            np.asarray(d["Tcs_axis"], dtype=np.float64),
            np.asarray(d["J"], dtype=np.float64),
            np.asarray(d["Vd"], dtype=np.float64),
            np.asarray(d["delta_V"], dtype=np.float64),
            np.asarray(d["phiE"], dtype=np.float64),
            np.asarray(d["phiC"], dtype=np.float64),
            safe,
        )
        return {key: d[key].copy() for key in d.files}


def _sample_inside_axes(data, n: int, seed: int = 622):
    rng = np.random.default_rng(seed)
    axes = [
        np.asarray(data["TE_axis"], dtype=np.float64),
        np.asarray(data["TC_axis"], dtype=np.float64),
        np.asarray(data["Vo_axis"], dtype=np.float64),
        np.asarray(data["Tcs_axis"], dtype=np.float64),
    ]
    return [rng.uniform(float(axis[0]), float(axis[-1]), size=n) for axis in axes]


def test_lookup_exact_grid_and_calc_path():
    te_solver.clear_emission_lookup()
    te_solver.set_emission_lookup_enabled(True)
    core_path = ROOT_DIR / "ThermoCalc" / "emission_database" / "chunks" / "core" / "core_0000.npz"
    data = _load_block(core_path, name="core_0000", priority=0)
    assert te_solver.emission_lookup_block_count() == 1

    idx = (1, 20, 30, 10)
    te = float(data["TE_axis"][idx[0]])
    tc = float(data["TC_axis"][idx[1]])
    vo = float(data["Vo_axis"][idx[2]])
    tcs = float(data["Tcs_axis"][idx[3]])
    lookup = te_solver.lookup_emission_point(te, tc, vo, tcs, 0.5)
    assert lookup["found"]
    for field in ("J", "Vd", "delta_V", "phiE", "phiC"):
        assert abs(float(lookup[field]) - float(data[field][idx])) <= 1.0e-12

    production_lookup = te_solver.calc_emission_point_production(te, tc, vo, tcs, 0.5)
    for field in ("J", "Vd", "delta_V", "phiE", "phiC"):
        assert abs(float(production_lookup[field]) - float(data[field][idx])) <= 1.0e-12

    te_solver.set_emission_lookup_enabled(False)
    production_analytic = te_solver.calc_emission_point_production(te, tc, vo, tcs, 0.5)
    assert abs(float(production_analytic["J"]) - float(data["J"][idx])) <= 1.0e-8


def test_optimized_accident_block_is_queryable():
    te_solver.clear_emission_lookup()
    te_solver.set_emission_lookup_enabled(True)
    path = ROOT_DIR / "ThermoCalc" / "emission_database" / "chunks" / "accident" / "accident_0010.optimized.npz"
    data = _load_block(path, name="accident_0010_optimized", priority=3)
    assert int(np.count_nonzero(data["unresolved_flag"])) == 0
    assert int(np.count_nonzero(data["lookup_safe_flag"])) == int(data["lookup_safe_flag"].size)
    assert int(np.count_nonzero(~np.isfinite(data["J"]))) == 0
    assert float(np.min(data["J"])) >= 0.0

    idx = np.argwhere(data["imputed_flag"])
    if idx.size == 0:
        idx = np.argwhere(data["zero_fill_flag"])
    assert idx.size > 0
    i, j, k, l = (int(x) for x in idx[0])
    lookup = te_solver.lookup_emission_point(
        float(data["TE_axis"][i]),
        float(data["TC_axis"][j]),
        float(data["Vo_axis"][k]),
        float(data["Tcs_axis"][l]),
        0.5,
    )
    assert lookup["found"]
    assert abs(float(lookup["J"]) - float(data["J"][i, j, k, l])) <= 1.0e-12


def benchmark_lookup_vs_analytic(n_points: int = 20000):
    te_solver.clear_emission_lookup()
    te_solver.set_emission_lookup_enabled(True)
    core_path = ROOT_DIR / "ThermoCalc" / "emission_database" / "chunks" / "core" / "core_0000.npz"
    data = _load_block(core_path, name="core_0000", priority=0)
    te, tc, vo, tcs = _sample_inside_axes(data, n_points)

    # Warm up both paths.
    te_solver.lookup_emission_points(te[:100], tc[:100], vo[:100], tcs[:100], 0.5)
    te_solver.set_emission_lookup_enabled(False)
    for i in range(20):
        te_solver.calc_emission_point_production(float(te[i]), float(tc[i]), float(vo[i]), float(tcs[i]), 0.5)

    te_solver.set_emission_lookup_enabled(True)
    t0 = time.perf_counter()
    result = te_solver.lookup_emission_points(te, tc, vo, tcs, 0.5)
    lookup_elapsed = time.perf_counter() - t0
    assert int(np.count_nonzero(result["found"])) == n_points

    te_solver.set_emission_lookup_enabled(False)
    t0 = time.perf_counter()
    for i in range(n_points):
        te_solver.calc_emission_point_production(float(te[i]), float(tc[i]), float(vo[i]), float(tcs[i]), 0.5)
    analytic_elapsed = time.perf_counter() - t0

    return {
        "points": n_points,
        "lookup_elapsed_s": lookup_elapsed,
        "analytic_elapsed_s": analytic_elapsed,
        "lookup_points_per_s": n_points / lookup_elapsed,
        "analytic_points_per_s": n_points / analytic_elapsed,
        "speedup": analytic_elapsed / lookup_elapsed,
    }


if __name__ == "__main__":
    test_lookup_exact_grid_and_calc_path()
    test_optimized_accident_block_is_queryable()
    stats = benchmark_lookup_vs_analytic()
    print("ThermoCalc lookup checks passed.")
    for key, value in stats.items():
        print(f"{key}: {value}")

