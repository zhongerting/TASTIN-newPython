import os
import sys
import time
import tempfile
import contextlib
import io
from argparse import Namespace
from pathlib import Path

import numpy as np


ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
BUILD_PYD_DIR = ROOT_DIR / "ThermoCalc" / "build_cp312" / "Release"
if str(BUILD_PYD_DIR) not in sys.path:
    sys.path.insert(0, str(BUILD_PYD_DIR))
TOOLS_DIR = ROOT_DIR / "ThermoCalc" / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import te_solver  # noqa: E402
from emission_database import cmd_export_runtime  # noqa: E402
from ThermoCalc.ThermoCalcWrapper import load_emission_lookup_database  # noqa: E402


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
        assert abs(float(lookup[field]) - float(data[field][idx])) <= 5.0e-6

    production_lookup = te_solver.calc_emission_point_production(te, tc, vo, tcs, 0.5)
    for field in ("J", "Vd", "delta_V", "phiE", "phiC"):
        assert abs(float(production_lookup[field]) - float(data[field][idx])) <= 5.0e-6

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
    assert abs(float(lookup["J"]) - float(data["J"][i, j, k, l])) <= 5.0e-6


def test_runtime_export_and_wrapper_loader():
    db_dir = ROOT_DIR / "ThermoCalc" / "emission_database"
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = Path(tmp) / "runtime"
        with contextlib.redirect_stdout(io.StringIO()):
            cmd_export_runtime(
                Namespace(
                    db_dir=db_dir,
                    out_dir=out_dir,
                    dtype="float32",
                    region=["core"],
                    limit_chunks=1,
                    zero_j_threshold=1.0e-3,
                    zero_compress=True,
                )
            )
        manifest_path = out_dir / "runtime_manifest.json"
        assert manifest_path.exists()
        runtime_path = out_dir / "core" / "core_0000.runtime.npz"
        assert runtime_path.exists()
        with np.load(runtime_path, allow_pickle=False) as data:
            assert data["J"].dtype == np.float32
            assert data["phiE"].dtype == np.float32
            assert data["phiC"].dtype == np.float32
            assert data["lookup_safe"].dtype == np.uint8
            assert data["zero_mask"].dtype == np.uint8
            idx = (1, 20, 30, 10)
            te_axis = np.asarray(data["TE_axis"], dtype=np.float64).copy()
            sample = {
                "TE": float(data["TE_axis"][idx[0]]),
                "TC": float(data["TC_axis"][idx[1]]),
                "Vo": float(data["Vo_axis"][idx[2]]),
                "Tcs": float(data["Tcs_axis"][idx[3]]),
                "J": float(data["J"][idx]),
                "Vd": float(data["Vd"][idx]),
                "delta_V": float(data["delta_V"][idx]),
                "phiE": float(data["phiE"][idx]),
                "phiC": float(data["phiC"][idx]),
            }

        te_solver.clear_emission_lookup()
        loaded = load_emission_lookup_database(str(out_dir), enable=True, force=True, regions=("core",))
        assert loaded == 1
        assert te_solver.emission_lookup_region_count() == 1
        lookup = te_solver.lookup_emission_point(sample["TE"], sample["TC"], sample["Vo"], sample["Tcs"], 0.5)
        assert lookup["found"]
        for field in ("J", "Vd", "delta_V", "phiE", "phiC"):
            assert abs(float(lookup[field]) - sample[field]) <= 5.0e-6

        te_gap = 0.5 * (
            float(te_axis[1])
            + float(te_axis[2])
        )
        gap_single = te_solver.lookup_emission_point(te_gap, sample["TC"], sample["Vo"], sample["Tcs"], 0.5)
        assert gap_single["found"]
        batch = te_solver.lookup_emission_points(
            np.asarray([sample["TE"], te_gap], dtype=np.float64),
            np.asarray([sample["TC"], sample["TC"]], dtype=np.float64),
            np.asarray([sample["Vo"], sample["Vo"]], dtype=np.float64),
            np.asarray([sample["Tcs"], sample["Tcs"]], dtype=np.float64),
            0.5,
        )
        assert tuple(batch["J"].strides) == (8,)
        assert tuple(batch["found"].strides) == (1,)
        assert int(np.count_nonzero(batch["found"])) == 2
        assert abs(float(batch["J"][1]) - float(gap_single["J"])) <= 5.0e-6


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
    test_runtime_export_and_wrapper_loader()
    stats = benchmark_lookup_vs_analytic()
    print("ThermoCalc lookup checks passed.")
    for key, value in stats.items():
        print(f"{key}: {value}")
