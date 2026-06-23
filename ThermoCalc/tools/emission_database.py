from __future__ import annotations

import argparse
import csv
import hashlib
import json
import struct
import subprocess
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from scan_emission_map import (
    BOOL_FIELDS,
    FLOAT_FIELDS,
    INT_FIELDS,
    cesium_pressure_from_tcs,
    import_te_solver,
    pressure_axis_from_args,
    save_npz,
    summarize,
    tcs_from_cesium_pressure,
)


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
DEFAULT_DB_DIR = ROOT / "emission_database"
DATABASE_VERSION = "emission-db-v1"
SOURCE_REGIONS = {
    "core": 0,
    "startup": 1,
    "high_power": 2,
    "accident": 3,
}
DEFAULT_PCS_AXIS = (0.02, 3.0, 61, "log")

REGIONS: dict[str, dict[str, Any]] = {
    "core": {
        "description": "normal high-density production region",
        "priority": 0,
        "te": (1300.0, 2150.0, 86, True),
        "tc": (700.0, 900.0, 41, True),
        "vo": (0.0, 3.5, 71, True),
        "pcs": DEFAULT_PCS_AXIS,
        "target_points_per_chunk": 250_000,
        "worker_group": "A",
    },
    "startup": {
        "description": "startup and low-temperature region",
        "priority": 1,
        "te": (700.0, 1300.0, 31, False),
        "tc": (500.0, 800.0, 31, True),
        "vo": (0.0, 3.5, 36, True),
        "pcs": DEFAULT_PCS_AXIS,
        "target_points_per_chunk": 250_000,
        "worker_group": "B",
    },
    "high_power": {
        "description": "high-power hot extension region",
        "priority": 2,
        "te": (2150.0, 2400.0, 26, True),
        "tc": (750.0, 1000.0, 26, True),
        "vo": (0.0, 3.5, 71, True),
        "pcs": DEFAULT_PCS_AXIS,
        "target_points_per_chunk": 250_000,
        "worker_group": "C",
        "drop_te_le": 2150.0,
    },
    "accident": {
        "description": "coarse accident and boundary-protection extension",
        "priority": 3,
        "te": (700.0, 2400.0, 86, True),
        "tc": (500.0, 1100.0, 61, True),
        "vo": (0.0, 3.5, 36, True),
        "pcs": DEFAULT_PCS_AXIS,
        "target_points_per_chunk": 300_000,
        "worker_group": "D",
    },
}

SMOKE_REGIONS: dict[str, dict[str, Any]] = {
    name: dict(
        spec,
        te=(spec["te"][0], spec["te"][1], 2, True),
        tc=(spec["tc"][0], spec["tc"][1], 2, True),
        vo=(0.0, 3.5, 2, True),
        pcs=(1.0, 1.2, 2, "linear"),
        target_points_per_chunk=8,
    )
    for name, spec in REGIONS.items()
}


def axis(start: float, stop: float, count: int, endpoint: bool = True) -> np.ndarray:
    if count <= 0:
        raise ValueError("axis count must be positive")
    return np.linspace(start, stop, count, endpoint=endpoint, dtype=np.float64)


def region_axes(region: dict[str, Any]) -> dict[str, np.ndarray]:
    te = axis(*region["te"])
    if "drop_te_le" in region:
        te = te[te > float(region["drop_te_le"])]
    pcs = pressure_axis_from_args(*region["pcs"])
    return {
        "TE_axis": te,
        "TC_axis": axis(*region["tc"]),
        "Vo_axis": axis(*region["vo"]),
        "Pcs_axis": pcs,
        "Tcs_axis": np.array([tcs_from_cesium_pressure(float(p)) for p in pcs], dtype=np.float64),
    }


def points_for_axes(axes: dict[str, np.ndarray]) -> int:
    return int(
        len(axes["TE_axis"])
        * len(axes["TC_axis"])
        * len(axes["Vo_axis"])
        * len(axes["Tcs_axis"])
    )


def git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=True,
        )
        return result.stdout.strip()
    except Exception:
        return None


def file_sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def pyd_info() -> dict[str, Any]:
    candidates = [
        ROOT / "build_cp312" / "Release" / "te_solver.cp312-win_amd64.pyd",
        ROOT / "te_solver.cp312-win_amd64.pyd",
    ]
    for path in candidates:
        if path.exists():
            return {
                "path": str(path.relative_to(REPO_ROOT)),
                "sha256": file_sha256(path),
                "size_bytes": path.stat().st_size,
                "mtime": path.stat().st_mtime,
            }
    return {"path": None, "sha256": None, "size_bytes": None, "mtime": None}


def chunk_te_ranges(n_te: int, points_per_te: int, target_points: int) -> list[tuple[int, int]]:
    if n_te <= 1:
        return [(0, n_te)]
    te_per_chunk = max(2, int(target_points // max(1, points_per_te)))
    ranges: list[tuple[int, int]] = []
    start = 0
    while start < n_te:
        stop = min(start + te_per_chunk, n_te)
        if stop < n_te:
            stop += 1  # include the next TE plane so interpolation spans chunk boundaries
        ranges.append((start, stop))
        if stop >= n_te:
            break
        start = stop - 1
    return ranges


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
    tmp.replace(path)


def load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def cmd_plan(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    axes_dir = db_dir / "axes"
    chunks: list[dict[str, Any]] = []
    regions_meta: dict[str, Any] = {}
    total_points = 0
    regions = SMOKE_REGIONS if args.preset == "smoke" else REGIONS

    for region_name, spec in regions.items():
        axes_data = region_axes(spec)
        axes_path = axes_dir / f"{region_name}_axes.npz"
        save_npz(axes_path, axes_data)

        n_points = points_for_axes(axes_data)
        total_points += n_points
        points_per_te = int(len(axes_data["TC_axis"]) * len(axes_data["Vo_axis"]) * len(axes_data["Tcs_axis"]))
        te_ranges = chunk_te_ranges(
            len(axes_data["TE_axis"]),
            points_per_te,
            int(spec["target_points_per_chunk"]),
        )
        regions_meta[region_name] = {
            "description": spec["description"],
            "priority": spec["priority"],
            "worker_group": spec["worker_group"],
            "axes_path": str(axes_path.relative_to(db_dir)),
            "point_count": n_points,
            "chunk_count": len(te_ranges),
            "axis_lengths": {name: int(len(values)) for name, values in axes_data.items()},
            "axis_min_max": {name: [float(values[0]), float(values[-1])] for name, values in axes_data.items()},
        }

        for chunk_index, (te_start, te_stop) in enumerate(te_ranges):
            chunk_id = f"{region_name}_{chunk_index:04d}"
            chunks.append(
                {
                    "chunk_id": chunk_id,
                    "region": region_name,
                    "worker_group": spec["worker_group"],
                    "priority": spec["priority"],
                    "te_start": int(te_start),
                    "te_stop": int(te_stop),
                    "point_count": int((te_stop - te_start) * points_per_te),
                    "output": f"chunks/{region_name}/{chunk_id}.npz",
                    "summary": f"chunks/{region_name}/{chunk_id}.summary.json",
                }
            )

    manifest = {
        "database_version": DATABASE_VERSION,
        "preset": args.preset,
        "created_at_unix": time.time(),
        "source_git_commit": git_commit(),
        "te_solver_pyd": pyd_info(),
        "diagnostic_interface": "te_solver.calc_emission_point",
        "d_gap": float(args.d_gap),
        "pressure_formula": "Pcs = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)",
        "field_groups": {
            "float_fields": list(FLOAT_FIELDS),
            "int_fields": list(INT_FIELDS),
            "bool_fields": list(BOOL_FIELDS),
            "derived_fields": [
                "valid_for_interpolation",
                "near_failed_region",
                "zero_emission_flag",
                "high_risk_flag",
                "source_region_id",
            ],
        },
        "source_region_ids": SOURCE_REGIONS,
        "regions": regions_meta,
        "chunk_count": len(chunks),
        "total_points": int(total_points),
        "region_priority": ["core", "startup", "high_power", "accident"],
        "overlap_policy": (
            "Regions are stored separately. Runtime lookup should prefer lower priority "
            "numbers: core, startup, high_power, then accident."
        ),
    }
    chunk_plan = {
        "database_version": DATABASE_VERSION,
        "db_dir": str(db_dir),
        "chunk_count": len(chunks),
        "total_points": int(total_points),
        "chunks": chunks,
    }

    write_json(db_dir / "manifest.json", manifest)
    write_json(db_dir / "chunk_plan.json", chunk_plan)
    print(json.dumps({"db_dir": str(db_dir), "chunk_count": len(chunks), "total_points": total_points}, indent=2))
    return 0


def chunk_arrays(shape: tuple[int, ...]) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for name in FLOAT_FIELDS:
        arrays[name] = np.full(shape, np.nan, dtype=np.float64)
    for name in INT_FIELDS:
        arrays[name] = np.full(shape, -1, dtype=np.int32)
    for name in ("converged", "finite_flag", "done"):
        arrays[name] = np.zeros(shape, dtype=bool)
    arrays["valid_for_interpolation"] = np.zeros(shape, dtype=bool)
    arrays["near_failed_region"] = np.zeros(shape, dtype=bool)
    arrays["zero_emission_flag"] = np.zeros(shape, dtype=bool)
    arrays["high_risk_flag"] = np.zeros(shape, dtype=bool)
    arrays["source_region_id"] = np.zeros(shape, dtype=np.int16)
    return arrays


def select_chunks(plan: dict[str, Any], args: argparse.Namespace) -> list[dict[str, Any]]:
    chunks = list(plan["chunks"])
    if args.region:
        chunks = [chunk for chunk in chunks if chunk["region"] == args.region]
    if args.chunk_id:
        wanted = set(args.chunk_id)
        chunks = [chunk for chunk in chunks if chunk["chunk_id"] in wanted]
    if args.worker_count > 1:
        chunks = [
            chunk
            for idx, chunk in enumerate(chunks)
            if idx % int(args.worker_count) == int(args.worker_index)
        ]
    return chunks


def load_axes(db_dir: Path, manifest: dict[str, Any], region: str) -> dict[str, np.ndarray]:
    axes_path = db_dir / manifest["regions"][region]["axes_path"]
    with np.load(axes_path, allow_pickle=False) as loaded:
        return {name: loaded[name] for name in loaded.files}


def compute_chunk(db_dir: Path, manifest: dict[str, Any], chunk: dict[str, Any], d_gap: float) -> dict[str, Any]:
    te_solver = import_te_solver()
    region = str(chunk["region"])
    axes_data = load_axes(db_dir, manifest, region)
    te_axis = axes_data["TE_axis"][int(chunk["te_start"]) : int(chunk["te_stop"])]
    tc_axis = axes_data["TC_axis"]
    vo_axis = axes_data["Vo_axis"]
    tcs_axis = axes_data["Tcs_axis"]
    pcs_axis = axes_data["Pcs_axis"]

    shape = (len(te_axis), len(tc_axis), len(vo_axis), len(tcs_axis))
    arrays = chunk_arrays(shape)
    arrays["TE_axis"] = te_axis
    arrays["TC_axis"] = tc_axis
    arrays["Vo_axis"] = vo_axis
    arrays["Tcs_axis"] = tcs_axis
    arrays["Pcs_axis"] = pcs_axis
    arrays["source_region_id"][:] = SOURCE_REGIONS[region]

    region_high_risk = region in {"high_power", "accident"}
    start = time.perf_counter()
    for idx in np.ndindex(shape):
        te = float(te_axis[idx[0]])
        tc = float(tc_axis[idx[1]])
        vo = float(vo_axis[idx[2]])
        tcs = float(tcs_axis[idx[3]])
        try:
            result = te_solver.calc_emission_point(te, tc, vo, tcs, d_gap)
            for name in FLOAT_FIELDS:
                arrays[name][idx] = float(result[name])
            for name in INT_FIELDS:
                if name != "error_code":
                    arrays[name][idx] = int(result[name])
            arrays["error_code"][idx] = 0
            arrays["converged"][idx] = bool(result["converged"])
            arrays["finite_flag"][idx] = bool(result["finite_flag"])
        except Exception:
            arrays["error_code"][idx] = 1
            arrays["regime"][idx] = -1
            arrays["converged"][idx] = False
            arrays["finite_flag"][idx] = False
        arrays["done"][idx] = True
        arrays["zero_emission_flag"][idx] = bool(abs(arrays["J"][idx]) < 1e-10)
        arrays["near_failed_region"][idx] = not bool(arrays["converged"][idx] and arrays["finite_flag"][idx])
        arrays["high_risk_flag"][idx] = bool(region_high_risk or arrays["near_failed_region"][idx])
        arrays["valid_for_interpolation"][idx] = bool(
            arrays["converged"][idx] and arrays["finite_flag"][idx] and not arrays["high_risk_flag"][idx]
        )

    output_path = db_dir / chunk["output"]
    save_npz(output_path, arrays)
    elapsed = time.perf_counter() - start
    summary = summarize(arrays, elapsed)
    summary.update(
        {
            "chunk_id": chunk["chunk_id"],
            "region": region,
            "output": chunk["output"],
            "te_start": chunk["te_start"],
            "te_stop": chunk["te_stop"],
            "zero_emission_points": int(np.count_nonzero(arrays["zero_emission_flag"])),
            "valid_for_interpolation_points": int(np.count_nonzero(arrays["valid_for_interpolation"])),
            "high_risk_points": int(np.count_nonzero(arrays["high_risk_flag"])),
        }
    )
    write_json(db_dir / chunk["summary"], summary)
    return summary


def cmd_worker(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    manifest = load_json(db_dir / "manifest.json")
    plan = load_json(db_dir / "chunk_plan.json")
    chunks = select_chunks(plan, args)
    if args.limit_chunks:
        chunks = chunks[: int(args.limit_chunks)]

    progress_path = db_dir / "summaries" / f"worker_{args.worker_id}_progress.json"
    completed = []
    failed = []
    start = time.perf_counter()
    for chunk in chunks:
        out = db_dir / chunk["output"]
        summary_path = db_dir / chunk["summary"]
        if out.exists() and summary_path.exists() and not args.force:
            completed.append({"chunk_id": chunk["chunk_id"], "status": "skipped"})
            continue
        try:
            summary = compute_chunk(db_dir, manifest, chunk, float(manifest["d_gap"]))
            completed.append({"chunk_id": chunk["chunk_id"], "status": "completed", "summary": summary})
        except Exception as exc:  # noqa: BLE001 - keep other chunks available.
            failed.append({"chunk_id": chunk["chunk_id"], "error": repr(exc)})
        write_json(
            progress_path,
            {
                "worker_id": args.worker_id,
                "updated_at_unix": time.time(),
                "completed_count": len(completed),
                "failed_count": len(failed),
                "completed": completed[-20:],
                "failed": failed,
            },
        )
        if args.max_runtime_s > 0 and (time.perf_counter() - start) >= args.max_runtime_s:
            break

    result = {
        "worker_id": args.worker_id,
        "selected_chunks": len(chunks),
        "completed_or_skipped": len(completed),
        "failed": len(failed),
        "elapsed_s": time.perf_counter() - start,
    }
    print(json.dumps(result, indent=2))
    return 0 if not failed else 1


def iter_chunk_summaries(db_dir: Path) -> Iterable[dict[str, Any]]:
    for path in sorted((db_dir / "chunks").glob("*/*.summary.json")):
        yield load_json(path)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="") as f:
            f.write("")
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def cmd_summarize(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    manifest = load_json(db_dir / "manifest.json")
    summaries = list(iter_chunk_summaries(db_dir))
    by_region: dict[str, dict[str, Any]] = {}
    global_totals = {
        "total_points": 0,
        "done_points": 0,
        "finite_points": 0,
        "converged_points": 0,
        "failed_or_nonfinite_points": 0,
        "zero_emission_points": 0,
        "valid_for_interpolation_points": 0,
        "high_risk_points": 0,
    }
    empty_totals = dict(global_totals)
    j_min = None
    j_max = None
    weighted_j_sum = 0.0
    weighted_j_count = 0

    for item in summaries:
        region = item["region"]
        bucket = by_region.setdefault(region, dict(empty_totals, region=region, chunk_count=0))
        bucket["chunk_count"] += 1
        for key in global_totals:
            value = int(item.get(key, 0))
            global_totals[key] += value
            bucket[key] += value
        if "J_min" in item:
            j_min = item["J_min"] if j_min is None else min(j_min, item["J_min"])
            j_max = item["J_max"] if j_max is None else max(j_max, item["J_max"])
            count = int(item.get("finite_points", 0))
            weighted_j_sum += float(item.get("J_mean", 0.0)) * count
            weighted_j_count += count

    global_summary = dict(global_totals)
    global_summary.update(
        {
            "database_version": manifest["database_version"],
            "chunk_count": len(summaries),
            "planned_chunk_count": manifest["chunk_count"],
            "planned_points": manifest["total_points"],
            "remaining_chunks": int(manifest["chunk_count"] - len(summaries)),
            "finite_rate_done": global_totals["finite_points"] / global_totals["done_points"]
            if global_totals["done_points"]
            else 0.0,
            "converged_rate_done": global_totals["converged_points"] / global_totals["done_points"]
            if global_totals["done_points"]
            else 0.0,
            "J_min": j_min,
            "J_max": j_max,
            "J_mean": weighted_j_sum / weighted_j_count if weighted_j_count else None,
        }
    )
    write_json(db_dir / "summaries" / "summary_global.json", global_summary)

    region_rows = []
    for row in by_region.values():
        row["finite_rate_done"] = row["finite_points"] / row["done_points"] if row["done_points"] else 0.0
        row["converged_rate_done"] = row["converged_points"] / row["done_points"] if row["done_points"] else 0.0
        region_rows.append(row)
    write_csv(db_dir / "summaries" / "summary_by_region.csv", sorted(region_rows, key=lambda r: r["region"]))

    if args.scan_chunks:
        axis_rows = aggregate_axes(db_dir)
        write_csv(db_dir / "summaries" / "summary_by_TE.csv", axis_rows["TE"])
        write_csv(db_dir / "summaries" / "summary_by_pressure.csv", axis_rows["Pcs"])
        write_csv(db_dir / "summaries" / "failure_boundary.csv", axis_rows["failure_boundary"])
        write_csv(db_dir / "summaries" / "zero_emission_boundary.csv", axis_rows["zero_emission_boundary"])

    print(json.dumps(global_summary, indent=2, sort_keys=True))
    return 0


def aggregate_axes(db_dir: Path) -> dict[str, list[dict[str, Any]]]:
    te_stats: dict[tuple[str, float], list[int]] = {}
    pcs_stats: dict[tuple[str, float], list[int]] = {}
    failure_rows = []
    zero_rows = []
    for chunk_path in sorted((db_dir / "chunks").glob("*/*.npz")):
        region = chunk_path.parent.name
        with np.load(chunk_path, allow_pickle=False) as data:
            bad = data["done"] & ~(data["finite_flag"] & data["converged"])
            zero = data["zero_emission_flag"]
            done = data["done"]
            for i, te in enumerate(data["TE_axis"]):
                key = (region, float(te))
                stats = te_stats.setdefault(key, [0, 0, 0])
                stats[0] += int(done[i, :, :, :].sum())
                stats[1] += int(bad[i, :, :, :].sum())
                stats[2] += int(zero[i, :, :, :].sum())
            for i, pcs in enumerate(data["Pcs_axis"]):
                key = (region, float(pcs))
                stats = pcs_stats.setdefault(key, [0, 0, 0])
                stats[0] += int(done[:, :, :, i].sum())
                stats[1] += int(bad[:, :, :, i].sum())
                stats[2] += int(zero[:, :, :, i].sum())
            if int(bad.sum()):
                failure_rows.append(
                    {
                        "region": region,
                        "chunk": chunk_path.stem,
                        "bad_points": int(bad.sum()),
                        "total_points": int(done.sum()),
                        "bad_rate": float(bad.sum() / done.sum()),
                        "TE_min": float(data["TE_axis"][0]),
                        "TE_max": float(data["TE_axis"][-1]),
                        "Pcs_min": float(data["Pcs_axis"][0]),
                        "Pcs_max": float(data["Pcs_axis"][-1]),
                    }
                )
            if int(zero.sum()):
                zero_rows.append(
                    {
                        "region": region,
                        "chunk": chunk_path.stem,
                        "zero_points": int(zero.sum()),
                        "total_points": int(done.sum()),
                        "zero_rate": float(zero.sum() / done.sum()),
                        "TE_min": float(data["TE_axis"][0]),
                        "TE_max": float(data["TE_axis"][-1]),
                    }
                )
    te_rows = [
        {
            "region": region,
            "TE": value,
            "done_points": stats[0],
            "bad_points": stats[1],
            "zero_emission_points": stats[2],
            "bad_rate": stats[1] / stats[0] if stats[0] else 0.0,
            "zero_rate": stats[2] / stats[0] if stats[0] else 0.0,
        }
        for (region, value), stats in te_stats.items()
    ]
    pcs_rows = [
        {
            "region": region,
            "Pcs": value,
            "Tcs": tcs_from_cesium_pressure(value),
            "done_points": stats[0],
            "bad_points": stats[1],
            "zero_emission_points": stats[2],
            "bad_rate": stats[1] / stats[0] if stats[0] else 0.0,
            "zero_rate": stats[2] / stats[0] if stats[0] else 0.0,
        }
        for (region, value), stats in pcs_stats.items()
    ]
    return {
        "TE": sorted(te_rows, key=lambda row: (row["region"], row["TE"])),
        "Pcs": sorted(pcs_rows, key=lambda row: (row["region"], row["Pcs"])),
        "failure_boundary": sorted(failure_rows, key=lambda row: row["bad_rate"], reverse=True),
        "zero_emission_boundary": sorted(zero_rows, key=lambda row: row["zero_rate"], reverse=True),
    }


def cmd_verify(args: argparse.Namespace) -> int:
    rng = np.random.default_rng(int(args.seed))
    db_dir = args.db_dir
    te_solver = import_te_solver()
    manifest = load_json(db_dir / "manifest.json")
    chunks = sorted((db_dir / "chunks").glob("*/*.npz"))
    if not chunks:
        raise RuntimeError("No chunk files found.")
    checked = []
    failures = []
    for chunk_path in rng.choice(chunks, size=min(int(args.samples), len(chunks)), replace=False):
        with np.load(chunk_path, allow_pickle=False) as data:
            shape = data["J"].shape
            idx = tuple(int(rng.integers(0, dim)) for dim in shape)
            te = float(data["TE_axis"][idx[0]])
            tc = float(data["TC_axis"][idx[1]])
            vo = float(data["Vo_axis"][idx[2]])
            tcs = float(data["Tcs_axis"][idx[3]])
            direct = te_solver.calc_emission_point(te, tc, vo, tcs, float(manifest["d_gap"]))
            ok = True
            diffs = {}
            for field in ("J", "Vd", "delta_V", "phiE", "phiC"):
                diff = abs(float(data[field][idx]) - float(direct[field]))
                diffs[field] = diff
                ok = ok and diff <= float(args.atol)
            item = {"chunk": str(chunk_path.relative_to(db_dir)), "idx": idx, "ok": ok, "diffs": diffs}
            checked.append(item)
            if not ok:
                failures.append(item)
    result = {"samples": len(checked), "failures": len(failures), "checked": checked, "failures_detail": failures}
    write_json(db_dir / "summaries" / "verification.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


def _mask_bbox(mask: np.ndarray, axes_data: dict[str, np.ndarray]) -> dict[str, list[float]] | None:
    if not bool(np.any(mask)):
        return None
    idx = np.argwhere(mask)
    axis_names = ("TE_axis", "TC_axis", "Vo_axis", "Pcs_axis")
    return {
        name.replace("_axis", ""): [
            float(axes_data[name][idx[:, axis_index]].min()),
            float(axes_data[name][idx[:, axis_index]].max()),
        ]
        for axis_index, name in enumerate(axis_names)
    }


def _top_axis_counts(mask: np.ndarray, axes_data: dict[str, np.ndarray], axis_name: str, limit: int = 10) -> list[dict[str, Any]]:
    if not bool(np.any(mask)):
        return []
    axis_order = {"TE_axis": 0, "TC_axis": 1, "Vo_axis": 2, "Pcs_axis": 3}
    axis_index = axis_order[axis_name]
    counts = np.sum(mask, axis=tuple(idx for idx in range(4) if idx != axis_index))
    order = np.argsort(counts)[::-1]
    rows = []
    for idx in order[:limit]:
        count = int(counts[idx])
        if count <= 0:
            break
        rows.append({"value": float(axes_data[axis_name][idx]), "count": count})
    return rows


def cmd_refine_risk(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    manifest = load_json(db_dir / "manifest.json")
    chunk_plan = load_json(db_dir / "chunk_plan.json")
    rows: list[dict[str, Any]] = []
    by_region: dict[str, dict[str, Any]] = {}

    zero_threshold = float(args.zero_j_threshold)
    for chunk in chunk_plan["chunks"]:
        chunk_path = db_dir / chunk["output"]
        if not chunk_path.exists():
            continue
        region = str(chunk["region"])
        axes_data = load_axes(db_dir, manifest, region)
        with np.load(chunk_path, allow_pickle=False) as data:
            done = np.asarray(data["done"], dtype=bool)
            finite = np.asarray(data["finite_flag"], dtype=bool)
            converged = np.asarray(data["converged"], dtype=bool)
            j = np.asarray(data["J"], dtype=np.float64)

            invalid = done & ~(finite & converged)
            zero_fill = invalid & finite & np.isfinite(j) & (np.abs(j) <= zero_threshold)
            analytic_fallback = invalid & ~zero_fill
            lookup_safe = (done & finite & converged) | zero_fill

            row = {
                "region": region,
                "chunk_id": chunk["chunk_id"],
                "total_points": int(done.sum()),
                "invalid_points": int(invalid.sum()),
                "zero_fill_points": int(zero_fill.sum()),
                "analytic_fallback_points": int(analytic_fallback.sum()),
                "lookup_safe_points": int(lookup_safe.sum()),
                "zero_fill_rate_invalid": float(zero_fill.sum() / invalid.sum()) if int(invalid.sum()) else 0.0,
                "analytic_fallback_rate_total": float(analytic_fallback.sum() / done.sum()) if int(done.sum()) else 0.0,
            }
            rows.append(row)

            bucket = by_region.setdefault(
                region,
                {
                    "region": region,
                    "total_points": 0,
                    "invalid_points": 0,
                    "zero_fill_points": 0,
                    "analytic_fallback_points": 0,
                    "lookup_safe_points": 0,
                    "fallback_bbox": None,
                    "zero_fill_bbox": None,
                    "top_TE_fallback": [],
                    "top_Vo_fallback": [],
                    "top_Pcs_fallback": [],
                },
            )
            for key in ("total_points", "invalid_points", "zero_fill_points", "analytic_fallback_points", "lookup_safe_points"):
                bucket[key] += int(row[key])

            chunk_axes = {
                "TE_axis": np.asarray(data["TE_axis"], dtype=np.float64),
                "TC_axis": np.asarray(data["TC_axis"], dtype=np.float64),
                "Vo_axis": np.asarray(data["Vo_axis"], dtype=np.float64),
                "Pcs_axis": np.asarray(data["Pcs_axis"], dtype=np.float64),
            }
            row["fallback_bbox"] = _mask_bbox(analytic_fallback, chunk_axes)
            row["zero_fill_bbox"] = _mask_bbox(zero_fill, chunk_axes)

            if args.write_flags:
                flags_path = chunk_path.with_suffix(".lookup_flags.npz")
                save_npz(
                    flags_path,
                    {
                        "lookup_safe_flag": lookup_safe,
                        "zero_fill_flag": zero_fill,
                        "analytic_fallback_flag": analytic_fallback,
                        "invalid_flag": invalid,
                        "zero_j_threshold": np.array([zero_threshold], dtype=np.float64),
                    },
                )

    global_summary = {
        "database_version": manifest["database_version"],
        "zero_j_threshold": zero_threshold,
        "total_points": int(sum(row["total_points"] for row in rows)),
        "invalid_points": int(sum(row["invalid_points"] for row in rows)),
        "zero_fill_points": int(sum(row["zero_fill_points"] for row in rows)),
        "analytic_fallback_points": int(sum(row["analytic_fallback_points"] for row in rows)),
        "lookup_safe_points": int(sum(row["lookup_safe_points"] for row in rows)),
        "chunk_count": len(rows),
    }
    global_summary["zero_fill_rate_invalid"] = (
        global_summary["zero_fill_points"] / global_summary["invalid_points"]
        if global_summary["invalid_points"]
        else 0.0
    )
    global_summary["analytic_fallback_rate_total"] = (
        global_summary["analytic_fallback_points"] / global_summary["total_points"]
        if global_summary["total_points"]
        else 0.0
    )
    global_summary["lookup_safe_rate_total"] = (
        global_summary["lookup_safe_points"] / global_summary["total_points"]
        if global_summary["total_points"]
        else 0.0
    )

    # Recompute concise region bboxes/tops from full region arrays for chunk-boundary correctness.
    for region, bucket in by_region.items():
        axes_data = load_axes(db_dir, manifest, region)
        shape = (
            len(axes_data["TE_axis"]),
            len(axes_data["TC_axis"]),
            len(axes_data["Vo_axis"]),
            len(axes_data["Pcs_axis"]),
        )
        fallback_full = np.zeros(shape, dtype=bool)
        zero_full = np.zeros(shape, dtype=bool)
        for chunk in chunk_plan["chunks"]:
            if chunk["region"] != region:
                continue
            flags_path = (db_dir / chunk["output"]).with_suffix(".lookup_flags.npz")
            if not flags_path.exists():
                continue
            with np.load(flags_path, allow_pickle=False) as flags:
                te_start = int(chunk["te_start"])
                te_stop = int(chunk["te_stop"])
                fallback_full[te_start:te_stop, :, :, :] = flags["analytic_fallback_flag"]
                zero_full[te_start:te_stop, :, :, :] = flags["zero_fill_flag"]
        bucket["fallback_bbox"] = _mask_bbox(fallback_full, axes_data)
        bucket["zero_fill_bbox"] = _mask_bbox(zero_full, axes_data)
        bucket["top_TE_fallback"] = _top_axis_counts(fallback_full, axes_data, "TE_axis")
        bucket["top_Vo_fallback"] = _top_axis_counts(fallback_full, axes_data, "Vo_axis")
        bucket["top_Pcs_fallback"] = _top_axis_counts(fallback_full, axes_data, "Pcs_axis")
        bucket["zero_fill_rate_invalid"] = (
            bucket["zero_fill_points"] / bucket["invalid_points"]
            if bucket["invalid_points"]
            else 0.0
        )
        bucket["analytic_fallback_rate_total"] = (
            bucket["analytic_fallback_points"] / bucket["total_points"]
            if bucket["total_points"]
            else 0.0
        )

    write_json(db_dir / "summaries" / "lookup_risk_refinement.json", {"global": global_summary, "regions": by_region})
    write_csv(db_dir / "summaries" / "lookup_risk_refinement_by_chunk.csv", rows)
    print(json.dumps({"global": global_summary, "regions": by_region}, indent=2, sort_keys=True))
    return 0


def _load_region_dense(db_dir: Path, manifest: dict[str, Any], plan: dict[str, Any], region: str) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    axes_data = load_axes(db_dir, manifest, region)
    shape = (
        len(axes_data["TE_axis"]),
        len(axes_data["TC_axis"]),
        len(axes_data["Vo_axis"]),
        len(axes_data["Pcs_axis"]),
    )
    arrays = chunk_arrays(shape)
    arrays["TE_axis"] = axes_data["TE_axis"]
    arrays["TC_axis"] = axes_data["TC_axis"]
    arrays["Vo_axis"] = axes_data["Vo_axis"]
    arrays["Tcs_axis"] = axes_data["Tcs_axis"]
    arrays["Pcs_axis"] = axes_data["Pcs_axis"]
    for chunk in plan["chunks"]:
        if chunk["region"] != region:
            continue
        chunk_path = db_dir / chunk["output"]
        if not chunk_path.exists():
            continue
        with np.load(chunk_path, allow_pickle=False) as data:
            sl = np.s_[int(chunk["te_start"]) : int(chunk["te_stop"]), :, :, :]
            for name in FLOAT_FIELDS:
                arrays[name][sl] = data[name]
            for name in INT_FIELDS:
                arrays[name][sl] = data[name]
            for name in ("converged", "finite_flag", "done"):
                arrays[name][sl] = data[name]
            for name in ("valid_for_interpolation", "near_failed_region", "zero_emission_flag", "high_risk_flag", "source_region_id"):
                if name in data.files:
                    arrays[name][sl] = data[name]
    return axes_data, arrays


def _nearest_valid_fill(
    field: np.ndarray,
    valid: np.ndarray,
    fill_mask: np.ndarray,
    *,
    log_space: bool = False,
    min_positive: float = 1e-300,
) -> tuple[np.ndarray, np.ndarray]:
    out = np.array(field, copy=True)
    if not bool(np.any(fill_mask)):
        return out, np.zeros(field.shape, dtype=bool)
    valid_values = field[valid]
    if valid_values.size == 0:
        return out, np.array(fill_mask, copy=True)

    if log_space:
        source = np.log(np.clip(field, min_positive, None))
    else:
        source = field

    try:
        from scipy import ndimage
    except Exception:
        # Fallback: global median/nearest representative if scipy is not present.
        fill_value = float(np.nanmedian(source[valid]))
        source_out = np.array(source, copy=True)
        source_out[fill_mask] = fill_value
        out[fill_mask] = np.exp(source_out[fill_mask]) if log_space else source_out[fill_mask]
        return out, np.zeros(field.shape, dtype=bool)

    invalid_for_distance = ~valid
    _, indices = ndimage.distance_transform_edt(invalid_for_distance, return_indices=True)
    nearest = tuple(index_array[fill_mask] for index_array in indices)
    source_out = np.array(source, copy=True)
    source_out[fill_mask] = source[nearest]
    out[fill_mask] = np.exp(source_out[fill_mask]) if log_space else source_out[fill_mask]
    return out, np.zeros(field.shape, dtype=bool)


def cmd_optimize_table(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    manifest = load_json(db_dir / "manifest.json")
    plan = load_json(db_dir / "chunk_plan.json")
    zero_threshold = float(args.zero_j_threshold)
    regions_requested = set(args.region or manifest["regions"].keys())
    region_summaries: dict[str, Any] = {}
    global_summary = {
        "total_points": 0,
        "raw_invalid_points": 0,
        "zero_fill_points": 0,
        "imputed_points": 0,
        "optimized_safe_points": 0,
        "unresolved_points": 0,
    }

    float_fields_to_optimize = ("J", "Vd", "delta_V", "phiE", "phiC")

    for region in manifest["regions"]:
        if region not in regions_requested:
            continue
        axes_data, arrays = _load_region_dense(db_dir, manifest, plan, region)
        done = arrays["done"]
        finite = arrays["finite_flag"]
        converged = arrays["converged"]
        j = arrays["J"]
        raw_valid = done & finite & converged
        raw_invalid = done & ~(finite & converged)
        zero_fill = raw_invalid & finite & np.isfinite(j) & (np.abs(j) <= zero_threshold)
        impute = raw_invalid & ~zero_fill
        optimized_safe = raw_valid | zero_fill | impute
        unresolved = done & ~optimized_safe

        optimized: dict[str, np.ndarray] = {
            "lookup_safe_flag": optimized_safe,
            "raw_valid_flag": raw_valid,
            "raw_invalid_flag": raw_invalid,
            "zero_fill_flag": zero_fill,
            "imputed_flag": impute,
            "unresolved_flag": unresolved,
            "zero_j_threshold": np.array([zero_threshold], dtype=np.float64),
            "optimization_version": np.array([1], dtype=np.int32),
        }
        for name in float_fields_to_optimize:
            optimized[name] = np.array(arrays[name], copy=True)
        for name in ("regime", "iteration_count"):
            optimized[name] = np.array(arrays[name], copy=True)

        optimized["J"][zero_fill] = 0.0
        optimized["regime"][zero_fill] = 0
        optimized["iteration_count"][zero_fill] = 0

        source_valid = raw_valid & np.isfinite(optimized["J"]) & (optimized["J"] >= 0.0)
        source_valid_nonzero_j = source_valid & (optimized["J"] > zero_threshold)
        unresolved_fields = []
        for name in float_fields_to_optimize:
            field_valid = source_valid_nonzero_j if name == "J" else source_valid
            field_valid = field_valid & np.isfinite(optimized[name])
            log_space = name == "J"
            optimized[name], unresolved_field = _nearest_valid_fill(
                optimized[name],
                field_valid,
                impute,
                log_space=log_space,
            )
            if bool(np.any(unresolved_field)):
                unresolved_fields.append(name)

        # Copy discrete fields from nearest raw-valid point.
        try:
            from scipy import ndimage

            discrete_source_valid = source_valid_nonzero_j if bool(np.any(source_valid_nonzero_j)) else source_valid
            _, indices = ndimage.distance_transform_edt(~discrete_source_valid, return_indices=True)
            nearest = tuple(index_array[impute] for index_array in indices)
            optimized["regime"][impute] = arrays["regime"][nearest]
            optimized["iteration_count"][impute] = -2
        except Exception:
            optimized["regime"][impute] = -2
            optimized["iteration_count"][impute] = -2

        for chunk in plan["chunks"]:
            if chunk["region"] != region:
                continue
            sl = np.s_[int(chunk["te_start"]) : int(chunk["te_stop"]), :, :, :]
            chunk_path = db_dir / chunk["output"]
            out_path = chunk_path.with_suffix(".optimized.npz")
            payload = {
                "TE_axis": arrays["TE_axis"][sl[0]],
                "TC_axis": arrays["TC_axis"],
                "Vo_axis": arrays["Vo_axis"],
                "Tcs_axis": arrays["Tcs_axis"],
                "Pcs_axis": arrays["Pcs_axis"],
            }
            for name, values in optimized.items():
                if isinstance(values, np.ndarray) and values.shape == done.shape:
                    payload[name] = values[sl]
                else:
                    payload[name] = values
            save_npz(out_path, payload)

        region_summary = {
            "region": region,
            "total_points": int(done.sum()),
            "raw_invalid_points": int(raw_invalid.sum()),
            "zero_fill_points": int(zero_fill.sum()),
            "imputed_points": int(impute.sum()),
            "optimized_safe_points": int(optimized_safe.sum()),
            "unresolved_points": int(unresolved.sum()),
            "zero_fill_bbox": _mask_bbox(zero_fill, axes_data),
            "imputed_bbox": _mask_bbox(impute, axes_data),
            "top_TE_imputed": _top_axis_counts(impute, axes_data, "TE_axis"),
            "top_Vo_imputed": _top_axis_counts(impute, axes_data, "Vo_axis"),
            "top_Pcs_imputed": _top_axis_counts(impute, axes_data, "Pcs_axis"),
            "unresolved_fields": unresolved_fields,
        }
        region_summary["optimized_safe_rate"] = (
            region_summary["optimized_safe_points"] / region_summary["total_points"]
            if region_summary["total_points"]
            else 0.0
        )
        region_summaries[region] = region_summary
        for key in global_summary:
            global_summary[key] += int(region_summary[key])

    global_summary["zero_j_threshold"] = zero_threshold
    global_summary["optimized_safe_rate"] = (
        global_summary["optimized_safe_points"] / global_summary["total_points"]
        if global_summary["total_points"]
        else 0.0
    )
    global_summary["unresolved_rate"] = (
        global_summary["unresolved_points"] / global_summary["total_points"]
        if global_summary["total_points"]
        else 0.0
    )
    result = {"global": global_summary, "regions": region_summaries}
    write_json(db_dir / "summaries" / "optimized_table_summary.json", result)
    write_csv(db_dir / "summaries" / "optimized_table_by_region.csv", list(region_summaries.values()))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


def _runtime_chunk_source_path(db_dir: Path, chunk: dict[str, Any]) -> Path:
    raw_path = db_dir / chunk["output"]
    optimized_path = raw_path.with_suffix(".optimized.npz")
    return optimized_path if optimized_path.exists() else raw_path


def _runtime_payload_from_source(
    source_path: Path,
    *,
    dtype: type[np.floating],
    zero_threshold: float,
    zero_compress: bool,
) -> dict[str, np.ndarray]:
    with np.load(source_path, allow_pickle=False) as data:
        if "lookup_safe_flag" in data.files:
            safe = np.asarray(data["lookup_safe_flag"], dtype=np.uint8)
        else:
            safe = np.asarray(data["done"] & data["finite_flag"] & data["converged"], dtype=np.uint8)

        j64 = np.asarray(data["J"], dtype=np.float64)
        zero_mask = safe.astype(bool) & np.isfinite(j64) & (np.abs(j64) <= zero_threshold)
        if "zero_fill_flag" in data.files:
            zero_mask |= np.asarray(data["zero_fill_flag"], dtype=bool)
        elif "zero_emission_flag" in data.files:
            zero_mask |= np.asarray(data["zero_emission_flag"], dtype=bool) & safe.astype(bool)

        payload = {
            "TE_axis": np.asarray(data["TE_axis"], dtype=dtype),
            "TC_axis": np.asarray(data["TC_axis"], dtype=dtype),
            "Vo_axis": np.asarray(data["Vo_axis"], dtype=dtype),
            "Tcs_axis": np.asarray(data["Tcs_axis"], dtype=dtype),
            "J": np.asarray(data["J"], dtype=dtype),
            "Vd": np.asarray(data["Vd"], dtype=dtype),
            "delta_V": np.asarray(data["delta_V"], dtype=dtype),
            "phiE": np.asarray(data["phiE"], dtype=dtype),
            "phiC": np.asarray(data["phiC"], dtype=dtype),
            "lookup_safe": safe.astype(np.uint8, copy=False),
            "zero_mask": zero_mask.astype(np.uint8, copy=False),
        }
        if zero_compress:
            payload["J"] = np.array(payload["J"], copy=True)
            payload["J"][zero_mask] = 0.0
        return payload


def _append_first_te_plane(payload: dict[str, np.ndarray], right_payload: dict[str, np.ndarray]) -> bool:
    if payload["TE_axis"].size == 0 or right_payload["TE_axis"].size == 0:
        return False
    if not float(payload["TE_axis"][-1]) < float(right_payload["TE_axis"][0]):
        return False
    for axis_name in ("TC_axis", "Vo_axis", "Tcs_axis"):
        if payload[axis_name].shape != right_payload[axis_name].shape:
            return False
        if not np.allclose(payload[axis_name], right_payload[axis_name], rtol=0.0, atol=1.0e-10):
            return False

    payload["TE_axis"] = np.concatenate((payload["TE_axis"], right_payload["TE_axis"][:1]))
    for field in ("J", "Vd", "delta_V", "phiE", "phiC", "lookup_safe", "zero_mask"):
        payload[field] = np.concatenate((payload[field], right_payload[field][:1, :, :, :]), axis=0)
    return True


def cmd_export_runtime(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    out_dir = args.out_dir
    manifest = load_json(db_dir / "manifest.json")
    plan = load_json(db_dir / "chunk_plan.json")
    regions_requested = set(args.region or manifest["regions"].keys())
    limit_chunks = int(args.limit_chunks)
    zero_threshold = float(args.zero_j_threshold)
    dtype = np.float32 if args.dtype == "float32" else np.float64
    chunks_out: list[dict[str, Any]] = []
    regions_out: dict[str, Any] = {}
    total_points = 0
    total_bytes = 0

    for region, meta in manifest["regions"].items():
        if region not in regions_requested:
            continue
        regions_out[region] = {
            "priority": int(meta["priority"]),
            "region_id": int(manifest.get("source_region_ids", SOURCE_REGIONS).get(region, SOURCE_REGIONS.get(region, 99))),
            "axis_min_max": meta.get("axis_min_max", {}),
            "axis_lengths": meta.get("axis_lengths", {}),
        }

    source_chunks = [chunk for chunk in plan["chunks"] if chunk["region"] in regions_requested]
    source_chunks_by_region: dict[str, list[dict[str, Any]]] = {}
    for chunk in source_chunks:
        source_chunks_by_region.setdefault(str(chunk["region"]), []).append(chunk)
    next_chunk_by_id: dict[str, dict[str, Any]] = {}
    for region_chunks in source_chunks_by_region.values():
        ordered = sorted(region_chunks, key=lambda item: (int(item["te_start"]), int(item["te_stop"])))
        by_start = {int(item["te_start"]): item for item in ordered}
        for chunk in ordered:
            next_chunk = by_start.get(int(chunk["te_stop"]))
            if next_chunk is not None:
                next_chunk_by_id[str(chunk["chunk_id"])] = next_chunk

    selected_chunks = list(source_chunks)
    if limit_chunks > 0:
        selected_chunks = selected_chunks[:limit_chunks]

    for chunk in selected_chunks:
        region = str(chunk["region"])
        source_path = _runtime_chunk_source_path(db_dir, chunk)
        if not source_path.exists():
            continue
        payload = _runtime_payload_from_source(
            source_path,
            dtype=dtype,
            zero_threshold=zero_threshold,
            zero_compress=bool(args.zero_compress),
        )
        stitched_right_boundary = False
        next_chunk = next_chunk_by_id.get(str(chunk["chunk_id"]))
        if next_chunk is not None:
            next_source_path = _runtime_chunk_source_path(db_dir, next_chunk)
            if next_source_path.exists():
                next_payload = _runtime_payload_from_source(
                    next_source_path,
                    dtype=dtype,
                    zero_threshold=zero_threshold,
                    zero_compress=bool(args.zero_compress),
                )
                stitched_right_boundary = _append_first_te_plane(payload, next_payload)

        chunk_id = str(chunk["chunk_id"])
        rel_output = Path(region) / f"{chunk_id}.runtime.npz"
        out_path = out_dir / rel_output
        save_npz(out_path, payload)
        size_bytes = out_path.stat().st_size
        total_bytes += size_bytes
        points = int(np.prod(payload["J"].shape))
        total_points += points
        chunk_meta = {
            "chunk_id": chunk_id,
            "region": region,
            "priority": int(regions_out[region]["priority"]),
            "region_id": int(regions_out[region]["region_id"]),
            "source": str(source_path.relative_to(db_dir)),
            "output": str(rel_output).replace("\\", "/"),
            "point_count": points,
            "size_bytes": size_bytes,
            "stitched_right_boundary": bool(stitched_right_boundary),
            "TE_min": float(payload["TE_axis"][0]),
            "TE_max": float(payload["TE_axis"][-1]),
            "TC_min": float(payload["TC_axis"][0]),
            "TC_max": float(payload["TC_axis"][-1]),
            "Vo_min": float(payload["Vo_axis"][0]),
            "Vo_max": float(payload["Vo_axis"][-1]),
            "Tcs_min": float(payload["Tcs_axis"][0]),
            "Tcs_max": float(payload["Tcs_axis"][-1]),
            "lookup_safe_points": int(np.count_nonzero(payload["lookup_safe"])),
            "zero_mask_points": int(np.count_nonzero(payload["zero_mask"])),
        }
        chunks_out.append(chunk_meta)

    runtime_manifest = {
        "runtime_database_version": "emission-runtime-db-v1",
        "source_database_version": manifest.get("database_version"),
        "created_at_unix": time.time(),
        "source_db_dir": str(db_dir),
        "d_gap": float(manifest.get("d_gap", 0.5)),
        "dtype": args.dtype,
        "zero_j_threshold": zero_threshold,
        "zero_compress": bool(args.zero_compress),
        "fields": ["J", "Vd", "delta_V", "phiE", "phiC", "lookup_safe", "zero_mask"],
        "regions": regions_out,
        "chunks": chunks_out,
        "chunk_count": len(chunks_out),
        "total_points": int(total_points),
        "total_size_bytes": int(total_bytes),
    }
    write_json(out_dir / "runtime_manifest.json", runtime_manifest)
    print(json.dumps(runtime_manifest, indent=2, sort_keys=True))
    return 0


def _pack_mask(mask: np.ndarray) -> np.ndarray:
    flat = np.asarray(mask, dtype=np.uint8).reshape(-1)
    return np.packbits(flat, bitorder="little")


def _dense_output_paths(out_dir: Path, region: str, fmt: str) -> list[Path]:
    paths: list[Path] = []
    if fmt in {"npz", "both"}:
        paths.append(out_dir / f"{region}.runtime.v2.npz")
    if fmt in {"binary", "both"}:
        paths.append(out_dir / f"{region}.runtime.v2.tedb")
    return paths


def _write_dense_binary(path: Path, region_name: str, meta: dict[str, Any], payload: dict[str, np.ndarray]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    point_count = int(payload["J"].size)
    bit_bytes = int(payload["lookup_safe_bits"].size)
    name_bytes = region_name.encode("utf-8")
    header = struct.pack(
        "<8sIIiidQQQQQQ",
        b"TEDBv2\0\0",
        1,
        len(name_bytes),
        int(meta["priority"]),
        int(meta["region_id"]),
        float(meta.get("d_gap", 0.5)),
        int(payload["TE_axis"].size),
        int(payload["TC_axis"].size),
        int(payload["Vo_axis"].size),
        int(payload["Tcs_axis"].size),
        point_count,
        bit_bytes,
    )
    with path.open("wb") as f:
        f.write(header)
        f.write(name_bytes)
        for name in ("TE_axis", "TC_axis", "Vo_axis", "Tcs_axis"):
            f.write(np.asarray(payload[name], dtype=np.float64).tobytes(order="C"))
        for name in ("J", "Vd", "delta_V", "phiE", "phiC"):
            f.write(np.asarray(payload[name], dtype=np.float32).reshape(-1).tobytes(order="C"))
        f.write(np.asarray(payload["lookup_safe_bits"], dtype=np.uint8).tobytes(order="C"))
        f.write(np.asarray(payload["zero_mask_bits"], dtype=np.uint8).tobytes(order="C"))
    return path.stat().st_size


def _make_dense_region_payload(
    db_dir: Path,
    manifest: dict[str, Any],
    plan: dict[str, Any],
    region: str,
    *,
    dtype: type[np.floating],
    zero_threshold: float,
    zero_compress: bool,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    axes_data = load_axes(db_dir, manifest, region)
    shape = (
        int(len(axes_data["TE_axis"])),
        int(len(axes_data["TC_axis"])),
        int(len(axes_data["Vo_axis"])),
        int(len(axes_data["Tcs_axis"])),
    )
    fields = {
        name: np.zeros(shape, dtype=dtype)
        for name in ("J", "Vd", "delta_V", "phiE", "phiC")
    }
    lookup_safe = np.zeros(shape, dtype=bool)
    zero_mask = np.zeros(shape, dtype=bool)
    filled = np.zeros(shape, dtype=bool)
    chunks_used = 0

    for chunk in plan["chunks"]:
        if str(chunk["region"]) != region:
            continue
        source_path = _runtime_chunk_source_path(db_dir, chunk)
        if not source_path.exists():
            continue
        chunk_payload = _runtime_payload_from_source(
            source_path,
            dtype=dtype,
            zero_threshold=zero_threshold,
            zero_compress=zero_compress,
        )
        te_start = int(chunk["te_start"])
        te_stop = int(chunk["te_stop"])
        sl = np.s_[te_start:te_stop, :, :, :]
        for name in fields:
            fields[name][sl] = np.asarray(chunk_payload[name], dtype=dtype)
        lookup_safe[sl] = np.asarray(chunk_payload["lookup_safe"], dtype=bool)
        zero_mask[sl] = np.asarray(chunk_payload["zero_mask"], dtype=bool)
        filled[sl] = True
        chunks_used += 1

    payload = {
        "TE_axis": np.asarray(axes_data["TE_axis"], dtype=np.float64),
        "TC_axis": np.asarray(axes_data["TC_axis"], dtype=np.float64),
        "Vo_axis": np.asarray(axes_data["Vo_axis"], dtype=np.float64),
        "Tcs_axis": np.asarray(axes_data["Tcs_axis"], dtype=np.float64),
        "J": fields["J"],
        "Vd": fields["Vd"],
        "delta_V": fields["delta_V"],
        "phiE": fields["phiE"],
        "phiC": fields["phiC"],
        "lookup_safe_bits": _pack_mask(lookup_safe),
        "zero_mask_bits": _pack_mask(zero_mask),
    }
    point_count = int(np.prod(shape))
    meta = {
        "region": region,
        "priority": int(manifest["regions"][region]["priority"]),
        "region_id": int(manifest.get("source_region_ids", SOURCE_REGIONS).get(region, SOURCE_REGIONS.get(region, 99))),
        "d_gap": float(manifest.get("d_gap", 0.5)),
        "shape": list(shape),
        "point_count": point_count,
        "bit_bytes": int(payload["lookup_safe_bits"].size),
        "chunks_used": chunks_used,
        "filled_points": int(np.count_nonzero(filled)),
        "lookup_safe_points": int(np.count_nonzero(lookup_safe)),
        "zero_mask_points": int(np.count_nonzero(zero_mask)),
        "axis_min_max": {
            "TE_axis": [float(payload["TE_axis"][0]), float(payload["TE_axis"][-1])],
            "TC_axis": [float(payload["TC_axis"][0]), float(payload["TC_axis"][-1])],
            "Vo_axis": [float(payload["Vo_axis"][0]), float(payload["Vo_axis"][-1])],
            "Tcs_axis": [float(payload["Tcs_axis"][0]), float(payload["Tcs_axis"][-1])],
        },
    }
    return payload, meta


def cmd_export_runtime_dense(args: argparse.Namespace) -> int:
    db_dir = args.db_dir
    out_dir = args.out_dir
    manifest = load_json(db_dir / "manifest.json")
    plan = load_json(db_dir / "chunk_plan.json")
    regions_requested = list(args.region or manifest["regions"].keys())
    dtype = np.float32 if args.dtype == "float32" else np.float64
    fmt = str(args.format)
    region_outputs: list[dict[str, Any]] = []
    total_size = 0
    total_points = 0

    for region in regions_requested:
        if region not in manifest["regions"]:
            raise ValueError(f"Unknown emission database region: {region}")
        payload, meta = _make_dense_region_payload(
            db_dir,
            manifest,
            plan,
            region,
            dtype=dtype,
            zero_threshold=float(args.zero_j_threshold),
            zero_compress=bool(args.zero_compress),
        )
        outputs: dict[str, str] = {}
        sizes: dict[str, int] = {}
        if fmt in {"npz", "both"}:
            npz_path = out_dir / f"{region}.runtime.v2.npz"
            save_npz(npz_path, payload)
            sizes["npz"] = npz_path.stat().st_size
            outputs["npz"] = npz_path.name
        if fmt in {"binary", "both"}:
            tedb_path = out_dir / f"{region}.runtime.v2.tedb"
            sizes["binary"] = _write_dense_binary(tedb_path, region, meta, payload)
            outputs["binary"] = tedb_path.name
        total_size += sum(sizes.values())
        total_points += int(meta["point_count"])
        region_outputs.append({**meta, "outputs": outputs, "size_bytes": sizes})

    runtime_manifest = {
        "runtime_database_version": "emission-runtime-dense-v2",
        "source_database_version": manifest.get("database_version"),
        "created_at_unix": time.time(),
        "source_db_dir": str(db_dir),
        "d_gap": float(manifest.get("d_gap", 0.5)),
        "dtype": args.dtype,
        "format": fmt,
        "bitorder": "little",
        "zero_j_threshold": float(args.zero_j_threshold),
        "zero_compress": bool(args.zero_compress),
        "fields": ["J", "Vd", "delta_V", "phiE", "phiC", "lookup_safe_bits", "zero_mask_bits"],
        "regions": {item["region"]: item for item in region_outputs},
        "region_count": len(region_outputs),
        "total_points": int(total_points),
        "total_size_bytes": int(total_size),
    }
    write_json(out_dir / "runtime_dense_manifest.json", runtime_manifest)
    print(json.dumps(runtime_manifest, indent=2, sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate and manage ThermoCalc emission database chunks.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_plan = sub.add_parser("plan")
    p_plan.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_plan.add_argument("--d-gap", type=float, default=0.5)
    p_plan.add_argument("--preset", choices=("full", "smoke"), default="full")
    p_plan.set_defaults(func=cmd_plan)

    p_worker = sub.add_parser("worker")
    p_worker.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_worker.add_argument("--worker-id", default="worker")
    p_worker.add_argument("--worker-index", type=int, default=0)
    p_worker.add_argument("--worker-count", type=int, default=1)
    p_worker.add_argument("--region", choices=tuple(REGIONS), default=None)
    p_worker.add_argument("--chunk-id", action="append", default=None)
    p_worker.add_argument("--limit-chunks", type=int, default=0)
    p_worker.add_argument("--max-runtime-s", type=float, default=0.0)
    p_worker.add_argument("--force", action="store_true")
    p_worker.set_defaults(func=cmd_worker)

    p_sum = sub.add_parser("summarize")
    p_sum.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_sum.add_argument("--scan-chunks", action="store_true")
    p_sum.set_defaults(func=cmd_summarize)

    p_verify = sub.add_parser("verify")
    p_verify.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_verify.add_argument("--samples", type=int, default=20)
    p_verify.add_argument("--seed", type=int, default=1)
    p_verify.add_argument("--atol", type=float, default=1e-12)
    p_verify.set_defaults(func=cmd_verify)

    p_refine = sub.add_parser("refine-risk")
    p_refine.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_refine.add_argument("--zero-j-threshold", type=float, default=1e-4)
    p_refine.add_argument("--write-flags", action="store_true")
    p_refine.set_defaults(func=cmd_refine_risk)

    p_opt = sub.add_parser("optimize-table")
    p_opt.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_opt.add_argument("--zero-j-threshold", type=float, default=1e-3)
    p_opt.add_argument("--region", action="append", default=None)
    p_opt.set_defaults(func=cmd_optimize_table)

    p_runtime = sub.add_parser("export-runtime")
    p_runtime.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_runtime.add_argument("--out-dir", type=Path, default=ROOT / "emission_runtime_db")
    p_runtime.add_argument("--dtype", choices=("float32", "float64"), default="float32")
    p_runtime.add_argument("--region", action="append", default=None)
    p_runtime.add_argument("--limit-chunks", type=int, default=0)
    p_runtime.add_argument("--zero-j-threshold", type=float, default=1e-3)
    p_runtime.add_argument("--zero-compress", action="store_true")
    p_runtime.set_defaults(func=cmd_export_runtime)

    p_dense = sub.add_parser("export-runtime-dense")
    p_dense.add_argument("--db-dir", type=Path, default=DEFAULT_DB_DIR)
    p_dense.add_argument("--out-dir", type=Path, default=ROOT / "emission_runtime_db_v2")
    p_dense.add_argument("--dtype", choices=("float32", "float64"), default="float32")
    p_dense.add_argument("--format", choices=("npz", "binary", "both"), default="npz")
    p_dense.add_argument("--region", action="append", default=None)
    p_dense.add_argument("--zero-j-threshold", type=float, default=1e-3)
    p_dense.add_argument("--zero-compress", action="store_true")
    p_dense.set_defaults(func=cmd_export_runtime_dense)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
