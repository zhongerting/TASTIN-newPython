"""Extract same-phase V14 radiator temperatures directly from NPZ checkpoints."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np


ORBITAL_PERIOD_S = 5668.144369
RING_WALL_RE = re.compile(r"^Solid_(Upper|Lower)_A([1-6])_[^/]+_Solid/T$")
HEAT_PIPE_RE = re.compile(
    r"^Solid_(Upper|Lower)_A([1-6])_[^/]+_RingHP_HP_node([0-2])_HP_inner/T$"
)
UPPER_MULTIPLIERS = ((8, 9, 9),) * 5 + ((8, 8, 8),)
LOWER_MULTIPLIERS = ((10, 10, 11),) * 6
RADIAL_FACES_M = np.array([0.0075, 0.0080, 0.0085, 0.0090])
AXIAL_WIDTHS_M = np.array([0.0605, 0.0415] + [0.47 / 12.0] * 12)
SECTION_COLUMNS = {"evaporator": slice(0, 1), "adiabatic": slice(1, 2), "condenser": slice(2, 14)}


def orbital_phase(global_time_s: float, phase_origin_s: float) -> float:
    return float((global_time_s - phase_origin_s) % ORBITAL_PERIOD_S)


def _stats(values: np.ndarray, weights: np.ndarray) -> dict:
    if values.size == 0 or values.shape != weights.shape or not np.all(np.isfinite(values)):
        raise ValueError("temperature values and weights must be finite, non-empty, and aligned")
    return {
        "min": float(values.min()),
        "volume_weighted_mean": float(np.average(values, weights=weights)),
        "max": float(values.max()),
        "cell_count": int(values.size),
    }


def _schema(npz: np.lib.npyio.NpzFile) -> dict[str, tuple[int, ...]]:
    return {key: tuple(npz[key].shape) for key in npz.files}


def _multiplier(level: str, segment: int, node: int) -> int:
    table = UPPER_MULTIPLIERS if level == "Upper" else LOWER_MULTIPLIERS
    return table[segment - 1][node]


def _temperature_groups(npz: np.lib.npyio.NpzFile) -> tuple[dict, dict]:
    ring = {"all": ([], []), "upper": ([], []), "lower": ([], [])}
    heat_pipe = {
        level: {section: ([], []) for section in SECTION_COLUMNS}
        for level in ("all", "upper", "lower")
    }
    ring_keys = hp_keys = 0
    cell_volumes = np.pi * np.diff(RADIAL_FACES_M**2)[:, None] * AXIAL_WIDTHS_M[None, :]

    for key in npz.files:
        match = RING_WALL_RE.fullmatch(key)
        if match:
            level = match.group(1).lower()
            values = np.asarray(npz[key], dtype=float).reshape(-1)
            if values.shape != (3,):
                raise ValueError(f"unexpected collector-ring wall shape for {key}: {values.shape}")
            for group in ("all", level):
                ring[group][0].append(values)
                ring[group][1].append(np.ones(values.shape))
            ring_keys += 1
            continue

        match = HEAT_PIPE_RE.fullmatch(key)
        if not match:
            continue
        level_name, segment_text, node_text = match.groups()
        level = level_name.lower()
        values = np.asarray(npz[key], dtype=float)
        if values.shape != (42,):
            raise ValueError(f"unexpected heat-pipe shape for {key}: {values.shape}")
        values = values.reshape(3, 14)
        multiplier = _multiplier(level_name, int(segment_text), int(node_text))
        for section, columns in SECTION_COLUMNS.items():
            section_values = values[:, columns].reshape(-1)
            section_weights = (cell_volumes[:, columns] * multiplier).reshape(-1)
            for group in ("all", level):
                heat_pipe[group][section][0].append(section_values)
                heat_pipe[group][section][1].append(section_weights)
        hp_keys += 1

    if ring_keys != 12 or hp_keys != 36:
        raise ValueError(f"unexpected V14 topology: {ring_keys} ring walls and {hp_keys} heat pipes")

    ring_result = {
        group: _stats(np.concatenate(items[0]), np.concatenate(items[1]))
        for group, items in ring.items()
    }
    hp_result = {
        group: {
            section: _stats(np.concatenate(items[0]), np.concatenate(items[1]))
            for section, items in sections.items()
        }
        for group, sections in heat_pipe.items()
    }
    return ring_result, hp_result


def extract_directory(directory: Path, phase_origin_s: float) -> dict:
    directory = Path(directory)
    checkpoints = sorted(directory.glob("*.npz"))
    if not checkpoints:
        raise FileNotFoundError(f"no NPZ checkpoints in {directory}")

    rows = []
    reference_schema = None
    for path in checkpoints:
        with np.load(path, allow_pickle=False) as npz:
            schema = _schema(npz)
            if reference_schema is None:
                reference_schema = schema
            elif schema != reference_schema:
                raise ValueError(f"checkpoint schema differs: {path.name}")
            if tuple(npz["System/global_time"].shape) != (1,):
                raise ValueError(f"invalid System/global_time shape in {path.name}")
            global_time_s = float(npz["System/global_time"][0])
            ring, heat_pipe = _temperature_groups(npz)
        phase_s = orbital_phase(global_time_s, phase_origin_s)
        rows.append({
            "file": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "global_time_s": global_time_s,
            "elapsed_from_phase_origin_s": global_time_s - phase_origin_s,
            "orbital_phase_s": phase_s,
            "orbital_phase_fraction": phase_s / ORBITAL_PERIOD_S,
            "collector_ring_wall_temperature_k": ring,
            "heat_pipe_temperature_k": heat_pipe,
        })

    rows.sort(key=lambda row: row["global_time_s"])
    return {
        "format_version": 1,
        "source_directory": str(directory.resolve()),
        "orbital_period_s": ORBITAL_PERIOD_S,
        "phase_origin_s": float(phase_origin_s),
        "checkpoint_count": len(rows),
        "schema_key_count": len(reference_schema),
        "temperature_weighting": {
            "collector_ring_wall": "equal V14 cylindrical cells",
            "heat_pipe": "cylindrical control-cell volume times physical heat-pipe multiplier",
        },
        "fin_temperature_k": {
            "available": False,
            "reason": "quasi-steady fin temperatures are not stored in these NPZ checkpoints",
        },
        "checkpoints": rows,
    }


def _default_phase_origin(directory: Path) -> float:
    config_path = directory / "run_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    return float(config["external_heat_time_origin_s"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint_directory", type=Path)
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--phase-origin-s", type=float)
    args = parser.parse_args()
    origin = args.phase_origin_s
    if origin is None:
        origin = _default_phase_origin(args.checkpoint_directory)
    result = extract_directory(args.checkpoint_directory, origin)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
