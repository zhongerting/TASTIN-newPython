"""Run the 15 requested V14/V15 pre-start temperature cases."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


TEMPERATURES_K = (340, 330, 320, 310, 300)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--max-parallel", type=int, default=6)
    args = parser.parse_args()
    if args.max_parallel < 1:
        parser.error("--max-parallel must be positive")

    runner = Path(__file__).with_name("run_prestart_cooldown.py")
    tasks = [
        (case, target, temperature)
        for temperature in TEMPERATURES_K
        for case, target in (("v14", "radiator"), ("v15", "radiator"), ("v15", "shield"))
    ]
    args.output_root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "duration_s": 12000.0,
        "stop_temperature_k": 260.55,
        "external_heat_period_s": 6552.0,
        "direct_external_heat_scale_factor": 1.0,
        "shield_outer_surface_absorptivity": 0.1,
        "max_dt_s": 0.1,
        "record_interval_s": 1.0,
        "restart_interval_s": 60.0,
        "max_parallel": args.max_parallel,
        "cases": [],
    }
    for case, target, temperature in tasks:
        name = f"{case}_{target}_{temperature}K"
        manifest["cases"].append({
            "name": name,
            "case": case,
            "external_heat_target": target,
            "initial_temperature_k": temperature,
            "output_dir": str(args.output_root / name),
            "status": "pending",
        })
    status_path = args.output_root / "batch_status.json"
    status_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    def run_one(item):
        case, target, temperature = item
        name = f"{case}_{target}_{temperature}K"
        output_dir = args.output_root / name
        output_dir.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable,
            str(runner),
            case,
            "--output-dir", str(output_dir),
            "--duration", "12000",
            "--max-dt", "0.1",
            "--record-interval", "1",
            "--restart-interval", "60",
            "--external-heat-target", target,
            "--stop-temperature", "260.55",
            "--initial-temperature", str(temperature),
        ]
        with (output_dir / "run.out").open("w", encoding="utf-8") as stdout, (
            output_dir / "run.err"
        ).open("w", encoding="utf-8") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        return name, result.returncode

    with ThreadPoolExecutor(max_workers=args.max_parallel) as executor:
        futures = {executor.submit(run_one, item): item for item in tasks}
        for future in as_completed(futures):
            name, returncode = future.result()
            for item in manifest["cases"]:
                if item["name"] == name:
                    item["status"] = "complete" if returncode == 0 else "failed"
                    item["returncode"] = returncode
                    break
            status_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            print(f"{name}: returncode={returncode}", flush=True)
    return 0 if all(item["status"] == "complete" for item in manifest["cases"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
