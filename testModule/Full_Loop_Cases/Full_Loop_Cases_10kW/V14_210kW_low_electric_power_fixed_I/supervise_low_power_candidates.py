"""Run and recover several isolated low-power candidates."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, TextIO


CASE_DIR = Path(__file__).resolve().parent
REPO_DIR = Path(__file__).resolve().parents[4]
RUNNER_MODULE = (
    "testModule.Full_Loop_Cases.Full_Loop_Cases_10kW."
    "V14_210kW_low_electric_power_fixed_I.run_v14_low_power_fixed_i"
)
PYTHON_EXE = Path(r"E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe")


@dataclass(frozen=True)
class Candidate:
    name: str
    output_dir: Path
    restart_in: Path
    duration_s: float
    dt_s: float = 0.05
    tec_update_interval_s: float = 0.5
    record_interval_s: float = 0.5
    hold_before_ramp_s: float = 20.0
    ramp_duration_s: float = 100.0
    ramp_shape: str = "cubic"
    final_power_w: float = 210000.0
    final_flow_kg_s: float = 2.46
    checkpoint_interval_s: float = 20.0
    staged_recording: bool = True

    @classmethod
    def from_dict(cls, raw: dict) -> "Candidate":
        values = dict(raw)
        for key in ("output_dir", "restart_in"):
            values[key] = Path(values[key])
        return cls(**values)


def build_command(candidate: Candidate, *, resume_from: Optional[Path] = None) -> list[str]:
    command = [
        str(PYTHON_EXE), "-m", RUNNER_MODULE,
        "--restart-in", str(candidate.restart_in),
        "--output-dir", str(candidate.output_dir),
        "--duration", str(candidate.duration_s),
        "--dt", str(candidate.dt_s),
        "--tec-update-interval", str(candidate.tec_update_interval_s),
        "--record-interval", str(candidate.record_interval_s),
        "--hold-before-ramp", str(candidate.hold_before_ramp_s),
        "--ramp-duration", str(candidate.ramp_duration_s),
        "--ramp-shape", candidate.ramp_shape,
        "--final-power", str(candidate.final_power_w),
        "--final-flow", str(candidate.final_flow_kg_s),
        "--checkpoint-interval", str(candidate.checkpoint_interval_s),
    ]
    if candidate.staged_recording:
        command.append("--staged-recording")
    if resume_from is not None:
        command.extend(("--resume-from", str(resume_from)))
    return command


def select_restart(output_dir: Path, *, exceptional_exit: bool) -> Optional[Path]:
    output = Path(output_dir)
    history = output / "history_control.csv"
    if not history.is_file():
        return None
    if exceptional_exit:
        failure = output / "failure.json"
        emergency = output / "emergency_restart.npz"
        if failure.is_file() and emergency.is_file():
            try:
                resumable = bool(json.loads(
                    failure.read_text(encoding="utf-8")).get("restart_resumable"))
            except (OSError, ValueError, TypeError):
                resumable = False
            if resumable:
                return emergency
    latest = output / "latest_restart.npz"
    return latest if latest.is_file() else None


def is_stalled(
    output_dir: Path, *, started_at_s: float, now_s: float,
    initialization_grace_s: float, stall_timeout_s: float,
    history_updated_at_s: Optional[float] = None,
) -> bool:
    history = Path(output_dir) / "history_control.csv"
    if not history.is_file():
        return now_s - started_at_s >= initialization_grace_s
    updated = history.stat().st_mtime if history_updated_at_s is None else history_updated_at_s
    return now_s - updated >= stall_timeout_s


def can_retry(retries: int, *, max_retries: int) -> bool:
    return int(retries) < int(max_retries)


def is_terminal_result(output_dir: Path) -> bool:
    output = Path(output_dir)
    return (output / "run_summary.json").is_file() or (output / "limit_trip.json").is_file()


@dataclass
class _Pending:
    candidate: Candidate
    retries: int = 0
    resume_from: Optional[Path] = None


@dataclass
class _Active:
    pending: _Pending
    process: subprocess.Popen
    started_at_s: float
    stdout: TextIO
    stderr: TextIO


def _start(item: _Pending) -> _Active:
    candidate = item.candidate
    candidate.output_dir.mkdir(parents=True, exist_ok=True)
    attempt = item.retries
    stdout = (candidate.output_dir / f"attempt_{attempt:02d}.stdout.log").open(
        "a", encoding="utf-8")
    stderr = (candidate.output_dir / f"attempt_{attempt:02d}.stderr.log").open(
        "a", encoding="utf-8")
    environment = os.environ.copy()
    for name in ("MKL_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS"):
        environment[name] = "1"
    try:
        process = subprocess.Popen(
            build_command(candidate, resume_from=item.resume_from),
            cwd=REPO_DIR, env=environment, stdout=stdout, stderr=stderr, text=True)
    except Exception:
        stdout.close()
        stderr.close()
        raise
    print(f"[supervisor] started {candidate.name} attempt={attempt} pid={process.pid}",
          flush=True)
    return _Active(item, process, time.time(), stdout, stderr)


def _stop(active: _Active) -> None:
    active.process.terminate()
    try:
        active.process.wait(timeout=10.0)
    except subprocess.TimeoutExpired:
        active.process.kill()
        active.process.wait(timeout=10.0)


def _close(active: _Active) -> None:
    active.stdout.close()
    active.stderr.close()


def supervise(
    candidates: Sequence[Candidate], *, max_parallel: int = 4,
    poll_interval_s: float = 30.0, initialization_grace_s: float = 900.0,
    stall_timeout_s: float = 1800.0, max_retries: int = 2,
) -> int:
    if not 1 <= int(max_parallel) <= 4:
        raise ValueError("max_parallel must be between 1 and 4")
    pending = [_Pending(candidate) for candidate in candidates]
    active: list[_Active] = []
    failures = 0
    while pending or active:
        while pending and len(active) < max_parallel:
            active.append(_start(pending.pop(0)))
        time.sleep(float(poll_interval_s))
        now = time.time()
        for job in list(active):
            returncode = job.process.poll()
            stalled = returncode is None and is_stalled(
                job.pending.candidate.output_dir,
                started_at_s=job.started_at_s, now_s=now,
                initialization_grace_s=initialization_grace_s,
                stall_timeout_s=stall_timeout_s,
            )
            if returncode is None and not stalled:
                continue
            if stalled:
                print(f"[supervisor] stalled {job.pending.candidate.name}", flush=True)
                _stop(job)
            _close(job)
            active.remove(job)
            output = job.pending.candidate.output_dir
            if is_terminal_result(output):
                print(f"[supervisor] terminal {job.pending.candidate.name}", flush=True)
                continue
            restart = select_restart(output, exceptional_exit=not stalled)
            retry_fresh = restart is None and not (output / "history_control.csv").exists()
            if can_retry(job.pending.retries, max_retries=max_retries) and (
                    restart is not None or retry_fresh):
                pending.append(_Pending(
                    job.pending.candidate, job.pending.retries + 1, restart))
                print(f"[supervisor] retrying {job.pending.candidate.name}", flush=True)
            else:
                failures += 1
                print(f"[supervisor] failed {job.pending.candidate.name}", flush=True)
    return 1 if failures else 0


def load_candidates(path: Path) -> list[Candidate]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return [Candidate.from_dict(item) for item in payload["candidates"]]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--max-parallel", type=int, default=4)
    parser.add_argument("--poll-interval", type=float, default=30.0)
    parser.add_argument("--initialization-grace", type=float, default=900.0)
    parser.add_argument("--stall-timeout", type=float, default=1800.0)
    parser.add_argument("--max-retries", type=int, default=2)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    candidates = load_candidates(args.candidates)
    if args.dry_run:
        for candidate in candidates:
            print(subprocess.list2cmdline(build_command(candidate)))
        return 0
    return supervise(
        candidates, max_parallel=args.max_parallel,
        poll_interval_s=args.poll_interval,
        initialization_grace_s=args.initialization_grace,
        stall_timeout_s=args.stall_timeout, max_retries=args.max_retries)


if __name__ == "__main__":
    raise SystemExit(main())
