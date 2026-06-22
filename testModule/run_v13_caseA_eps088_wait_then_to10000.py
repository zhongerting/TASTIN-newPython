import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def is_process_running(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        import psutil  # type: ignore

        return psutil.pid_exists(pid)
    except Exception:
        if os.name == "nt":
            result = subprocess.run(
                ["tasklist", "/FI", f"PID eq {pid}"],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                check=False,
            )
            return str(pid) in result.stdout
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Wait for the active V13 eps088 run, then continue to about 10000 s."
    )
    parser.add_argument("--wait-pid", type=int, required=True)
    parser.add_argument(
        "--restart-in",
        default=(
            "testModule/v13_caseA_closed_loop_eps088_long3000s/"
            "v13_caseA_closed_loop_eps088_long3000s_latest_restart.npz"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="testModule/v13_caseA_closed_loop_eps088_to10000s",
    )
    parser.add_argument("--duration", type=float, default=4500.0)
    parser.add_argument("--poll-interval", type=float, default=60.0)
    args = parser.parse_args()

    while is_process_running(args.wait_pid):
        time.sleep(max(5.0, float(args.poll_interval)))

    python = sys.executable
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        python,
        "testModule/run_v13_caseA_closed_loop.py",
        "--restart-in",
        args.restart_in,
        "--tube-emissivity",
        "0.88",
        "--fin-emissivity",
        "0.88",
        "--duration",
        str(float(args.duration)),
        "--record-interval",
        "500",
        "--restart-interval",
        "500",
        "--max-dt",
        "0.5",
        "--output-dir",
        str(output_dir),
        "--case-prefix",
        "v13_caseA_closed_loop_eps088_to10000s",
    ]
    with (output_dir / "run.out").open("ab", buffering=0) as stdout, (
        output_dir / "run.err"
    ).open("ab", buffering=0) as stderr:
        return subprocess.call(cmd, cwd=os.getcwd(), stdout=stdout, stderr=stderr)


if __name__ == "__main__":
    raise SystemExit(main())
