#!/usr/bin/env python3
"""
Runs RunSkimTrigBit8.py once for every run number listed in
F18_in_RunList.dat, one run at a time: all of a run's skim jobs (started in
parallel by RunSkimTrigBit8.py) must finish before the next run starts.

Usage:
    python RunAllRuns.py
"""

import os
import sys
import subprocess

# ── Configure paths here ──────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUNLIST_FILE = os.path.join(SCRIPT_DIR, "F18_in_RunList.dat")
RUN_SKIM_SCRIPT = os.path.join(SCRIPT_DIR, "RunSkimTrigBit8.py")
# ─────────────────────────────────────────────────────────────────────────────


def read_run_numbers(runlist_file):
    runs = []
    with open(runlist_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            runs.append(int(line.split()[0]))
    return runs


def main():
    if not os.path.isfile(RUNLIST_FILE):
        print(f"Error: run list not found at '{RUNLIST_FILE}'.")
        sys.exit(1)

    if not os.path.isfile(RUN_SKIM_SCRIPT):
        print(f"Error: RunSkimTrigBit8.py not found at '{RUN_SKIM_SCRIPT}'.")
        sys.exit(1)

    runs = read_run_numbers(RUNLIST_FILE)
    print(f"Found {len(runs)} run(s) in {RUNLIST_FILE}.\n")

    n_failed = 0

    for i, run in enumerate(runs, start=1):
        print(f"=== [{i}/{len(runs)}] Run {run} ===")
        result = subprocess.run(["python3", RUN_SKIM_SCRIPT, str(run)], cwd=SCRIPT_DIR)
        if result.returncode != 0:
            n_failed += 1
            print(f"WARNING: run {run} exited with code {result.returncode}")

    print(f"\nAll {len(runs)} run(s) processed.")
    if n_failed > 0:
        print(f"{n_failed} run(s) failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
