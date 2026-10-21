#!/usr/bin/env python3
"""
For a given run number, find all rec_clas_<RUN>.evio.*.hipo files for that run
under the RG-A fall2018 pass2 recon directory and run bin/SkimTrigBit8.exe on
each of them, up to MAX_PARALLEL_JOBS at a time. As soon as one job finishes,
the next queued file is started, keeping up to MAX_PARALLEL_JOBS running.

Usage:
    python RunSkimTrigBit8.py <RUN>

Example:
    python RunSkimTrigBit8.py 5163
"""

import os
import sys
import glob
import time
import subprocess

# ── Configure paths here ──────────────────────────────────────────────────────
RECON_DIR = "/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/dst/recon"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SKIM_EXE = os.path.join(SCRIPT_DIR, "bin", "SkimTrigBit8.exe")
# ─────────────────────────────────────────────────────────────────────────────

MAX_PARALLEL_JOBS = 18
POLL_INTERVAL_SEC = 5


def main():
    if len(sys.argv) != 2 or not sys.argv[1].isdigit():
        print(f"Usage: python {sys.argv[0]} <RUN>")
        print("  RUN must be an integer (e.g. 5163)")
        sys.exit(1)

    run6 = f"{int(sys.argv[1]):06d}"

    if not os.path.isfile(SKIM_EXE):
        print(f"Error: SkimTrigBit8.exe not found at '{SKIM_EXE}'.")
        sys.exit(1)

    run_dir = os.path.join(RECON_DIR, run6)
    if not os.path.isdir(run_dir):
        print(f"Error: '{run_dir}' is not a valid directory.")
        sys.exit(1)

    pattern = os.path.join(run_dir, f"rec_clas_{run6}.evio.*.hipo")
    input_files = sorted(glob.glob(pattern))

    if not input_files:
        print(f"No files found matching: {pattern}")
        sys.exit(1)

    print(f"Found {len(input_files)} file(s) for run {run6}.\n")

    running = {} # inputFile -> Popen
    queue = list(input_files)
    n_failed = 0

    while queue or running:

        # Start new jobs while there is room and files left to process
        while queue and len(running) < MAX_PARALLEL_JOBS:
            input_file = queue.pop(0)
            print(f"Starting: {os.path.basename(input_file)}  ({len(running) + 1}/{MAX_PARALLEL_JOBS} running)")
            running[input_file] = subprocess.Popen([SKIM_EXE, input_file], cwd=SCRIPT_DIR)

        # Check for jobs that finished, freeing up slots for the next ones
        finished = [f for f, proc in running.items() if proc.poll() is not None]

        for input_file in finished:
            proc = running.pop(input_file)
            if proc.returncode != 0:
                n_failed += 1
                print(f"WARNING: {os.path.basename(input_file)} exited with code {proc.returncode}")
            else:
                print(f"Finished: {os.path.basename(input_file)}")

        if not finished and (queue or running):
            time.sleep(POLL_INTERVAL_SEC)

    print(f"\nAll {len(input_files)} file(s) processed for run {run6}.")
    if n_failed > 0:
        print(f"{n_failed} job(s) failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
