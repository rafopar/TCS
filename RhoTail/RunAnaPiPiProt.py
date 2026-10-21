#!/usr/bin/env python3
"""
Runs AnaPiPiProt.exe on every skimmed hipo file found for each run listed in
F18_in_RunList.dat, up to MAX_PARALLEL_JOBS at a time. As soon as one job
finishes, the next queued file is started, keeping up to MAX_PARALLEL_JOBS
running. Unlike RunAllRuns.py/RunSkimTrigBit8.py, jobs for different runs do
not need to wait on one another -- AnaPiPiProt.exe jobs are independent of
each other, so all files across all runs are pooled into a single queue.

Skim files may be laid out either flat in OutSkims/ (e.g. a merged
PionPairSkim_5343_All.hipo) or nested under a per-run subdirectory (e.g.
OutSkims/005343/PionPairSkim_005343_00055_00059.hipo, one per original input
file), and the run number in the file name may or may not be zero-padded.
All of these are searched for each run.

Usage:
    python RunAnaPiPiProt.py
"""

import os
import sys
import glob
import time
import subprocess

# ── Configure paths here ──────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUNLIST_FILE = os.path.join(SCRIPT_DIR, "F18_in_RunList.dat")
ANA_EXE = os.path.join(SCRIPT_DIR, "bin", "AnaPiPiProt.exe")
OUTSKIMS_DIR = os.path.join(SCRIPT_DIR, "OutSkims")
HISTS_DIR = os.path.join(SCRIPT_DIR, "Hists")
MERGED_HISTS_FILE = "Hists_PionPairSkim_All.root"
# ─────────────────────────────────────────────────────────────────────────────

MAX_PARALLEL_JOBS = 18
POLL_INTERVAL_SEC = 5


def read_run_numbers(runlist_file):
    runs = []
    with open(runlist_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            runs.append(int(line.split()[0]))
    return runs


def find_skim_files(run):
    run6 = f"{run:06d}"
    patterns = [
        os.path.join(OUTSKIMS_DIR, f"PionPairSkim_{run}_*.hipo"),
        os.path.join(OUTSKIMS_DIR, f"PionPairSkim_{run6}_*.hipo"),
        os.path.join(OUTSKIMS_DIR, str(run), "PionPairSkim_*.hipo"),
        os.path.join(OUTSKIMS_DIR, run6, "PionPairSkim_*.hipo"),
    ]
    files = set()
    for pattern in patterns:
        files.update(glob.glob(pattern))
    return sorted(files)


def merge_histograms():
    if not os.path.isdir(HISTS_DIR):
        print(f"No Hists directory at '{HISTS_DIR}', skipping hadd.")
        return

    pattern = os.path.join(HISTS_DIR, "Hists_PionPairSkim_*_All.root")
    hist_files = sorted(glob.glob(pattern))

    if not hist_files:
        print(f"No histogram files matching {pattern}, skipping hadd.")
        return

    print(f"\nMerging {len(hist_files)} histogram file(s) with hadd...")
    cmd = ["hadd", "-f", MERGED_HISTS_FILE] + [os.path.basename(f) for f in hist_files]
    result = subprocess.run(cmd, cwd=HISTS_DIR)
    if result.returncode != 0:
        print(f"WARNING: hadd exited with code {result.returncode}")
    else:
        print(f"Merged histograms written to {os.path.join(HISTS_DIR, MERGED_HISTS_FILE)}")


def main():
    if not os.path.isfile(RUNLIST_FILE):
        print(f"Error: run list not found at '{RUNLIST_FILE}'.")
        sys.exit(1)

    if not os.path.isfile(ANA_EXE):
        print(f"Error: AnaPiPiProt.exe not found at '{ANA_EXE}'.")
        sys.exit(1)

    runs = read_run_numbers(RUNLIST_FILE)

    # AnaPiPiProt.exe names its output purely from the input file's basename,
    # so two different paths sharing a basename (e.g. a stray duplicate copy
    # of a run's merged skim) would race to write the same output file if
    # both were queued. Keep only the first path seen per basename.
    input_files = []
    seen_basenames = set()
    for run in runs:
        for skim_file in find_skim_files(run):
            base = os.path.basename(skim_file)
            if base in seen_basenames:
                print(f"Skipping duplicate (same basename already queued): {skim_file}")
                continue
            seen_basenames.add(base)
            input_files.append(skim_file)

    if not input_files:
        print(f"No skim files found under {OUTSKIMS_DIR} for the runs in {RUNLIST_FILE}")
        sys.exit(1)

    print(f"Found {len(input_files)} skim file(s) across {len(runs)} run(s).\n")

    running = {}  # inputFile -> Popen
    queue = list(input_files)
    n_failed = 0

    while queue or running:

        # Start new jobs while there is room and files left to process
        while queue and len(running) < MAX_PARALLEL_JOBS:
            input_file = queue.pop(0)
            print(f"Starting: {os.path.basename(input_file)}  ({len(running) + 1}/{MAX_PARALLEL_JOBS} running)")
            running[input_file] = subprocess.Popen([ANA_EXE, input_file], cwd=SCRIPT_DIR)

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

    print(f"\nAll {len(input_files)} file(s) processed.")
    if n_failed > 0:
        print(f"{n_failed} job(s) failed.")

    merge_histograms()

    if n_failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
