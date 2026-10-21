#!/usr/bin/env python3
"""
Runs AnaPiPiProt.exe on every merged MC hipo file for one job, i.e. every
MCin/Job_<jobID>/MC_Rho_pipi_<jobID>_<index>.hipo file, up to MAX_PARALLEL_JOBS
at a time. As soon as one job finishes, the next queued file is started.

AnaPiPiProt.exe recognizes the MC_Rho_pipi_<jobID>_<index>.hipo naming
convention and writes each file's histograms to
Hists/MC/Job_<jobID>/Hists_PionPair_<index>.root.

Once every file has been processed successfully, the per-index histogram files
are hadd'ed into Hists/MC/Job_<jobID>/Hists_PionPair_All.root and then deleted,
leaving only the merged file.

Usage:
    python RunAnaPiPiProt_MC.py <jobID>

Example:
    python RunAnaPiPiProt_MC.py 11833
"""

import os
import sys
import glob
import time
import subprocess

# ── Configure paths here ──────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ANA_EXE = os.path.join(SCRIPT_DIR, "bin", "AnaPiPiProt.exe")
MCIN_DIR = os.path.join(SCRIPT_DIR, "MCin")
# ─────────────────────────────────────────────────────────────────────────────

MAX_PARALLEL_JOBS = 18
POLL_INTERVAL_SEC = 5
MERGED_HISTS_NAME = "Hists_PionPair_All.root"


def merge_and_clean(job_dir):
    pattern = os.path.join(job_dir, "Hists_PionPair_*.root")
    hist_files = sorted(
        f for f in glob.glob(pattern)
        if os.path.basename(f) != MERGED_HISTS_NAME
    )

    if not hist_files:
        print(f"No histogram files matching {pattern}, skipping hadd.")
        return

    print(f"\nMerging {len(hist_files)} histogram file(s) with hadd...")
    cmd = ["hadd", "-f", MERGED_HISTS_NAME] + [os.path.basename(f) for f in hist_files]
    result = subprocess.run(cmd, cwd=job_dir)
    if result.returncode != 0:
        print(f"WARNING: hadd exited with code {result.returncode}; leaving input histogram files in place.")
        return

    merged_file = os.path.join(job_dir, MERGED_HISTS_NAME)
    print(f"Merged histograms written to {merged_file}")

    for f in hist_files:
        os.remove(f)
    print(f"Removed {len(hist_files)} input histogram file(s).")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <jobID>")
        sys.exit(1)

    job_id = sys.argv[1]

    if not os.path.isfile(ANA_EXE):
        print(f"Error: AnaPiPiProt.exe not found at '{ANA_EXE}'.")
        sys.exit(1)

    job_dir = os.path.join(MCIN_DIR, f"Job_{job_id}")
    if not os.path.isdir(job_dir):
        print(f"Error: job directory not found at '{job_dir}'.")
        sys.exit(1)

    pattern = os.path.join(job_dir, f"MC_Rho_pipi_{job_id}_*.hipo")
    input_files = sorted(glob.glob(pattern))

    if not input_files:
        print(f"No MC hipo files matching {pattern}")
        sys.exit(1)

    print(f"Found {len(input_files)} MC hipo file(s) for job {job_id}.\n")

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
        print(f"{n_failed} job(s) failed. Skipping merge/cleanup.")
        sys.exit(1)

    hists_dir = os.path.join(SCRIPT_DIR, "Hists", "MC", f"Job_{job_id}")
    merge_and_clean(hists_dir)


if __name__ == "__main__":
    main()
