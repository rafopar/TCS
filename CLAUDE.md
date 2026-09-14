# TCS

Codes and documents for the Timelike Compton Scattering (TCS) experiment analysis.

## Build system

The project is built with CMake (`cmake_minimum_required` 3.16), styled after the sibling
`../uRWellTestProto` project.

```
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<install dir> ..
make
make install
```

Top-level `CMakeLists.txt` finds ROOT, LZ4, Hipo4, and CLAS12AnaTools, then adds
`Eb_Change` and `RhoTail` as subdirectories, each with its own `CMakeLists.txt`.

External dependencies are located via custom find modules in `cmake_modules/`
(`FindLZ4.cmake`, `FindHipo4.cmake`, `FindCLAS12AnaTools.cmake`), which hardcode the
paths already used by `Eb_Change/compile_anaTCS.sh`:
- Hipo4: `/home/rafopar/work/git/hipo` (headers in `hipo4/`, libs in `lib/`)
- LZ4: header in `/home/rafopar/work/git/hipo/lz4/lib`, lib in `/home/rafopar/work/git/hipo/lib`
- CLAS12AnaTools: `/home/rafopar/work/git/clas12AnaTools` (headers in `include/`, lib in `lib/`)

These paths are specific to this machine/user, not portable — update the find modules if
the dependency locations change.

**Known caveat #1:** `libclas12AnaTools.so` was built with an unresolved `TSpectrum` symbol
(same failure reproduces with the original `compile_anaTCS.sh`, so it's a pre-existing
issue in that sibling project, not something introduced here). `Eb_Change/CMakeLists.txt`
works around it by explicitly linking `ROOT::Spectrum` into `anaTCS`.

**Known caveat #2 — hipo4 ABI mismatch:** `clas12AnaTools` vendors its own copy of hipo4
(`clas12AnaTools/hipo/hipo4`), whose `hipo::bank`/`structure` layout has diverged from the
standalone `/home/rafopar/work/git/hipo` checkout used elsewhere in this project (different
struct members/offsets — confirmed by diffing `bank.h` between the two). `RecEvent`/`RecPart`
and `MCEvent`/`MCPart` (clas12AnaTools's higher-level particle wrappers, used by
`RhoTail/AnaPiPiProt.cc` for reconstructed and generated particles respectively) embed
`hipo::bank` objects **by value**, so building against the wrong hipo4 corrupts memory at
runtime (segfault/`std::length_error` from `hipo::bank`'s internals — not a link error, so it
isn't caught until you actually run it; this reproduced exactly as predicted the first time
`MCEvent` was smoke-tested against a real file, compiled against the standalone hipo4 headers
by mistake). Any target using `RecEvent`/`RecPart`/`MCEvent`/`MCPart` must build against
clas12AnaTools's own vendored `hipo4`/`libHipo4.a` instead of the standalone one — see how
`RhoTail/CMakeLists.txt` gives the `AnaPiPiProt` target its own `target_include_directories`/
`target_link_libraries` pointing at `clas12AnaTools/hipo/hipo4` and `clas12AnaTools/lib/libHipo4.a`,
separate from the plain `Hipo4_INCLUDE_DIRS`/`Hipo4_LIBRARY` used by `SkimTrigBit8`/`Skim_empipi`.
Targets that only use the raw hipo4 reader/writer/bank API directly (no clas12AnaTools
`RecEvent`/`RecPart`/`MCEvent`/`MCPart`) are unaffected and should keep using the standalone hipo4.

`MCPart` (clas12AnaTools, `include/MCPart.h`) is `RecPart`'s generated-particle counterpart:
pid, momentum components/magnitude, theta/phi, vertex position/time, and mass, with **no**
detector response (there isn't one for a generated particle). `MCEvent` loads it from
`MC::Particle` (fields: `pid, px, py, pz, vx, vy, vz, vt` — no `p`/angles/mass, and unlike
`RecPart`'s φ, no CLAS12-sector `+30°` offset is applied, since that convention is specific to
reconstructed-track sectors and meaningless for generator truth). Mass isn't a
`MC::Particle` field, so it's looked up from pid via ROOT's `TDatabasePDG` (needs `ROOT::EG`
linked — not pulled in by plain `ROOT_LIBRARIES`, same as the `ROOT::Spectrum` workaround
above); this covers the standard leptons/mesons/baryons this experiment generates (e⁻, γ,
π⁰/π±, K±, p, n, Λ⁰, ...) but not nuclear/ion PDG codes (e.g. deuteron), which come back as
mass 0 with a one-time warning per unrecognized pid. Constructing `MCEvent` on a real-data
file (no `MC::Particle` bank actually written) is safe and intentional — `MC::Particle`'s
schema is still present in the dictionary (`hasSchema` true even with 0 rows written), so
`getRows()` just comes back 0 and `Particles()` is empty rather than erroring.

Also, the *installed* `libclas12AnaTools.so` at the flat `clas12AnaTools/lib/` location (the
one the find module and old `compile_anaTCS.sh` actually link against) can go stale relative
to `clas12AnaTools`'s own `include/` headers if that sibling repo has uncommitted/unbuilt
changes — its real CMake install prefix is `/home/rafopar/work/builds/CLAS12AnaTools/`, a
*third* location distinct from both the flat `lib/` dir and its in-tree `build/` dir. If a
build fails with `undefined reference` to a clas12AnaTools symbol that exists in the headers,
check whether `clas12AnaTools/build/` has a newer `.so` that never got copied/reinstalled to
the flat `lib/` location (this happened once for `RecEvent`; fixed by rebuilding
`clas12AnaTools` and copying the refreshed `.so` into the flat dir, including a matching
`SONAME` symlink since the two locations use different casing:
`libclas12AnaTools.so` vs. `libCLAS12AnaTools.so`). Repeat this rebuild-and-copy step any time
`clas12AnaTools/include/RecPart.h`, `RecEvent.{h,cc}`, `MCPart.h`, or `MCEvent.{h,cc}` change
(e.g. this is how `REC::ForwardTagger` support — `RecPart::FTCal()`/`FTHodo()`/`FTTRK()`,
detector IDs 10/11/13 — got added, since `AnaPiPiProt.cc` now uses `FTCal()` for
forward-tagger-electron vertex time).

### Install layout

Each subdirectory installs into its own mirrored path under `CMAKE_INSTALL_PREFIX`, e.g.
`<prefix>/RhoTail/bin/`, `<prefix>/Eb_Change/bin/`. Only executables (and ROOT macros with
no `main()`, following the same convention `uRWellTestProto` uses for non-compiled scripts —
e.g. `DrawPlots.cc`, `RhoFitGUI.cc`) go in `bin/` or the subdirectory root; data/output
directories are installed as siblings of `bin/`, not inside it (e.g. `<prefix>/RhoTail/OutSkims/`,
`<prefix>/RhoTail/MCin/`, `<prefix>/RhoTail/Hists/`, `<prefix>/RhoTail/Figs/`,
`<prefix>/Eb_Change/Figs/`, `<prefix>/Eb_Change/Data/`).

The user's real deployed install lives at `/local/work/builds/TCS` (also reachable as
`/work/clas12/rafopar/builds/TCS`), separate from this git checkout — that's where real skim
output and merged histograms actually accumulate, not under this source tree.

## Directory layout

- `Eb_Change/` — beam-energy-dependence studies. `anaTCS.cc` builds to `anaTCS.exe`;
  `DrawPlots.cc` is a ROOT macro (no `main()`, run as `root DrawPlots.cc`, not compiled).
- `RhoTail/` — trigger/PID skimming, π⁺π⁻p analysis, and interactive mass-fitting tools for
  the rho-tail analysis. See below.
- `cmake_modules/` — custom `Find*.cmake` modules for LZ4, Hipo4, CLAS12AnaTools.
- `miscCodes/`, `nbproject/`, `Docs/` — misc scripts, NetBeans project metadata, docs.

## RhoTail/

Pipeline, in order: `SkimTrigBit8.exe` (or `Skim_empipi.exe`) skims raw recon files →
merged per-run skim file → `AnaPiPiProt.exe` selects π⁺π⁻p events and fills histograms →
`hadd`-merged into one histogram file → `RhoFitGUI.cc` fits the π⁻π⁺ invariant mass
interactively. `CompareHists.cc` overlays two histogram files (e.g. data vs. MC) for QA at
any stage. `RunSkimTrigBit8.py` / `RunAllRuns.py` / `RunAnaPiPiProt.py` (data) and
`RunAnaPiPiProt_MC.py` (simulation) drive the skim/analysis steps in parallel batches.

### SkimTrigBit8.cc

Skims hipo files, keeping only events that satisfy **both**:
- trigger bit 8 (counting from bit 0) is set in `RUN::config.trigger` (an `int64`/`L`
  field — must be read with `bank::getLong`, not `getInt`), **and**
- at least one `REC::Particle` row with `pid == 211` and at least one with `pid == -211`.

Usage: `SkimTrigBit8.exe <input1.hipo> [input2.hipo ...]` — each input file is skimmed
independently into its own output file.

The run number and file index range are deduced from the input file's basename via regex
on the pattern `<run>.evio.<low>-<high>.hipo`, e.g.
`rec_clas_005163.evio.00055-00059.hipo` → run `005163`, low `00055`, high `00059`.

Output is written to `OutSkims/<run>/PionPairSkim_<run>_<low>_<high>.hipo` (directory and
run-number segment zero-padded to match the input filename; created relative to the working
directory the executable is run from), preserving the full original event content (whole
dictionary copied to the writer, `writer.addEvent(event)` on the untouched event object —
not a rebuilt subset of banks).

### Skim_empipi.cc

Same skim-file mechanics as `SkimTrigBit8.cc`, but no trigger-bit requirement — instead
requires an e⁻ (`pid==11`), a π⁻ (`pid==-211`) and a π⁺ (`pid==211`) each with
`abs(status) ∈ [2000,4000)` (forward-detector status cut). Output prefix `empipiSkim_`.

### RunSkimTrigBit8.py / RunAllRuns.py

`RunSkimTrigBit8.py <run>`: for one run, globs `rec_clas_<run6>.evio.*.hipo` under
`/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/dst/recon/<run6>/` and runs
`bin/SkimTrigBit8.exe` on each file with a sliding-window pool capped at 18 concurrent jobs
(refills a slot as soon as one finishes, rather than waiting for a whole batch).

`RunAllRuns.py`: reads run numbers from `F18_in_RunList.dat` and calls `RunSkimTrigBit8.py
<run>` once per run, **sequentially** — each call already blocks until that run's whole
parallel pool finishes, so no extra synchronization is needed to guarantee one run's jobs
complete before the next run starts.

### AnaPiPiProt.cc

Analyzes a (skimmed) hipo file and selects events with exactly one π⁻, one π⁺, and one
proton, using `RecEvent`/`RecPart` (clas12AnaTools) — see the hipo4 ABI caveat above for the
build-time gotcha this implies. Base selection/kinematics follow
`~/work/git/DDVCS/elDDVCS/AnaElDDVCS.cc`'s pattern (proton `chi2pid` cut, `vt` from the
highest-priority valid scintillator response with beta from the particle's mass); this has
since grown substantially and now also applies per-species `vt`/`vz` window cuts, computes
`Mx2` (missing mass² off the beam+target−π⁻π⁺p system), `Q2`, `Eg`, and `tM`, and separately
identifies a forward-tagged scattered electron candidate (`pid==11 && abs(status) < 2000`,
i.e. *not* the main forward-detector status range used for the pions) via `RecPart::FTCal()`
(added to clas12AnaTools specifically for this — see caveat #2) — its vertex time is computed
from the FTCal cluster position/time and the proton's vz/vt, and events with a valid in-time
FT electron get their own parallel set of `_HasFT` histograms. `h_Minv_pippim_Mx2_Q2_cuts` is
the histogram normally used for fitting (`Mx2` and `Q2` cuts both applied). Check the file
directly for current cut values/histogram names — it changes often.

Also fills `h_MC_Minv_pippim_All`/`_RecSel`/`_Mx2Cut`: the MC-truth (generated, via `MCEvent` —
see caveat #2) π⁺π⁻ invariant mass, filled from the *same* per-event value at three points of
the reconstructed-event selection cascade (every event regardless of any `RecPart` found, once
the π⁻/π⁺/proton selection passes, once that selection also passes the `Mx2` cut) — a
reconstruction/selection-efficiency-vs-generated-mass shape, not three different mass
calculations. Real data leaves all three empty (no `MC::Particle` rows to find a generated
pair in).

Two output layouts, chosen by input filename:
- Normal data: `Hists/Hists_<input basename>.root`.
- MC, recognized by the naming convention `MC_Rho_pipi_<jobID>_<index>.hipo` (e.g.
  `MC_Rho_pipi_11833_11.hipo`): `Hists/MC/Job_<jobID>/Hists_PionPair_<index>.root` instead.

Usage: `AnaPiPiProt.exe <input.hipo>`.

### RunAnaPiPiProt.py / RunAnaPiPiProt_MC.py

`RunAnaPiPiProt.py`: runs `bin/AnaPiPiProt.exe` over every **data** skim file for every run in
`F18_in_RunList.dat`, pooled into a **single** flat queue across all runs (capped at 18
concurrent — unlike the skim step, `AnaPiPiProt.exe` jobs don't depend on each other, so
there's no need to finish one run before starting the next). After all jobs finish, `hadd -f`s
every `Hists/Hists_PionPairSkim_*_All.root` into `Hists/Hists_PionPairSkim_All.root` (run from
inside `Hists/`, so plain filenames work).

Skim files in practice are **not** always laid out the way `SkimTrigBit8.cc`'s own naming
convention would suggest — the real deployed `OutSkims/` has flat, **unpadded**-run-number
merged files like `OutSkims/PionPairSkim_5343_All.hipo` (from a separate manual `hadd` step
the user runs), not nested zero-padded per-file skims. `find_skim_files()` therefore searches
all four combinations (flat/nested × padded/unpadded run number) per run, and results are
deduplicated by output **basename** (not full path) before queueing — the real `OutSkims/`
was found to contain one straight-up duplicate (run 5163, both a symlink into a different
tree and a separate regular file with identical content under a stray subdirectory), which
would otherwise have raced to write the same output filename.

`RunAnaPiPiProt_MC.py <jobID>`: the MC counterpart. Globs `MCin/Job_<jobID>/MC_Rho_pipi_<jobID>_*.hipo`
(same 18-way parallel pool), then `hadd`s the resulting `Hists/MC/Job_<jobID>/Hists_PionPair_*.root`
files into `Hists_PionPair_All.root` **and deletes the per-index inputs afterward** (only on a
successful merge) — unlike the data driver, which leaves per-run inputs in place.

### CompareHists.cc

QA tool: overlays every same-named `TH1`-derived histogram found in two root files (e.g. a
data merge vs. an MC merge) onto one multi-page pdf — 1D histograms normalized to unit-max
and overlaid directly (file1 red, file2 blue); 2D histograms get a 2×2 page with the raw
`COLZ` 2D histogram from each file on the left and normalized X/Y projections overlaid on the
right. Histograms present in file1 but missing from file2 are skipped with a warning, not
silently dropped.

Titles/axis labels and cut-value reference lines aren't read from the histograms themselves
(`AnaPiPiProt.cc` creates every histogram with an empty title) — `CompareHists.cc` keeps its
own `name -> "Title;XTitle;YTitle"` lookup (`HistTitles`) and per-name cut-value lookups
(`VCuts`/`HCuts`), hand-matched against `AnaPiPiProt.cc`'s current selection logic and cut
constants, and falls back to the histogram's own name if it isn't in the table. Cuts are drawn
as dashed black lines: vertical for a 1D histogram's own cut(s) or a 2D histogram's X-variable
cut(s) (`VCuts`), horizontal for a 2D histogram's Y-variable cut(s) (`HCuts`) — and a 2D
histogram's X/Y projections each inherit whichever cut set matches their own axis (e.g. the
Y-projection gets `HCuts` drawn as *vertical* lines, since that projection's X axis *is* the
2D's Y variable). Since none of this is derived from the histograms, **it goes stale silently
if `AnaPiPiProt.cc`'s histogram definitions or cut constants change without a matching update
here** — worth checking both files together.

Usage: `CompareHists.exe <file1.root> <file2.root> <keyword1> <keyword2>` → writes
`./Figs/RhoTail_Comparisons_<keyword1>_<keyword2>.pdf`.

### RhoFitGUI.cc

Interactive ROOT GUI (`root -l RhoFitGUI.cc` or `RhoFitGUI.cc+`; no `main()`, not compiled by
CMake — installed as a plain file like `DrawPlots.cc`) to fit the π⁻π⁺ invariant mass with the
sum of three relativistic Breit-Wigners: ρ(770) (P-wave, L=1), f2(1270) (D-wave, L=2), and
f0(1500) (S-wave, L=0), using the Jackson mass-dependent-width form. Lets the user pick which
histogram to fit from a combo box (defaults to `h_Minv_pippim_Mx2_Q2_cuts` if present), set
separate fit/draw ranges, fix or range-limit any of the 9 (+1 optional resolution σ) shape
parameters, optionally convolve with a Gaussian (via `TF1Convolution`) for mass resolution,
and save the plot as pdf/png/root. Draws the total (blue), ρ (red), f2 (green), and f0
(magenta) components separately with a legend, `SetNpx(4500)` for smoothness.

Two non-obvious bugs were found and fixed here, both worth knowing if this file is touched again:
- **FFT-convolution range bug:** `TF1Convolution`'s FFT grid is fixed at construction time and
  is inherently periodic over that range. The total fit function used to be built over the
  narrow *fit* range and then widened via `SetRange()` for display — which only changes what
  range the outer `TF1` reports, not the underlying FFT grid, so evaluating past the original
  edge wrapped around and picked up the function's behavior from the *other* end of the range
  (visible as ρ's sharp near-threshold rise "leaking" into the high-mass tail). Fixed by always
  building over the wide *draw* range from the start, narrowing to the fit range only
  temporarily (via `SetRange()`, safe — doesn't touch the FFT grid) for the actual `Fit()` call.
- **Parameter-bounds-vs-Minuit gotcha:** the "value" widgets get overwritten with each fit's
  result, so a parameter that saturates at its `SetParLimits` boundary (which happens
  routinely — Minuit correctly reports it "(limited)" right at the edge, not a bug by itself)
  becomes the *starting* value for the next fit. Starting a bounded Minuit fit exactly on its
  boundary is numerically unstable (the internal arcsin-based transform has an infinite
  derivative there). Fixed by nudging the starting value strictly inside `(min, max)` before
  every fit, and — since this couldn't be reliably reproduced in synthetic tests, so the exact
  trigger condition on real data is unconfirmed — by also adding an explicit post-fit check
  that flags in the status line if any free parameter's fitted value ends up outside its
  declared range, rather than silently plotting an inconsistent curve.

## Testing notes

This machine has no X server or Xvfb, so `RhoFitGUI.cc`'s actual widget rendering/interaction
has never been visually verified in this environment — only that it compiles cleanly via
ACLiC. Its underlying physics (RBW lineshapes, component-sum consistency, convolution
behavior) has been independently verified numerically in standalone (non-GUI) scripts each
time it changed. When testing any of these tools with real or synthetic data, prefer
`/volatile/clas12/rafopar/tmp/` for scratch files of non-trivial size — see
[[feedback-tmp-scratch-location]] in memory for why (filling the default session scratchpad
once broke the Bash tool entirely, host-wide).
