# Developer notes

Details for the people who work on the code (Judah and Claude). The short, user-facing
introduction is in [README.md](README.md); blind-search specifics are in
[blind_search/DEVELOPER_NOTES.md](blind_search/DEVELOPER_NOTES.md). Decisions and results,
with dates, live in [PLAN.md](PLAN.md) and [QUESTIONS.md](QUESTIONS.md); working rules in
`CLAUDE.md`.

## Blind-search quick start in detail (moved from the README)

`blind_search/` finds companions blindly in Gaia DR4 epoch astrometry: it adds planets one
at a time while the data demand it, with a detection threshold built from simulated
datasets. It needs only Python with numpy, scipy, pandas and matplotlib (the full
environment below also covers the census and RV scripts; a TeX installation is needed only
for `blind_search/make_flowchart.py`). The algorithm and modules are described in
[blind_search/README.md](blind_search/README.md).

### 1. Install

Clone the repository (https://github.com/jvanzand/dr4_sensitivity). Then create and activate the conda environment (Python 3.14) from the repository root:

```bash
conda env create -f environment.yml
conda activate gaia_dr4_sensitivity
```

Check the installation (about 1.5 minutes; every line should print `PASS` and the last
line `all passed`):

```bash
python -I tests/test_blind_search.py --quick
```

### 2. Input data

The loader `data.load_gaia_dr4` reads Gaia epoch astrometry in the archive CSV layout, one
row per CCD crossing, with the columns `transit_id, epoch, scan_pos_angle,
centroid_pos_al, centroid_pos_error_al, ipd_error_al, parallax_factor_al,
used_by_agis_al, outlier_flag`. It drops unused rows and outliers, then reduces the
crossings to one measurement per transit (median centroid, bootstrap error, as in the
Octofitter DR4 tutorial). Two pre-release examples are included:
`data_conversion/GaiaDR4_data_Gaia4.csv` (Gaia-4, 93 transits) and
`data_conversion/GaiaDR4_data_HD114762.csv`.

### 3. Run a search

Save this as `my_search.py` in the repository root:

```python
import sys
sys.path.insert(0, 'blind_search')      # make the blind_search modules importable
import data
import model
import search


def main():
    channel = data.load_gaia_dr4('data_conversion/GaiaDR4_data_Gaia4.csv')
    stack = model.Stack([channel])      # the data plus the linear model
    result = search.run_search(
        stack,
        M=50,                 # noise-sampled datasets per threshold (default)
        fap=1e-3,             # false-alarm probability of each decision
        max_companions=4,     # stop after this many detections
        out_dir='my_results', # one summary figure per step is written here
        label='gaia4')
    for k, (P, e, M0) in enumerate(result['known'], start=1):
        print(f'companion {k}: P = {P:.1f} d, e = {e:.2f}, M0 = {M0:.2f} rad')


if __name__ == '__main__':             # required: the search uses worker processes
    main()
```

```bash
python -I my_search.py
```

The search uses all cores but one and took about 10 s per step on the 12-core machine
this was written on (each step draws M noise-sampled datasets). It prints one line per step,
for example (Gaia-4, M = 50):

```
step 1: highest peak 181.3 at 562.3 d vs threshold 24.6 -> DETECTED  [9 s]
step 2: highest peak -13.6 at 191.5 d vs threshold 11.4 -> not detected  [11 s]
companion 1: P = 577.4 d, e = 0.40, M0 = 5.75 rad
```

`M`, the number of noise-sampled datasets behind each threshold (datasets sampled from
the fitted model, so any signal in them arises from noise only), defaults to 50. On 2000 injections,
50 and 100 noise-sampled datasets gave the same completeness map (rms difference 0.005, maximum 0.029) and
changed the detection decision for only 15 injections, because the threshold of any one
dataset scatters by about 2 in Delta-BIC between the two settings. More noise-sampled datasets make single
thresholds smoother and cost proportionally more time.

Each line reads: at this step the highest peak of the Delta-BIC periodogram (a measure of
how much better the data fit with one more companion; larger is better) was X, at trial
period P; the detection threshold built from the noise-sampled datasets was T; the companion
is added only if X exceeds T. The search stops at the first step with no detection, so
the number of "DETECTED" lines is the number of companions found. The final periods,
eccentricities and mean anomalies are the maximum-likelihood refit of all detected
companions.

### 4. Read the step figures

`my_results/gaia4_step1.png`, `gaia4_step2.png`, ... (the repository keeps this pair for Gaia-4 as `blind_search/figures_search/gaia4_real_step1.png` and `gaia4_real_step2.png`) show one step each: the real
periodogram (blue) with 20 of the simulated ones (grey, their highest peaks as dots) and
the threshold (dashed orange line) on top; how the threshold is built (histogram of the
simulated highest peaks with an extreme-value fit, thresholds for false-alarm
probabilities 0.01 and 0.001) at bottom left; and the residuals to the model before and
after adding the detected companion at bottom right. The step with no detection shows the
highest remaining peak below its threshold.

### 5. Measure what the search can find (injection/recovery map)

To learn which companions a search of these data (the Gaia-4 transit times, scan angles
and errors) would detect, inject planets into noise and count recoveries. Each injection
has its own noise and its own simulated threshold, so the run is slow: about 33 s of one
core per injection with 50 noise-sampled datasets. Test the pieces, run 2000 injections (about 95 minutes
on 11 workers; the run can be stopped and resumed), then draw the map:

```bash
python -I blind_search/injection_map.py --test
python -I blind_search/injection_map.py run --n 2000 --M 50 --out blind_search/results/my_map.jsonl
python -I blind_search/injection_map.py plot --out blind_search/results/my_map.jsonl
```

The map (`blind_search/figures_search/my_map.png`) shows each injected planet (mass 0.5-20
Jupiter masses, semi-major axis 0.5-5 AU, circular orbits, host 0.644 solar masses at
13.58 mas) as recovered or not, with smoothed completeness contours at 0.1, 0.5 and 0.9.
"Recovered" means detected above the threshold with the fitted period and mass within 25%
of the injected values. To add a variant of the injections (another host, range or
detection rule), edit the constants at the top of `blind_search/injection_map.py`.

### Limits to keep in mind

The simulated noise is white (the real noise may not be), the map injects one circular
companion per star, and an undetected companion hidden by a poor fit of an earlier one is
a known failure mode. See the end of [blind_search/README.md](blind_search/README.md) and
`QUESTIONS.md` for the full list.

- [PLAN.md](PLAN.md): the science question, data-format notes, analysis steps, workflow,
  roadblocks, extras, and who does what.
- [QUESTIONS.md](QUESTIONS.md): open decisions, items to remember, and settled decisions.
- `sample/lit_search/`: the literature search for brown-dwarf and other RV companions
  ([README.md](sample/lit_search/README.md) has the executive summary and leads).
- `sample/`: the census of RV companions from the NEA and survey catalogs
  ([CENSUS.md](sample/CENSUS.md), `census_table.csv`, `alpha_vs_period.png`,
  `snr_vs_period_by_distance.png`, `provenance.png`) and the scripts that build it
  (`raw_lit/<source>/build.py`, `literature.py`, `census.py`, `census_plot.py`,
  `provenance_plot.py`).
- `environment.yml` and `requirements.txt`: the project's Python environment, the conda
  environment `gaia_dr4_sensitivity` (Python 3.14; set up 2026-11-04). Every Python
  script runs in it: `conda activate gaia_dr4_sensitivity`, then run scripts from this
  folder with `python -I`. Recreate it with `conda env create -f environment.yml`.
  (Before 2026-11-04 the census ran in Python 3.7; rebuilding it in 3.14 changed nothing
  beyond floating-point rounding, at most 1e-13 relative.)
- `matplotlibrc`: Judah's plotting style, loaded by the plotting scripts.
- `Project.toml`/`Manifest.toml`: the project's Julia environment, with Octofitter v9.0.0
  (upgraded 2026-11-04). Run Julia with `--project=.` from this folder. Loading it creates
  `.CondaPkg/` (about 100 MB), a Python environment that a dependency's PythonCall
  extension installs.
- `tests/test_scan_angle.jl`: checks that the converted DR4 files and Octofitter agree on
  the scan-angle unit (`julia --project=. tests/test_scan_angle.jl`).
- `blind_search/`: a standalone, Octofitter-free blind companion search for DR4 epoch
  astrometry (Delta-BIC periodogram, simulated detection threshold, injection/recovery map,
  reliability tests); [blind_search/README.md](blind_search/README.md) describes the algorithm,
  modules and results. `tests/test_blind_search.py` tests it.
- `rv_data/`: RV collection (PLAN.md, Step 2). The Jump exports (`jump_raw/`, `jump_rvs.csv`) are not in
  the public repository (kept on Judah's machine, listed in `.gitignore`); `rvs_merged.csv` and the system tables
  carry the Jump-sourced rows (access flagged public) with their provenance. `build_targets.py` (identities and
  aliases), `jump_queries.sql` (run by Judah in Jump's SQL Explorer; exports kept in
  `jump_raw/`), `ingest_jump.py`,
  `pull_dace.py`, and `merge_rvs.py` (one table per system in `rv_data/systems/`), with
  the duplicate rules in `rv_data/DUPLICATE_RULES.md`.
- `tests/test_duplicate_rules.py`: one synthetic case per duplicate rule
  (`python -I tests/test_duplicate_rules.py`).
- `sample/catalog.py`: the stage-1 snapshot catalog and stage-2 fit set built from the census
  (`stage1_catalog.csv`, `stage1_systems.csv`, `resonance_candidates.csv`,
  `stage1_resonance.png`, `missing_errors.csv`); the census now carries error columns and
  precision flags (PLAN.md, Step 1b).
- `sample/system_names.py`, `system_names.csv`: readable names for systems (hierarchy HD, GJ,
  TOI, Kepler, K2, WASP, HAT-P, HIP, then SIMBAD's main id), used for figures and result
  folders; `fits_rv/kepler.c` is the compiled Kepler solver.
- `fits_rv/`, `configs/`: the RV companion search (`companion_search.py`, `--selftest`) and the
  Octofitter v9 fit pipeline (`pipeline.jl`, `run_fit.jl`, `--test`), with settings in
  `configs/pilot.toml` (PLAN.md, Step 3).
- `data_conversion/`: the 12-source DR4 pre-release epoch astrometry, RVs for Gaia-4 and
  HD 114762, and Judah's `data_converter.py`, which converts both to Octofitter inputs
  (following the Octofitter DR4 pre-release tutorial). Since 2026-11-04 it writes
  `scan_pos_angle` in degrees, for Octofitter v9.

The folder is a git repository (remote: https://github.com/jvanzand/dr4_sensitivity, public since 2026-11-13). `.gitignore` leaves out the
arXiv source bundles behind the literature tables (~250 MB, kept on Judah's machine):
the per-source `build.py` scripts need them, but everything downstream, from
`literature.py --cached` on, runs from a clone without them (checked: a fresh clone
rebuilds `census_table.csv` byte for byte).

Rebuild the census from the cached downloads with one command (from this folder, in the
`gaia_dr4_sensitivity` environment):

```bash
conda activate gaia_dr4_sensitivity
python -I sample/build_all.py
```

Add `--refresh` to re-query the NEA, Gaia, VizieR, and SIMBAD (a new dated snapshot), or
`--download` to fetch any missing literature source files. See
[CENSUS.md](sample/CENSUS.md) and the Reproducibility section of
[lit_search/README.md](sample/lit_search/README.md).

Key date: Gaia DR4 is scheduled for **2 December 2026** (66 months of data, epoch
astrometry for all sources).

Tools on this machine (checked 2026-10-23):
- Julia 1.12.7. This folder's environment has Octofitter **v9.0.0**,
  OctofitterRadialVelocity v9.0.0, and Pigeons v0.4.11 (upgraded 2026-11-04). The default
  environment, which `RV_DR4_fit_demo` uses, still has v8.2.4.
- RadVel 1.3.8, orbitize! 1.16.1, orvara.
- `../GaiaForecasts` (Lammers & Winn 2025, arXiv:2511.04673): mock DR4/DR5 catalogs and a
  DR4 simulation setup built on `gaiamock` (El-Badry et al.).
- `../rvsearch`: iterative RV planet search, a candidate for blind companion counts.
