# Gaia DR4 sensitivity from known RV giant planets (planning folder)

Working title: *Did Gaia find the planets we think it should? An empirical DR4 sensitivity
map from radial-velocity giant planets.*

This folder holds the plan and the DR4 pre-release template. It becomes a repo once the
plan is firm.

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
- `rv_data/`: RV collection (PLAN.md, Step 2). `build_targets.py` (identities and
  aliases), `jump_queries.sql` (run by Judah in Jump's SQL Explorer; exports kept in
  `jump_raw/`), `ingest_jump.py`,
  `pull_dace.py`, and `merge_rvs.py` (one table per system in `rv_data/systems/`), with
  the duplicate rules in `rv_data/DUPLICATE_RULES.md`.
- `tests/test_duplicate_rules.py`: one synthetic case per duplicate rule
  (`python -I tests/test_duplicate_rules.py`).
- `data_conversion/`: the 12-source DR4 pre-release epoch astrometry, RVs for Gaia-4 and
  HD 114762, and Judah's `data_converter.py`, which converts both to Octofitter inputs
  (following the Octofitter DR4 pre-release tutorial). Since 2026-11-04 it writes
  `scan_pos_angle` in degrees, for Octofitter v9.

The folder is a git repository (local only, since 2026-11-04). `.gitignore` leaves out the
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
