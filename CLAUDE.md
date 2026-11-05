# Project: empirical Gaia DR4 sensitivity map from RV companions

Test how reliably a Gaia DR4-only fit recovers known RV companions, using RV(+astrometry)
reference fits as the answer. Gaia DR4 is released 2 December 2026.

## Where things are
- `PLAN.md`: science framing, analysis Steps 0-6 and their decisions, data-format notes
  (Octofitter v9, time standard, DR4 transit reduction).
- `QUESTIONS.md`: open questions, and settled decisions with dates. Treat "Settled" as
  decided; don't reopen them unless Judah does.
- `sample/`: the census (`CENSUS.md`, `census_table.csv`, figures, `build_all.py`).
- `rv_data/`: RV collection. `build_targets.py` -> Jump queries (`jump_queries.sql`, run by
  Judah; exports in `jump_raw/`) -> `ingest_jump.py`, `pull_dace.py` -> `merge_rvs.py`
  (one table per system in `systems/`). Duplicate rules: `DUPLICATE_RULES.md`.
- `data_conversion/`: DR4 pre-release data and converter. `tests/`: scan-angle (Julia) and
  duplicate-rule (Python) tests.
- Judah's demo fits (Gaia-4, HD 114762, Octofitter v8):
  `/Users/judahvz/academics/applications/Postdocs_2026/RV_DR4_fit_demo`.

## Working rules
- **Python:** the conda environment `gaia_dr4_sensitivity` (Python 3.14; `environment.yml`).
  Run scripts from the project root as `python -I ...` (e.g.
  `~/opt/anaconda3/envs/gaia_dr4_sensitivity/bin/python -I sample/build_all.py`).
- **Julia:** 1.12, `julia --project=.` (Octofitter and OctofitterRadialVelocity v9.0.0).
  v9 takes `scan_pos_angle` in degrees; `tests/test_scan_angle.jl` checks it.
- Bash commands that use the network need the sandbox disabled. ADS and exoplanet.eu are
  unreachable; use arXiv and the mirrors noted in `sample/lit_search/README.md`.
- Be token-efficient: do small tasks directly; use fresh, narrowly briefed subagents only
  for work that needs no history.
- Record every new decision in `PLAN.md` and `QUESTIONS.md` (dated). Never invent values.
- Every script needs a local test mode. Code runs locally first; production jobs go to
  the Sirius cluster (Slurm), which Judah runs: code locally -> push to GitHub -> pull on
  Sirius. The repo is local only for now (no remote yet).
- Plots: load the project `matplotlibrc` with `matplotlib.rc_file`, large figures, and
  every figure self-explanatory (title, a legend entry for every marker, region, line).
- Commit only when Judah asks.

## Conventions already decided
- Times: BJD_TCB everywhere (barycentric). Jump HIRES/APF and DACE HARPS times are
  verified BJD_TDB before conversion.
- Reference solution per companion: newest paper per system unless another has a 2x
  better period precision; period wins over SMA (Kepler SMA if > 20% off). Feng+2022 and
  Gaia-only solutions are excluded.
- RVs: DACE (public data only) and Jump; no KPF RVs (not yet stable). Each dataset gets
  its own offset and jitter. Provenance per RV: database and study.
- Our fits use RVs plus the astrometry Octofitter provides (Hipparcos IAD, HGCA, G23H); the
  DR3/DR4 overlap is acceptable. Literature values never set priors or starting points.

## Status (2026-11-04)
- Done: census (2667 companions, 189 plausibly recoverable); RV pipeline built and run on
  the five pilot systems (TOI-1694, HD 111232, HD 191939, TOI-1736, HD 43691: 1026 RVs in
  `rv_data/systems/`); Octofitter v9 upgrade with the scan-angle test.
- Pending from the original brief: literature M sin i uncertainties; the stage-1 catalog
  and stage-2 fit set (`SNR_MIN`, PLAN.md Step 1b); the alias list for the full fit set
  (`build_targets.py`), then the full Jump export and DACE pull.
- **Next step:** the Octofitter fit structure (PLAN.md Step 3): a configurable pipeline,
  generalized from the demo, that builds one model per system from `rv_data/systems/` and
  `targets.csv`, runs locally on the pilot systems first, and later scales to many
  systems in parallel on Sirius.
- Open decisions for the pilot: which astrometry set (Hipparcos IAD + HGCA, or G23H,
  without double counting); priors (stellar mass and parallax; the census has no
  stellar-mass errors yet); the companion-addition criterion (QUESTIONS.md, Open 1);
  sampler and convergence settings.
