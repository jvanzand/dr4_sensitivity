# Project: empirical Gaia DR4 sensitivity map from RV companions

Test how reliably a Gaia DR4-only fit recovers known RV companions, using RV(+astrometry)
reference fits as the answer. Gaia DR4 is released 2 December 2026.

## Where things are
- `PLAN.md`: science framing, analysis Steps 0-6 and their decisions, data-format notes
  (Octofitter v9, time standard, DR4 transit reduction). Its Section 4 "Overview"
  summarizes the project flowchart (2026-11-05): three lanes (our pipeline: census -> RV
  fits -> synthetic DR4 -> fits; real DR4 -> DR4-only fits; the Epochalypse simulated
  catalog -> DR4-only fits) feeding four comparisons. The main result is our RV fits vs.
  our real-DR4 fits. Three checks: our RV fits vs. the literature (RV fitting); Epochalypse
  truths vs. our fits to its data (DR4 fitting); our synthetic-DR4 fits vs. our real-DR4
  fits (whether the simulations are realistic).
- `QUESTIONS.md`: open questions, and settled decisions with dates. Treat "Settled" as
  decided; don't reopen them unless Judah does.
- `sample/`: the census (`CENSUS.md`, `census_table.csv`, figures, `build_all.py`).
  `catalog.py` builds the stage-1 snapshot catalog and stage-2 fit set from it
  (`stage1_catalog.csv`, `stage1_systems.csv`, `resonance_candidates.csv`).
- `sample/system_names.py` -> `sample/system_names.csv`: readable system names (HD111232,
  TOI1736) for figures and result folders; hierarchy HD > GJ > TOI > Kepler > K2 > WASP >
  HAT-P > HIP, else SIMBAD's main id (coverage 88% / 12%). The Gaia DR3 source_id stays the
  key in every file. `fits_rv/kepler.c` is the compiled Kepler solver (`kepler_c.py`).
- `fits_rv/`: `companion_search.py` (RV-only periodogram + max-likelihood + Delta-BIC
  search for how many companions the RVs support; `--selftest`), and the Octofitter v9
  pipeline `pipeline.jl` / `run_fit.jl` (NUTS fits with convergence checks; `--test`).
  `configs/pilot.toml` holds the priors, sampler, and convergence settings.
- `blind_search/` (started 2026-11-10; PLAN.md "Blind-search tool"): standalone, Octofitter-free blind
  search for DR4 astrometry (+RVs). Modules: `orbit.py`, `data.py` (channels, `select`), `model.py` (stacked
  linear model), `periodogram.py`, `efap.py` (rvsearch port, comparison rule), `simthreshold.py` (simulated
  per-step detection threshold; PLAN.md 2026-11-11), `fit.py` (refit of (P, e, M0) with free jitter),
  `search.py` (`run_search`, the iterative search; M = 50 noise-sampled datasets by default, Judah 2026-11-13), `injection.py`
  (physical companions, Gaia-4 b), `plots.py` (`plot_step_sim`, the step figure), `injection_map.py`
  (injection/recovery map on the Gaia-4 epochs, `--method sim|efap`; QUESTIONS.md 'Plan 2026-11-12'),
  `reliability_tests.py` (noise-only calibration, detection, multi-companion and real-data tests; QUESTIONS.md
  2026-11-13), `make_flowchart.py` (draws `search_flowchart.png`; update its box texts when the code changes).
  Tests: `tests/test_blind_search.py` (`--quick` skips the noise calibration). `blind_search/README.md` describes
  the tool for users and `blind_search/DEVELOPER_NOTES.md` holds the details; `blind_search/FIGURES.md` lists the figure code. Docstrings are wrapped at 80 columns (Judah,
  2026-11-13; keep new ones that way).
  **Core figures only (Judah, 2026-11-13):** the step figure, the injection map and the flowchart. Testing and
  comparison figures were removed (new ones may be added later and pruned again). Removed 2026-11-13 as unused,
  with copies in `../gaia_dr4_sensitivity_blind_search_backup_2026-11-13.tar.gz`: `threshold.py`, the red-noise
  scripts, `make_comparison_figures.py`, `make_figures.py`, `make_sim_figures.py`, `make_method_diagram.py`,
  `make_slide_figures.py`, `make_threshold_comparison_figure.py`, `run_injection_tests.py`,
  `run_threshold_comparison.py`, `injection_map_compare.py`, `injection_map_compare_methods.py`.
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
- **No new tests without asking (Judah, 2026-11-12).** Run the experiments he asked for; check with him before introducing any other test, variant or comparison (e.g. extra fake-noise models), even as a diagnostic.
- **Define terms when you use them (Judah, 2026-11-11).** The first time a code name appears in a reply or note, say
  what kind of thing it is (function, class, module = a `.py` file, folder, script, image file) and where it lives
  (e.g. "`plot_step`, a function in `blind_search/plots.py` that draws one search step's figure"). Define any other jargon
  at first use too (e.g. "stride 4"), and add new terms to the glossary in `QUESTIONS.md`.
- **READMEs stay short and human-readable (Judah, 2026-11-13).** `README.md` and `blind_search/README.md` are for users:
  succinct, with images, no developer detail. Developer detail goes to `DEVELOPER_NOTES.md` and
  `blind_search/DEVELOPER_NOTES.md`; update those, not the READMEs, when recording details.
- **Terminology (Judah, 2026-11-13).** A *model+noise dataset* is a dataset built as a known model plus noise (e.g. the
  zero-planet model plus white noise at the formal errors; stands in for real data in tests: reliability test T1 uses
  N = 45, the injection map 2000). A *noise-sampled dataset* is one of the M datasets sampled from the model fitted to a
  model+noise (or real) dataset plus noise, so any signal in it arises from noise only (M = 50 per simulated threshold).
  Do not call them "fakes", "fake datasets" or "simulated datasets" in new text; convert old text when touching it.
  M and N are different numbers (N x M noise-sampled datasets in total). Defined also in `blind_search/README.md` and the
  QUESTIONS.md glossary.
- **Deliverable format (Judah, 2026-11-11):** a diagram or flowchart is a real `.png` image drawn by a script, not a
  Markdown file with ASCII art. Diagrams use mathematical notation typeset
  with LaTeX (matplotlib `text.usetex`, TeX is installed at `/Library/TeX/texbin`; add it to PATH), and separate each box's
  steps/decisions (left) from the background needed to understand them (grey panel, right).

## Conventions already decided
- Times: BJD_TCB everywhere (barycentric). Jump HIRES/APF and DACE HARPS times are
  verified BJD_TDB before conversion.
- Reference solution per companion: newest paper per system unless another has a 2x
  better period precision; period wins over SMA (Kepler SMA if > 20% off). Feng+2022 and
  Gaia-only solutions are excluded.
- RVs: DACE (public data only) and Jump; no KPF RVs (not yet stable). Each dataset gets
  its own offset and jitter. Provenance per RV: database and study.
- Our fits use RVs plus the astrometry Octofitter provides (Hipparcos IAD, HGCA, G23H); the
  DR3/DR4 overlap is acceptable. Literature values set no priors on the DR4 or sensitivity
  results, with one optional exception decided 2026-11-06: the N literature companions can enter
  the RV search with Gaussian priors around their quoted values (see "Companions"); since
  2026-11-09 this is off by default (`--literature`).

- Fits (decided 2026-11-06; PLAN.md Step 3): G23H astrometry (it contains Hipparcos IAD
  and the HGCA); blind priors from `configs/pilot.toml`: P log-uniform (0.05 d to 200 yr for
  the literature companions, 1 d floor for extra ones), q = M_p/M_* log-uniform
  1e-5 to 1, e uniform 0 to 0.99, isotropic i, per-dataset RV jitter log-uniform 0.1 to
  30 m/s (watch it: ballooned jitter masks poor fits; the fit flags it). Offsets are
  marginalized analytically. NUTS (`octofit`), 4 chains, convergence checked after
  each fit (R-hat, ESS, divergences, tree depth, E-BFMI), rerun with doubled iterations
  on failure. Pigeons and a RadVel-style dynamic loop are later options.
- Companions (decided 2026-11-06, revised 2026-11-09; `fits_rv/companion_search.py`). Default
  since 2026-11-09 (Judah): a fully blind search (`--literature` off), the setting closest to
  what Gaia alone sees; `fits_rv/recovery.py` then tabulates how many literature companions the
  RVs alone recover (period within 10%). Results of the literature-prior mode (the option below)
  are kept in `fits_rv/results_WithLitPlanets/` (run with the older 3x period cap). In that
  mode: a fully blind search cannot find every literature companion (e.g., transiting planets with small K),
  so the N literature companions are in the RV model from the start, with Gaussian priors on
  ln P, ln K, and e around the quoted values (sd = published error, never narrower than 1% in
  P, 10% in K, 0.1 in e; 20% in P and K where none is known; K width from the mass error since
  the catalog carries no K errors; K derived from M sin i where missing; e = 0 assumed where
  unreported; phase and omega free). They are refit as maximum a posteriori; BIC uses the data
  lnL at that fit. Each one's Delta-BIC for dropping it is recorded: small means the prior,
  not the RVs, holds it. Then extra companions are added blindly one at a time if
  Delta-BIC >= 20 (RV count N; no MCMC). The trend (order 0, linear, quadratic) is re-tested at
  every step: first simplified if its top term is worth < 20, then the options are add a
  companion, add a companion and lower the trend order, or raise the trend order; trend moves
  are made before the extras. A blind companion needs period <= 2x the RV baseline (set from 3x on 2026-11-09; longer
  signals can only be a trend) and
  >= 2x the mean RV spacing (span/(N-1); decided 2026-11-09: sparser sampling cannot
  constrain a shorter period, and the fit tends to place a spurious high-e signal there) and
  e <= 0.7 (decided 2026-11-09; literature companions keep e <= 0.99 and their prior).
  Flowchart comparison B is a sanity check, not a test of recovering such planets.
- RV outliers (decided 2026-11-09; `rv_data/merge_rvs.py`): per system, center each dataset on
  its median, take the pooled MAD sigma (1.4826 x median |deviation|), and mark an RV `isjunk`
  if it is more than 10 of these from the median (`rv_outliers.csv`, `rv_outlier_stats.csv`).
  RVs with errors >= 1e5 m/s are placeholders and are marked first.
- Errors (decided 2026-11-06): fill every missing error where possible: from the NEA, from
  the source papers (`e_K_ms`, `e_msini_mj`, `e_e`, `e_a_au` columns of the orbit tables),
  or by derivation (SMA error from P and M_*; P error from the SMA error; M sin i error
  from K, P, e, M_*, flagged `msini_err_from_K`), and say so wherever the catalog is
  shown. Only if no route works is the error omitted: the companion then passes the
  3-sigma cut ("innocent until proven guilty"), is flagged `precision_unverified` and
  listed in `sample/missing_errors.csv` (concerning: neither the quantity's error nor a
  related one's was published). Stellar-mass errors: published where available, else a
  stand-in fractional error (median of stars of similar mass that have one;
  `sample/st_mass_err_standin.csv`), labeled in `st_mass_err_source`.
- Resonant systems (2026-11-06): resonant companions are poorly described by Keplerians, and
  adopting literature orbits would defeat the purpose. The period-ratio proxy
  (`sample/catalog.py`) is too crude to drop systems, so it is a second-tier check, used
  after a system gives cause for concern (e.g., non-convergence). For now nothing is
  dropped: flagged companions are marked (`resonance_member`), counted, and plotted
  (`sample/stage1_resonance.png`). The intended rule, available as `--resonance-policy
  exclude`: resonant companions stay in the model but are not scored as recoveries, and a
  system whose companions are all resonant is not fit.

## Status (2026-11-06)
- Done: census (2667 companions, with error columns, literature errors from the source
  papers, and precision flags; 149 companions still lack an orbit or mass error);
  RV pipeline on the five pilot systems; Octofitter v9 with the scan-angle test; stage-1
  catalog (1989 systems, 2503 companions) and stage-2 fit set (686 systems, 888
  companions at `SNR_MIN` 0.1; resonance flagged, not dropped: 70 stage-1 companions in
  34 systems, 22 in 11 fit-set systems); `companion_search.py` (self-test passes);
  `run_fit.jl --test` runs end to end RV-only. Kepler solver: about half of the
  likelihood time; C or Julia is only 1.7-3x faster than numpy, so no change for now.
- Pilot companion searches (outputs in `fits_rv/results/<display name>/`: `companion_search.json`,
  `.log`, and for HD 111232 and TOI-1736 `trend_vs_planet.png`): the literature periods are
  recovered blindly; HD 111232 prefers a ~22,000 d companion (coverage 0.5) over a linear
  or quadratic trend; TOI-1736 prefers a linear trend plus the 7.08 d companion; HD 43691
  and the MIKE data of HD 111232 have jitter at or near the 30 m/s bound (check); HD 191939
  is to be revisited (forced companions 5-6 have negative Delta-BIC).
- Not yet done: G23H runs (catalog approved and downloaded to data/g23h/, git-ignored; the G23H
  branch of `build_model` is written but untested); the first real pilot fits; the
  Pigeons and dynamic-loop options; the literature companions' M sin i and K errors (only
  period errors are carried: optional, per source paper); the alias list for the full fit
  set (`build_targets.py`), then the full Jump export and DACE pull.
- **Next step:** run the pilot systems through `companion_search.py` and `run_fit.jl`
  (RV-only first), review the jitter flags and convergence checks, then switch on G23H.
