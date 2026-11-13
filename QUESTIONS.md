# Decisions for Judah

## Open

1. **Resonance, second tier.** The period-ratio proxy (within 1% of 2:1, 3:2, 4:3, 5:4,
   5:3, 5:2, or 3:1; `sample/resonance_candidates.csv`) flags 70 stage-1 companions in 34
   systems (22 in 11 fit-set systems; 15 systems are all-flagged) and over-flags (it
   includes HD 191939). Nothing is dropped; assess resonance later, for systems that
   give cause for concern (e.g., a Keplerian fit that will not converge).
2. **VizieR fill-in for RVs:** decide once the DACE + Jump coverage audit shows how many
   systems lack RVs.
3. **Criterion for comparing the number of companions Gaia prefers** with the RV fit
   (the RV-side criterion is settled, 2026-11-06). Bayes factors are one option.
4. **G23H download** (~14 GB catalog): approved 2026-11-06; files in data/g23h/ (git-ignored).

## To remember (later refinements)

- **Pathological solutions** (noted, not corrected; listed in `sample/census_output.txt`):
  - 2× period-precision switches whose period disagrees with the newest paper's at > 3σ
    (41 of 161 on 2026-10-30, e.g., ε Eri b, GJ 3512 c, HIP 14810 c, HD 47536 b).
  - Published SMAs replaced by the Kepler value (21). For 11 Stassun+2017 solutions the
    stellar mass, not the SMA, looks wrong; Kepler-34 b and Kepler-35 b are circumbinary.

- **Number of companions in the fits.** Superseded on 2026-10-27: RV fits add companions
  iteratively (criterion settled 2026-11-06), and Gaia refits are compared on the number of
  companions preferred as well as the orbital parameters (PLAN.md, Section 2.6). Formal
  model comparison (Bayes factors) is one way to make the Gaia-side criterion concrete.
- **Extra companions in updated RV data.** Where RVs have been added since the published
  fit, the data may be sensitive to further companions. Expected to be rare; flag cases
  where the residuals suggest it.
- **Transit-level errors.** The demo already collapses CCDs to one row per transit, so CCDs
  are never treated as independent. Its transit errors are a single CCD's error, about 2×
  the Octofitter tutorial's bootstrap-median errors. Default to the tutorial's reduction and
  test the demo's as the conservative case in the pilot (PLAN.md, Section 3).
- **Companions still lacking an error (149 of 2667; `sample/missing_errors.csv`).** After
  the NEA, derivations, and the source papers: 112 NEA rows (mostly transit+RV planets and
  RV+astrometry fits whose papers give no mass or period error) and ~35 literature rows
  (SB9, exoplanet.eu, Stevenson+2023, and zero-error K in Troup+2016). They pass the
  cut and are flagged. This is concerning for the NEA rows: look for other NEA solutions
  or the papers. Two marvels_bd M sin i errors come from the text, not a table (Fleming+2012
  ±3.0 M_J, Jiang+2013 ±2.9 M_J); the SB9 M sin i has no published error.
- **Dynamic NUTS loop (later, if the pilot timings call for it).** A custom loop on
  AdvancedHMC that keeps the adapted metric, samples in blocks, checks RadVel-style
  criteria between blocks (R-hat, autocorrelation factor, Tz), stops early, and prints
  real-time progress to the terminal for single-system diagnostics (a table of many
  systems for parallel runs is a possible extension).
- **Trend vs. planet ordering** in the companion search is handled by letting both compete
  at every step; check cases where a signal could be either.
- **Look-elsewhere effect** on the fixed Delta-BIC = 10 threshold: the self-test found 0
  false additions in 10 noise-only datasets; a larger calibration on realistic cadences is
  worth doing before trusting counts on sparse systems.
- **Feng+2022 exclusion (2026-10-30).** 28 NEA planets whose only NEA solution is
  Feng+2022, and Feng+2022's 114 transcribed literature rows, are out of the census while
  those fits are considered suspect. Revisit if that view changes, or fit their RVs
  ourselves. (Of the 94 planets in `sample/no_rv_only_solution.csv`, which lacked an
  RV-only solution, 28 are now in the census.)
- **Literature follow-ups:** prioritized leads, known issues (mixed time standards,
  transcription completeness, unverified citations), and the excluded Feng+2022 systems
  are in `sample/lit_search/README.md`.
- **Bulk data-quality checks.** Duplicates and similar issues are handled per system for now;
  Judah will add automated checks for bulk runs.

## Settled (2026-11-06): fit setup, companion search, errors, resonances

Details in PLAN.md (Step 1b, Step 3).
- **Astrometry set: G23H** (Hipparcos IAD, HGCA, Gaia DR2/DR3; Thompson et al.), not
  Hipparcos + HGCA separately, so nothing is double counted.
- **Priors:** the G23H defaults, with P log-uniform from 0.05 d (known companions; the
  census minimum is 0.07 d) or 1 d (extra companions) to 200 yr; q log-uniform 1e-5 to 1;
  e uniform 0 to **0.99**; i isotropic; per-dataset jitter log-uniform 0.1 to 30 m/s,
  watched in every fit; RV offsets marginalized analytically. Stellar mass: Gaussian,
  from the census mass and its error. No extra short-period census cut is needed:
  `SNR_MIN` already removes undetectable short-period companions from the fit set, while
  the RV model still carries them.
- **Stellar-mass errors:** published errors where available (NEA solution, NEA composite,
  Gaia FLAME bounds), else a stand-in fractional error: the median among stars in the
  same mass bin that have one (about 5%, 9% for 1.6-3 M_sun), labeled in the table.
- **Sampler:** NUTS first (`octofit`); Pigeons as needed.
- **Convergence:** checked after each fit (R-hat, bulk/tail ESS, independent draws per
  chain, divergences, tree depth, E-BFMI), with an automatic rerun at doubled iterations
  on failure. The dynamic loop is the later option (see "To remember").
- **Companions:** start with the literature's N (blind orbits), then add one at a time if
  Delta-BIC >= 20 from a maximum-likelihood fit to the RVs alone, with N in the BIC the RV
  count. The trend competes with the extra companion at every step and is re-tested after
  each addition. The periodogram is a fast generalized Lomb-Scargle for ranking periods;
  the threshold is fixed, not calibrated to each system's noise (a problem with
  `rvsearch`'s). Extra companions outside the census cuts are found by the RVs, not claimed
  from the literature.
- **Errors:** a missing or zero error never excludes a companion from the 3σ cut
  ("innocent until proven guilty"). Companions that pass only for that reason are flagged
  (`precision_unverified`). SMA errors are derived from the period and stellar-mass errors
  where the published SMA is not kept or has no error. M sin i errors are derived from K
  (with P, e, M_* terms) when only a K error exists, and flagged (`msini_err_from_K`,
  `mass_err_source`); the catalog's system table carries a system-level flag. The orbit
  precision test uses P/σ_P, or the published SMA error when P has none.
- **Resonance** (revised the same day): the period-ratio check is a crude, second-tier
  proxy, so no system is dropped for now (Open 1). The intended rule, once resonance is
  confirmed: resonant companions stay in the model but are not scored, and a system with
  only resonant companions is not fit. Flagged companions are marked, counted, and
  plotted.
- **Missing errors are filled wherever possible** (NEA, source papers, derivation); only if
  no route works is the error omitted and the companion flagged (it still passes the
  cut). Derivations include the period error from the SMA error and, for transiting
  planets, the mass error from the K error.
- **Period cap: 2x the RV baseline** (3x until 2026-11-09, when Judah tightened it). A companion
  needs at least half of its orbit covered; longer-period signals can only be a trend. The trend allows quadratic curvature.
- **Trend ladder:** at every step the trend is first simplified if its top term is worth
  < Delta-BIC 20, then the options are add a companion, add a companion and lower the
  trend order, or raise the trend order.
- **Kepler solver:** profiled at 40-58% of the likelihood time; C and Julia are only
  1.7-3x faster than numpy (0.04-0.06 ms vs. 0.09-0.2 ms for 410 points), so it is not
  the bottleneck (the numerical-gradient fits are). Revisit for large runs.
- **Stage 1 built** (1989 systems, 2503 companions after the literature errors went in):
  systems with at least one precise companion. Stage 2 (`SNR_MIN` 0.1): 686 systems, 888
  companions.

## Settled (2026-11-05): project structure and Epochalypse

- **Four comparisons** structure the project (PLAN.md, Section 4, "Overview"): our RV fits
  vs. real DR4 fits (main result); our RV fits vs. the literature (RV-fitting
  consistency); Epochalypse truths vs. our fits to its data (DR4-fitting consistency);
  our synthetic-DR4 fits vs. our real-DR4 fits (realism check).
- **Epochalypse** (Yahalomi et al. 2026, arXiv:2611.06915) is a validation set for our
  DR4 fitter, not a substitute for Step 4: it injects random companions, not our RV
  orbits, and doesn't test known planets.

## Settled (2026-11-04): RV collection and refits

- **Our fits use RVs plus the astrometry Octofitter provides** (Hipparcos IAD, HGCA, G23H)
  where available. The DR3 overlap with DR4 is acceptable: the question is how reliable a
  DR4-only fit is, which needs the best available answer. It would matter only for a
  "DR4 vs. everything but DR4" comparison, which is not planned.
- **RV sources:** DACE is the primary hub, then Jump (Judah's SQL query), then the RV
  catalogs bundled with Octofitter. VizieR is deferred (Open, 2).
- **Two-tier RV provenance:** each RV records the database it came from and, where known,
  the study that published it.
- **Octofitter v9** upgrade approved, and done the same day: v9.0.0 in this folder's
  environment, with the scan-angle unit test passing (PLAN.md, Section 3).
- **Pilot systems:** TOI-1694 (Jump), HD 111232, HD 191939 (multi-planet), TOI-1736
  (subgiant host), HD 43691 (= HIP 30057; giant host, 1.38 M_sun).
- **No KPF RVs** for now: the instrument isn't yet stable over the long term.
- **Public DACE data only** (no DACE account needed) are acceptable for the fits.
- **Duplicate RVs** follow the standard rules in `rv_data/DUPLICATE_RULES.md`, including
  the size tie-break (a larger dataset wins among otherwise equal ones, since one longer
  series needs fewer offsets).
- **Python:** one conda environment, `gaia_dr4_sensitivity`, on a current Python (3.14);
  `environment.yml` recreates it.
- **Jump and DACE overlap:** where both hold the same exposures, keep Jump's CPS
  reduction for Keck and APF data, and DACE's for the instruments it hosts.
- **Duplicate star records** are merged by system (Gaia DR3 ID) in our pipeline, with
  positions checked against SIMBAD (PLAN.md, Step 2).

## Settled (2026-10-27): study sample and pipeline inputs

Details in PLAN.md, Step 1b and Section 3.
- **Sample unit:** a system. Its companions are of interest, but a system is fit as a whole.
- **Precision cut:** a companion is precisely characterized if SMA (or P) and M sin i (or M)
  are both measured at ≥ 3σ. The M sin i requirement means the mass comes from RVs, so
  planets found purely by transits don't qualify. The same rule applies to every discovery
  method.
- **Two stages of cuts.**
  - *Stage 1, the snapshot catalog:* every census system with at least one precisely
    characterized companion, with all its companions, chosen on orbit precision alone (no
    S/N). Fixed before DR4 together with the census snapshot; changes are rare and
    recorded. Later restrictions are flags for specific fits.
  - *Stage 2, the fit set:* catalog systems in which any companion, precise or not, has
    forecast S/N ≥ `SNR_MIN`. `SNR_MIN` is a top-level, user-settable parameter (default
    0.1), tuned to balance comprehensiveness against compute time. A system with a precise
    low-S/N planet and an imprecise high-S/N planet is fit.
- **No separate host-quality cuts** (brightness, RUWE, distance, partial orbits): the S/N
  and precision cuts decide. Precise companions that Gaia misses are useful results.
- **Whole systems are modeled:** marginal detections and trends are included in systems with
  at least one precise companion. RV fits add companions iteratively; Gaia refits are
  compared on the number of companions and on the orbital parameters.
- **RV availability** is a later cut (Steps 2-3). Systems without RVs we can refit drop out
  there, flagged rather than removed.
- **RVs:** fresh retrieval for every system, from DACE and from Jump. Judah runs the Jump
  SQL query that cross-matches the sample's names; the pipeline supplies the names.
- **Fitting code:** generalize Judah's `RV_DR4_fit_demo` (joint RV + DR4 fits of Gaia-4 and
  HD 114762). No known-answer priors or starting values.
- **Astrometry simulation** is its own sub-project; Octofitter combines it with the RVs.
- **Time standard:** all time stamps in TCB, barycentric (Gaia's `obs_time_bary_corr`
  applied). This replaces the BJD_TDB standard for periastron times below.
- **DR3-DR4 source IDs:** use Gaia's published correspondence.

## Settled (2026-10-23)

- **Period limit for "plausibly recoverable" (2026-10-26):** P ≤ 5.5 yr, the DR4 baseline
  (189 companions as of 2026-10-30), for estimation. Lammers & Winn's 4-yr limit for
  accurate periods (156 companions) stays marked in the figures for reference.
- **Sample cuts (revised 2026-10-26):** M sin i ≤ 80 M_J and SMA ≤ 5 AU, with no lower
  bounds and no brightness limit; any star with an RV companion is eligible.
- **Keep all companions** (2026-10-26), including low-precision APOGEE and MARVELS orbits
  and partial orbits.
- **Distant hosts:** sensitivity beyond ~100 pc is expected to be essentially zero, so no
  extra effort goes into extending the distance axis.
- **Periastron times** are standardized to BJD_TDB during ingestion (superseded 2026-10-27:
  TCB, barycentric; the census conversion is still BJD_TDB).
- **Substellar companions** come from RV survey catalogs as well as the NEA. An M vs. M sin i
  mismatch is acceptable: the other orbital parameters are compared directly, and masses
  separately. Luminous-companion contamination is considered minimal in RV samples.
- **Eligible solutions (amended 2026-10-30; was RV-only):** census values may come from
  RV-only or RV+astrometry solutions (Hipparcos-Gaia, Gaia DR2/DR3). "Ground truth" means
  the preferred solution for a system, chosen by the rule below; RV+astrometry fits are
  often better but are not preferred automatically (clarified 2026-11-04). The recovery
  stage uses only DR4, so Gaia-informed inputs are acceptable. Exception: Feng+2022 is
  excluded, since its fits are considered suspect. Every solution must use RVs, so
  Gaia-only orbits (the NEA's Gaia Collaboration 2023 solutions) are excluded too.
- **Choice of NEA solution (2026-10-30; replaces the archive default; confirmed
  2026-11-04):** newest paper per system (covering the most planets), then any solution
  with less than half the period error replaces it. No extra guard for now; tension
  cases are noted (see "To remember").
- **Period over SMA (2026-11-04):** period is read directly from RV and astrometric data,
  while SMA needs the host mass. Where they disagree by > 20%, SMA comes from Kepler's
  law with the period. A newer, more precise period wins even if its solution's SMA is
  inconsistent and an older solution's is not.
- **Published sources only** (arXiv is fine); no numbers from talks.
- **RV fits:** an initial "blind" search with weakly informative priors (broad uniform
  bounds), number of planets fixed to the literature count (revised 2026-10-27: whole
  system, companions added iteratively). Disagreements with the
  literature are investigated; heavy ones are flagged.
- **Astrometry fits:** far more blind, with wide priors; literature companion count for now
  (revised 2026-10-27: compare the number of companions Gaia prefers with the RV fit).
- **Literature values never set priors** in either fit; they are post-hoc checks only.
- **Octofitter:** pinned to v8.2.4 at first; upgraded to v9 on 2026-11-04 (see above).
- **Detection criterion:** provisional (mass posterior excludes zero at a significance set by
  simulated no-planet fits, and the period agrees with the RV period); see "To remember."
- Detections are defined by our pipeline; DPAC's solutions are a secondary comparison.
- Control stars and the DPAC comparison are extras for later.
- "Gaia-only" means DR4 epoch astrometry alone.
- Sensitivity is fit in (α/σ, P) and projected onto distance.
- RVs come from DACE plus Judah's team database (after permissions).
- Claude works locally only; every script needs a local test mode. Judah runs production
  jobs on the remote cluster.
- DR4 pre-release data and the conversion template are in `data_conversion/`.
- Timing is not a constraint; build the scaffolding first and sharpen novelty later.
- Judah handles contacting collaborators once the initial checks look promising.

## Settled 2026-11-06 (review of the first companion searches)
- The automated search need not recover every literature companion: it gets far less
  attention per system than a paper does. What matters is whether Gaia recovers what the RVs
  show; the RV(+astrometry) refit is the intermediate comparison (flowchart comparison B/A).
- Gallery of the pilot searches: `fits_rv/plot_gallery.py` -> `fits_rv/results/search_gallery.png`.
- The draft project flowchart is `sample/project_flowchart.png` (`sample/plot_flowchart.py`),
  reconstructed from the PLAN.md Overview because the original sketch was not in the repo.
- More optimizer starts for the extra-companion step wait until Judah has checked the BIC
  (`fits_rv/companion_search.py`, `n_params` and `bic`).
- Claude asks Judah before using any subagent (usage limit).
- Open: DACE reports a SOPHIE RV of 999001.8 m/s with error 999000 m/s (HD 43691) that is not
  flagged `isjunk`; it has no weight in the fit but counts toward N and wrecks plots. Proposed:
  drop rows with error >= 1e5 m/s in `merge_rvs.py` and `load_system`. Awaiting Judah's OK.

## Settled 2026-11-06 (outliers, threshold, figures)
- Placeholder RVs (error >= 1e5 m/s) are marked `isjunk` in `rv_data/merge_rvs.py`
  (HD 43691: SOPHIE 999001.8 +/- 999000 m/s; N 89 -> 88). A general outlier rule is not yet
  decided. Preliminary pass (median, STD of all RVs, worst |RV - median|/STD): TOI-1694
  -4.0, 25.5, 1.7; HD 191939 -2.5, 3042, 3.0; TOI-1736 32.3, 128.9, 2.6; HD 111232
  104464, 31449, 3.3; HD 43691 -28880, 108133, 9.5 (the placeholder). Pooled statistics mix
  datasets with different zero points, and the placeholder inflates the STD so it sits below
  a 10-STD cut. Candidate rule: center each dataset on its median, use a MAD-based sigma.
- The Delta-BIC adoption threshold is 20 (was 10); `DELTA_BIC` in `companion_search.py`.
  Pilot searches rerun.
- Figure DPI: Claude's 60-80 was arbitrary, kept low to keep images light to view; Judah set
  the flowchart to 250. Gallery now 120.
- Pilot searches rerun at Delta-BIC 20 (2026-11-06): TOI-1694 and TOI-1736 now stop at the
  literature companions (extras at dBIC 2.4 and 15.9 rejected); HD 111232 keeps 3 (19.7,
  1168, ~17,100 d; the long one moved from 22,754 d); HD 43691 unchanged (N = 88). HD 191939's
  forced fit changed a lot (now 0.087, 0.589, 8.88, 28.6, 101.7, 3631 d; before 0.14, 8.88,
  28.6, 102.2, 764, 1941 d): the multi-start optimizer is not reproducible there, which the
  planned extra starts should address once Judah has checked the BIC.

## Decided 2026-11-06 (literature priors; details open)
- A fully blind search will not find every literature companion (transiting planets with
  small K, for example). The N literature companions go into the RV fit from the start,
  not forced to be found blindly, with Gaussian priors around the quoted values; values are
  refit under those priors. Extra companions are still found blindly (Delta-BIC >= 20).
  This supersedes "literature values never set priors" for the RV search only.
- Comparison B (RV fits vs. literature fits, flowchart) is a sanity check, not a scientific
  result; the pipeline is not being tested on recovering transiting planets.
- Open (proposed to Judah): which parameters get priors (P, K, e; phase and omega free?);
  widths (published errors with floors, and a default where none is published); BIC uses
  the data log-likelihood at the maximum-a-posteriori fit; flag companions whose posterior
  stays at the prior ("prior-dominated") so they are not read as recoveries.
- Gallery: "found" now means the final model's companions, marked as added by the search
  (passed Delta-BIC), forced with Delta-BIC >= threshold, or forced with Delta-BIC below it
  (hollow); legends are placed away from data (`place_legend` in `plot_search.py`).

## Settled 2026-11-06 (later)
- BIC code (`n_params`, `bic`, and the likelihood) checked by Judah: it counts the model
  parameters properly. Extra optimizer starts are approved (to be done with the literature-
  prior base model).
- G23H (CANFAR DOIs 10.11570/26.0002, 14.5 GB, and the mandatory DR2-transit sidecar
  26.0016, 314 MB) downloaded to `data/g23h/`; `/data/` is git-ignored; `configs/pilot.toml`
  and `pipeline.jl` read the local files (untested until a G23H run).
- `trend_vs_planet.png` now exists for all five pilot systems (`plot_search.py --compare-from`).

## Settled 2026-11-09 (literature priors implemented; MAD rule)
- Prior widths accepted as proposed: ln P sd = max(published error / P, 1%), or 20% if no error;
  ln K sd = max(relative error, 10%), or 20%; e sd = 0.1; phase and omega free. The catalog has
  no K or e errors, so the K width uses the mass error as a proxy (flagged in the JSON notes)
  and e uses 0.1; K is derived from M sin i where missing (flagged); e = 0 is assumed where
  unreported (flagged). Caveat: some period errors are derived from SMA errors (HD 191939 c, d),
  looser than the transit periods really are.
- Implemented in `companion_search.py` (`literature_priors`, `fit_literature`,
  `polish_literature`, `support_check`): base model from the quoted values, periodogram-window
  and phase starts, coordinate polish, trend seeded by least squares, and a from-scratch
  rebuild for trend raises (a plain refit from the no-trend solution got stuck: TOI-1736 BIC
  2749 vs 2332). Extra-companion search has more starts (8 peaks x 3 eccentricity starts, 4
  joint refits). `--no-literature` gives the old blind search.
- MAD outlier rule implemented in `rv_data/merge_rvs.py` (see CLAUDE.md); the five pilots have
  no outliers beyond the placeholder row.
- Pilot searches rerun with literature priors (2026-11-09): TOI-1736 (2 + linear trend) and
  HD 43691 (1) are clean; HD 191939 keeps all 6 literature companions, two held only by their
  priors (g, 288.6 d, K 0.8 m/s, and d, 38.35 d, K 0.7 m/s; dBIC for dropping -19 and -17), no
  extra, no trend. Concerns, awaiting Judah: HD 111232 adds five blind companions plus a linear
  trend, including a spiky high-e long-period companion whose periastron falls in a data gap
  (model reaches ~6000 m/s near MJD 51,400); TOI-1694 adds a 4.095 d companion at dBIC 22.0
  (28 RVs, close to the 3.77 d literature period).
- Sampling floor for blind companions (decided 2026-11-09, Judah): a blind companion needs
  P >= 2 x the mean RV spacing, span/(N-1) (Nyquist-style; literature companions exempt), because
  a shorter period is not constrained by the sampling and the fits placed spurious, often
  high-e, signals there (`NYQUIST_FACTOR`, `RVData.p_min_blind`; self-test 7). Floors: TOI-1694
  90 d, HD 43691 140 d, HD 111232 165 d, HD 191939 30 d, TOI-1736 11 d. Rerun result: no blind
  companions for TOI-1694 (the 4.095 d one is gone), HD 43691, HD 191939, TOI-1736. HD 111232
  still adds a 408 d companion (e 0.99, dBIC 48.7) and a ~9,370 d one (e 0.986, K 4,260 m/s,
  model spikes to ~8,200 m/s in a data gap) after a linear trend; the floor does not touch these.
  Open nuance: the mean spacing is set by season gaps in clumped series (HD 111232, HD 191939,
  HD 43691 have dense runs with median gaps of 1-17 d), so a dense-subset criterion (some
  window of >= ~3 cycles with >= ~10 RVs at spacing <= P/2) would be looser. Not adopted.
- Eccentricity cap for blind companions (decided 2026-11-09, Judah): e <= 0.7 (`E_BLIND_MAX`;
  a soft wall in the fit plus projection onto the cap; literature companions keep 0.99 and
  their prior; self-test 8). Only HD 111232 rerun (the other four had no blind companions, so
  the cap only matters for their rejected candidates). HD 111232 now: literature 1168 d
  (K 156 m/s) plus a 27,300 d companion (K 224 m/s, e 0, 0.45 of an orbit covered, near the
  3x-baseline cap, so it is largely a curvature stand-in) and a 935 d companion (K 4.0 m/s,
  e 0.65, dBIC 26.5, marginal). No trend. BIC 1160 vs 1069 before the cap: the e ~ 0.99
  companions were buying 90 BIC units, but with a model that spikes to ~8,000 m/s in a data
  gap. Concern: the 935 d signal is only 6 over the threshold of 20 (look-elsewhere effect).

## Settled 2026-11-09 (second batch of test systems)
- Judah asked for five more test systems, run with the current decisions (priors, DELTA_BIC 20,
  sampling floor 2x mean spacing, blind e <= 0.7). Judah's framing: the RV search need not
  recover the literature exactly (it is a consistency check), and for the astrometric fits,
  spurious extra signals have to be handled separately; what we judge there is whether the
  literature signals are recovered alongside any spurious extras.
- Chosen: HD 39091 (pi Men), HD 89707, HD 147513, HD 10180, HD 168443 (see PLAN.md Step 0).
  RVs: DACE public data only (`build_targets.py` now lists 13 systems: the 5 pilots, these 5, and
  70 Vir, HD 94340, HD 160691 (mu Ara) as unsearched alternates).
- Data fix: `merge_rvs.py` lacked instrument families for HIRES-PRE04/POST04 (Rosenthal+2021
  reprocessing of the same HIRES spectra as Butler+2017), APF, HARPS (publication series),
  HAMILTON, UCLES, CES, so HD 168443's HIRES exposures were counted twice (138 RVs). Added the
  families; the five pilot tables are byte-identical, 70 Vir lost 67 copies, HD 168443 138.
  Leftover small HIRES remainders of HD 168443 (14, 2, 3 RVs) keep their own offset/jitter.
- Results (final models; blind = added by the search):
  - pi Men: 3 literature companions all RV-supported (dBIC for dropping 17000, 272, 221); linear
    trend; one blind companion (72.7 d, K 0.6 m/s, e 0.68, dBIC 64.7), suspicious (K below the
    noise, e just under the 0.7 cap, likely stellar activity).
  - HD 89707: b recovered (e 0.94, K 4.9 km/s); one blind companion (1006 d, K 23 m/s, e 0.4,
    dBIC 24.2, marginal).
  - HD 147513: b recovered (P 540 vs 528 d); a blind 3102 d companion (K 78 m/s, e 0.58,
    dBIC 144) replaces the trend: needs review (the CORAVEL points are placeholders-quality).
  - HD 10180: all six literature companions supported (dBIC for dropping 61 to 617); a linear
    trend (dBIC 26.4, marginal); no blind companions (floor 51 d excludes the 7.8 d candidate).
  - HD 168443: b and c recovered; linear trend plus a blind 17,336 d companion (K 7 m/s, e 0,
    dBIC 36), a curvature stand-in.

## Settled 2026-11-09 (blind default, 2x period cap, blind rerun of all ten systems)
- Literature priors are now optional: `companion_search.py --literature` (default off = fully
  blind, closest to what Gaia alone sees). `--results-dir` picks the output folder. Judah renamed
  the earlier folder `fits_rv/results_WithLitPlanets/` (run with literature priors; the older
  3x period cap for most, not regenerated).
- Blind period range: floor = 2 x mean RV spacing (min 1 d), cap = 2 x RV baseline (min with 200
  yr), e <= 0.7. A signal longer than the cap can only be a trend/curvature term.
- `fits_rv/recovery.py` scores recovery: a literature companion is recovered if a model
  companion has period within 10% (`RECOVERY_TOL`), one to one. Literature companions below the
  floor or above the cap are labeled "not searchable", not counted as misses.
- Blind results (`fits_rv/results/`, galleries `search_gallery.png`, `search_gallery_batch2.png`,
  table `recovery.csv`, `recovery_extras.csv`): 25 literature companions; 15 in the searchable
  range; 10 of those 15 recovered; 38 extra model companions (22 of them in HD 168443).
  Not recovered although searchable: HD 10180 g, h; HD 191939 d, g, f.
- Known failure modes seen in the blind run (not yet addressed): (1) HD 168443 hit the loop
  limit (3 x `max-extra` = 24 steps; `max-extra` limits iterations, not companions) because its
  58 d K = 477 m/s companion is below the 59 d floor, so the search tiles 60-640 d with 22
  e <= 0.7 signals to mimic it; (2) HD 89707 b (e 0.9, K 4.2 km/s) is fit as e 0.7, K 2.6 km/s
  plus a 293 d K 683 m/s extra; (3) HD 43691 b (37 d) is below the 140 d floor and six
  high-e signals fill in. Decision pending on whether to cap extras, tighten the high-e extras
  (e.g., require dBIC jump), or accept as-is for the astrometry stage.

## Open 2026-11-09 (how to run the blind search on Gaia astrometry)
- Judah's question after the blind RV rerun: rebuild an RVSearch-like search (per-system BIC
  threshold) in Octofitter for astrometry, RVs, and both, or accept the literature fits as the
  reference and run a blind Gaia search with an arbitrary Delta-BIC threshold (known to produce
  spurious signals, so less informative)?
- Of the 38 extras in the blind RV run, only 7 have Delta-BIC < 60 (HD 10180 x2, HD 39091 x2,
  HD 43691 x3) and 31 have Delta-BIC >= 60 (22 in HD 168443), so a different threshold alone
  would not remove most of them. (Correction, same day: the claim that a per-system, noise-based
  threshold could not remove them rested on a wrong description of RVSearch; see below.)
- Figures changed 2026-11-09: `sample/stage1_resonance.png` legend now gives companions and
  systems per group (legend moved below the axes); the flowchart's census box now reads its
  counts from `stage1_catalog.csv` / `stage1_systems.csv` (1989 / 2503 -> 686 / 888).

### Update 2026-11-09 (RVSearch's threshold, and what the blind test is for)
- RVSearch's threshold (checked against `rvsearch/periodogram.py` and `search.py` on GitHub; the
  paper's Sec. V.1 "Planet Search", read via ar5iv): there is NO residual shuffling or noise
  injection. `Periodogram.eFAP()` takes the Delta-BIC values of the real data's periodogram over all
  trial periods, keeps the 50th to 95th percentile, histograms them (10 bins), fits a straight line
  to log10(counts) (an exponential tail with decay rate A), and sets
  threshold = ln(FAP / N_periods) / (-A) + (lowest kept Delta-BIC). FAP = 0.001 by default.
  `run_search` recomputes it for every periodogram (every added planet); a detection needs
  best Delta-BIC > threshold. A system whose periodogram has power at many periods has a slow
  decay (small A), hence a high threshold. Not yet tested on our systems.
- Judah's decision: the period floor (2 x mean spacing) and the e <= 0.7 cap are artifacts of our
  RV test and cannot be in the actual (Gaia) run. The 2x-baseline period cap was not mentioned
  (kept; to confirm).
- Question clarification: a search with a FIXED number of companions answers only "given N known
  companions, does DR4 alone recover their parameters?". The target use (a Gaia user who does not
  know N) needs a true blind search; the fixed-N fit is a secondary diagnostic (if even fixed-N
  fails, a blind search cannot succeed).

## Proposed 2026-11-10 (blind astrometric search: Octofitter scaffold; not yet approved in detail)
- Judah (2026-11-10): a blind search is needed to answer the question (a Gaia user who does not know
  N); build it by scaffolding a search around Octofitter, inspired by RVSearch, unless a better plan.
- Proposal: three layers. (A) a fast linear periodogram of the DR4 along-scan data (for fixed period,
  eccentricity and phase the model eta = dRA*sin(psi) + dDec*cos(psi) + plx*f_al + orbit is linear
  in the 5 reference parameters and the 4 Thiele-Innes constants), written new; (B) RVSearch's eFAP
  threshold (exponential-tail fit to the periodogram, FAP 0.001) validated on noise-only mocks;
  (C) Octofitter for the maximum-a-posteriori refit of all companions after each detection and for
  final posteriors (NUTS only at the end). Fixed-N fits become a mode of the same code.
- Validation before DR4 (2026-12-02): unit test (periodogram chi^2 equals Octofitter's likelihood),
  noise-only false-alarm rate, injection-recovery, then the ten test systems; freeze all settings
  before touching real DR4.
- Open decisions: BIC parameters per companion (5 vs 7), eccentricity handling in the periodogram,
  period grid and whether the 2x-baseline cap stays, astrometric jitter term, mock generator
  (Octofitter's `generate_from_params`, mentioned in its source, vs gaiamock) and scanning-law source.

### Update 2026-11-10 (blind-search plan, edits from Judah)
- Scope: first build and validate only the blind search software (no 686-system run yet).
- Phase 1a (real pre-release DR4 epochs, real scan schedules): Gaia-4 b (11.8 M_J, 571 d, forecast
  S/N 6.1) and HD 114762 (83.9 d, S/N 1.85, near threshold). Phase 1b (simulated DR4, needs the
  scanning-law forecast, PLAN Step 4): astrometry-amenable systems that already have RV tables and
  reference fits: pi Men b (10.5 M_J, 2089 d, S/N 19.7; inner c, d should NOT be detected),
  HD 111232 b (7.2 M_J, 1169 d, S/N 8.7), HD 168443 c (17.8 M_J, 1750 d, S/N 20.2; b at 58 d,
  S/N 0.94, is a stress test), HD 89707 b (54 M_J, e 0.9, S/N 21), optionally TOI-1736 c (S/N 2.4).
  Forecast S/N is the census estimate (`snr` in stage1_catalog.csv).
- State object: RVSearch's `add_planet` (checked in search.py) builds a NEW radvel Parameters and
  posterior from copies of the old values plus defaults for planet n+1 (with `vary` flags), and
  replaces `self.post`. Octofitter's `LogDensityModel` has fixed dimension, so the equivalent is an
  immutable `SearchState` (values by name, fixed/free flags, trend order, jitter, lnL, BIC,
  periodogram, threshold, fit status) from which a fresh Octofitter model is built each step; every
  state is saved. Fixed parameters are constants in the model's variables block. Recompile cost per
  rebuilt model is unmeasured (M0).

### Update 2026-11-10 (blind search is Octofitter-free; basis; data selection)
- Decision (Judah, agreed): the search code does not import Octofitter. Periodogram, least squares,
  joint refits and model selection are standalone; its only output for Octofitter is the final model
  (number of companions, trend order, jitter, starting values) written as a readable Julia model
  script that runs on vanilla Octofitter v9. The test suite (not the tool) uses Octofitter: V1
  equivalence of the likelihood, and an independent mock generator, so a convention error cannot
  hide by appearing in both the generator and the fitter.
- Parameterization: for fixed (P, e, mean anomaly at a reference epoch) both the along-scan
  astrometric model and the RV model are linear in the remaining parameters (Thiele-Innes constants
  A, B, F, G; K cos(omega), K sin(omega)), so omega, Omega, i, a never enter the nonlinear fit. One
  internal basis (ln P, e, M0) plus one conversion function to Octofitter's physical parameters at
  export (round-trip tested). No radvel-style basis switching.
- Data channels: any subset of datasets (DR4 epochs, Hipparcos IAD epochs, RVs per instrument) can be
  selected; astrometric datasets share the companion's orbit but each has its own reference-motion
  block; RV datasets get their own offsets and jitter. HGCA/G23H proper-motion constraints are not
  epoch data and are deferred.

### Update 2026-11-10 (priors, plots, G23H; Judah)
- Tests use Octofitter (V1 likelihood equivalence, independent mock generator): agreed.
- Priors for the final Octofitter fit: the classic uninformative set. P log-uniform; e via
  sqrt(e) cos(omega) and sqrt(e) sin(omega) uniform on the unit disk (this is uniform in e and uniform
  in omega); mass or amplitude log-uniform; isotropic inclination (Sine); Omega and M0 uniform. The
  search stage is maximum likelihood plus BIC, so priors enter it only as parameter bounds.
- Plots (requirement for the periodogram search): produce informative plots like those of the
  preliminary RV code (`fits_rv/plot_search.py`, `plot_gallery.py`): per search step, the periodogram
  (Delta-BIC against period, with the period floor and cap, the eFAP threshold, the peak, and any
  literature period where known), and the data with the model and residuals. Also draw the eFAP
  tail fit: the histogram of the periodogram's Delta-BIC values (log counts) with the 50th-95th
  percentile range used, the fitted straight line, and the resulting threshold, so the false-alarm
  level can be checked by eye at every step. Every figure self-explanatory (title, legend entries).
- G23H: the catalog (data/g23h) holds measurements (catalog proper motions and covariances,
  Hipparcos IAD, unit-weight-error terms, per-transit information), not orbit posteriors, as far as
  the Octofitter documentation shows (the 14 GB file's columns were not inspected). No pre-computed
  companion posterior to reuse; the pre-DR4 constraint on companions has to come from fitting those
  data. In Octofitter v9 the G23H likelihood is in `src/legacy/`. Its proper-motion terms are not
  epoch data, so they do not fit the linear periodogram; deferred as already noted.
- G23H file inspected 2026-11-10 (via Octofitter's Arrow; the project Python env lacks pyarrow): 18,061,961 stars x 114
  scalar columns, no array columns. Contents: catalog proper motions and epochs (Hip, HG, DR2, DR32, DR3) with errors and
  correlations; Gaia DR3 quality terms (RUWE, excess noise, matched transits); noise-term estimates with 1-sigma
  (sig_AL, sig_att_radec, sig_cal); Gaia RV-variability statistics; and `*_msms` columns (best-fit single-star mass,
  metallicity, binary mass ratio q and secondary mass m2, flux ratio, each with an error where given). These are point
  estimates with errors, not sample chains, and none is an orbit posterior. The sidecar file holds only DR2 matched
  observation counts. Judah says full posteriors exist in G23H: ask which columns he means (open).
- eFAP (rvsearch, last lines): fit log10(count) = a x + b to a 10-bin histogram of Delta-BIC from the median to the 95th
  percentile (x = Delta-BIC - median). A = -a ln10 is the decay rate, B = 10^b the amplitude (B is never used afterwards).
  Threshold = median + ln(N_periods/FAP)/A; fap_min (FAP of the best peak) = N_periods exp(-A (max - median)).
  Detection (max > threshold) is equivalent to fap_min < FAP. See the 2026-11-10 chat for the derivation.

### Update 2026-11-10 (blind search: first code, and a finding on eFAP; Judah to decide)
- Written: `blind_search/{orbit,data,model,efap,periodogram,plots}.py` and `tests/test_blind_search.py`
  (23 checks, all pass; Octofitter-free; V1 vs Octofitter not yet done). Circular periodogram first (e = 0,
  M0 = 0 absorbed by the amplitudes); new-companion BIC cost k_new = 1 + 4 (astrometry) + 2 (RV).
  Speed: 564 periods x 93 astrometric transits in 0.06 s; 9221 periods x 200 RVs in 1.2 s.
- Finding: with plain eFAP, a clear astrometric signal is NOT detected. On the Gaia-4 scan geometry a
  circular companion with Delta-BIC ~ 60 (S/N ~ 9) has threshold 127; at Delta-BIC ~ 600 the threshold is
  ~800-1300. The strong signal's own sidelobes (uneven, clumped scan sampling; 93 points; 4 free
  amplitudes) populate the 50th-95th percentile range the tail is fitted to, flattening the slope (A ~ 0.01-0.1
  instead of ~0.5). Noise-only: 0 false alarms in 100 realizations (the threshold errs high). The RV
  tests pass with plain eFAP at the injected strengths (K ~ 11 m/s, 200 RVs).
  Prototype of one fix (not adopted): compute the threshold from the periodogram of the residuals after
  removing the best peak (prewhitening), then compare the ORIGINAL maximum with it. Circular injections:
  Delta-BIC 580 and 40 vs thresholds 47 and 32 (detected); S/N ~ 5 still not detected. An e = 0.4
  injection needs the peak removed with a full (P, e, M0) refit, not a circular one (peak 33.5 vs threshold 65
  with a circular removal). Options to decide: (a) plain eFAP (misses strong signals); (b) prewhitened
  threshold with the full refit of the best peak (deviates from RVSearch; keeps the property that a system
  with genuine noise power at many periods gets a high threshold); (c) fixed Delta-BIC floor as well.
- Update 2026-11-10 (later): Kepler solver: `orbit.kepler_E` now uses the project's C solver (`fits_rv/kepler_c.py`, numpy
  fallback; agreement 5e-14). Measured: C 2.8 ms vs numpy 3.7 ms per 20,000 solves; the circular periodogram needs no
  solve. Plots: eFAP fit bins are coloured bars; new `plot_injection` (original vs injected data, baseline removed).
- Threshold calibration on Gaia-4 scan geometry (circular injection at P = 300 d, 30 noise draws per strength, period
  within 3%): plain eFAP detects 0% at S/N 5.5-16.6 (noise-only: 0 false alarms in 300). Prewhitened (circular removal of
  the best peak): detects 33% at S/N 5.5, 50% at 6.9, 83% at 8.2, 100% at >= 11, but 4 false alarms in 300 noise-only
  periodograms (1.3% vs the nominal 0.1%): too liberal, because the removed peak takes some noise power with it and
  lowers the threshold. So the prewhitened rule needs calibration before use. Alternative on the table: a Monte Carlo
  null (simulate noise-only data at the real epochs and scan angles, with the fitted jitter, and take the 1 - FAP
  quantile of the maximum Delta-BIC); about 0.05 s per realization for astrometry, calibrated by construction for white
  noise but blind to red noise unless residuals are resampled.
- Update 2026-11-10 (Gaia-4 b injection; Judah's requests): `data.load_gaia_dr4` now defaults to the tutorial reduction
  (PLAN.md decision; median error 0.0376 mas, matching the 0.038 quoted there; `reduction='demo'` gives 0.082).
  `blind_search/injection.py` injects physical companions; GAIA4B = published P 571.3 d, e 0.338, omega 180.3 deg,
  K 368 m/s, 11.8 M_J, host 0.644 M_sun (catalog), parallax 13.58 mas (73.6 pc, the demo fit's), Tp from the catalog;
  inclination 115.3 deg derived from K and the true mass (i > 90 deg branch); Omega 2.96 rad = median of Judah's demo joint
  fit (the catalog has none). Result: astrometric semi-major axis 273 uas (catalog 272), a = 1.17 AU. Sky orientation
  unverified until V1. Injected into simulated noise at the real epochs and errors: S/N 37.7, peak 562 d (571 true;
  grid and resolution), dBIC 1383 vs plain eFAP threshold 2175 (missed).
  REAL Gaia-4 pre-release data (no injection): peak 562.3 d, tutorial errors dBIC 1188 vs plain threshold 1591 (missed),
  prewhitened threshold 84 (detected); demo errors dBIC 183 vs 302 (missed) vs 11 (detected). So plain eFAP fails on the
  one clear real detection we have. Figures in `blind_search/figures/` (gaia4b_injection*.png, gaia4_real_*.png).

### Update 2026-11-10 (prewhitened threshold implemented and evaluated; Judah's request)
- `blind_search/threshold.py`: `refit_peak` (variable projection over P, e, M0; Nelder-Mead from 6 starts; P within 10%,
  e <= 0.7; 0.4 s per call) and `prewhitened_efap` (eFAP on the periodogram with the refitted best peak in the model;
  the ORIGINAL maximum is compared with that threshold). Gaia-4 b injection refits to P 569.4 d, e 0.27 (true 571.3, 0.338).
- Injection test redone as asked: data = realistic reference motion (baseline-only fit to the real Gaia-4 data) + noise +
  ONE injected Gaia-4 b; the true reference motion is subtracted for plotting; data minus model is pure noise
  (chi2 80 for 93 transits; `gaia4b_injection_data.png`, a closed ellipse in the sky plane).
- Pitfalls, measured (Gaia-4 geometry, tutorial errors):
  * White noise, 300 draws, FA nominal 0.1%: plain eFAP 0 false alarms; prewhitened with circular removal 7 (2.3%);
    with the full refit 12 (4.0%). Median thresholds 21, 4.7, 2.5 dBIC. The removed peak takes noise power with it, so
    the residual periodogram runs low. White-noise maxima are small anyway (max over 300 draws 13.8 dBIC), so a floor
    dBIC >= 20 removes these false alarms: rule to adopt = detect if max > max(20, prewhitened threshold).
  * Two companions (Gaia-4 b S/N 38 + circular 200 d companion), 20 draws: plain eFAP detects neither; prewhitened
    detects Gaia-4 b 20/20 and the second companion (after refitting the first) 16/20 at S/N 8, 20/20 at S/N 15.
  * Correlated noise (6 sinusoids of 0.5-1.0 x mean error, periods 30-600 d; strictly periodic, so arguably signals):
    100 draws, median max dBIC 65. False alarms: plain 0 (median threshold 269), prewhitened 29 (median 96), floor 20
    alone 96, white-noise Monte Carlo (max of 300) 99. The prewhitened threshold adapts partly, not fully; plain eFAP
    adapts fully but is blind to strong real signals. A stochastic (GP-like) noise test is still to do.
  * Other: depends on the candidate it tests (circular); an imperfect removal leaves leakage (conservative); the refit can
    fail (optimizer) or hit the e cap; it characterizes the residual periodogram for the Keplerian model class, not the noise.
- Why RVSearch rarely needs it (my emulation of its eFAP on RV-like geometries, 20 draws each; RVSearch itself not
  run): leakage floor is about (k/N) S^2 for a signal of power S^2 in N points with k amplitudes. Plain eFAP detected
  20/20 at S/N 10 and 37 with 300 RVs over 3000 d (thresholds 25 and 160 vs peaks 91 and 1366), but at S/N 37 with 93
  random epochs and 4 amplitudes 19/20 (borderline), and 0/20 for the real clumped Gaia-4 epochs (at S/N 10 and 37).
  Gaia-4: N = 93, k = 4, clumped epochs, 564 trial periods vs tens of thousands for RVs.

### Update 2026-11-10 (prewhitening not needed: the jitter was the difference; Judah's challenge)
- Judah asked why RVSearch needs no prewhitening. Checked in the RVSearch source (`per_bic`, `_fit_period`): at every trial
  period `maxlike` refits all parameters that are free, including each instrument's jitter (only P, secosw, sesinw are
  fixed), and the baseline model refits jitter too. My periodogram held the jitter fixed. With fixed jitter the
  likelihood gain at a wrong period grows as (amplitudes/N) x signal power ("leakage": power at trial periods other than
  the true one, because the uneven sampling lets a sinusoid at the wrong period correlate with the real signal), which
  swamps eFAP. With jitter refitted, the null model absorbs a strong signal into its jitter, the gain at wrong periods is
  bounded by about the number of fitted amplitudes, and plain eFAP works.
- Result (Gaia-4 geometry, tutorial errors): injected Gaia-4 b (S/N 38): peak dBIC 212 vs eFAP threshold 76 (fixed jitter:
  1383 vs 2175, missed). REAL Gaia-4: 181 vs 76 detected (demo errors, 0.082 mas: 69 vs 78, not detected). Weaker
  injections: S/N 11 detected (61 vs 59), S/N 8 not (27 vs 62). White noise: 0 false alarms in 100 (test) and 40
  draws. So prewhitening (threshold.py, archived in the scratchpad as threshold_prewhitening.py, not in the repo) is
  dropped; `periodogram.delta_bic(..., jit=None)` is the default.
- Real RV test systems, plain eFAP with free jitter, first step (RV-only, circular; period <= 2 x baseline): detected
  3 of 10 (HD 168443: 161 vs 95; HD 191939: 452 vs 341; TOI-1736: 604 vs 51). Missed: HD 10180 (46 vs 150), HD 111232
  (242 vs 442, the peak is the right period), HD 147513 (235 vs 1972, wrong period), HD 39091 = pi Men (2290 vs 7169; peak
  at P/8 of b: a high-e orbit), HD 43691, HD 89707 (30 vs 98, right period), TOI-1694. So eFAP at step 1 is strict on
  real multi-planet / high-e RV systems; not adopted for the RV search (it uses Delta-BIC >= 20), and RVSearch itself was
  not run (not installed; running it needs a pip install, which needs Judah's OK).
- Two injected companions (Gaia-4 b + a physical 8 M_J, 200 d, e 0.1, i 70 deg companion, 92 uas, S/N 11.8 given
  Gaia-4 b), `two_companions_*.png`: step 1 peak 562 d, dBIC 136 vs threshold 119 (detected); step 2 (step-1 companion
  refitted and in the model) peak 202 d, dBIC 73 vs threshold 86: NOT formally detected, but the peak is at the right
  period; a forced refit gives P 201.0 d, e 0.25 (injected 200 d, 0.1); Gaia-4 b refits to 571.9 d, e 0.22.
- Figures redone (`make_figures.py`; stale ones removed): sky-plane points now sit at the model position plus the
  along-scan residual; the 1D panel shows the model at each transit (east sin psi + north cos psi, so no single smooth
  curve passes through it) with the bound sqrt(east^2 + north^2); chi2/N is stated as chi2 divided by N data points.

### Update 2026-11-11 (plot conventions; multi-companion sensitivity)
- Plot conventions (Judah): "SNR" not "signal-to-noise"; reduced chi2 = chi2/dof with dof = N - free parameters
  (`Stack.n_params`: baseline + per companion 3 nonlinear + 4 linear astrometric amplitudes; jitter not counted);
  legends inside the axes; periodogram x-range = searched range (no grey bars); "astrometric SMA";
  step figures show residuals to every model (n-1 and n companions). "eFAP" is RVSearch's name for the empirical
  false-alarm-probability method (Howard & Fulton 2016); plots say "FAP threshold" for the line and "eFAP" for the method.
- Literature on 2-planet astrometric signatures: no paper found that plots the 2D sky-plane signature of a 2-planet system
  as a figure (searched twice). Related: Casertano+2008 (A&A 482, 699), double-blind Gaia test: > 70% of 2-planet systems
  with well-separated periods (0.2-9 yr), astrometric SNR 2-50 and e <= 0.6 correctly identified; pi Men (arXiv:2007.08549)
  uses the reflex motion of the outer planet. Perryman & Schulze-Hartung 2011 is cited (arXiv:1411.1173) for the
  growing complexity of multi-planet signals (not read).
- Multi-companion sensitivity (Gaia-4 geometry, tutorial errors; Gaia-4 b SNR 38 + illustrative companion 2 (200 d) and 3
  (90 d, 6 M_J); SNR k|previous = sqrt(delta chi2) of adding companion k to the true earlier ones; 8 noise draws per cell;
  period within 3%; plain eFAP, jitter refitted; script in the scratchpad `exp7.py`):
  * Companion 2, step 2 detected: SNR 6: 0/8 (peak 15, threshold 58); 9: 0/8 (39, 82); 12: 2/8 (62, 93); 16: 4/8 (94, 96).
    Companion 2 ALONE (no Gaia-4 b) at the same SNR: 0/8, 0/8, 0/8, 5/8 (thresholds 50, 64, 79, 91). Step 1 detects Gaia-4 b
    8/8 up to SNR2 = 12 but only 5/8 at SNR2 = 16.
  * Three companions (companion 2 at SNR 16): step 1 / 2 / 3 detections for SNR3 = 6: 3/3/0; 9: 1/1/0; 12: 1/0/0; 16: 0/0/0
    of 8; step-3 thresholds 47, 68, 77, 152 vs peaks 10, 24, 38, 63.
  * Reading: the data support the second and third companions (peak dBIC 39 at SNR 9, well above the white-noise maximum
    of 13.8 in 300 draws), but the plain eFAP threshold rises with the strength and number of signals (a signal's own
    leakage inflates its threshold, even with jitter refitted; noise only gives about 32), so formal detection needs
    SNR >~ 14-16 for each, and the first companion can stop being detected when others are present. A threshold that
    does not scale with signal power (dBIC >= 20, 0 false alarms in 300 white-noise draws; or a white-noise Monte Carlo
    null) would detect the second companion from SNR ~7-9 on, at the cost of no adaptation to non-white noise.
    Decision needed on the threshold.

### Update 2026-11-11 (Judah: keep the name eFAP; threshold rule still to be decided)
- Keep "eFAP" (RVSearch origin). NEXT TIME the plots are regenerated: label the threshold line "eFAP threshold"
  (not "FAP threshold"); no regeneration needed now.
- SNR in the sensitivity runs was varied by scaling companion 2's (and 3's) Thiele-Innes amplitudes by one factor at fixed
  P, e, orientation and phase: equivalent to changing the companion's mass (about linearly; or the distance) at fixed orbit.
  Companion 2 nominal 8 M_J = SNR 11.8, so SNR 6, 9, 12, 16 correspond to about 4.1, 6.1, 8.1, 10.8 M_J at 73.6 pc.
- RVSearch itself supports both rules: `run_search(fixed_threshold=None)` uses eFAP; a number replaces it with a fixed
  Delta-BIC threshold (search.py lines 443-447). Open: report both "eFAP-detected" and "Delta-BIC >= 20" in the sensitivity runs.

### Update 2026-11-11 (red-noise trial redone with the free-jitter periodogram; `blind_search/red_noise_test.py`)
- CORRECTION: the earlier red-noise numbers (floor 20: 96/100 false alarms; plain eFAP 0; prewhitened 29) came from the
  old fixed-jitter periodogram and no longer apply. Redone (100 draws per kind, Gaia-4 epochs and tutorial errors, jitter
  refitted, no companion): white noise: median highest peak -9.5, median eFAP threshold 17.9, false alarms eFAP 0,
  dBIC >= 20 0. Six sinusoids (30-600 d): median peak 19.8, median threshold 89.2, false alarms eFAP 0, dBIC >= 20 50.
  Ornstein-Uhlenbeck noise (60 d, rms 1.5 x mean error): median peak 17.3, median threshold 87.5, eFAP 0, dBIC >= 20 41.
  max(20, eFAP): 0 in all three. eFAP raises its threshold about 5-fold in correlated noise and never fires; a fixed 20 fires
  on 40-50% of correlated-noise draws. Figures: red_noise_summary.png, red_noise_example_{sinusoids,OU}.png.

### Update 2026-11-11 (are the errors realistic; is the multi-companion limit the threshold or the data?)
- The "tutorial" errors are the real pre-release Gaia-4 per-transit errors reduced by the tutorial method (median 0.038 mas).
  Check on the real data after a one-companion refit (formal errors): Gaia-4 reduced chi2 1.03 (chi2 83, dof 81), ML jitter
  0.010 mas; HD 114762 reduced chi2 0.73 (N 63, median error 0.057 mas). Demo reduction (single-CCD errors, median 0.082 mas):
  Gaia-4 reduced chi2 1.50, jitter 0.043 mas, total per-transit error 0.093 mas, so the tutorial reduction describes the real
  scatter well for these two stars. Within-transit CCD correlations are not modelled by either reduction (open).
- Real-data residuals are white-like: step-2 periodogram after the one-companion refit peaks at dBIC -14 (Gaia-4, eFAP
  threshold 7) and -11 (HD 114762, threshold 8), inside the white-noise range (largest noise peak 13.8 in 300 draws).
  HD 114762 b is detected at step 1 (84 d, dBIC 199 vs eFAP threshold 172).
- Does the pipeline lose information? Two injected companions (Gaia-4 b + 200 d), 8 draws per SNR: median dBIC of the TRUE
  two-companion model over the true one-companion model (7 parameters penalised) 29 / 58 / 91 for SNR2|1 = 9 / 12 / 16;
  periodogram peak with companion 1 at its true parameters 38 / 66 / 100 (eFAP thresholds 65 / 86 / 103); with companion 1
  refitted from the data 37 / 65 / 98 (thresholds 69 / 78 / 82). The periodogram and the refit recover the available
  information; the eFAP threshold is what exceeds the peak.
- Checks proposed: (1) inject companions into the REAL residuals of Gaia-4 (after removing the fitted Gaia-4 b), many
  realizations by random sign flips; (2) per-step parametric-bootstrap threshold (simulate n-1 refitted companions + white
  noise, refit, take the max of the step-n periodogram); (3) injection-recovery grids with more draws, period ratios, and
  both rules; (4) Epochalypse two-companion systems with known truth; (5) more real stars.

## Update 2026-11-11 (code audit after the fixed-jitter bug; flowchart image: blind_search/search_flowchart.png, drawn by blind_search/make_flowchart.py)

Checks run today (scratchpad audit1-5.py), all passed:
- `fit_free_jitter` lnL (4 iterations) equals 200 iterations and a brute-force profile likelihood (exact weighted least
  squares at each jitter, 1D maximisation over jitter) to 1e-4 at 5 periods on real Gaia-4, white noise and Gaia-4 b
  injected data; the dBIC error from the 4-iteration warm start is below 0.001.
- Null calibration: the trial columns are full rank (trace of the fitted projector 4.000; smallest/largest singular value
  0.03 minimum); fixed-jitter white-noise mean 2 dlnL = 4.06 +- 0.055 over 400 draws (4 expected); free-jitter mean is
  about 0.2 lower on identical draws (the null absorbs some noise into jitter), as expected.
- `efap.py` gives the same threshold and fap_min as RVSearch's own `eFAP` source (to 6 digits) on three periodograms.
- Multi-companion bookkeeping: in 16 two-companion step-1 peaks (SNR2|1 = 12, 16) the peak was Gaia-4 b every time; at
  SNR2|1 = 16 the peak (100-132) sits at the eFAP threshold (83-124), i.e. the misses are threshold self-masking, not
  bookkeeping.
Not verified (gaps): orientation/sign conventions against Octofitter (V1); injection and recovery share orbit.py;
the loop is not a tested function (no search.py); joint astrometry+RV path (k_new = 7) untested end to end;
`plot_step` residual panels use a zero-jitter fit while the search uses free jitter (small for Gaia-4);
`plot_injection(recovered=...)` slices amplitudes with stride 4 (wrong if RV channels are present).
2D plots: the origin (0, 0) is the system barycentre (the reference-motion track after subtracting the true reference
motion), not the star; the star is the moving point and, for e > 0, the barycentre is at the orbit's focus.

## Settled 2026-11-11 (period degeneracies; format of the flowchart; glossary)
- No special handling near 1 yr (parallax) or at long periods (proper motion): let those signals be fit, judge their
  reliability later. This replaces the "not special-cased" gap noted above.
- The flowchart is a PNG (`blind_search/search_flowchart.png`), drawn by `blind_search/make_flowchart.py`; the Markdown
  version was deleted.
- Terms are defined when first used (CLAUDE.md working rule). Glossary of terms used in the 2026-11-11 audit:
  - module: a Python file (e.g. `blind_search/plots.py`); function: a named block of code inside a module; class: a
    function-like template that holds data and methods (e.g. `Stack` in `model.py`); script: a file run directly
    (e.g. `make_figures.py`); scratchpad: a temporary folder outside the repo for experiments (`exp7.py` lives there).
  - `plot_step` (function, `blind_search/plots.py`): draws one search step's figure (periodogram, eFAP tail fit,
    residual panels). `plot_injection` (same file): draws the injection test (2D sky panels, 1D along-scan panel).
  - stride 4: a list of numbers is stored flat, one companion's four amplitudes (A, B, F, G) after another. The code finds
    companion k's block by jumping 4 entries per companion (`coef[k_base + 4*k : k_base + 4*k + 4]`), the "stride". If an
    RV channel is present each companion has 6 entries (4 astrometric + 2 RV), so a stride of 4 would read the wrong
    numbers. Not fixed yet because `plot_injection` is only used for astrometry.
  - null model: the model without the companion being tested. trial model: the null model plus one companion at a trial
    period. Delta-BIC: 2 x (gain in ln-likelihood) minus a penalty of (number of new parameters) x ln N.
  - projector trace: for a least-squares fit, the sum of the fitted-value weights; equals the number of independent
    columns, and is the average chi-squared gain expected from pure noise (4 here).
  - profile likelihood: the likelihood after maximising over every other parameter at each value of the one of interest
    (used as the brute-force check on the jitter fit).
  - Thiele-Innes constants A, B, F, G: four linear numbers (mas) that fix the star's orbit size and orientation on the sky
    once P, e and M0 are known. barycentre: the system's centre of mass. jitter: extra error added in quadrature to each
    transit's error. self-masking: the eFAP threshold rising when strong signals fill the periodogram.
    V1: the planned equivalence test of this code's model against Octofitter's.

## Update 2026-11-11 (flowchart redesigned)
`blind_search/search_flowchart.png` now uses LaTeX-typeset mathematical notation instead of code names (code pointers are
small grey tags), and each step has two parts: the flowchart box (what it does, inputs and outputs, decision) and a grey
background panel to its right (definitions and reasons). Redraw with `python -I blind_search/make_flowchart.py` (needs
`/Library/TeX/texbin` on PATH; about 30 s).

## Update 2026-11-11 (does the jitter handling match RVSearch? read-only comparison)
Sources read: `rvsearch_src/{periodogram,search,utils}.py` (scratchpad copy of RVSearch) and radvel 1.3.8 `likelihood.py`,
`fitting.py` (conda env `radvel`; RVSearch's `utils.py` mentions radvel >= 1.4, not checked).
Matches:
- Jitter is a free parameter at every trial period in the trial model (`Periodogram.per_bic` -> `utils.maxlike` fits every
  parameter with `vary=True`: gamma, jitter, tc, k; only per, secosw, sesinw are fixed). `initialize_default_pars` never
  fixes the jitter; each trial starts from the previous best-fit jitter (`default_pdict`).
- The null model's BIC (`basebic`) comes from a full fit with jitter free, computed once; the null does not depend on the
  period, so one fit equals refitting at every period.
- Likelihood: `loglike_jitter` = -1/2 sum[r^2/(sigma^2 + s^2)] - sum ln sqrt(2 pi (sigma^2 + s^2)), identical to ours.
  Negative jitter is flipped (only s^2 enters), equivalent to our s^2 >= 0.
- Jitter is counted as a varied parameter in both models, so it cancels in dBIC (as in ours, where it is counted in neither).
Differences:
1. Penalty: RVSearch's trial model adds k and tc (period fixed), so dBIC = ... - 2 ln N (circular; its `floor` is -2 ln N).
   Ours subtracts 5 ln N for astrometry (period + 4 amplitudes), i.e. one ln N (4.5 for N = 93) more than RVSearch's convention
   would give (4 ln N). A constant shift does not change eFAP detections (threshold shifts equally) but does change a fixed
   dBIC >= 20 rule by 4.5. Undecided.
2. Known companions: in RVSearch's trial fits the known planets' per, tc, k, secosw, sesinw also stay free (set by
   `fit_orbit`), so they are re-optimised at every trial period. Ours holds their (P, e, M0) fixed and refits only their
   amplitudes. Expected effect: small; not measured.
3. Optimiser: RVSearch uses Powell (a derivative-free local optimiser) with up to three restarts when dBIC < floor - 1; ours
   solves the linear parameters exactly and the jitter by maximum-likelihood iteration (verified against a brute-force
   profile likelihood), so ours cannot stop at a local optimum at a given period.

## Update 2026-11-11 (benchmark: known companions fixed vs re-optimised at every trial period, RVSearch convention)
Added `refit_known=True` (default False) to `periodogram.delta_bic` (function in `blind_search/periodogram.py`): the known
companions' (P, e, M0) are re-optimised with Nelder-Mead at every trial period, and the null model likewise; tests added
(`test_refit_known`, all pass). Paired benchmark on Gaia-4 geometry (same data, same refits; only the periodogram differs;
scratchpad `bench_known.py`, 65 tasks, 11 cores, 44 min wall):
- Real Gaia-4, step 2: peak -13.6 (fixed) vs -12.7 (free); eFAP threshold 6.7 vs 6.5. No difference.
- Gaia-4 b injected alone, 12 draws, step 2: false alarms 0/12 both; median threshold 10.7 (fixed) vs 18.8 (free).
- Two companions (Gaia-4 b + 200 d), 12 draws per SNR, step-2 hits (fixed vs free): SNR2|1 = 6: 0 vs 0; 9: 0 vs 0; 12: 3 vs 1;
  16: 5 vs 4. Paired medians: peak +1.1 / +1.8 / +2.7 / +4.9, threshold +9.6 / +8.9 / +15.7 / +21.1, so (peak - threshold)
  worsens at every SNR (fixed -41.6 / -36.8 / -18.1 / -0.7; free -53.7 / -45.8 / -35.2 / -12.9).
- Step 3 (|K| = 2, SNR3 = 12, 4 draws): hits 0/4 both; median peak 36 vs 56, threshold 75 vs 138.
- Time per periodogram (564 periods): step 2, 1.3 s fixed vs 315 s free (about 235x); step 3, 1 s vs 2153 s (about 2000x).
Reading: re-optimising the known companions raises the peaks a little but raises the eFAP threshold more, so it does not
help detection and is far slower. Cause of the threshold rise not measured (likely the whole periodogram body rises because
the known companions soak up leaked signal at every period). Default stays fixed-known pending Judah's decision.
Penalty convention (5 ln N vs RVSearch-style 4 ln N for astrometry) still undecided; explained 2026-11-11: RVSearch's
circular RV trial has 2 new numbers (K, phase) and P is a fixed grid value; the astrometric circular trial has 4 new
numbers (A, B, F, G: ellipse axes, orientation, phase) and we additionally counted P.

## Update 2026-11-11 (diagnosis: why eFAP misses companions, even a single one; scratchpad diag_threshold.py)
Medians over draws (16 per row; 150 for noise only); Gaia-4 geometry, white noise at the formal errors:
| case | peak | eFAP threshold | body p50 / p95 | decay rate A |
| noise only, step 2 (Gaia-4 b refit in the model) | -11.0 | 14.9 | -19.7 / -15.3 | 0.385 |
| two companions, SNR2|1 = 9, step 2 | 37.6 | 77.3 | -16.0 / -0.3 | 0.143 |
| two companions, SNR2|1 = 12, step 2 | 61.9 | 88.3 | -15.4 / 3.4 | 0.128 |
| companion 2 ALONE, SNR 9, step 1 | 40.6 | 67.1 | -17.0 / -3.3 | 0.158 |
| companion 2 ALONE, SNR 12, step 1 | 66.8 | 79.7 | -16.4 / 0.0 | 0.138 |
| companion 2 ALONE, SNR 16, step 1 | 99.1 | 95.0 | -15.9 / 2.9 | 0.119 |
- Empirical noise-only step-2 highest peak over 150 draws: median -11.0, 99th percentile 1.3, maximum 2.5 (eFAP threshold
  there: 14.9, 0 false alarms). Companion 2 peaks at the right period (within 3%) in 16/16 draws at SNR 9 and 12, 16/16 above
  the empirical 99% noise level.
- Companion 1 refitted from the data vs at its true parameters gives the same threshold (77.3 vs 78.1; 88.3 vs 90.4): imperfect
  removal of companion 1 is NOT the cause. A single companion alone at SNR 12 is detected by eFAP in only 1/16 draws although
  its peak (67) is at the right period every time: the failure is not specific to multiple planets.
- Mechanism: a real signal lifts the upper part of the periodogram body (p95 from -15 to about 0..3), which flattens the
  fitted tail (A from 0.385 to 0.12-0.16); threshold = median + ln(N_P/FAP)/A with ln(564/0.001) = 13.2, so a 3x smaller A gives
  a 3x larger distance above the median. (That the lift is leakage of the signal into the body is inferred, not measured
  directly.)
- Caveat for any calibrated or fixed threshold: this is white noise only; correlated noise gave 41-50% false alarms for a plain
  dBIC >= 20 in the red-noise trial, which eFAP suppresses.

## Update 2026-11-11 (how a real signal raises the eFAP threshold although the fit uses only the 50th-95th percentiles; diag_skirt.py, diag_leak.py)
The peak itself is not the cause: only 6-7 periods lie within 1/T of the true frequency with dBIC > 2.5, and removing all
periods within 2/T changes the threshold by at most 3 (SNR 9: 67.1 -> 63.5; 12: 79.7 -> 81.8; 16: 95.0 -> 97.5).
The cause is leakage at distant periods. The 95th-percentile cut is by rank: it removes only the top 5% = 28 periods however
large their values. A noise-free SNR-12 signal (pure leakage) puts 61 of 526 periods farther than 3/T from the signal at
dBIC > 0 (strongest near 38 d, dBIC 47, and 20.7 d, dBIC 25, for a 200 d signal), so more than 28 leakage periods remain in
the fitted range. Mean dBIC at rank percentiles 50/60/70/80/90/95, noise only: -18.8, -18.1, -17.3, -16.2, -14.5, -12.9;
SNR-12 signal + noise: -16.4, -14.8, -12.7, -9.8, -4.6, +0.6. The fitted range widens from about 6 to about 17 in dBIC, the
histogram bins widen with it, and A falls. Keeping only the 50th-80th percentile instead lowers the threshold only slightly
(SNR 9: 67 -> 72, i.e. no help; 12: 80 -> 71; 16: 95 -> 88), so tuning the cut is not a fix.

## Glossary additions 2026-11-11 (simulated-threshold experiments)
- reference motion: the star's own motion without companions, 5 numbers (position 2, proper motion 2, parallax 1). In the injection
  tests: the baseline-only (no companion) fit to the real Gaia-4 data, used as the fixed true values.
- injected Gaia-4 b parameters: the catalogue values (Stefansson+2025: P 571.3 d, e 0.338, K 368 m/s, 11.8 M_J, host 0.644 M_sun,
  parallax 13.58 mas; inclination derived from K and mass; Omega 2.96 rad from Judah's demo joint fit), taken as the truth of the
  simulated world. Sky orientation unverified until V1.
- "noise-only" dataset (imprecise): a fully synthetic dataset with no companion beyond Gaia-4 b: reference motion + injected
  Gaia-4 b + Gaussian noise at the formal errors, on the real epochs, scan angles and errors. No real measured positions are used.
  Prefer the name "dataset with no second companion".
- "one fit" behind the 200 fakes: one null-model fit per starting dataset (Gaia-4 b refitted + jitter); the 200 fakes are drawn
  around it and each refits Gaia-4 b from its own data. Stage 1 used 20 starting datasets (20 fits, 4000 fakes).

## Terminology fix 2026-11-11 (Judah)
- Do not call the model with the already-found companions the "null model": that sounds like no planets. Say "k-planet model"
  (reference motion + k companions + jitter); at step 2 it is the "1-planet model". The code function is still named `null_model`
  in simthreshold.py (docstrings updated).
- "Drawn around the fit" is replaced by: at each real transit, take the fitted model's predicted along-scan position and add an
  independent Gaussian random number of standard deviation sqrt(err^2 + s^2); a fake dataset is one such set of 93 numbers, and
  each of the 200 fakes uses new random numbers.

## Results 2026-11-11 (simulated per-step threshold, M = 200, Gumbel tail; Judah away; scratchpad stage1/stage2/diag_metric1/stage3)
Plan and metrics: PLAN.md 2026-11-11. Module `blind_search/simthreshold.py` (tests pass). Time: about 35 s per step-2 threshold at M = 200 on 11 cores.
- Metric 3 PASS: 10 datasets with no second companion, simulated 95th/99th percentiles median -4.5/-1.3 vs reference -2.8/+1.3
  (differences -1.7/-2.5, scatter about 1.0); Gumbel FAP 1e-3 median 8.8 vs reference 10.6. The simulated thresholds sit 1-2.5 low
  (slightly anti-conservative); the reference itself is noisy (99th percentile from 1-2 of 150 draws). Fakes skip the step-1 periodogram
  and refit from the fitted values (single polish): possible cause, not tested.
- Metric 4 PASS: bootstrap sd of the Gumbel FAP 1e-3 threshold 1.1 at M = 200; M = 100 differs from M = 200 by a median 0.6, M = 50 by 1.2,
  so M = 100 suffices (half the cost).
- Metric 5 PASS: 0 of 20 datasets (Gaia-4 b injected, no second companion) exceed their simulated threshold at FAP 0.01 or 1e-3 (max real
  peak 1.7 vs threshold 4.6 for that dataset).
- Metric 6 PASS: real Gaia-4: step 1 peak 181 (562 d) vs threshold 16.1, step 2 peak -13.6 vs 13.8 (88% of simulated maxima are >= it);
  real HD 114762: step 1 peak 199 (84.1 d) vs 18.6, step 2 peak -11.4 vs 7.3 (62% of simulated maxima >= it).
- Metric 2 PASS: lone companion (step 1) SNR 9 and 12: 16/16 and 16/16 detected at FAP 1e-3 (eFAP: 0/16 and 1/16); second companion at
  step 2 SNR 9, 12, 16: 16/16, 16/16, 8/8 (eFAP step 2: 0, 3, 4 of 12 at SNR 9, 12, 16); the peak is at 200 d every time; step 1 finds
  Gaia-4 b 16/16 at every SNR (eFAP 5/8 at SNR 16).
- Metric 1 FAIL as run, PASS with a variant: with the fakes' noise at sqrt(err^2 + s^2), s = the 1-planet model's fitted jitter, the
  threshold with an unmodelled companion is 17 (SNR 9, 12, 16: 17.2, 17.0, 16.8) vs 8.9 with none (+8; criterion +5). The rise does not
  grow with SNR. Cause: that jitter (0.029-0.041 mas against a formal median 0.038) absorbs part of the unmodelled signal. With formal-error noise only: 10.1 at
  SNR 9 and 12 vs 8.6 with none (+1.5, pass). Caveat: formal-only assumes no real extra jitter (fine for the tutorial reduction: real Gaia-4
  ML jitter 0.010 mas, reduced chi2 1.03; not for the demo reduction, jitter 0.043 mas). Undecided; code default is noise='fitted'.
- Sensitivity floor (6 draws): SNR 6: peak median 14.5; detected 2/6 (fitted-jitter noise) or 3/6 (formal noise); eFAP gave 0/12.
- Metric 7 (red noise, 10 draws per kind, no companion, step 1): simulated white-noise threshold false alarms: sinusoids 5/10 at FAP 1e-3
  (10/10 at 1e-2), OU 2/10 (6/10); eFAP 0/10; dBIC >= 20: 5/10 and 2/10. As expected the simulated white-noise threshold gives no protection
  against correlated noise (same as a fixed floor). The two real stars checked have white-like residuals.
Open: noise setting for the fakes (fitted vs formal); correlated-noise protection (e.g. fakes with block-bootstrapped real residuals, or keep eFAP
as a cross-check for flagging); more draws at SNR 6-8; whether the fakes' slight low bias (1-2.5) matters; update the flowchart Step B once decided.

## Check 2026-11-12 (cursory): Octofitter "transdimensional" planet-count modelling vs the iterative blind search
Read locally: Octofitter v9.0.0 docs and source (`~/.julia/packages/Octofitter/AyX9v`), plus the G23H paper (arXiv 2602.00235).
- What exists: `docs/src/index.md` lists "multiple planets (zero or more); marginalize over the number of planets (transdimensional
  models)". The only implementation shown (docs/src/g23h-example.md) is a spike-and-slab construction, not reversible-jump: a system-level
  `n_planets ~ truncated(NegativeBinomial(), upper = num_planets)`; each planet has `planet_present = n_planets >= planet_i` and
  `mass = mass' * planet_present` (mass 0 when absent; its other parameters then sample the prior and are inert). The maximum number of planets is
  fixed. Discrete variables need Pigeons (non-reversible parallel tempering, slice-sampler explorer; HMC cannot move them and the variational
  reference is switched off). Outputs: posterior of n_planets, Bayes factors from the indicator's posterior mean (p/(1-p)), and stepping-stone evidence
  ratios between fixed-count models (docs/src/rv-multi-planet.md). docs/src/completeness.md shows completeness maps using such an indicator
  (`planet_present ~ Bernoulli(0.5)`, detection if BF > 3 or 10), with a Gaia DR4 example (examples/completeness_dr4, Slurm-scale).
- Use in papers (cursory, web search found little): G23H (Sec. III.7) includes one-companion and zero-companion models in one Pigeons run with a binary
  indicator, about 30 min on a 16-core node for the combined astrometry; its multi-companion runs (14 Her: 2- and 3-companion fits, two-planet RV +
  astrometry fit) use a fixed number of companions, and no injection-recovery test of an inferred planet COUNT was found. No paper found that infers the
  number of planets with the `n_planets` NegativeBinomial variable. The Octofitter paper (arXiv 2402.01971) describes multi-planet systems as added
  planets in the model.
- Assessment: not a replacement for the blind search as the first-pass population search (cost, period multimodality, unvalidated sensitivity and
  false-inclusion rates, prior-dependent Bayes factors); a sensible final confirmation and posterior stage on selected systems. Untested here.

## Update 2026-11-12 (comparison figures: simulated threshold vs eFAP)
Figures in `blind_search/figures_compare/` (drawn by `blind_search/make_comparison_figures.py`): `compare_detection.png`,
`compare_peak_vs_threshold.png`, `compare_falsealarm_real_cost.png`. The scratch folder was cleared between sessions, so the per-draw files of the first
batch are gone; the figures use the counts and medians recorded above (hardcoded in the script with their provenance) plus a new lone-companion run
saved in `blind_search/results/solo_extra.jsonl` (M = 100, same draws for eFAP and the simulated threshold).
Lone companion, step 1, detections (eFAP / simulated fitted-jitter noise / simulated formal-error noise): SNR 4: 0/9, 0/9, 1/9; SNR 5: 0/12, 4/12, 6/12;
SNR 6: 0/12, 6/12, 7/12; SNR 7.5: 0/12, 10/12, 12/12; SNR 9: 0/16, 16/16, -; SNR 12: 1/16, 16/16, -; SNR 16: 8/16 for eFAP (simulated not run). The simulated
threshold is flat near 15-17 (fitted noise) or 11.3 (formal noise) at every SNR; eFAP's threshold rises 40, 42, 48, 59, 67, 80, 95 for SNR 4 ... 16.
Half-detection SNR: about 6 (simulated) vs about 16 (eFAP). eFAP's advantage: no false alarms on correlated noise (0/10 on both red-noise kinds, where
the simulated threshold gives 5/10 and 2/10) and speed (0.3 s vs about 35 s per step).

## Update 2026-11-12 (summary figures for the simulated threshold, real Gaia-4)
`blind_search/figures_sim/gaia4_real_step1_sim.png` and `gaia4_real_step2_sim.png`, drawn by `blind_search/make_sim_figures.py` (function `plot_step_sim` in
`blind_search/plots.py`, the analogue of `plot_step` for eFAP; `simthreshold.simulated_threshold(..., keep=N)` now also returns N example fake periodograms).
Layout: top, the real periodogram with 20 of the 200 fake periodograms (grey) and the highest peak of each (grey dots), the threshold and the result line
(DETECTED / NOT DETECTED, how many fakes reach the real peak); bottom left, how the threshold is built (histogram of the 200 highest peaks, Gumbel fit,
thresholds at FAP 0.01 and 0.001, real peak); bottom right, residuals. Results: step 1 peak 181.3 (562 d) vs threshold 16.1 (0 of 200 fakes reach it), fitted jitter
0.118 mas; step 2 peak -13.6 vs threshold 13.8 (175 of 200 fakes reach it), fitted jitter 0.010 mas. Takes about 2 minutes.

## Results 2026-11-12 (iterative search with the simulated threshold on injected systems; blind_search/search.py, run_injection_tests.py)
Figures: `blind_search/figures_search/` (`gaia4_plus_200d_step{1,2,3}.png`, `hd114762_plus_2_step{1,2}.png`, `method_diagram.png`); summaries in `blind_search/results/`.
- Gaia-4 + 200 d companion (synthetic data on the real epochs/errors; Gaia-4 b SNR 37.8, the 8 M_J 200 d companion SNR 12.2 given Gaia-4 b), M = 200: step 1 peak 136.4 at 562.3 d vs
  threshold 16.5 (Gaia-4 b found); step 2 peak 72.7 at 201.7 d vs 13.5 (200 d companion found); step 3 peak -11.5 vs 5.7 (stop). Both found, no spurious. Ten more noise draws (M = 100):
  10/10 found both, 0 spurious.
- HD 114762 (REAL DR4 pre-release data, which contain HD 114762 b, plus two synthetic planets added to the measured positions; host 0.80 M_sun (Stassun+2017; archive range 0.80-1.07),
  parallax 26.2 mas (Octofitter docs); orientation angles arbitrary): planet 1 = 5 M_J, 400 d, e 0.2, i 75 deg, astrometric SMA 154 uas, SNR 11.4 alone (10.8 given b); planet 2 = 1 M_J, 1095 d, e 0.1,
  i 70 deg, 60 uas, SNR 3.6 alone (2.7 given b). Main run M = 200: step 1 peak 201.9 at 84.1 d vs threshold 22.2 (HD 114762 b); step 2 peak 14.3 at 51.5 d vs 17.4: not detected. Planet 1 and planet 2
  both missed. Six random-orientation draws (M = 100; planet-1 SNR 13.6, 12.7, 8.6, 8.6, 7.7, 10.4; planet-2 SNR 3.6-5.4): only HD 114762 b found in 5 of 6; draw 0 also detected a peak at 387.2 d, 3.2% from the
  injected 400 d (outside the 3% matching tolerance, so scored spurious; probably planet 1).
- Diagnosis (main HD dataset): the step-1 refit of b is biased by the unmodelled planet 1 (refit e 0.318, M0 6.17 vs 0.344, 5.92 when b is fitted to the real data alone). With that b fixed, adding the
  true planet-1 orbit gains chi2 = 31 (noise-free injected signal alone: 117; real data alone with this b: 47, i.e. destructive interference), and the fixed-known step-2 periodogram has dBIC 1.2 at 400 d.
  Re-optimising b at every period (`refit_known=True`, 96 s) gives a highest peak 18.8 at 401.0 d (vs 14.3 at 51.5 d with b fixed); a joint refit of b + planet 1 started at 400 d fits well
  (reduced chi2 0.86; planet 1 at 405.0 d with e at the 0.7 cap). The threshold for the free-known statistic was not computed (its fake datasets would cost about 100 s each). Conclusion: the greedy search
  fails here because the strongest companion's step-1 fit absorbs distortions from a weaker unmodelled signal; this is not a thresholding problem. The 400 d period is close to 1 yr (365 d) but the 4200 d beat
  period exceeds the 1798 d baseline; parallax absorbs essentially none of the signal (noise-free gain 130 without the parallax column vs 129 with it).
- Options, none implemented: (1) joint-refit the top K periodogram peaks (all known companions free) and use the maximum refit dBIC as the statistic, with the same procedure in the fakes (about 20-40 s per
  fake); (2) a final joint backward refit after each detection; (3) accept the sequential limit and report it. 
- Method diagram: both thresholds are FAP thresholds (FAP = 1e-3); the difference is where the noise distribution comes from: eFAP = empirical (tail of this periodogram, union bound); simulated = model-based
  (distribution of the highest peak of fake datasets, Gumbel fit). Measured costs in the diagram: eFAP about 0.3 s per step; simulated about 35 s per step at M = 200 on 11 cores (about 2 s per fake).

## Results 2026-11-12 (threshold rules on real Gaia-4 + one injected planet; run_threshold_comparison.py, make_threshold_comparison_figure.py)
Setup: REAL Gaia-4 data (contains Gaia-4 b) plus one synthetic planet added to the measured positions: P = 1428 d (2.5 x Gaia-4 b), e 0.2, 9.8 M_J (main case), astrometric SMA 420 uas,
SNR 56.6 = 1.5 x the SNR of Gaia-4 b injected with the catalogue parameters (37.8); orientation arbitrary (main) or random (10 draws, mass re-solved each time for SNR 56.6: 7.0-11.2 M_J).
Rules (iterative search, search.py): eFAP (FAP 1e-3), simulated threshold (FAP 1e-3; M = 200 main, 100 in the draws), fixed dBIC >= 20. A detection matches a truth period when the frequencies agree
within half a resolution element of the 1805 d baseline (a fixed 3-10% tolerance mislabelled the injected planet's 1586 d step-1 peak as spurious; with the new rule the 387.2 d detection in HD 114762 draw 0 also matches the
injected 400 d).
- Main case: all three rules find both companions, no spurious: step 1 (injected planet) peak 93.3 at 1586 d vs thresholds eFAP 90.5 / simulated 13.7 / floor 20; step 2 (Gaia-4 b) peak 129.8 at 562.3 d vs 102.6 / 15.6 / 20;
  step 3 peak -11.0 vs 39.9 / 8.3 / 20 (stop). eFAP passes the step-1 peak by only 3%.
- Ten random-orientation draws: simulated 10/10 and floor 10/10 find both companions with 0 spurious; eFAP finds NOTHING in 10/10 (the search stops at step 1, so not even Gaia-4 b is found): step-1 peak 80.6-111.1
  against eFAP thresholds 105.9-179.0 (peak/threshold 0.46-0.96); the eFAP decay rate A falls to 0.07-0.11 (0.125 in the main case). Simulated step-1 thresholds 11-17 in all cases.
- Over the 11 cases: both found by eFAP 1/11, simulated 11/11, floor 11/11; spurious detections 0 in all.
- Reading: two strong signals (SNR 35 and 57) lift the periodogram body through leakage and flatten the eFAP tail, so eFAP's threshold passes the peaks (self-masking). This scenario does NOT discriminate the simulated
  threshold from the fixed floor of 20: both detect everything. Their known differences (floor does not adapt to non-white noise or to the number of trial periods; simulated threshold is calibrated for white noise at the chosen FAP,
  costs about 35 s per step) need other tests: low SNR (floor 20 vs simulated near the SNR 5-7 floor), correlated noise, and the biased-first-fit case (HD 114762).
- Files: figures_search/gaia4_real_inj_{efap,sim}_step{1,2,3}.png, gaia4_real_inj_summary.png; results/gaia4_real_inj_{main,repeats,step1_efap,step1_sim}.json.

## Settled 2026-11-12 (noise model of the simulated threshold)
Judah: keep Gaussian white noise for the fakes (sqrt(err^2 + s^2) per transit, fitted jitter s). Later, account for the harmonics in the Gaia noise noted by the Gaia team (not yet specified here; a source is needed).
The comparison tests so far support: the simulated threshold (and a fixed floor of 20) miss obviously detectable planets much less often than eFAP (self-masking); they do not yet show an advantage of the
simulated threshold over the fixed floor, and eFAP's advantages (no false alarms on correlated noise in the 10-draw red-noise trial; 0.3 s per step) are unaffected.

## Results 2026-11-12 (why the simulated threshold does not rise with correlated noise; red_noise_fakes_test.py)
Answer: the 200 fakes are made with WHITE Gaussian noise (the decided noise model), so they cannot contain the correlated signals however strong they are in the real data. eFAP reads its threshold from the data's own
periodogram, which does contain them (empirical vs model-based FAP). The earlier red-noise trial used the same correlated-noise datasets for every rule; only the fakes differed.
Test (Gaia-4 epochs/formal errors, NO companion, 20 datasets per kind, step 1, FAP 1e-3, M = 100): median highest peak / median thresholds eFAP | white fakes | block-bootstrap fakes | oracle fakes; false alarms of 20 each:
- sinusoids: peak 20.2; 86.8 | 15.8 | 31.9 | 89.4; false alarms eFAP 0, white 12, block 2, oracle 0.
- OU (60 d): peak 12.5; 76.2 | 15.9 | 40.4 | 80.7; false alarms eFAP 0, white 7, block 1, oracle 0.
- Block-bootstrap fakes: moving blocks of 8 consecutive transits of the baseline-fit residuals (simthreshold noise='block'); oracle fakes: white noise plus a fresh realization of the true red-noise process (red_noise_test.OracleNoise).
- Why the white fakes barely move: the fitted jitter only inflates the white level (0.044-0.050 mas vs 0.038 formal), and with free jitter the highest-peak distribution of white noise is nearly scale-free (median fake maximum about -7).
- Consequence: a model-based threshold is only as good as its noise model. Fakes that carry the correlation raise the threshold as expected (oracle slightly above eFAP, since it knows the process); a data-driven block bootstrap
  gets most of the way (2/20 and 1/20 false alarms versus nominal 0.1%) and is a candidate for the correlated-noise extension (with the Gaia harmonics). Any threshold that rises with correlated noise also loses sensitivity to real
  planets in such data (here thresholds 30-90 vs about 16).
Files: results/red_noise_fakes.json, figures_search/red_noise_fakes_summary.png.

## Settled 2026-11-12 (white-noise fakes for now; merits of a noise model kept for later)
Judah leans toward the simulation-based threshold with white-noise fakes for now. Merits of a noise model, for a possible later implementation:
- Calibration under correlated noise: white fakes give false alarms when the data carry correlated or periodic systematics (red-noise trial, no companion, 20 datasets per kind, FAP 1e-3: sinusoidal noise 12/20, Ornstein-Uhlenbeck 7/20 false alarms
  with white fakes; eFAP 0/20 on both). Fakes that contain the correlation raised the threshold in the same test (block-bootstrap fakes: 2/20 and 1/20 false alarms; oracle fakes with the true process: 0/20 and 0/20).
- It is the natural place to handle the Gaia noise harmonics (mask or model the known periods), and to make the stop rule honest when stochastic correlated noise is present.
- Costs: any threshold that rises with correlated noise also loses sensitivity to real planets in such data (thresholds 30-90 against about 16 for white fakes); a process model must be right or data-driven (block bootstrap, which needs a block length);
  strictly periodic systematics cannot be told from planets by any threshold and need vetting against other data.
- Candidates seen so far (not adopted): moving-block bootstrap of the fitted model's residuals (`noise='block'` in simthreshold.py); a parametric noise process; masking known harmonic periods.

## Results 2026-11-12 (eFAP vs simulated threshold on real HD 114762 + the Gaia-4 injected planet; run_threshold_comparison.py hd114762)
Setup: REAL HD 114762 data (contains HD 114762 b, 83.9 d; 63 transits) plus the same planet as the Gaia-4 main case (1428.25 d, e 0.2, 9.8 M_J, same orientation angles) around a 0.80 M_sun host at 26.2 mas:
astrometric SMA 702 uas, SNR 50.3 (Gaia-4 case: 420 uas, 56.6). Rules: eFAP and simulated threshold (FAP 1e-3; M = 200 main, 100 in the draws); detection matched to truth (83.9 d, 1428.25 d) by half a resolution element in frequency.
Ten random-orientation draws re-solve the mass (8.3-12.3 M_J) to keep SNR 50.3.
- Main case: eFAP finds nothing: step-1 highest peak 117.0 at 83.5 d (HD 114762 b) vs eFAP threshold 203.6. Simulated finds both: step 1 peak 117.0 at 83.5 d vs 21.4; step 2 peak 76.9 at 1386.4 d (injected 1428 d) vs 20.1;
  step 3 peak -13.4 vs 2.7 (stop). No spurious detections.
- Over the 11 cases: both companions found by eFAP 3/11 (draws 4, 7, 9; margins of 1-11% over its step-1 threshold), simulated 11/11; spurious detections 0 for both. Step-1 eFAP thresholds 137-204 against peaks 117-160;
  simulated step-1 thresholds 18-25.
- Same pattern as Gaia-4 (eFAP 1/11, simulated 11/11): the eFAP threshold is lifted by the strong signals present. Files: figures_search/hd114762_real_inj_{efap,sim}_step*.png, hd114762_real_inj_summary.png;
  results/hd114762_real_inj_{main,repeats,step1}.json.

## Update 2026-11-12 (HD 114762 figures regenerated with the known signals)
Every HD 114762 comparison figure now lists the known signals in the title (short names) and draws each as its own labelled line: HD 114762 b (real, fitted to the real data alone): P = 83.9 d, e = 0.34,
a* = 1863 uas, SNR = 186 (square root of the chi-squared gain at formal errors); injected planet: P = 1428 d, e = 0.2, m = 9.8 M_J, a* = 702 uas, i = 65 deg, SNR = 50. HD 114762 b is therefore a much stronger
signal than the injected planet (the periodogram peak understates this because the jitter is refitted: highest step-1 peak 117). The eFAP step figure now has its legend below the axes like the simulated one.
"Main" = the injected planet at the same orientation angles as the Gaia-4 case, run at M = 200 with per-step figures; the 10 "draws" use random orientation angles (mass re-solved for the same SNR), M = 100, no figures.

## Accepted limitation 2026-11-12 (to reconsider later): planets absorbed by bad first fits
Judah: for now live with the possibility that additional planets are "absorbed" by bad initial fits of strong signals. Mechanism (measured on HD 114762 + a 5 M_J, 400 d planet, SNR 11): the step-1 refit of the dominant companion shifts to absorb
part of the weaker unmodelled signal (b: e 0.318, M0 6.17 vs 0.344, 5.92 when fitted alone); the step-2 periodogram holds b at those values, the misfit and the weaker planet partly cancel (chi2 gain 31 instead of 117, dBIC 1.2 at 400 d), and
the planet is missed; the joint refit of all companions only runs after a detection, so it never gets the chance (the search already polishes all companions jointly after each detection, so a "backward refit" would add nothing).
Candidate fix kept for later: at each step joint-refit the top K periodogram peaks (all known companions free) and use the best refit dBIC as the statistic, with the identical procedure in the fakes (estimated 6-20x slower than the current
simulated threshold, about 35 s per step at M = 200). Re-optimising the known companions at every trial period (`refit_known=True`) found the planet's peak (18.8 at 401 d) but is about 235x slower and its threshold is not known.

## Plan 2026-11-12 (injection/recovery map on Gaia-4 epochs and errors; scaffold written, nothing run; blind_search/injection_map.py)
Decisions (Judah): mass log-uniform 0.5-20 M_J; physical semi-major axis log-uniform 0.5-5 AU (kept to show the long-period edge: the search cap of 3609 d is reached at a = 3.7 AU, the baseline at 2.2 AU); e ~ Beta(0.867, 3.03) (capped at 0.95); angles
uniform; cos i uniform; host mass and parallax fixed at Gaia-4's (0.644 M_sun, 13.58 mas); data = injected planet + white noise at the formal errors on the Gaia-4 epochs (no real positions); the detection threshold is SIMULATED FOR EACH INJECTION
(each injection is a different system with its own noise and sampling); no full fits: period = periodogram peak, mass derived from the fitted amplitudes at the peak (Thiele-Innes -> astrometric semi-major axis, Kepler's third law, fixed host);
recovered = detected AND period and mass within 25%. SMA means the physical a in AU. Earlier SNR estimate for the original ranges (1-50 M_J, 0.5-5 AU): percentiles 5/25/50/75/95 = 3.5/8.9/24/62/172; transition near SNR 5-10 at 1-10 M_J.
Implementation: one injection per process, fakes serial inside (about 35 core-s per injection at M = 100, so about 3 s effective on 11 cores: 2000 injections about 1.8 h, 10 injections about 40 s); results appended to results/injection_map.jsonl (resumable);
map = Gaussian-kernel-smoothed recovered fraction in (log a, log m) with all injections drawn (bins are an option). Self-test of the pieces (draws, mass inversion round trip, kernel machinery) passes; no search has been run on injected data.
Expected caveats: the circular template under-estimates the size of eccentric orbits (leading-order about 10-20% at e = 0.5-0.7), so the mass criterion can fail at high e; periods beyond the cap cannot be recovered.

## TO REVISIT (2026-11-12): temporary assumptions of the injection map
- CIRCULAR ORBITS ONLY for now (Judah): injections use e = 0 (`ECC_MODE = 'circular'` in blind_search/injection_map.py). Restore e ~ Beta(0.867, 3.03) (Kipping 2013, capped at 0.95) by setting `ECC_MODE = 'beta'`; the circular template then
  under-estimates eccentric orbits (about 10-20% at e = 0.5-0.7), which will affect the mass criterion.
- Period grid: 20 d to the longest period in the injection range (5087 d: a = 5 AU, 0.5 M_J), so every injected period is searched (Judah). The grid is the usual one, uniform in frequency (spacing 1/(2 pi T)), so it is coarse in period at long periods
  (points near 2170, 2680, 3510, 5090 d).
- Map: kernel-smoothed recovered fraction for now; binning to be considered later.

## Result 2026-11-12 (10-system end-to-end check of the injection map; results/injection_map_check10.jsonl; Gaia-4 epochs, circular orbits, M = 100 fakes per injection)
Recovered 3/10, detected 6/10. Misses: SNR 1.4 (0.5-0.6 M_J at 0.6-0.9 AU) undetectable; two near-threshold non-detections (SNR 8.5 at 167 d, peak 11.9 vs threshold 12.6; SNR 6.0 at 4401 d, 11.9 vs 13.1);
three detections failing the 25% criteria at long period (injected P 1466, 3626, 3589 d; recovered P 2167, 2680, 2680 d; masses off by 30-80%): beyond about the baseline (1805 d) the period is poorly constrained and the peak sits on a coarse grid point.
Time: 81 s wall for 10 injections, about 76 s each per core in parallel (about twice the 35 s estimated from solo timings), so 2000 injections would take about 3.8 h and 1000 about 1.9 h on 11 cores at M = 100.
Search thresholds across the 10 datasets 10.1-18.8. The plot code runs; contours are masked where fewer than 3 effective injections. The baseline period corresponds to a = 2.5 AU (an earlier note said 2.2 AU).

## Decision 2026-11-12 (injection map: every detected signal is optimised as usual; refit window; circular refit)
Judah: "no full fits" meant no Octofitter/MCMC; each signal that passes the threshold is optimised as usual. Implemented as `fit.refit_circular` (blind_search/fit.py): a scan of 400 log-spaced periods plus a 1-D polish, with the linear amplitudes and the jitter refitted at
every period, over a window of a factor 1.5 either way around the periodogram peak. The window must reach the neighbouring grid periods: the grid (565 periods, uniform in frequency) has adjacent period ratios up to 1.45 (3510 -> 5087 d) and above 1.2
only beyond 2680 d, so +/-40% would fall just short and 50% (factor 1.5) is used. The mass and the 25% criteria use the refitted period and amplitudes.
- TO REVISIT: the refit is CIRCULAR (e = 0, Judah) and must be generalised to eccentric orbits when the injections become eccentric (the nonlinear refit `fit.add_and_refit` already handles e, with a +/-10% period window that would need widening).
- Check of the refit (noise-free circular 3600 d signal, SNR 13.8): recovered exactly from starts at 2680 d and 5087 d. With noise the likelihood is nearly flat in period for a 1.26-cycle orbit (lnL 189.91 at 3600 d vs 189.95 at the refit optimum 4469 d),
  so even a perfect optimiser can return a period more than 25% off: this is a data limitation at long period, not a search failure.
- Duration estimates (M = 100 fakes per injection, 11 workers, about 77 s per injection per core): 10 injections about 1.4 min; 100 about 13 min; 1000 about 2.0 h; 2000 about 3.9 h (M = 50 halves these).

## Note 2026-11-12 (fakes per injection for the injection map)
Recommendation: M = 50 for a first full pass (2000 injections about 2 h instead of 3.9 h at M = 100). Basis: the earlier convergence check (20 datasets, step 1): the Gumbel threshold differs from its M = 200 value by a median 1.2 dBIC at M = 50
and 0.6 at M = 100, with a bootstrap sd of 1.1 at M = 200; thresholds are 10-19, and the per-injection error is random, so it averages out in a map built from thousands of injections. Not tested: the effect of M on the map itself, or bias at M = 50.
Docs: the Stack class and all its methods in blind_search/model.py now have full docstrings (no code change).

## Note 2026-11-12 (model.py docstrings)
All docstrings in blind_search/model.py are wrapped to at most 80 columns; the code is unchanged (checked by comparing syntax trees with docstrings removed, against the original code as listed when the file was written, and before/after the rewrap; model.py is not in git yet, so
there is no committed baseline). Code lines over 80 columns (10) were left alone.

## Result 2026-11-13 (injection map, 2000 injections, M = 50 fakes per threshold; results/injection_map_M50.jsonl, figures_search/injection_map_M50.png)
Gaia-4 epochs/errors, circular orbits, mass 0.5-20 M_J, a 0.5-5 AU, period grid to 5087 d, detected peaks refitted (circular, window factor 1.5), recovered = detected and P and mass within 25%. The run took about 95 min (median 33 s per injection per core).
- Detected 1351/2000 (68%), recovered 1092/2000 (55%). Simulated thresholds: median 14.9 (5-95%: 9.6-20.2).
- By injected SNR (detected / recovered): <3: 0% / 0%; 3-5: 2% / 0%; 5-10: 63% / 34%; 10-30: 100% / 82%; 30-100: 100% / 93%; >100: 100% / 100%.
- By period: P < 1805 d (baseline) recovered 61% (detected 66%); P >= 1805 d recovered 40% (detected 72%).
- Of the 1351 detected: 259 fail the 25% criteria (period only 16, mass only 98, both 145); for detected injections with SNR > 30, 93% are recovered and every failure has P > baseline.
- Map: the 50% completeness contour runs at about 3-4 M_J at 0.5 AU, falls to about 1.5-2 M_J near 2-2.5 AU (best sensitivity), and rises again toward 5 AU, where periods exceed the baseline and period/mass recovery degrades (50% at about 5-10 M_J at 5 AU).

## Result 2026-11-13 (injection map: M = 50 versus M = 100 fakes per threshold, the same 2000 injections; results/injection_map_M100.jsonl, results/injection_map_compare.json, figures_search/injection_map_M100.png and injection_map_compare.png)
Identical injections, data noise and periodogram peaks in both runs; the first 50 fakes of the M = 100 run are the 50 fakes of the M = 50 run, so only the simulated threshold differs.
- Per-injection thresholds: M = 50 minus M = 100 has mean -0.21 and standard deviation 1.93 dBIC (median absolute difference 1.28; correlation 0.80); means 14.92 (M = 50) and 15.13 (M = 100).
- Outcomes: detected 1351 (M = 50) vs 1348 (M = 100), 15 injections change detection; recovered 1092 vs 1092 (54.6% overall), 8 injections change outcome.
- Completeness maps (kernel 0.12 dex): rms difference 0.005, maximum 0.029 over the grid; the 50% contours coincide.
- Cost: median 33 s per injection per core at M = 50 against 63 s at M = 100; the full 2000 injections took about 95 min and about 3.3 h on 11 workers.
- Conclusion: M = 50 gives the same map as M = 100 for 2000 injections; the per-injection threshold noise (about 2 dBIC between the two) averages out.

## Result 2026-11-13 (injection map: eFAP versus simulated threshold, the same 2000 injections; results/injection_map_efap.jsonl, results/injection_map_compare_methods.json, figures_search/injection_map_compare_methods.png and injection_map_efap.png)
Run with `injection_map.py run --method efap` (same seed, so same parameters, noise and periodogram peaks; checked: peaks, peak periods and SNRs are identical to the M = 50 run). eFAP threshold from each periodogram at FAP 1e-3; the refit and recovery rule are unchanged, so only the detection decision differs. The run took 225 s (median 1.0 s per injection) against about 95 min for M = 50.
- Detected 813 (eFAP) vs 1351 (simulated); recovered 725 vs 1092 of 2000. Every injection eFAP detects, the simulated threshold also detects (0 eFAP-only); the 538 simulated-only detections have median SNR 10.7.
- 50% detection at SNR 16.3 (eFAP) vs 6.4 (simulated); 50% recovery at SNR 17.5 vs 7.7.
- Thresholds: simulated median 14.9 (sd 3.2, flat in SNR); eFAP median 82 (sd 48) and it rises with SNR, tracking the highest peak (self-masking: a strong signal's leakage lifts the tail fit). eFAP thresholds are well above the simulated ones even at low SNR (about 30-100 against about 15); not investigated.
- Recovered given detected: 89% (eFAP) vs 81% (simulated); the extra simulated detections are marginal ones whose period or mass misses the 25% window more often.
- Completeness difference (simulated minus eFAP): mean +0.18, rms 0.26, maximum +0.74. The 50% contour sits about 1 M_J higher (about 3 M_J against about 1.5 M_J near 2-2.5 AU; about 10 against about 4 M_J at 0.5 AU); the difference is largest along the transition band, zero at high mass and at the faint end.
- Caveat: the data are white noise by construction (no correlated noise), so these maps measure completeness only, not false alarms; eFAP's advantage on correlated noise (earlier test) is not tested here. Lone companion only: eFAP's multi-planet failure is not part of this map.

## Note 2026-11-13 (slide figures for the eFAP and simulated threshold; blind_search/make_slide_figures.py -> blind_search/figures_slides/)
Eight bare panels, recomputed (not cropped) from Gaia-4 real data + the injected planet (run_threshold_comparison main case); eFAP step 2 reproduces the original figure (peak 129.8 at 562.3 d, threshold 102.6) and the simulated step 1 its thresholds (13.7 / 6.4 at M = 200, seed 500).
- Error bars ARE used by the periodogram: every fit is weighted by 1/sqrt(err_i^2 + s^2) (formal error plus a jitter s refitted at each trial period; model.Stack.sigma). The time-series panels show the formal errors; the simulated panels the error the fakes use (fitted 0-companion jitter s = 0.227 mas).
- The eFAP step-2 periodogram is of the data with the step-1 companion (the injected planet, fitted at 1745 d) already removed, whereas the time-series panel (as requested) has only the baseline subtracted, so it shows both planets summed. Flagged to the user; an alternative panel with the step-1 companion also subtracted was not made.
- Simulated panels: panel 1 = the 0-companion model at the real epochs (zero) with error bars sqrt(err^2 + s^2); panel 2 = fake #0 of the step-1 figure (identical to the first grey periodogram, checked); panel 3 its periodogram; panels 4-5 the grey periodograms and histogram without the real-data curve, horizontal lines, texts or arrow.

## Plan 2026-11-13 (cleanup of blind_search/: remove unused code, readable docstrings, README, reliability tests)
Requested by Judah 2026-11-13. `blind_search/` is untracked in git, so deletions cannot be undone from history: step 0 is a backup tarball
(`../gaia_dr4_rv_sensitivity_blind_search_backup_2026-11-13.tar.gz`, outside the repo; Judah may delete it once satisfied).
1. **Delete unused code.** Rule: a module goes if nothing imports it and it produces no figure or result still cited in these notes.
   - Whole files: `threshold.py` (prewhitened eFAP; nobody imports it; the idea was dropped when the simulated threshold replaced eFAP);
     `red_noise_test.py`, `red_noise_fakes_test.py`, `make_red_noise_fakes_figure.py` (the correlated-noise experiments; settled 2026-11-12
     "white-noise fakes for now"; their results stay in `results/` and `figures_search/red_noise_fakes_summary.png`);
     `make_comparison_figures.py` (plots numbers typed in from experiment logs that no longer exist; superseded by the 2000-injection comparison).
   - Functions: `data.synthetic_astro`, `model.bic`, and, in `simthreshold.py`, the block-bootstrap and callable noise options
     (`block_resample`, `noise='block'`, `block=`), which only the red-noise experiments used.
   - Kept although borderline (Judah decides): `make_sim_figures.py`, `make_method_diagram.py`, `make_flowchart.py`, `data.load_rv_project`
     (the RV channel loader, not yet exercised by any test).
   - After each deletion: grep for dangling imports; run `tests/test_blind_search.py --quick` and every kept script's `--test`.
2. **Readability.** Every docstring (module, function, class) wrapped at 80 columns with a script (`blind_search/wrap_docstrings.py`, itself
   removed afterwards) that rewraps plain prose paragraphs only and leaves lists/tables for hand edits; checked by comparing the code's AST with
   docstrings stripped before and after (so no code changes), then the tests. Comments and code lines are not reflowed. New
   `blind_search/README.md` (what the tool does, the algorithm in plain words, module map, how to run, results, limitations) and a pointer from
   the top-level `README.md`; `CLAUDE.md`'s module list updated.
3. **Reliability tests** (new script `blind_search/reliability_tests.py`, with `--test`; Judah asked for "a handful"; these are my choices):
   - T1 rule out noise: 300 datasets of pure white noise at the Gaia-4 epochs and formal errors, each through the step-1 decision with the
     simulated threshold (M = 50). Expected false alarms: 3 of 300 at FAP 1e-2, 0.3 at 1e-3. About 15 min.
   - T2 detect signals: the 2000-injection map already holds detection versus SNR (no new compute: tabulated). Plus a regression check that the
     cleaned code reproduces the stored result of the first 20 injections exactly (guards against cleanup changing behaviour).
   - T3 iterative search, synthetic Gaia-4 epochs, M = 100: (a) no planet: nothing found; (b) Gaia-4 b alone: found, then stops;
     (c) Gaia-4 b + the 200 d companion: both found, then stops; (d) a faint companion (SNR about 3): correctly not detected.
   - T4 the real Gaia-4 pre-release data: Gaia-4 b found, nothing more.
   Each prints PASS/FAIL against stated criteria; the whole set is about 25 min on 11 workers. No other variants (e.g. eFAP on the null
   datasets, eccentric injections) unless Judah asks.
4. Record results here and in PLAN.md; commit nothing (Judah asks).

## Decision 2026-11-13 (Judah: core figures only, M = 50 default, smaller reliability tests)
- Figures: keep only the core set (the step figure `plots.plot_step_sim`, the injection map `injection_map.py plot`, the flowchart `make_flowchart.py`); everything else deleted (list in `blind_search/FIGURES.md`; copies in `../gaia_dr4_sensitivity_blind_search_backup_2026-11-13.tar.gz`, taken before the deletions). `search.run_search` lost the options only the removed scripts used (fixed-floor rule, `truth_info`, `title_prefix`, `resolution_tol`); it draws step figures only for the simulated threshold. Result files in `results/` were kept.
- Default number of fakes per threshold: M = 50 everywhere (`run_search`, `simulated_threshold`, `injection_map.py run --M`), because the 2000-injection comparison showed no sensitivity (completeness rms difference 0.005 against M = 100).
- Reliability tests redesigned before the first full run (Judah asked whether 300 datasets were needed): T1 uses 100 noise-only datasets (each with its own M = 50 threshold) and checks both exceedance counts at FAP 0.1 and 0.01 and the uniformity of the p-values; FAP 1e-3 is not testable at that size and is no longer claimed. T2b reruns 5 injections (was 20), T3 uses 3 draws at M = 50 (was 5 at M = 100).
- The flowchart still described the eFAP threshold and a scratch script (`exp7.py`); its Step B, decision box and notes now describe the simulated threshold (re-rendered `search_flowchart.png`).
- Cleanup of 2026-11-13: also removed `threshold.py`, the red-noise scripts, `make_comparison_figures.py`, `data.synthetic_astro`, `model.Stack.bic` and the block-bootstrap/callable noise options of `simthreshold.py`; all docstrings wrapped at 80 columns (AST-checked: no code change); `blind_search/README.md` added and a quick-start section in the top-level `README.md` (checked by running the example on Gaia-4).

## Result 2026-11-13 (reliability tests; blind_search/reliability_tests.py, results/reliability_tests.json and .log; 492 s on 11 workers; 12 of 12 checks passed; final run with N = 45 model+noise datasets)
- T1, N = 45 model+noise datasets (zero-planet model plus white noise at the Gaia-4 epochs and formal errors), each with M = 50 noise-sampled datasets: 2 above the FAP 0.1 threshold (4.5 expected; 95% interval 0.5-15%), 0 above FAP 0.01 (0.45 expected; weak), KS test of the p-values: statistic 0.17, p = 0.13. No sign of miscalibration, but weak. An earlier run with N = 100 (the 45 were its first 45: same seeds, so not independent) gave 5 above FAP 0.1 (10 expected) and KS p = 0.61, 737 s in all. Both fall below the expected rate: a hint of a slightly conservative threshold, not significant. FAP 1e-3 is not testable at this size.
- T2a, detection versus SNR from the 2000-injection map (M = 50): 0% below SNR 2, 0.4% at 2-4, 9.6% at 4-6, 65% at 6-8, 98% at 8-12, 100% at 12 and above (1020 of 1020). Recovered (period and mass within 25%): 32% at 6-8, 65% at 8-12, 84%, 86% and 95% in the bins 12-20, 20-40 and above 40. About 11% of detections above SNR 12 fail the recovery window; the cause has not been investigated (not part of the request).
- T2b, regression: the first 5 injections reproduce the stored map results exactly with the cleaned code, so the cleanup did not change behaviour.
- T3, synthetic Gaia-4 epochs, 3 noise draws each: no companion, nothing found; Gaia-4 b alone, found then stop; Gaia-4 b + the 200 d companion, both found then stop; a faint (SNR 3) companion, not detected; no spurious detection anywhere.
- T4, real Gaia-4 data: Gaia-4 b found (peak 181.3 against threshold 14.3), step 2 peak -13.6 against 16.9, nothing more.
- Time: T1 took 259 s of the 492 s (45 model+noise datasets x 50 noise-sampled datasets).

## Glossary and decision 2026-11-13 (terminology: model+noise and noise-sampled datasets; T1 uses N = 45)
Judah's terms (also in `CLAUDE.md` "Terminology" and `blind_search/README.md`):
- **Model+noise dataset**: a dataset built as a known model plus noise (here the zero-planet model plus white noise at the formal errors, or a model with injected planets plus noise). It stands in for real data in tests. Replaces my earlier "noise-only datasets" / "stand-in real dataset" (T1) and "injection" datasets.
- **Noise-sampled dataset**: one of the M datasets sampled from the model fitted to a model+noise (or real) dataset plus noise, so any signal in it arises from noise only. Replaces "fake", "fakes", "fake dataset" and "simulated dataset". M = 50 per simulated threshold (default).
- **N** (model+noise datasets in a test) and **M** (noise-sampled datasets per threshold) are different numbers; a test uses N x M noise-sampled datasets in total.
- Decision: T1 of `reliability_tests.py` uses N = 45 model+noise datasets (was 100; "similar to the 50, slightly off for easier distinction"; Judah thought 100 overkill). Consequence: expected false alarms 4.5 at FAP 0.1 and 0.45 at FAP 0.01, so the FAP 0.01 exceedance check is weak and the KS test of the p-values carries T1. Older entries in these notes (and result files such as `results/*.jsonl` keys `M`, `n_fake_ge`) still use "fake"/"simulated dataset"; the code text, figure labels, READMEs and flowchart were converted on 2026-11-13.
