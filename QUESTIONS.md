# Decisions for Judah

## Open

1. **Criterion for adding companions** in the RV fits (and for comparing the number of
   companions Gaia prefers): to be defined.
2. **VizieR fill-in for RVs:** decide once the DACE + Jump coverage audit shows how many
   systems lack RVs.

## To remember (later refinements)

- **Pathological solutions** (noted, not corrected; listed in `sample/census_output.txt`):
  - 2× period-precision switches whose period disagrees with the newest paper's at > 3σ
    (41 of 161 on 2026-10-30, e.g., ε Eri b, GJ 3512 c, HIP 14810 c, HD 47536 b).
  - Published SMAs replaced by the Kepler value (21). For 11 Stassun+2017 solutions the
    stellar mass, not the SMA, looks wrong; Kepler-34 b and Kepler-35 b are circumbinary.

- **Number of companions in the fits.** Superseded on 2026-10-27: RV fits add companions
  iteratively, and Gaia refits are compared on the number of companions preferred as well as
  the orbital parameters (PLAN.md, Section 2.6). Formal model comparison (Bayes factors) is
  one way to make that criterion concrete.
- **Extra companions in updated RV data.** Where RVs have been added since the published
  fit, the data may be sensitive to further companions. Expected to be rare; flag cases
  where the residuals suggest it.
- **Transit-level errors.** The demo already collapses CCDs to one row per transit, so CCDs
  are never treated as independent. Its transit errors are a single CCD's error, about 2×
  the Octofitter tutorial's bootstrap-median errors. Default to the tutorial's reduction and
  test the demo's as the conservative case in the pilot (PLAN.md, Section 3).
- **Census uncertainties.** The 3σ precision cut needs errors on P or SMA and on M sin i or M.
  The NEA download has none, and the literature tables carry only `e_P_d`. Add them before
  building the sample (PLAN.md, Step 1b).
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
