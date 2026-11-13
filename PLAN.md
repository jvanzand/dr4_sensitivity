# Plan: an empirical Gaia DR4 sensitivity map from RV giant planets

Drafted 2026-10-23; revised the same day after Judah's first two reviews, and on 2026-10-27
with the study-sample definition (Step 1b), the data-format decisions (Section 3), and the
pipeline starting points (Steps 0, 2-4). Citations marked
*(verify)* were written from memory and need checking on ADS or arXiv.

## 1. Science question

Forecasts of the Gaia planet yield come from injection-recovery tests on simulated data.
Lammers & Winn (2025, arXiv:2511.04673) predict roughly 7,500 DR4 detections, with peak
sensitivity to 3-13 M_J planets at 2-5 au around 0.4-1.3 M_sun stars within ~500 pc; their
simulation assumes a per-transit error of ~54 µas for G < 14, worsening by 10^(0.2(G-14))
for fainter stars (`GaiaForecasts/Python scripts/DR4_exoplanet_catalog_fits.py`). Real DR4
data will differ from those simulations in ways injections can't capture: calibration of
bright stars, outliers, and choices in data processing.

The DR3 test was small and mixed. The Gaia Collaboration (Arenou et al. 2023 *(verify)*)
found consistent periods and eccentricities for 11 known RV planets, while Winn (2022,
*(verify)*) found agreement for three systems and disagreement for four. A quick search
(2026-10-23) found no published full-sample DR4 test against the known RV giants.

**Question:** Does DR4 detect the RV giant planets it should, and how does the empirical
detection probability depend on astrometric signal, period, host brightness, and distance?
The RV planets act as a "known-answer" sample, so the sensitivity map needs no injections
and includes every real-data effect.

## 2. Framing

1. **RVs give M sin i, so the predicted signal is a lower bound.** The astrometric signature
   is α = (M_p/M_*)(SMA/d), and RVs fix it only up to 1/sin i. Whether Gaia "should" detect a
   planet is therefore a probability, computed by marginalizing over inclination.
2. **Sensitivity is natively a function of (α/σ(G), P/T_baseline).** Distance enters through
   α ∝ 1/d and through G. Fit the detection model in (α/σ, P), then project it onto distance
   for representative planet/host pairs.
3. **Detection means detection by our pipeline**: our own fit to the DR4 epoch astrometry
   prefers a companion and recovers the RV orbit. DPAC's published orbital solutions are a
   secondary comparison (Section 7).
4. **"Gaia-only" means DR4 epoch astrometry alone.** Hipparcos (HGCA) can be a secondary run.
5. **What "the same answer" means.** Astrometry and RVs share P, e, T_p, and ω (astrometry
   alone has the (ω, Ω) → (ω+π, Ω+π) ambiguity, so compare ω mod 180°). A Gaia-only fit also
   gives α and i, hence M_p sin i (with M_* and parallax), which is directly comparable to
   the RV value; equivalently, the Gaia-only fit predicts K.
6. **Literature values never set priors.** No published orbital parameters enter the priors
   or initialization of either fit; literature values appear only afterward, as checks.
   - *RV fits* are an initial "blind" search with weakly informative priors (broad uniform
     bounds on period, K, eccentricity, and so on). No prior or starting value comes from a
     known answer. **The whole system is modeled** (revised 2026-10-27): every companion the
     literature reports, including marginal detections and trends, not only the precisely
     measured ones that put the system in the sample. Companions are added iteratively by a
     stated criterion (Delta-BIC >= 20 on RV maximum-likelihood fits; Step 3), so the RV fit also says how many companions the RVs
     support. With the same data and similar software we should get roughly the same answer
     as the literature, so any disagreement is investigated, and heavy disagreements are
     flagged.
   - *Astrometry fits* are far more blind: wide priors throughout, and inclination
     isotropic. Disagreements with the RV fits may then reflect the different data. Each
     system is refit with Gaia, and the comparison asks both whether the **same number of
     companions** is preferred and whether the orbital parameters agree.
   - Inputs that are not planet results are allowed: stellar mass from an isochrone or
     catalog fit, and the data themselves. The DR4 parallax prior should be wide rather
     than taken from DR3, since DR3 and DR4 share data.
7. **The comparison region.** The sample has only upper bounds (SMA ≤ 5 AU; Step 1), so it
   spans regions of both good and poor Gaia sensitivity. P at 5 AU is ~11 yr for a
   1 M_sun host (~16 yr at 0.5 M_sun), about twice the 5.5-yr baseline. Beyond
   P ≈ T_baseline, Gaia sees mostly curvature, part of which the fitted proper motion
   absorbs. The census reports P/T_baseline for every planet; P ≤ 4 yr marks the
   cleanest comparisons (the limit of reliable period recovery in Lammers & Winn's DR4
   simulations).
8. **Published sources only.** Noise models, thresholds, and comparisons come from
   published or arXiv papers, not talks.

## 3. Data format (from the pre-release template in `data_conversion/`)

Judah's joint RV + DR4 demo fits of Gaia-4 and HD 114762 are the starting point for the
pipeline: `/Users/judahvz/academics/applications/Postdocs_2026/RV_DR4_fit_demo`
(`data_conversion/data_converter.py`, `model_definition_<target>.jl`, `model_run.jl`,
`model_plot.jl`; chains in `results/`). It runs on Octofitter v8.2.4. A copy of the
converter is in this folder's `data_conversion/`; the two copies differ only in the
hand-picked RV subset. Decisions and findings (reviewed 2026-10-27):

- **Time standard: TCB, barycentric** (decided 2026-10-27). Every time stamp (DR4 epochs,
  RV epochs, Tp) is converted to one standard. The raw `obs_time_tcb` is "effective observing
  time as TCB" at Gaia, in ns since JD 2455197.5. The VOTable's separate
  `obs_time_bary_corr` column (up to ±6 min for both demo stars) is not applied by the
  demo converter, nor mentioned in the Octofitter tutorial, so the demo's DR4 epochs are
  not barycentric. RV epochs are BJD (usually TDB, about 19 s from TCB). Both offsets are
  small for these periods but should be removed for consistency. The census currently
  stores Tp as BJD_TDB (`Tp_bjd_tdb`) and will need converting.
- **Collapsing CCD rows to transits.** The converter writes one row per CCD (Gaia-4: 1077
  rows from 109 transits). The demo then reduces them to one row per transit in Julia
  (`median_row` in `model_definition_<target>.jl`), so the likelihood never treats the
  ~9 CCDs of a transit as independent. That settles the overweighting concern; what
  remains is the choice of transit position and error. The demo and the current Octofitter
  tutorial (https://sefffal.github.io/Octofitter.jl/stable/gaia-dr4-prerelease/, written
  for v9) differ:
  - *Tutorial (`reduce_transits`):* keep `used_by_agis_al` CCDs; position = median
    `centroid_pos_al`; error = standard deviation of 256 bootstrap medians (RNG seeded by
    `source_id`), floored at median(`ipd_error_al`)/√n; single-CCD transits use the IPD
    error; scan angle and parallax factor from the first CCD. It notes that correlated
    within-transit errors are not modeled.
  - *Demo:* keep `used_by_agis_al` CCDs; sort by `centroid_pos_error_al` and take the row
    of the CCD with the median *error* (mean of the two middle ones for even n). The
    position is that CCD's centroid, not the median centroid, and the error is a single
    CCD's error, with no √n reduction.
  - *Size of the difference:* for the demo stars the demo's transit errors are 2.1-2.2×
    the tutorial's (median 0.082 vs. 0.038 mas for Gaia-4; 0.124 vs. 0.058 for
    HD 114762), so Gaia gets ~4.5× less weight. Positions differ by a median 0.05-0.07 mas
    (up to 1.9 mas). The demo is effectively the fully-correlated limit.
  - *Plan:* adopt the tutorial's reduction by default; run the demo's as a conservative
    sensitivity test in the pilot.
- **Octofitter version: v9** (decided 2026-10-27; upgraded 2026-11-04 to Octofitter and
  OctofitterRadialVelocity v9.0.0, Pigeons v0.4.11). The demo still runs on v8.2.4 in the
  default environment. Migration guide:
  https://sefffal.github.io/Octofitter.jl/stable/v9-migration/. Changes that touch this
  pipeline:
  - `scan_pos_angle` is in **degrees** in v9 (radians in v8.2.4, whose likelihood calls
    `sincos(scan_pos_angle)`). Done 2026-11-04: `data_converter.py` sets the unit in
    `SCAN_ANGLE_UNIT` (degrees), the two converted files were regenerated, and
    `tests/test_scan_angle.jl` checks that Octofitter reproduces independently computed
    abscissae for a known single-star model with zero χ², and rejects radians or swapped
    sin/cos (χ² > 10³). It passes.
  - `GaiaDR4AstromObs` takes no `gaia_id` and has no `gaia_sol`; the DR3 row comes from
    `gaia_dr3_solution(; gaia_id)`.
  - RV likelihoods: `StarAbsoluteRVObs` → `RadialVelocityObs(tab; target=A, ref=Barycentre)`;
    `offset` and `jitter` are no longer added automatically.
  - Model structure: `Planet` → `Body`, the host star is a body (`bodies=[A, b, ...]`),
    `basis=`, `M_pri`, and system `M` are gone (each body has a `mass` in M_sun),
    `θ_at_epoch_to_tperi` is removed, and `initialize!` overrides nest under `bodies=`.
  - Chains: save with `Octofitter.savehdf5(fname, model, chain)`; load with
    `Octofitter.loadchain(fname; model)`.
- **Reference positions.** `centroid_pos_al` is measured from each transit's reference point
  (`ra0`, `dec0`), which the converter drops. v8.2.4 measures offsets from the DR3 solution,
  queried online by `gaia_id`; the broad `ra_offset_mas`/`dec_offset_mas` priors absorb
  the difference. The pipeline should cache DR3 solutions (as `_gaia_dr3_final/` does) and
  keep DR3 and DR4 source IDs separately. Gaia is expected to publish a DR3-DR4 ID
  correspondence.
- **Outliers.** `used_by_agis_al == False` CCDs are dropped: 23.5% of rows for Gaia-4 and
  36.5% for HD 114762. Worth checking what drives it (whole transits vs. single CCDs,
  bright-star gating). Masked CCDs (no centroid) are dropped too: 73 of 1150 for Gaia-4 and
  21 of 900 for HD 114762, none of them otherwise usable.
- **Not issues for the pipeline:** duplicated rows and the different raw RV formats in the
  demo files. Data are retrieved fresh for every system (Steps 2 and 5). The demo's
  known-answer priors and starting values (e.g., Gaia-4's `a ~ Uniform(0, 2)`) and its
  per-target configuration by editing code are not carried over (Section 2.6, Section 5).

## 4. Analysis steps

### Step 0: pilot on a handful of systems

Pilot systems chosen by Judah (2026-11-04): TOI-1694 (Jump data), HD 111232, HD 191939
(multi-planet), TOI-1736 (subgiant host), and HD 43691 (= HIP 30057; giant host,
1.38 M_sun).

Second batch of five (chosen 2026-11-09 from the stage-2 fit set, public DACE RVs only, so
no Jump export is needed): pi Men (HD 39091; three companions including an 11 M_J outer one,
max astrometric SNR 19.7), HD 89707 (54 M_J brown dwarf, e ~ 0.9, SNR 21), HD 147513
(bright host, one 1.2 M_J planet, SNR 0.98, near threshold), HD 10180 (six low-signal
companions, SNR 0.25), and HD 168443 (7.9 M_J + 18 M_J, SNR 20, mixed HIRES/APF/HARPS/CORALIE
data). Alternates with RV tables built but not searched: 70 Vir, HD 94340, mu Ara.

Choose systems that span the outcomes, not the best cases. Start with the two pre-release
systems that have converted data (Gaia-4 and HD 114762), whose real DR4 epoch astrometry
lets the pilot test the real format instead of only simulations. Judah has already run joint
RV + DR4 fits of both (`RV_DR4_fit_demo`, Section 3); the pilot generalizes that code
into a configurable pipeline on Octofitter v9. Then add simulated DR4
data for a few RV systems that span brightness, multiplicity, and expected signal, chosen
from the census ([sample/CENSUS.md](sample/CENSUS.md)). Include a recoverable planet, one
near threshold, a bright-host planet that misses the threshold, a multi-planet system,
and a low-signal planet (expected non-detection).

Pilot goals: get the full chain running end to end locally, measure per-system run time to
plan cluster jobs, and settle the data-format issues in Section 3 (time standard, transit
reduction, scan-angle unit test).

### Step 1: sample selection and census (done 2026-10-23; see [sample/CENSUS.md](sample/CENSUS.md))

- Every NEA planet with RV data and **M sin i ≤ 80 M_J, SMA ≤ 5 AU** (SMA: semi-major axis); no lower bounds and
  no brightness limit, so the sample covers regions where Gaia has poor sensitivity.
- Brown dwarfs beyond the NEA's 30 M_J limit come from 41 RV survey papers and
  catalogs (CORALIE, HARPS, SOPHIE, Keck/Lick, MARVELS, APOGEE, giant-star surveys,
  exoplanet.eu, and others), audited against their sources and deduplicated against the
  NEA. See [sample/lit_search/README.md](sample/lit_search/README.md).
- M sin i (or mass), P, and SMA come from one NEA solution per planet, **RV-only or
  RV+astrometry** (amended 2026-10-30; previously RV-only). The preferred solution may be
  any published form, including fits that use Hipparcos-Gaia or Gaia DR2/DR3 astrometry,
  because the recovery stage uses only DR4, but every solution must use RVs.
  RV+astrometry fits are eligible, not preferred automatically: the rule below decides. **Feng+2022
  is excluded** (NEA solutions and the literature table), since its fits are considered
  suspect, and so are Gaia-only NEA solutions (Gaia Collaboration 2023).
- **Which NEA solution (decided 2026-10-30; replaces the archive default):** per system,
  the newest paper covering the most planets (uncovered planets take their own newest
  solution); then, per planet, any other solution with less than half the period error
  replaces it (the newest such, repeated). Where SMA and period disagree by > 20%, the
  period wins and SMA comes from Kepler's law (2026-11-04). Details and outliers in
  [sample/CENSUS.md](sample/CENSUS.md).
- Cross-match to Gaia DR3 (and later DR4) source IDs, with a positional fallback and
  flags for suspect matches; record G, parallax, and RUWE.
- Forecast α_min/σ_FOV(G) with the `gaiamock` error curve (which rises for bright stars)
  and flag planets as plausibly recoverable: S/N ≥ 1.5 (Lammers & Winn's DR4 threshold)
  and P ≤ 5.5 yr (the DR4 baseline; settled 2026-10-26).

Result (2026-10-30, new solution rule and NEA snapshot): 2667 companions around 2136
hosts, including 338 brown dwarfs. 189 are plausibly recoverable and 171 more are near
threshold (S/N 0.5-1.5). (Earlier on 2026-10-30, archive defaults: 2627 companions, 192
recoverable.) The "recoverable" flag is descriptive only; no companion is
cut by it. (2026-10-26, RV-only: 2616 companions, 184 recoverable.)

### Step 1b: the study sample (rules decided 2026-10-27; built 2026-11-06)

The goal is to test Gaia against RV orbits that are trusted. The rules:

- **The sample unit is a system.** Companions are of interest, but a system must be fit as a
  whole.
- **Precisely characterized companion:** both SMA (or P) and M sin i (or M) measured at
  **≥ 3σ** (value/uncertainty ≥ 3). The M sin i requirement means the companion's mass comes
  from RVs, so a planet found purely by transits (no RV mass) does not qualify. This applies
  the same way to every discovery method: a transiting planet qualifies only with an RV
  mass at ≥ 3σ.
Two stages of cuts (decided 2026-10-27):

- **Stage 1, the snapshot catalog (orbit precision only).** Every census system with at
  least one precisely characterized companion (within the census cuts, M sin i ≤ 80 M_J and
  SMA ≤ 5 AU), with all of its companions. **S/N plays no part in this stage.** The catalog
  is chosen once and then fixed before DR4 (2 December 2026), together with the census
  snapshot. It may still undergo occasional deliberate changes, recorded with their reasons.
  Later restrictions for specific fits are applied as flags; the catalog's systems and
  companions stay the same.
- **Stage 2, the fit set (S/N only).** The snapshot-catalog systems in which **any**
  companion, precise or not, has forecast S/N (α_min/σ_FOV) ≥ **`SNR_MIN`**. `SNR_MIN` is a
  top-level, user-settable parameter, default **0.1**. It will be tuned as a compromise
  between comprehensiveness and compute time. Example: a system with one precise but
  low-S/N planet and one imprecise but high-S/N planet is in the catalog (thanks to the
  first) and is fit (thanks to the second). The aim is to refit as many systems with a
  known ground truth as possible. So forecast S/N is needed for every companion in a
  catalog system, including marginal ones. For scale: the current census has 752 companions
  around 695 hosts at S/N ≥ 0.1, before any precision cut.
- **No separate host-quality cuts.** Bright hosts, RUWE or non-single-star flags, distance,
  and partial orbits (P > 5.5 yr) are not cut separately. The S/N and precision cuts
  decide. A precisely measured companion with detectable forecast S/N that Gaia misses is
  a useful result (RVs see it; Gaia doesn't).
- **The whole system is modeled.** Systems with at least one precisely characterized
  companion also carry their marginal detections and trends into the fits (Section 2.6),
  including companions outside the census cuts (e.g., beyond 5 AU).
- **RV availability is a later cut.** Every system needs RVs we can refit. Systems without
  them drop out at Step 2 or 3 (insufficient RVs, or a/M not constrained by our own fit);
  they are flagged, not removed from the snapshot catalog.
- **Source IDs:** use Gaia's DR3-DR4 correspondence once it is published.

**Built 2026-11-06** (`sample/catalog.py`; results in `sample/stage1_catalog.csv`,
`stage1_systems.csv`, and `resonance_candidates.csv`):
- *Errors.* (Correction: the cached NEA download does carry P, SMA, and mass errors; the
  census just didn't pass them on.) `census.py` now carries them, adds the NEA's
  eccentricity, K, and stellar-mass errors (new cache files, so the 2026-10-30 snapshot is
  unchanged) and Gaia FLAME mass bounds, and records how each error was obtained. The
  source papers' errors on P, K, M sin i, e, and SMA were added to the 40 literature
  orbit tables (new optional columns `e_K_ms`, `e_msini_mj`, `e_e`, `e_a_au`; CLS I and
  Kiefer+2019 through `raw_lit/<source>/errors.py`), only where a paper prints them.
- *Filling missing errors, in this order* (2026-11-06): the published error; a derivation
  from a related published error; the source paper. SMA error from the period and
  stellar-mass errors (σ_a/a = sqrt((2σ_P/P)² + (σ_M/M)²)/3); period error from the SMA error
  (σ_P/P = 1.5 σ_a/a); M sin i (or a transiting planet's mass) from the K error, with
  the relative errors of K, P (1/3 power), M_* (2/3 power), and e's term, in quadrature,
  each where known, flagged `msini_err_from_K`. Stellar-mass errors: published, else a
  stand-in fractional error (the median among stars in the same mass bin with an error,
  ~5%; `sample/st_mass_err_standin.csv`), labeled in `st_mass_err_source`. If no route
  works (the paper gives neither the quantity's error nor a related one's), the error is
  omitted, the companion passes the cut ("innocent until proven guilty"), and it is
  flagged `precision_unverified` and listed in `sample/missing_errors.csv`: 149
  companions (112 NEA, 37 literature), which is concerning. The orbit test uses P/σ_P (else
  the published SMA error).
- *Result.* Of 2667 census companions, 2435 pass the precision cut (105 only because an
  error is missing). Stage 1: 1989 systems, 2503 companions. Stage 2 at `SNR_MIN` = 0.1:
  686 systems, 888 companions.
- *Resonance* (2026-11-06): resonant companions are poorly described by independent
  Keplerians, and adopting literature orbits would defeat the purpose of refitting. But
  the proxy available (two census companions within 1% of 2:1, 3:2, 4:3, 5:4, 5:3, 5:2,
  or 3:1) is crude (it flags HD 191939), so it is a second-tier check, to be consulted
  after a system gives cause for concern (e.g., a Keplerian fit that will not converge).
  Nothing is dropped now. Flagged companions are marked (`resonance_member`), counted (70
  stage-1 companions in 34 systems; 22 in 11 fit-set systems), and plotted
  (`sample/stage1_resonance.png`). The intended rule, `--resonance-policy exclude`:
  resonant companions stay in the model, so the non-resonant ones fit well, but are not
  scored as recoveries; a system whose companions are all resonant is not fit (15 stage-1
  systems, 6 in the fit set).

### Step 2: RV data assembly

- Retrieve fresh RVs for every system; nothing is reused from the demo's files. Collect
  all available RVs from bulk sources rather than hunting for each paper's data
  (decided 2026-11-04), in this order of priority:
  1. **DACE** (`dace-query`), the primary hub.
  2. **Jump (Judah's team database):** Judah runs an SQL query that cross-matches the
     sample's host names with Jump's and pulls RVs for the matches. The pipeline supplies
     the name list (identifiers and aliases from the census) and ingests the output.
     The queries are in [rv_data/jump_queries.sql](rv_data/jump_queries.sql), written
     against Jump's schema page, to run directly in Jump's SQL Explorer; exports go to
     `rv_data/jump_raw/` (2026-11-04). Besides HIRES, APF, and KPF, Jump holds
     NEID RVs, literature RVs (with bibcodes), and imported datasets (DACE, HARPS-RVBank,
     OHP, VizieR) that record their source, time scale, corrections, and bibcode, so
     both provenance tiers come from Jump directly.
     - *First run (2026-11-04, pilot):* 815 RVs for four of the five stars (HD 111232,
       a far-southern star, isn't in Jump). HIRES and APF BJDs are BJD_TDB (HIRES to
       < 0.01 s, APF to ~2 s). The observations Jump has without RVs are iodine-free
       templates and recon spectra (checked with a one-off query). Jump splits HD 191939 and TOI-1736
       over duplicate star records (e.g., T001339, TIC269701147) that hold only
       templates. Jump's imports for HD 191939 hold only HARPS-N (Orell-Miquel+2023),
       not DACE's SOPHIE RVs, so the DACE pull is still needed.
  - **KPF RVs are left out** (decided 2026-11-04): the instrument isn't yet stable over
    the long term. (The SQL Explorer shows no KPF RVs anyway.)
- **Identity and merging, for every system** (decided 2026-11-04):
  1. *One key per system:* the host's Gaia DR3 source_id from the census.
     `rv_data/build_targets.py` resolves it in SIMBAD and writes `targets.csv` (SIMBAD
     position, proper motion, and aliases) and the VALUES list for the Jump queries.
  2. *Matching:* each source (Jump, DACE, Octofitter's catalogs) is queried with every
     alias; any record matched by any alias belongs to that system, so duplicate records
     of one star in a database are merged. A record matched to two systems stops the
     ingest. Each matched record's position is checked against SIMBAD's (30" plus 30 yr
     of proper motion); misfits are flagged for a look, which matters most for close
     binaries whose components have separate records.
  3. *One table per system:* every RV with its source database, study, source record
     name, and dataset (instrument and reduction, split at upgrades), times in BJD_TCB.
     Each dataset gets its own offset and jitter.
  4. *Duplicates:* the same exposure in two datasets (same instrument, < 60 s apart) is
     kept once: Jump's CPS pipeline first, then DACE, then other imports (the dataset
     Jump marks preferred first), then Octofitter's catalogs. Dropped copies are logged.
  5. *Flags, not cuts:* junk and quality flags, unknown time scales, and programs that
     may need permission are columns.
  4b. *Standard duplicate rules* (decided 2026-11-04; full text in
     [rv_data/DUPLICATE_RULES.md](rv_data/DUPLICATE_RULES.md), tested in
     `tests/test_duplicate_rules.py`): copies only within an instrument family and never
     within one dataset; datasets ranked by source, then reduction before publication,
     then size; a copy is recognized by a shared raw-frame ID, else by time (one to one,
     nearest first, within 2 min for exposure-level times or 15 min when a published
     series is involved, since Butler+2017's HIRES times differ from Jump's by up to
     ~9 min for the same exposures); each pair of datasets is checked for a constant
     offset, and pairs far from it are kept and flagged.
  Scripts, run in order (conda environment `gaia_dr4_sensitivity`, see `environment.yml`):
  `build_targets.py` (identities), Judah's Jump queries, `ingest_jump.py`,
  `pull_dace.py` (public DACE, cached in `rv_data/dace_raw/`), `merge_rvs.py` (one table
  per system in `rv_data/systems/`, plus `rv_coverage.csv`).
  - *Pilot result (2026-11-04):* 1080 RVs in, 54 copies dropped, 1026 kept for the five
    systems. Copies: Butler+2017's HIRES RVs of HD 43691 (18, offset 27-31 m/s, scatter
    4-5 m/s against Jump's), da Silva+2007's SOPHIE and ELODIE RVs (36, against DACE's
    and OHP's reductions). DACE has nothing public for TOI-1694 or TOI-1736; HD 111232
    (not in Jump) has 149 DACE RVs from five instruments. DACE's HARPS BJDs are BJD_TDB
    (checked to < 1 s); its uploaded series (SOPHIE, ELODIE, publications) carry no raw
    file times, so their time scale is taken as BJD_TDB and flagged.
  - *DACE lessons:* DACE resolves some names itself (HD with or without a space;
    TOI-1339 = HD 191939) but not HIP numbers, and Durchmusterung names time out, so
    each system is queried under several names and the results combined. DACE labels
    some archive reductions in its publication fields ("HR DRS-3.1", "TACOS"); these
    count as reductions, not studies. Only public data are visible without a DACE
    account (`~/.dacerc`).
  3. **RV catalogs bundled with Octofitter v9:** HARPS RVBank (Trifonov+2020), ESO HARPS
     DR1, HIRES (Butler+2017), Lick (Fischer+2014), CES, and UVES M dwarfs.
  - Where sources hold the same spectra, the higher-priority source wins for that
    instrument. VizieR tables are deferred until the coverage audit shows the need.
- **Two-tier provenance:** every RV records its database and, where known, the study
  that published it. DACE supplies the study per RV (`pub_reference`, `pub_bibcode`);
  its own reductions carry none.
- **DACE check of the pilot systems (2026-11-04, public mode, `dace-query` 3.0.2):**
  - HD 111232: 149 RVs, from CORAVEL (Baranne+1979), CORALIE (Mayor+2004), MIKE
    (Minniti+2009), and HARPS (DACE reductions, already split at the 2015 upgrade).
  - HD 43691: 111 RVs, from ELODIE, SOPHIE, and HIRES (Butler+2017). The same spectra
    appear more than once: 22 ELODIE RVs from da Silva+2007 and 24 from DACE's reduction
    over the same dates, and da Silva+2007's 14 SOPHIE RVs inside DACE's SOPHIE series.
    Deduplicate by instrument and epoch, keeping one reduction.
  - HD 191939: 5 SOPHIE RVs only; its orbits rest on Keck/APF data, so Jump is needed.
  - TOI-1694 and TOI-1736: nothing under their TOI names. DACE lookups need the alias
    list (identifiers by catalog), as Jump does. Per system, also record instruments, epochs, and whether the data
  match those behind the published solution.
- Split instruments at known upgrades (e.g., HIRES 2004, HARPS 2015, CORALIE 2007 and
  2014), each with its own offset and jitter; convert times to TCB.
- Before fitting, report coverage per system: number of RVs, instruments, and baseline
  relative to each companion's period.

### Step 3: "blind" ground-truth fits (RVs plus available astrometry)

- Build on the demo scripts (`RV_DR4_fit_demo`, Section 3), iterating on them as needed.
  Fit the RVs together with the astrometry Octofitter provides, where available:
  Hipparcos IAD, HGCA, and G23H (decided 2026-11-04). The overlap of DR2/DR3 with DR4 is
  acceptable, since the question is how reliable a DR4-only fit is, judged against the
  best available answer.
- Fit each system in Octofitter with weakly informative priors (broad uniform bounds;
  Section 2.6) and the same orbit parameterization the astrometry fits use; cross-check a
  subset with RadVel.
- Model the whole system: precise companions, marginal detections, and trends, with
  companions added iteratively by a stated criterion (Section 2.6). Per-instrument offsets
  and jitter.
- **Settings decided 2026-11-06** (code in `fits_rv/`, configuration in
  `configs/pilot.toml`):
  - *Astrometry:* G23H (it already includes Hipparcos IAD, the HGCA, and Gaia DR2/DR3
    proper motions; ~14 GB catalog, first use needs approval). `freeze_epochs` is a config
    switch (fast and approximate vs. production).
  - *Priors:* P log-uniform, 0.05 d to 200 yr for the literature companions (below the
    census minimum of 0.07 d, so hot companions are modeled) and a 1 d floor for extra
    ones; q = M_p/M_* log-uniform 1e-5 to 1; e uniform 0 to 0.99; i isotropic (`Sine()`),
    Ω and M0 uniform; stellar mass Gaussian from the census mass and its error; per-dataset
    jitter log-uniform 0.1 to 30 m/s; offsets marginalized analytically
    (`MarginalizedRVObs`). RV-only fits use edge-on orbits, so `mass` is M sin i. Keep an
    eye on jitter: ballooned jitter often masks a poor match. The fit flags a dataset whose
    median jitter exceeds 3x its formal errors or presses the upper bound.
  - *Sampler:* NUTS (`octofit`), 4 independent chains. Pigeons is for posteriors that turn
    out to be multimodal. A dynamic, RadVel-style loop is a later option (QUESTIONS.md):
    `octofit` has fixed iterations and no resume, so stopping early needs a custom loop on
    AdvancedHMC. For now every fit is checked afterward (R-hat <= 1.01; bulk and tail
    ESS >= 400; >= 40 independent draws per chain, the analogue of RadVel's
    autocorrelation factor; divergences <= 1%; tree-depth saturation <= 10%; E-BFMI >= 0.3)
    and rerun with doubled iterations (up to twice) on failure.
  - *Number of companions* (`fits_rv/companion_search.py`): RVs only, maximum likelihood,
    no MCMC. Start with the N literature companions, found blindly (periodogram peak of the
    residuals, then a joint refit); then repeatedly test one more and adopt it if
    Delta-BIC >= 20 (5 parameters per Keplerian, N = RV count, per-dataset offsets and
    jitter in every model). The trend is a polynomial of order 0, 1 (linear), or 2
    (quadratic). At every step it is first simplified if its top term no longer earns
    Delta-BIC 20; then the candidates are add a companion, add a companion while lowering
    the trend order (a companion that replaces a trend), and raise the trend order; the
    best one with Delta-BIC >= 20 is adopted. So a signal that could be either is not
    decided by test order. A blind companion's period is capped at 2x the RV baseline (3x
    until 2026-11-09; at least half an orbit covered), so longer-period signals can only be
    a trend or curvature. The search is blind by default (literature companions are an
    optional toggle, decided 2026-11-09); `fits_rv/recovery.py` scores how many literature
    companions it recovers. A blind
    companion's period must also be at least 2x the mean RV spacing, span/(N-1) (decided
    2026-11-09), since sparser sampling cannot constrain a shorter period, and its
    eccentricity is capped at 0.7 (decided 2026-11-09; literature companions keep 0.99 and
    their prior), because e ~ 0.99 long-period "companions" with periastron in a data gap
    mimic curvature (HD 111232). While the N
    literature companions are being found, trend moves are free (not counted toward N) and
    come first at each step once one companion exists; without them, TOI-1736's curvature
    was fit as a spiky high-eccentricity "second companion" (found in the 2026-11-06
    pilot runs).
  - *Periodogram:* a generalized Lomb-Scargle (linear in the amplitudes and offsets,
    fast) used only to rank trial periods, on a grid uniform in frequency with spacing
    1/(5 x baseline): a peak's width in frequency is about 1/baseline whatever its period,
    so five points per width keep the sampling loss small, and a uniform-in-period grid
    would waste points at long periods and miss peaks at short ones. Grid sizes: 0.4-1.2 s
    per periodogram for a 1 d floor (9,500 to 61,000 trial periods for the pilots),
    6-24 s for the 0.05 d floor of the known companions. The threshold is fixed rather
    than calibrated to each system's noise.
  - *Self-test* (`--selftest`): recovers injected planets, a linear trend, and a quadratic
    trend without a long-period companion; no false additions in 10 noise-only datasets.
  - *Kepler solver:* was 40-58% of the likelihood time. Switched (2026-11-06) to a C
    implementation, `fits_rv/kepler.c` (Halley iteration from a series start, folded to
    M in [0, pi], cos and sin of the true anomaly from E without an arctangent, one fused
    pass over the times). It matches the numpy solver to ~1e-13 and makes the whole
    likelihood 1.2x (1 companion) to 1.8x (6 companions) faster; a Julia version of the
    solver alone ran at the same speed as the C one, and compile flags (-O3, -march=native,
    -ffast-math, sincos) changed nothing, since the time is in libm's sin and cos. It is
    compiled on first use (`kepler_c.py`), with a numpy fallback.
  - *Resonance:* a second-tier check (Step 1b); nothing is dropped now.
- Only after fitting, compare to the published solution (P, K, e, M sin i). Disagreement
  calls for investigation, since the same data and similar software should give roughly
  the same answer. Assign a ground-truth tier:
  - *Gold:* orbit well covered; agrees with the literature.
  - *Silver:* agrees with caveats (partial phase coverage, different data, extra RVs since
    publication that hint at another companion).
  - *Bronze (flagged):* heavy disagreement, or P exceeds the RV baseline. Investigate these
    individually; they don't feed the sensitivity map unweighted.

### Step 4: simulated DR4 astrometry (pre-DR4)

This is its own sub-project. Once simulated astrometry exists, combining it with the RVs is
handled by Octofitter.

- Generate DR4-like data with `gaiamock` (`data_release='dr4'`, as in `GaiaForecasts`), in the
  same format the converter produces, so one loader serves simulated and real data.
- For each real system, inject the RV-fit orbit(s) with many inclination and Ω draws, using
  the real host's position, parallax, proper motion, and G.
- Add realism step by step: a calibration-floor jitter term, outliers, and per-CCD
  correlations if Section 3 shows they matter.
- Fit each realization blind with Octofitter, then compare the preferred number of
  companions and the orbital parameters with the RV fit. Cross-check with `gaiamock`'s own
  fitter.
- Detection criterion (provisional): the companion's mass posterior excludes zero at a set
  significance, and its period agrees with the RV period. Calibrate the significance level
  on fits to injected no-planet data (the simulated false-alarm rate). Formal 0-vs.-N model
  comparison comes later.
- Output, per system: predicted detection probability and predicted posterior widths on
  P, e, and M_p.

**About multimodality.** With wide, blind period priors, astrometry-only posteriors often
have several separated modes: period aliases from the scanning pattern, P vs. 2P
solutions, and an eccentric single planet that resembles a 2:1 pair (arXiv:2512.02007).
This doesn't require multiple fits per system. It requires one fit with a sampler that can
jump between modes (Pigeons parallel tempering), which costs more per system. The pilot
will measure how much.

### Step 5: real DR4 (after 2 December 2026)

- Convert the epoch astrometry for every sample star with the Step 0 loader.
- Repeat the Step 4 fits on real data.
- Compare each detection's P, e, ω, T_p, and M_p sin i to the RV truth. For each
  non-detection, the astrometry gives an upper limit on α, which constrains i.
- Diagnose failures: bright-star calibration, data processing, or a wrong RV "truth."

### Step 6: the empirical sensitivity map

- Model each system's detection outcome as a Bernoulli draw with probability
  f(α/σ(G), P/T), with α marginalized over each system's inclination posterior.
- Use a low-dimensional f (e.g., a logistic threshold in log α/σ with a period roll-off):
  a few hundred systems can support a few parameters, not a fine grid.
- Compare the empirical surface to the simulated one (Step 4) and to Lammers & Winn's.
  The difference measures how much real data degrade the forecasts.
- Project onto distance for representative planet/host pairs.

### Overview: the three lanes and the four comparisons (sketched 2026-11-05)

Judah's flowchart of the project, a candidate paper figure (a rough draft so far; a paper
version would be made in matplotlib). Three parallel lanes:

1. **Our pipeline:** census -> RV fits (+ pre-DR4 astrometry; Step 3) -> synthetic DR4
   with our RV orbits injected at our hosts (Step 4) -> fits to the synthetic DR4, which
   give predicted detection probabilities and posterior widths per system.
2. **Real DR4** (2 December 2026): Gaia DR4 epoch astrometry -> DR4-only fits (Step 5).
3. **Epochalypse** (Yahalomi et al. 2026, arXiv:2611.06915): their simulated DR4-like
   epoch astrometry -> our DR4-only fits to it. This lane needs only the DR4 fitter, so it
   can run before DR4, in parallel with the RV fits.

The four comparisons, each stating both sides:
- **Our RV fits vs. our fits to real DR4 data.** The main result: does DR4 alone recover
  the known RV companions? (Steps 5-6.)
- **Our RV fits vs. published literature fits.** A consistency check of the RV fitting:
  do we reproduce published orbits? Agreement shows consistency with published practice,
  not correctness (the same RVs, mostly); disagreements are cases to investigate.
- **Epochalypse injected orbits vs. our fits to its data.** A consistency check of the DR4
  fitting against exact truth. Limits: per-epoch data (real DR4 is per CCD); no RVs; their
  reported errors deliberately omit an injected noise term, so a fit without jitter should
  come out overconfident. Check their scan-angle units and Ω convention first, and pick a
  subset matched to our hosts and companions.
- **Our fits to synthetic DR4 vs. our fits to real DR4.** A realism check: do simulated
  detection odds match what DR4 actually finds? Per system it is one noisy outcome, but
  across the sample it is a calibration test (in bins, the summed predicted probability
  vs. the number detected; predicted vs. real posterior widths; simulated vs. real
  residual scatter). It tests our simulations and, by extension, injection-based
  forecasts such as Lammers & Winn's and Epochalypse.

Left out of the sketch, but part of the plan: the false-alarm threshold that defines
"detection," set on companion-free data (our own zero-companion injections, or
Epochalypse's no-companion systems); and DPAC's DR4 orbital solutions as a secondary
comparison (Section 7).

**About Epochalypse.** Yahalomi et al. inject random companions (log-uniform in mass and
SMA, no occurrence model) onto all ~5.04 million main-sequence G23H stars within 500 pc
and brighter than G = 16, in zero-, one-, and two-companion versions, and recover them
with a kepmodel periodogram. They don't model or test known planets: the data are
simulated from scratch, so a known host simply gets random companions. So they don't do
our Step 4 (which injects our real RV orbits), and nobody has yet published a test of DR4
against known RV companions. What may be reusable: their per-star noise model (Kiefer et
al. terms from G23H, checked against the DR4 prerelease to ~10%), their public injection
code (a possible replacement or cross-check for `gaiamock`), and their long-period S/N
penalty, S/N_tot = S/N_1 √N / [1 + (P/T_baseline)²], which is relevant to `SNR_MIN`.
Their prediction that companions beyond the baseline are often detected with wrong
periods could be tested with our long-period RV companions (a possible extension; it
would mean relaxing the P ≤ 5.5 yr cut). Several authors (W. Thompson, Lammers, Winn) are the
possible collaborators named in Section 8.

## 5. Workflow: local first, then the cluster

Claude works only on this machine. Every script must run locally on a toy configuration
(one or a few systems, short chains) and produce a sensible result before Judah moves it
to the remote cluster for production runs. In practice:

- One config file per run (sample, priors, sampler settings) with a `--test` mode that runs
  in minutes.
- No hard-coded paths, so the same code runs locally and on the cluster.
- Pinned environments (Julia `Project.toml`/`Manifest.toml`, a Python environment file).
- Outputs in a fixed layout, so cluster results can be copied back and analyzed locally.

## 6. Anticipated roadblocks

- **Bright hosts.** Many RV hosts have G ≈ 4-8, where Gaia's calibration is hardest. In the
  census, host brightness alone keeps 32 high-signal planets below threshold. Whether real
  DR4 data follow the bright-star error curve is itself a result.
- **Transit-level errors and outlier fraction** (Section 3) could change the effective
  precision, and so the sensitivity, by a large factor: the demo's and the tutorial's
  transit reductions already differ by ~2× in error. Settle these in the pilot.
- **Octofitter versions.** The DR4 likelihood changed between v8.2.4 and v9.0.0 (scan-angle
  units, model structure). Upgrade to v9 before fitting in earnest (Section 3), with a
  unit test for the scan angle, and validate on injected data before trusting it.
- **Other fitters.** orbitize! and orvara likely don't ingest DR4 epoch astrometry yet
  (orvara works from the HGCA). `gaiamock`'s fitter is the practical astrometry cross-check;
  RadVel stays the RV cross-check.
- **Run time.** Blind, wide-prior astrometry fits are multimodal (Step 4), so production runs
  belong on the cluster. The local test mode keeps this manageable.
- **Small numbers.** If only tens of systems in the baseline region have α/σ near threshold,
  the empirical map will be coarse. Report it as such.
- **Novelty.** Others are preparing similar work. The plan is to build the scaffolding first
  and sharpen the differentiation later.
- **Minor: unseen stellar companions.** RV samples are generally vetted for binaries, and wide
  binaries shouldn't affect orbits at the separations considered here.

## 7. Extras (when the time comes)

- **Joint DR4 + RV refits.** Once both data sets are in hand, refit every system jointly. For
  most systems the gain is incremental, but the result is a uniform set of true masses and,
  for multi-planet systems, mutual inclinations (connects to the CLS true-mass catalog).
- **Inclination distribution of RV-detected planets.** From the joint fits, measure the
  population-level distribution of i. It should be close to isotropic (uniform in cos i,
  i.e., p(i) ∝ sin i) but skewed slightly toward edge-on, because RV detectability scales
  with K ∝ sin i. Measuring that skew is itself a test of the RV selection function.
- **Comparison to DPAC's orbital solutions.** Whether DPAC published an orbital or
  acceleration solution for each host, compared with our pipeline's detections.
- **Control stars.** RV-quiet stars as an empirical false-alarm test. This needs an RV
  completeness map to say confidently that a star has no giant, which is hard, and
  non-detections are inconsistently published.
- **Gaia + Hipparcos runs** for the long-period edge.

## 8. Who does what

**Claude can do without input:**
- NEA queries, Gaia cross-match, the baseline-region census, and the sample table.
- DACE download scripts and a general DR4/RV loader built from the template.
- Blind RV fits (Octofitter and RadVel), post-hoc literature comparisons, and tier
  assignments, with a list of ambiguous cases for review.
- `gaiamock` simulations, local test-mode astrometry fits, and the detection-probability
  model.
- Pipeline structure, figures (with Judah's matplotlibrc), and running notes.

**Needs Judah:**
- The decisions in [QUESTIONS.md](QUESTIONS.md).
- The Jump SQL query for team-database RVs (the pipeline provides the name list).
- Moving tested code to the cluster and running production jobs.
- Data-quality checks for bulk runs.
- Judgment calls on bronze-tier systems and on Step 5 failures.
- Contacting collaborators (e.g., Will Thompson, Lammers, Winn) once the initial checks
  show the project is promising, and authorship decisions.

## 9. Proposed repo layout (once firm)

```
gaia_dr4_rv_sensitivity/
  data_conversion/  pre-release data and the conversion template (existing)
  sample/           NEA query, Gaia cross-match, baseline census
  rv_data/          per-system RVs and provenance files
  fits_rv/          blind Octofitter and RadVel RV fits
  sim_dr4/          gaiamock realizations
  fits_gaia/        astrometry-only fits (simulated, then real)
  sensitivity/      detection-probability model and maps
  configs/          run configs, each with a local test mode
  paper/            figures and draft
  matplotlibrc
```

**Draft flowchart file (2026-11-06).** `sample/project_flowchart.png`, made by
`sample/plot_flowchart.py` from the text above (reconstruction; replace with the original
sketch if it differs). Comparisons are labeled A (RV vs. real DR4), B (RV vs. literature),
C (Epochalypse truth vs. our fits), D (synthetic vs. real DR4).

**Blind-search tool (decided 2026-11-10; details in QUESTIONS.md, "Proposed/Update 2026-11-10").**
A standalone, Octofitter-free Python search (`blind_search/`): linear periodogram, joint refit,
eFAP threshold and BIC model selection; it writes a Julia model script for vanilla Octofitter v9.
Tests use Octofitter. Final-fit priors: log-uniform P, uniform sqrt(e)cos(omega) and sqrt(e)sin(omega),
isotropic inclination. Data channels are selectable (DR4 epochs, Hipparcos IAD, RVs per instrument).
Required outputs include per-step periodogram plots and the eFAP log-histogram tail fit.

Decision 2026-11-11: no special handling of trial periods near 1 yr (degenerate with parallax) or at P >~ the baseline
(degenerate with proper motion). The search lets those signals be fit; whether they are reliable is judged later.
Algorithm flowchart: `blind_search/search_flowchart.png` (drawn by `blind_search/make_flowchart.py`).

## 2026-11-11: simulated per-step threshold (replacing eFAP) - approved plan and success metrics
Idea (Judah approved): at step k, fit the null model (reference motion + the k-1 companions found + jitter) to the real data,
simulate M fake datasets from it (real epochs/scan angles, Gaussian noise at sqrt(err^2 + s^2)), run the same pipeline on each
(refit known companions, Step A periodogram), record each highest peak, and take the threshold from those maxima. Start with
M = 200 and a Gumbel (extreme-value) tail fit to the maxima to reach FAP = 1e-3; also report the direct 0.99 quantile.
Module: `blind_search/simthreshold.py`. Judah was away for about 2 hours; the plan below was written to run unattended.
Success metrics (Judah: 1, 2, 5, 6 matter most; do not agonise over 3, 4, 7):
 1. Independence: step-2 threshold with an unmodelled companion (SNR 9, 12, 16) within 5 dBIC of the noise-only threshold.
 2. Detection: lone companion (step 1) and second companion (step 2) at SNR 9 and 12: at least 14 of 16 draws.
 3. Accuracy: simulated 95th/99th percentiles for 10 noise-only datasets within 3 dBIC of the 150-draw reference (-2.8, 1.3),
    scatter between datasets under 3.
 4. Convergence: bootstrap uncertainty of the threshold under 2 dBIC at M = 200 (compare M = 50, 100, 200).
 5. No false extra companions: at most 1 of 20 noise-only datasets (Gaia-4 b injected) exceed the simulated threshold at FAP 0.01.
 6. Real data: Gaia-4 and HD 114762: known planet detected at step 1 and nothing at step 2.
 7. Red noise (sinusoids, OU): record false alarms; not pass/fail.
Stages: 0 build + time one threshold; 1 metrics 3, 4, 6 (about 20 min); 2 metrics 1, 2, 5 (about 35 min); 3 metric 7 (about 5 min).
If a stage runs more than 1.5x its estimate, cut sample sizes (e.g. 8 draws) and say so.
On failure (try in this order, one change at a time, record each result):
 - 3 or 5 fail (threshold off or too many false alarms): (a) simulate noise with formal errors only instead of the fitted jitter;
   (b) M = 400; (c) direct 0.99 quantile instead of the Gumbel extrapolation; (d) skip the refit in the simulations (also a speed test).
 - 4 fails: M = 400, report the direct quantile at a looser FAP.
 - 1 fails (threshold rises with the signal): the inflated fitted jitter is the suspect; try the formal-error noise or a jitter
   fixed at the noise-only value; report both.
 - 2 fails: split misses into peak below threshold versus wrong period; report at FAP 0.01 as well as 1e-3.
 - 6 fails: compare the real step-2 peak with the simulated maxima (its percentile); check whether the real residuals are white.
 - If 1, 2, 5, 6 pass: write up, update the flowchart (Step B becomes the simulation), keep eFAP and fixed dBIC >= 20 as options,
   and stop. Do not commit. Do not spawn subagents. Stop and report if a bug is found.

Decision 2026-11-12: the simulated threshold keeps Gaussian white noise for the fake datasets; Gaia noise harmonics are handled later. Comparison results: QUESTIONS.md "Results 2026-11-12".

Note 2026-11-12: the simulated threshold keeps white-noise fakes for now; the merits of a noise model for the fakes (correlated noise, Gaia harmonics; costs in sensitivity) are recorded in QUESTIONS.md "Settled 2026-11-12 (white-noise fakes for now ...)".

## 2026-11-13: blind_search cleanup and reliability tests
Unused code and non-core figures removed, docstrings wrapped at 80 columns, `blind_search/README.md` and a quick-start in the top-level README written, M = 50 made the default, and `blind_search/reliability_tests.py` added (12 of 12 checks pass: noise-only calibration with 45 model+noise datasets, detection versus SNR, multi-companion search, real Gaia-4; terms "model+noise" and "noise-sampled" datasets defined in CLAUDE.md). Details and numbers in QUESTIONS.md ("Decision 2026-11-13", "Result 2026-11-13").
