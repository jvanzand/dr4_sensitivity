# Aggregate catalogs of RV companions (search of 2026-10-24)

Goal: find catalogs that supply many RV companions at once (priority: M sin i = 13–80 M_J, a ≤ 5 AU),
especially those missing from the NASA Exoplanet Archive (NEA). Clean tables are in
`sample/raw_lit/<short_id>/<short_id>_orbits.csv` (standard 14-column schema). Nothing here has been merged into
the census.

## Conventions used in every table

- **Range filter.** Rows have M sin i ≤ 80 M_J (SB9: ≤ 100 M_J, because the primary masses are rough) and
  P ≤ 4000 d, or a ≤ 5 AU when P is missing. Exceptions are noted per source.
- **Host names.** Every host was resolved through CDS Sesame (SIMBAD). The notes give `SIMBAD main_id=... (otype ...)`.
  When a source name did not resolve, I tried a cleaned variant (spacing, Greek letters, dropping a trailing " A",
  exoplanet.eu alternate names). If that worked, `host` was replaced and the notes say `host renamed from ...`.
  10 exoplanet.eu hosts never resolved (flagged `host NOT resolved`). All other tables resolve 100%.
- **NEA overlap.** I flagged overlap by matching the host's SIMBAD HD/HIP/Gaia DR3 identifiers against NEA
  `pscomppars` (hd_name, hip_name, gaia_dr3_id, hostname; queried 2026-10-24), then requiring an NEA planet on that
  host whose period is within 10%. The notes say `companion YES (<NEA name>, mass)`, `host yes, companion no`, or `no`.
  Rows that also appear in the project's `literature_companions.csv` (CLS I, Kiefer+2019) carry a note.
- **fit_type.** "RV-only" means the orbit comes only from RVs, as far as the source states. "RV+astrometry" means a
  joint fit with Hipparcos/Gaia data, so the orbital elements are not independent of Gaia.
- **Times.** Tp is copied as published. SB9 and EOD give JD/HJD-like epochs, not necessarily BJD (see notes).

## Summary

| short_id | Source | Status | Rows | 13–80 M_J | 30–80 M_J | Not in NEA, 13–80 (30–80) | fit_type | K given |
|---|---|---|---|---|---|---|---|---|
| exoplanet_eu | exoplanet.eu, via OEC nightly mirror | CSV ready (no K/Tp) | 229 | 155 | 86 | 118 (85) | assumed RV-only | no |
| stevenson_2023 | Stevenson+2023 BD-desert catalogue | CSV ready (no K/ω/Tp) | 153 | 105 | 66 | 82 (65) | 74 RV-only / 79 RV+astrometry | no |
| sb9 | SB9 (VizieR B/sb9), SB1 orbits | CSV ready, rough masses | 62 | 51 | 46 | 51 (46) | RV-only | yes |
| xiao_2023 | Xiao+2023, orvara RV+HGCA fits | CSV ready | 76 | 22 | 15 | 18 (15) | RV+astrometry | no |
| feng_2022 | Feng+2022 (J/ApJS/262/21) | CSV ready (list only) | 114 | 20 | 7 | 10 (7) | RV+astrometry | yes |
| gaia_dr3_val | Arenou+2023 table + Holl+2023 refsols | CSV ready (P, e only) | 35 | 10 | 9 | 10 (9) | RV-only (literature) | no |
| exoplanets_org | Exoplanet Orbit Database (Han+2014) | CSV ready, not useful | 34 | 16 | 0 | 1 (0) | RV-only | yes |
| oec | Open Exoplanet Catalogue (main) | not useful | — | — | — | — | — | — |
| dace | DACE public exoplanet tables | not useful | — | — | — | — | — | — |

The rows overlap heavily across sources. Unique SIMBAD hosts with a 13–80 M_J companion not in the NEA, summed over
these seven tables: about 210 (about 166 with M sin i > 30 M_J). This count still includes transiting BDs and
cataclysmic-variable/WD systems from exoplanet.eu, and some objects that the southern/northern agents also cover
(Wilson+2016, Dalal+2021, Kiefer+2019/2021, Grieves+2017, Ma & Ge 2014, Sahlmann+2011). Cross-source overlap of
not-in-NEA BD hosts: Stevenson ∩ exoplanet.eu = 36, Stevenson ∩ SB9 = 12, Stevenson ∩ Xiao = 12,
Stevenson ∩ Gaia-DR3-validation = 10.

## Per-source notes

### exoplanet.eu (exoplanet_eu)
- **What.** The Encyclopaedia of Exoplanetary Systems lists objects up to ~60 M_J (plus uncertainties); the largest
  value in the dump is 78 M_J.
- **Access.** exoplanet.eu and its TAP (voparis-tap-planeto.obspm.fr) still time out from this machine (retried
  2026-10-24, 15 s timeouts). I used a mirror instead. The Open Exoplanet Catalogue's `oec_continuity` job downloads
  `http://exoplanet.eu/catalog/csv/` nightly, converts it to OEC XML, and publishes
  `github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems-exoplaneteu.xml.gz`. The last commit was
  2026-10-24 07:02 UTC, and per-planet `lastupdate` stamps run to 2026-10. Raw file:
  `raw_lit/oec/systems-exoplaneteu.xml.gz`. I read the conversion script
  (`oec_continuity/generate_exoplaneteu_systems.py`) but did not run it.
- **Mirror caveats.** The mirror keeps P, e, ω, a, i, mass, star mass, discovery method and alternate names. It drops
  K, Tp, mass_detection_type and the separate M sin i column. Its `mass` field is exoplanet.eu's true mass when one
  exists, otherwise M sin i. The notes flag rows with a listed inclination ("may be a TRUE mass").
- **Selection.** I kept discovery method RV with mass 13–80 M_J, any RV planet below 13 M_J not matched to the NEA,
  and transit-discovered objects with mass 13–80 M_J not in the NEA. Transiting BD masses come from RV follow-up,
  so I included them, flagged by `discovery method=transit` (44 rows). Then I dropped 11 rows below 13 M_J that the
  SIMBAD-ID match placed in the NEA.
- **Counts.** 229 rows: 155 at 13–80 M_J (111 RV-discovered, 44 transit), and 118 of the BDs are not in the NEA.
  There are 74 RV planets below 13 M_J not in the NEA. Many are retracted/controversial planets or are absent from
  the NEA for a reason, so check them before use.
- **Contamination.** Among ≥13 M_J rows: 9 CV*, 8 WD*, 2 WD? hosts and 28 with P < 1 d (BD donors in CVs and
  WD+BD binaries such as ZTF J1840-1742, SDSS J1557+0916, Gaia 19cwm). These are useless for Gaia astrometry. Filter
  on the `otype` in the notes. A few rows have no period, only a (for example NSC J0148-0357, [BHB2007]-1).
- **Status.** CSV ready, but it is a pointer list. For K/Tp, get the discovery papers or exoplanet.eu once it is
  reachable again.

### Stevenson et al. 2023 brown-dwarf-desert catalogue (stevenson_2023), a new aggregate source
- **What.** arXiv:2310.02695, "Combing the Brown Dwarf Desert with Gaia DR3" (MNRAS, accepted 2023). It catalogues
  214 BD companions with P < 10^4 d, each with method, mass, M sin i, P, e, M_star, Teff, logg, [Fe/H] and the DOI of
  the source paper.
- **Access.** Machine-readable at `github.com/adam-stevenson/brown-dwarf-desert`. I used the latest version,
  `BD_catalogue_20-11-23.csv`, which has the 2023-11-20 additions. The publication version `BD_catalogue.csv` is also
  saved, plus the LaTeX source in `src/`. The paper says the catalogue is also on ORDO (doi:10.21954/ou.rd.24156393).
- **Selection.** I dropped `Gaia binary_masses` (no RV orbit) and `imaging`, and kept the rest when M sin i (or the
  transit true mass) ≤ 80 M_J and P ≤ 4000 d.
- **fit_type.** `rv` and `transit` map to RV-only. `rv+H`, `rv+HG`, `rv+HGCA`, `rv+DR3` and `rv+GASTON` map to
  RV+astrometry, with M sin i from RV and the true mass in the notes.
- **Counts.** 153 rows: 105 at 13–80 M_J, 82 not in the NEA (65 above 30 M_J).
- **Caveats.**
  - No K, ω, Tp or a; each row's DOI points to the orbit paper.
  - The catalogue contains a duplicate (TYC 3130-00160-1 and TYC 3130-160-1).
  - `rv+DR3` entries use Gaia DR3 orbits, so they are not independent of Gaia.
  - Frozen since 2023-11.
- **Status.** CSV ready. This is the best single BD list found.

### SB9 (sb9)
- **What.** The 9th Catalogue of Spectroscopic Binary Orbits (Pourbaix+2004). The VizieR copy B/sb9 has 4077 systems
  and 5099 orbits, last refreshed 2024-04-22 according to its ReadMe history. The SB9 home page sb9.astro.ulb.ac.be
  timed out.
- **Raw files** (in `raw_lit/sb9/`): main/orbits/alias via VizieR TAP (`tap_*.csv`) plus FTP copies;
  `EEM_dwarf_UBVIJHK_colors_Teff.txt` (Mamajek dwarf table, v2022.04.16); Gaia DR3 FLAME masses and Gaia DR3
  parallaxes for the candidates (`gaia_dr3_*_for_candidates.csv`, from VizieR I/355).
- **Method for M sin i** (rough):
  1. I used single-lined orbits only (no K2) with P ≤ 4000 d.
  2. Mass function: f(m) = 1.0361e-7 (1−e²)^1.5 K1³ P [M_sun], with K1 in km/s and P in days.
  3. I solved (M2 sin i)³/(M1 + M2 sin i)² = f(m) for M2 sin i. This is the usual minimum-mass approximation, which
     replaces sin³i in the denominator term.
  4. Prefilter: M2 sin i < 200 M_J for M1 = 1 M_sun (515 systems).
- **Primary mass M1**, in priority order:
  - (a) Gaia DR3 FLAME mass (`Mass-Flame`), using the Gaia DR3 ID from SB9's alias list. 35 of 62 rows; FLAME
    assumes a single star.
  - (b) Otherwise, the Mamajek dwarf mass for SB9's Sp1. If Sp1 has no subtype I assumed 5. Class IV used the dwarf
    mass, which likely underestimates M1. Rows were rejected in three cases:
    - Gaia DR3 GSP-Phot logg < 3.8.
    - Gaia absolute G (or V) more than 1.5 mag brighter than the dwarf value.
    - Sp1 has no luminosity class and no parallax check was possible.
  - (c) Giants or supergiants without FLAME, and unparsable or blank spectral types, were excluded. That is 242
    orbits, listed with reasons and the M1 = 1 M_sun M sin i in `sb9_excluded_no_primary_mass.csv`; 25 of them have
    M sin i < 60 M_J at M1 = 1 M_sun.
- **Output.** a is from Kepler's third law with M1 + M2 sin i. When several SB9 orbits exist for a system, the table
  keeps the most recent reference (then the highest grade) and lists the others in the notes.
- **Counts.** 62 rows: 51 with M sin i ≤ 80 M_J, 46 at 30–80 M_J, plus 11 at 80–100 M_J. Only 1 host is in the NEA
  at all, and no companion matches. Frequent orbit references (by bibcode): 2011A&A...525A..95S (Sahlmann+2011, also
  covered by the southern agent), 2018A&A...619A..81H, 2002AJ....124.1144L, 2002ApJS..141..503N, 2009AJ....137.3743G
  (NGC 188 binaries), and several papers in The Observatory ("Obs" bibcodes).
- **Caveats.**
  - The masses are rough: ±30% or worse from M1 alone.
  - SB9 orbits are often old, sometimes with assumed (fixed) e or P (flagged in notes), and T0 is in JD as published.
  - HIP 64426 = HD 114762 is included (13 M_J minimum mass; known to be stellar from Gaia).
  - The VizieR files are pipe-delimited, not in the fixed-width layout the ReadMe describes, so I used the TAP CSV.

### Xiao et al. 2023 (xiao_2023), a new aggregate source
- **What.** arXiv:2303.12409, accepted by RAA. orvara fits of literature RVs plus HGCA for 115 RV companions.
- **Access.** Table "Posteriors of RV Companions" parsed from the arXiv LaTeX source (two parts). It gives Mp, i, a, e,
  P, Ω, ω, Tp and M sin i, but no K.
- **Counts.** 76 rows in range, 22 at 13–80 M_J (15 above 30 M_J), 18 not in the NEA. Several of the paper's "BD
  candidates" have true masses that are stellar; the true mass is in the notes.
- **fit_type.** RV+astrometry. ω is in orvara's convention, and i has a prograde/retrograde ambiguity.

### Feng et al. 2022 (feng_2022)
- **What.** VizieR J/ApJS/262/21 `table3.dat`: 230 companions, RV+Hipparcos-Gaia fits (the paper's 167 cold giants
  and 68 others).
- **Selection.** True mass mc ≤ 80 M_J (mc sin i for P < 1000 d) and P ≤ 4000 d or a ≤ 5 AU: 114 rows. Another 9
  companions with mc > 80 M_J and P ≤ 4000 d were excluded as stellar.
- **Mass column.** For fits with an inclination I computed msini_mj = mc·sin i; the true mass is in the notes.
- **Counts.** 20 at 13–80 M_J (7 above 30 M_J); 10 of those are not in the NEA (7 above 30 M_J): GJ 864 b,
  HD 122562 b, HD 127506 c, HD 203771 b, HD 211847 b, HD 39213 b, HD 48679 b, HIP 22203 b (HD 30246),
  HIP 26196 b (HD 244957), HIP 48953 B (BD-08 2813). M sin i is mc·sin i, so some companions with BD true masses
  (for example HD 105811 b, mc = 55.5 M_J at low i) fall below 13 M_J here.
- **Use.** As a system list only; every orbit is an RV+astrometry joint fit.

### Gaia DR3 validation papers (gaia_dr3_val)
- **Arenou+2023** (Gaia Collaboration, A&A 674, A34). Not on VizieR (J/A+A/674/A34 gives 404), so I parsed the table
  "Known substellar companions..." from the arXiv:2206.05595 LaTeX source.
  - 19 objects: 9 planets and 10 BDs.
  - Literature values: M sin i, P_lit and e_lit. Gaia values: true mass, i, P_Gaia, e_Gaia. No K, ω or Tp.
  - The 10 BDs come from Dalal+2021, Wilson+2016, Kiefer+2019 and Sahlmann+2011, all of which the survey agents
    cover.
- **Holl+2023** (A&A 674, A10; VizieR J/A+A/674/A10).
  - `refsols.dat` lists 204 Gaia DR3 orbits with reference data: reference P/e, Gaia P/e, a pseudo-mass index
    assuming 1 M_sun, and the reference bibcode.
  - I kept pseudo-mass index 0 (< 20 M_J) or 1 (20–120 M_J) with non-Gaia references: 16 more rows (P and e only).
  - Several references are "Barbato et al. (submitted)", "Unger et al. (submitted)", or "DACE re-analysis of public
    data", which point to newer CORALIE/HARPS orbits.
  - The 12 RV time series published with the paper are saved in `holl2023/rv/`. They are new CORALIE data for
    BD-18 113, CD-46 10046, HD 134237, HIP 9095 and others.
- **Status.** CSV ready, 35 rows. Mainly useful as a list of RV companions with Gaia DR3 orbits already known.

### Exoplanet Orbit Database / exoplanets.org (exoplanets_org)
- **Access.** The site still serves `http://exoplanets.org/csv-files/exoplanets.csv` (HTTP only; HTTPS fails). Saved
  in `raw_lit/exoplanets_org/`.
- **How current.** Frozen: 5747 rows, latest DATE 2018.
- **Content.** Only 16 objects have M sin i > 13 M_J, and the maximum is 22.6 M_J. EOD effectively stops around
  24 M_J.
- **Selection.** RV-discovered planets not in the NEA, plus all objects ≥ 13 M_J: 34 rows. 15 of these are in the
  NEA, including 15 of the 16 BD-mass ones. The 19 not in the NEA include retracted or disputed planets (for example
  HD 114762 b, alpha Cen B b, HD 131664 b at 18 M_J, M 67 SAND 364 b, 30 Ari B b), so don't add them without
  checking.
- **Status.** Not useful for BDs; K/ω/T0 are available if needed.

### Open Exoplanet Catalogue (oec)
- **Contents.** The main catalogue `systems.xml.gz` has 5370 planets. Since 2020 it is maintained mostly by importing
  the NEA, and its latest per-planet update is 2023. Only 79 objects exceed 13 M_J and only 4 RV objects exceed
  30 M_J.
- **Status.** Not useful, except that its sister file `systems-exoplaneteu.xml.gz` is the exoplanet.eu mirror used
  above. `systems-exoplanetarchive.xml.gz` (an NEA mirror) is also saved.

### DACE (dace)
- **Access.** Public mode works without login over plain HTTP. The endpoint is
  `https://exo-webapp.obsuksprd2.unige.ch/search/{planets|atreides}?limit=10000&filters={}&sort={}`; I found it by
  reading the `dace-query` 3.0.2 wheel from PyPI. I did not install or run it. I re-implemented its JSON flattening.
- **Contents.** The "planets" catalogue (PlanetS) has 958 rows, mostly transiting planets (930 transit, 27 RV), with a
  maximum mass of 28 M_J. "atreides" has 2 rows. Saved in `raw_lit/dace/api/`.
- **Status.** Not useful for this purpose. DACE's public RV time series (spectroscopy module) remain a lead for
  re-fitting orbits.

## Leads
1. **exoplanet.eu directly**: K, Tp, `mass_sini` vs `mass`, and `mass_detection_type` are missing from the mirror.
   Retry `https://exoplanet.eu/catalog/csv/` or the TAP `exoplanet.epn_core` table from another network. The OEC
   mirror proves the CSV endpoint is alive from elsewhere.
2. **Stevenson BD catalogue follow-up**: for each of the 82 not-in-NEA BDs, the `source` DOI gives the orbit paper.
   - Most common DOIs: 10.1093/mnras/stu134 (24), 10.1051/0004-6361/201527581 (11), 10.3847/1538-3881/ab9b84 (11),
     10.1051/0004-6361/202243782 (9; Gaia DR3/Arenou+2023), 10.1093/mnras/stx334 (8),
     10.1088/1674-4527/accb7e (6; Xiao+2023), 10.1051/0004-6361/201935113 (4), 10.1051/0004-6361/202141145 (3).
   - I have not checked these DOIs, but they appear to be Ma & Ge 2014, Wilson+2016, a 2020 AJ transiting-BD paper,
     Grieves+2017, Kiefer+2019 and Dalal+2021. Most are covered by the southern/northern agents.
   - The methods for these 82 are rv 25, transit 24, rv+DR3 19, rv+HGCA 11, rv+H 2 and rv+GASTON 1.
3. **Barbato et al. 2023, CORALIE XIX** (arXiv:2303.16717), "Brown dwarfs and stellar companions unveiled by radial
   velocity". It corresponds to the "Barbato et al. (submitted)" references in Holl+2023 and belongs with the
   southern (CORALIE) agent.
4. **Recent compilations to check:**
   - Ren et al. 2026 (arXiv:2610.16662, AJ accepted): a sample of 159 BDs with absolute masses.
   - Cui et al. 2026 (arXiv:2603.01808, PNAS): 55 RV+astrometry companions around FGK stars.
   - Barbato et al. 2026, GAPS LXXV (arXiv:2606.29001): HARPS-N validation of Gaia substellar candidates.
   - Jaros et al. 2026, PTPS/JOTA (arXiv:2608.30546): 16 SB1s.
   I verified that each exists via the arXiv API but did not download their tables.
5. **SB9 refinement**: the 242 excluded SB9 orbits (`sb9_excluded_no_primary_mass.csv`) need primary masses, for
   example from TIC v8 (VizieR IV/39) or StarHorse. 25 have M sin i < 60 M_J at M1 = 1 M_sun. The newer SB9 at
   sb9.astro.ulb.ac.be was unreachable.
6. **Gaia DR3 substellar candidates**: Arenou+2023's 1843 BD and 72 planet candidates from the astrometric mass
   catalogue, and its SB1 substellar candidates, are Gaia-only and not RV ground truth. Skip them except as a
   cross-check.

## Data-quality concerns
- Heavy cross-source duplication. De-duplicate on `SIMBAD main_id` (in notes) plus period.
- Feng, Xiao, and the `rv+*` rows in Stevenson are RV+astrometry joint fits that include Hipparcos/Gaia DR2/EDR3 or
  DR3 information. Using them as "ground truth" for Gaia DR4 is partly circular.
- SB9 masses are rough (M1 from FLAME or spectral type), and its orbits are older and heterogeneous.
- exoplanet.eu's mass is ambiguous between M sin i and true mass, and its list includes CV/WD systems and objects
  without periods.
- NEA overlap uses ID matching plus a 10% period tolerance. Misses are possible when the NEA period differs (aliases)
  or when the host has no HD/HIP/Gaia DR3 ID in the NEA.

Scripts (scratch, not in the repo):
`/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/agg/scripts/`
(`*_build.py` per source, `finalize.py` for the SIMBAD/NEA pass, `resolve.py` for Sesame).
