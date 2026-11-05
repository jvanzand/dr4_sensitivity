# Southern RV surveys and brown-dwarf-desert compilations: literature search notes

Search date: 2026-10-24. Scope: published RV orbits of substellar companions (M sin i up to ~80 MJ,
priority 13-50 MJ; 5-13 MJ welcome) with a <= 5 AU, from southern surveys (CORALIE, HARPS, FEROS,
AAT/UCLES) and BD-desert compilations. Already covered elsewhere and skipped: NASA Exoplanet Archive
(<= 30 MJ), CLS I (Rosenthal+2021), Kiefer+2019 SOPHIE.

## How the CSVs were made

- Every table was read from the arXiv LaTeX source (or VizieR / GitHub where noted) with small
  Python parsers run with `python -I` from the session scratchpad (not saved in the project).
  Downloads sit in `sample/raw_lit/<short_id>/` (`arxiv_<id>.src` = raw e-print; `arxiv_src/` =
  extracted source).
- Columns follow the requested schema. Asymmetric errors on P were averaged into `e_P_d`.
- Each CSV lists every companion the paper reports with M sin i >= ~5 MJ (and <= 80 MJ where the
  paper also lists stars). Out-of-range rows stay in the CSV, flagged at the start of `notes` with
  `OUT OF TARGET RANGE (a>5 AU)`, `OUT OF TARGET RANGE (Msini>80 MJ)`, `LIKELY OUT OF TARGET RANGE
  (P>4000 d, a not given)` or `BELOW 5 MJ`. Filter on these flags.
- "In range" counts below mean 5 <= M sin i <= 80 MJ and a <= 5 AU (or P <= 4000 d when a is not given).
- Where a paper gives a true mass (Hipparcos, Gaia, or HGCA astrometry), it goes in `notes`. The
  `msini_mj` value always comes from the RV fit.
- Tp conversions: all are additive offsets stated in the paper, except two noted below (Barbato
  lambda0 epoch; Frensch rBJD) that were inferred and checked.

## Per-paper notes

### 1. Sahlmann et al. 2011, A&A 525, A95 (arXiv:1009.5991). Search for brown-dwarf companions of stars
- Data: LaTeX Tables 3 and 4 (Keplerian elements; M2 sin i, a_rel), Table 2 (M1), Table 7 (Hipparcos
  true masses), and Table 8 (spectroscopic elements of literature targets). VizieR J/A+A/525/A95 has
  only the RVs (not downloaded).
- CSV: `sahlmann_2011/sahlmann_2011_orbits.csv`, 20 rows: 15 survey stars (CORALIE plus HIP 103019 from
  HARPS) and 5 literature orbits that Sahlmann compiled (HD 131664, HD 174457, HD 190228, HD 168443 c,
  HD 191760; flagged "not a new fit"). All RV-only. T0 is in BJD-2400000, converted.
- In range: 15 (11 survey stars and 4 literature orbits). HD 74014, HD 167665, HD 4747 and HD 211847 have
  a > 5 AU. HD 190228 has M sin i 3.6 MJ.
- Notes: HD 43848 has two RV solutions. I list the adopted low-e solution and give the high-e one in
  `notes`. Hipparcos shows that HD 17289, HD 30501 (90 MJ), HD 43848, HD 53680, HD 164427A and HIP 103019
  are stellar. HD 17289's companion biases the CCF.
- Status: CSV ready.

### 2. Unger et al. 2023, A&A 680, A16 (arXiv:2310.02758). BD desert with precision RVs and Gaia DR3 orbits
- Data: LaTeX Table 3 (RV-only MCMC), Table 4 (joint Gaia DR3 + RV MCMC, "robust" and "challenging"
  sections) and Table 1 (stellar parameters). VizieR has only the RVs.
- CSV: `unger_2023/unger_2023_orbits.csv`, 61 rows: 31 `RV-only` rows (Table 3) and 30 `RV+astrometry`
  rows (Table 4; true mass, inclination and Omega in `notes`, "CHALLENGING CASE" flagged).
- In range (RV-only): 29 companions (HIP 66074 and HR 810 are below 5 MJ). The joint fits show 13
  companions above 80 MJ in true mass (stellar).
- Tp in the paper is relative to BJD 2457389.0 (stated for Table 4). I assumed the same for Table 3 and
  checked it: HD 3277's Tp matches Sahlmann+2011 to 0.001 in phase. Converted.
- The sample is not purely southern. It includes RVs of HD 81040, BD+29 1539, HD 132406 and others from
  SOPHIE, ELODIE and HIRES. Many orbits refit literature data.
- HD 68638 appears as HD 68638A in Tables 1 and 4. HD 111232 (multi-planet) is excluded in the paper
  ("Unger et al. in prep").
- Overlap with existing `literature_companions.csv`: BD+29 1539, HD 114762, HD 48679, HD 82460.
- Status: CSV ready.

### 3. Barbato et al. 2023, A&A 674, A114 (arXiv:2303.16717). CORALIE survey XIX: BDs and stellar companions (RV + Hipparcos-Gaia PMa)
- Not on the assigned list. Found as a CORALIE follow-on to Sahlmann+2011. Its paper table is truncated,
  but the full table is at VizieR J/A+A/674/A114 (downloaded: `barbato_2023/vizier/ReadMe`,
  `table1.dat`, `table2.dat`).
- CSV: `barbato_2023/barbato_2023_orbits.csv`, 29 rows (all companions with RV-only M sin i <= 80 MJ out
  of 221). The rows use the RV-only columns. The orvara RV+PMa true mass, P, e, a and i go in `notes`
  (the PMa solution gives no K or Tp).
- In range: 19. Nine have a > 5 AU, and HD 196885 b (2 MJ planet) is flagged. New relative to
  Sahlmann/Unger: HD 20916, HD 28454, HD 30774, HD 78746, HD 94340, HD 153284, HD 184860A, HD 217580,
  HD 219709. Most have stellar PMa true masses (for example HD 20916 160 MJ, HD 94340 409 MJ).
- Unusual Tp: Barbato gives the mean longitude lambda0, not Tp, and the paper does not state the
  reference epoch. I took BJD 2455500 (DACE convention) and confirmed it empirically. The computed Tp
  matches 16 independent Tp values from Sahlmann+2011 and Unger+2023 to within 0.005 in phase for
  well-sampled orbits (worse only for HD 4747, HD 43848 and HD 89707, whose orbits differ).
  Tp_bjd = 2455500 + (omega - lambda0)/360 * P; lambda0 is kept in `notes`. K is converted from km/s.
- Status: CSV ready (Tp derived as described).

### 4. Rickman et al. 2019, A&A 625, A71 (arXiv:1904.01573). CORALIE XVIII: massive planets and BDs beyond 5 AU
- Data: LaTeX Tables 3-4. P is given in years; I converted it with 365.25 d/yr. Tp is in BJD-2450000.
- CSV: `rickman_2019/rickman_2019_orbits.csv`, 4 rows (HD 13724 b 26.8 MJ, HD 92987 b 16.9 MJ,
  HD 181234 b 8.4 MJ, HD 98649 b 6.8 MJ), all with a = 6.6-12.4 AU.
- In range: 0. All rows are flagged as out of range.
- Status: CSV ready, nothing in range.

### 5. Segransan et al. 2010, A&A 511, A45 (arXiv:0908.1479). CORALIE XVI
- Data: LaTeX Tables 2-3. Only HD 147018 c (6.56 MJ, P=1008 d, a=1.92 AU) is >= 5 MJ. The others
  (HD 147018 b, HD 171238 b, HD 204313 b) have 2.1-4.1 MJ.
- CSV: `segransan_2010/segransan_2010_orbits.csv`, 1 row, RV-only, in range.
- Unusual: Table 2 gives T0 = 55301 +/- 22 in "JD* = JD-2450000", which is inconsistent (probably meant
  5301, i.e. 2455301). I left Tp blank. HD 147018 c is very likely already in the NASA Exoplanet
  Archive.
- Status: CSV ready (Tp blank).

### 6. Kiefer et al. 2021, A&A 645, A7 (arXiv:2009.14164). True masses of RV exoplanets with Gaia DR1 (GASTON)
- Data: LaTeX Table 3 (29 detection-sample companions) and appendix longtable (227 non-detection
  companions). The RV elements come from exoplanets.org, so they are compiled values, not new fits. The
  GASTON true masses (Table 4) are in `notes`.
- CSV: `kiefer_2021/kiefer_2021_orbits.csv`, 34 rows with M sin i >= 5 MJ. In range: 33
  (31 hosts; HD 75784 b has a = 6.5 AU).
- Almost all are <= 30 MJ exoplanet-catalogue objects that the NASA Exoplanet Archive likely covers.
  Overlap with existing CLS rows: HD 114762, HD 156279, HD 168443, HD 178911 B, HD 222582, HD 33636,
  HD 74156. The CSV adds little that is new. Its value is GASTON mass constraints (e.g., GASTON gives 30 Ari B b 147 MJ,
  HD 16760 b 292 MJ, HIP 65891 b 312 MJ, HD 5388 b 87 MJ).
- Status: CSV ready (compiled elements).

### 7. Ma & Ge 2014, MNRAS 439, 2781 (confirmed via Crossref, DOI 10.1093/mnras/stu134, 2026-10-26; arXiv:1303.6442). Statistical properties of BD companions
- Data: deluxetable Table 1 (the arXiv version is the Feb 2013 preprint; I did not check the published
  table against it). It gives Mc or M sin i, P, e, M*, Teff, log g, [Fe/H] and reference numbers. It
  gives no K, omega or Tp.
- CSV: `ma_ge_2014/ma_ge_2014_orbits.csv`, 65 rows. Reference numbers are mapped to citation keys in
  `notes`. Where Ma & Ge list only a true mass (astrometric or transiting), msini is blank and Mc is in
  `notes`. Includes transiting BDs (CoRoT, KELT-1, WASP-30, KOI) and the MARVELS BDs.
- In range: 58 rows (57 hosts). All need full orbits from the original papers.
- Status: compilation only. Use it as a pointer list. It is fully superseded by Stevenson+2023 (below).

### 8. Grieves et al. 2017, MNRAS 467, 4264 (confirmed via Crossref, DOI 10.1093/mnras/stx334, 2026-10-26; arXiv:1702.01784). MARVELS BD candidates
- Northern (SDSS-III MARVELS). I parsed Tables 4-5 (12 companions, 10 with M sin i <= 80 MJ).
- Directory conflict: another search agent working in the same `grieves_2017/` directory (its download
  is `src_1702.01784*`) overwrote `grieves_2017_orbits.csv` after me (10 rows, without the two stellar
  MARVELS-14b/15b). Its values match my parse, so I left its file in place. My downloads in that
  directory are `arxiv_1702.01784.src` and `arxiv_src/`.
- In range: 10 BD candidates (MARVELS-8b to 13b and 16b to 19b). The MARVELS RVs have large rms (up to
  430 m/s). Treat these as candidates.
- The paper points to Wilson+2016 Table 7, a BD candidate compilation (see Leads).

### 9. Jenkins et al. 2009, MNRAS 398, 911 (confirmed via Crossref, DOI 10.1111/j.1365-2966.2009.15097.x, 2026-10-26; arXiv:0905.2985). First CHEPS results: BD in the desert (HD 191760)
- Data: LaTeX Table 5. HD 191760 b (38.17 MJ, P=505.65 d, HARPS CHEPS) and HD 154672 b (5.02 MJ,
  Magellan plus HARPS). T0 is labelled JD-2,450,000 and assumed to be periastron.
- CSV: `jenkins_2009/jenkins_2009_orbits.csv`, 2 rows, RV-only, both in range. HD 191760 M* = 1.28
  (this paper's value). Sahlmann+2011 lists the same orbit.
- Status: CSV ready.

### 10. Jenkins et al. 2015, MNRAS 453, 1439-1457 (from the Barbato+2023 .bib; arXiv:1507.04749). Spectroscopic binaries from the AAPS
- Found through Barbato's references. AAT/UCLES and HARPS. Tables 2(i)-(ii) give P, omega, e, Tp, K and
  a for 16 binaries, with no M sin i column. M sin i appears only in the per-star text.
- CSV: `jenkins_2015/jenkins_2015_orbits.csv`, 2 rows: HD 39213 b (M sin i 0.07 Msun = 73 MJ,
  P=1309 d, a=2.3 AU; in range; authors flag a poor fit) and HD 145825 b (0.06 Msun, a=6.5 AU; out of
  range). M sin i was converted from Msun at 1047.57 MJ/Msun. The other 14 are stellar (M sin i >= 0.1 Msun
  per text) and are not listed.
- Status: CSV ready.

### 11. Parc et al. 2026, accepted in A&A (arXiv:2610.25152, posted 2026-10-18; volume not yet assigned). Companions of hot Jupiters: two decades of CORALIE
- Found while checking Ségransan's recent CORALIE papers. Appendix table (label tab:orbital_params)
  lists the outer companions to WASP hot Jupiters (kima fits).
- CSV: `parc_2026/parc_2026_orbits.csv`, 7 rows with M sin i >= 5 MJ. In range: 5 (WASP-8 c 8.3 MJ,
  WASP-34 c 5.3, WASP-53 c 16.7, WASP-81 c 55.6, WASP-132 c 5.6). WASP-68 c and WASP-99 c have a > 7 AU.
  WASP-22 c (4.9 MJ) is not listed.
- omega is converted from radians. The paper gives the mean anomaly M0 at an unstated epoch, so Tp is
  blank (M0 in `notes`). M* is only in the online table. WASP-34 c drops to 4.6-4.9 MJ in its
  alternative models.
- Just accepted, so values may still change in proofs. Hosts are faint (V ~ 10-12).
- Status: CSV ready (no Tp or M*).

### 12. Stevenson et al. 2023, MNRAS 526, 5155 (confirmed via Crossref, DOI 10.1093/mnras/stad3041, 2026-10-26; arXiv:2310.02695). Combing the brown dwarf desert with Gaia DR3
- The BD-desert compilation with the widest coverage I found: 214 BDs and candidates with
  P < 10^4 d, 13-80 MJ. It is machine-readable on GitHub
  (https://github.com/adam-stevenson/brown-dwarf-desert; also ORDO doi:10.21954/ou.rd.24156393).
  I downloaded `BD_catalogue.csv` (the version used in the paper), `BD_catalogue_20-11-23.csv` (the
  latest) and `README.md` to `stevenson_2023/github/`.
- CSV: `stevenson_2023/stevenson_2023_orbits.csv`, from the 20-11-23 version, 193 rows (excluding the
  19 "Gaia binary_masses" and 3 imaging entries). It gives P, e, M sin i, true mass, M*, method and
  source DOI. It gives no K, omega or Tp. `fit_type` is RV-only for methods rv and transit and
  RV+astrometry otherwise (P and e may come from joint fits).
- In range: 149 rows (145 hosts) before deduplication against the other sources. 76 are rv or transit
  only. It covers Ma & Ge 2014, Wilson+2016, Kiefer+2019 and many transiting BDs.
- Host names: TYC leading zeros removed, CoRoT, LHS and NGC names fixed. The original name is in
  `notes` where changed.
- Status: CSV ready as a master pointer list. Full orbits must come from the source DOIs.

### 13. Frensch et al. 2023, A&A 675, A173 (arXiv:2306.04420). HARPS XLVII: long-period planets, one BD, binaries
- Found in an arXiv search for southern HARPS BD papers. Table 5 (table:brown_dwarfs) gives four BD
  candidates with HARPS RV orbits plus orvara true masses. The planets are all < 2.5 MJ.
- CSV: `frensch_2023/frensch_2023_orbits.csv`, 4 rows. In range: 1 (HD 56380 B, M sin i 33.2 MJ,
  a=4.17 AU, but orvara true mass 375 MJ, so stellar). HD 62364 b (12.7 MJ, true mass 18.8 MJ) is at
  6.15 AU, just outside the cut.
- Tp is given in "rBJD", which the paper does not define. Values such as 55738 only make sense as
  BJD-2400000 (DACE convention), so I converted that way and noted it.
- Status: CSV ready.

### Not processed
- Wilson et al. 2016 (SOPHIE IX, arXiv:1602.02749): another agent owns `wilson_2016/`. I briefly put a
  duplicate e-print there by mistake, then removed it. Its Table 7 compilation is covered by
  Stevenson+2023.
- Feng et al. 2022 (ApJS 262, 21; arXiv:2208.12720): `feng_2022/` exists from another agent.

## Totals (my CSVs)
- 57 unique hosts have full RV-only orbits (P, K, e, omega, and Tp for most) in range, from Sahlmann,
  Unger, Barbato, Jenkins 2009, Parc, Frensch and Ségransan. 45 of these have M sin i of 13-80 MJ.
- Many of these are stellar in true mass (Hipparcos, Gaia DR3 or PMa). For a Gaia DR4 sensitivity
  test they remain valid RV ground truth, but their larger astrometric signals make them easy
  detections.
- The compilations (Stevenson, Ma & Ge, Kiefer 2021) add about 130 more hosts with P, e and M sin i
  only.

## Data-quality concerns
- High-e orbits disagree between papers. HD 89707: K = 4190 (Sahlmann), 5110 (Unger) and 5302 m/s
  (Barbato), with e = 0.90-0.95. HD 43848: Sahlmann's adopted e = 0.70 against Barbato's 0.92 (Sahlmann's
  rejected high-e branch). HD 4747: P = 11593 d against 12197 d. Prefer the most recent fit with the
  longest baseline, but flag these.
- The same orbits appear in several CSVs (e.g., HD 3277, HD 52756, HD 17289 in Sahlmann, Unger,
  Barbato and Ma & Ge). Deduplicate by host and prefer the newest RV-only fit.
- Tp derived or assumed: Barbato (lambda0 epoch inferred, verified), Frensch (rBJD assumed
  BJD-2400000), Jenkins 2009/2015 (JD-2450000, may be HJD or JD rather than BJD), Grieves (JD).
  Parc and Ségransan have Tp blank.
- Kiefer 2021 and Ma & Ge rows are compiled values, not independent fits.
- Host names were normalized for SIMBAD, but not resolved against SIMBAD. Spot-check component
  suffixes (HD 164427A, HD 184860A, HD 33473A, HD 68638) and cluster or survey names (NGC 2423 3,
  NGC 4349 127, KOI-205.01, GPX-TF8A-1859).

## Leads
1. **Feng et al. 2022, ApJS 262, 21 (arXiv:2208.12720), "3D selection of 167 substellar companions".**
   RV+HGCA joint fits for many AAT, HARPS and UCLES targets, with full orbits. These are RV+astrometry
   fits, so flag them, but Feng also lists RV elements. Another agent has started `feng_2022/`; check its
   coverage rather than repeat it.
2. **Xiao et al. 2023, RAA 23, 055022 (arXiv:2303.12409), "Masses of a sample of RV exoplanets with
   astrometric measurements".** orvara masses for RV companions; a cross-check on true masses.
3. **Stevenson+2023 source papers without full orbits here.** The 76 rv/transit-only rows each cite a
   DOI. The most useful southern ones to transcribe: Díaz+2012 (HD 14651, HD 22781, HD 30246,
   HD 92320; `diaz_2012/` exists from another agent), Jones+ EXPRESS (HIP 67537, HIP 97233;
   `express_jones/` exists), Wittenmyer+2009 AAPS (HD 91669), Moutou+2011 (HD 217786), and the CORALIE
   HD 62364/HD 162020 papers.
4. **Barbato+2023 full catalogue.** It has 221 RV+PMa solutions. Companions with RV-only M sin i up to
   ~100 MJ near the cut could be added by relaxing the 80 MJ limit. The RVs are in `rvdata.dat` on VizieR
   (not downloaded), if refits are wanted.
5. **Unger et al. in prep (HD 111232) and Rickman et al. 2022 (arXiv:2209.12957)**: CORALIE follow-ups
   cited by Unger and Barbato. Rickman 2022 objects are > 10 AU, so low priority.
6. **Halbwachs et al. 2000, A&A 355, 581, "Exploring the brown dwarf desert with Hipparcos"**: CORAVEL
   and ELODIE BD candidates with Hipparcos masses (e.g. HIP 21832). Older orbits; no arXiv source.
7. **Sozzetti & Desidera 2010** (cited by Ma & Ge as sozzetti10): Hipparcos upper limits for BD
   candidates (HD 43848, HD 131664). Useful only for true-mass cross-checks.
8. **Parc+2026 online tables**: M* and full sample (147 WASP hosts) are at CDS after publication.
   Re-check WASP-81 c (55.6 MJ) and WASP-53 c (16.7 MJ) then.
9. **AAPS and Pan-Pacific Planet Search.** No dedicated southern BD compilation turned up beyond
   Jenkins+2015. The AAPS "Jupiter analogs" and "cool Jupiters" papers (Wittenmyer+) may list massive
   planets above 5 MJ. Not verified.
10. **Gaia DR3 validation (Holl+2023, Gaia Collaboration Arenou+2023)**: lists the RV-known companions
    with NSS orbits. A `gaia_dr3_val/` directory already exists from another agent.
