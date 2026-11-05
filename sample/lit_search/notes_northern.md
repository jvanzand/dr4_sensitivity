# Northern / US RV surveys: substellar-companion orbit tables

Search date: 2026-10-24. Scope: M sin i up to ~80 MJ (priority 13-50 MJ), a <= 5 AU, RV orbits.
Each source has raw downloads (arXiv LaTeX source `src_<arXiv id>/`, plus VizieR files where used)
and a clean table `sample/raw_lit/<short_id>/<short_id>_orbits.csv` in the standard columns.

Conventions used in every CSV:
- `Tp_bjd` holds the paper's periastron epoch converted to a full Julian date (offsets such as
  BJD-2450000 or MJD were added back). Many papers give JD or HJD, not BJD; the `notes` column says
  which. When a paper gives only a conjunction time (Tc) or a mean anomaly, see the row's notes.
- Arguments of periastron given in radians were converted to degrees.
- `e_P_d` is the symmetric error, or the mean of the +/- errors when the paper gives asymmetric ones.
- Objects slightly outside the target range (a = 5-6.5 AU, or M sin i = 80-85 MJ) are kept and flagged
  in `notes`; objects clearly outside are left out and listed below.
- The CSVs were built by short parsing scripts (kept in the session scratchpad, not in the project)
  that read the LaTeX tables directly; single-object tables were transcribed by hand from the LaTeX.
- All fits listed are RV-only (`fit_type`), unless noted. Where a paper also ran a separate Hipparcos or
  Gaia analysis that did not change the RV parameters, the RV row is still "RV-only" and the
  astrometric mass is quoted in `notes`.

Overlap checks done: the hosts were matched by name against CLS I (`cls_2021/planets.dat`), Kiefer+2019
(`kiefer_2019/tableb*.dat`) and `sample/raw/nea_pscomppars_rv.csv` (NEA). Results are in each section.
Gaia-ID matching was not done.

---

## SOPHIE

### Wilson et al. 2016, A&A 588, A144 (arXiv:1602.02749), "SOPHIE IX: Populating the brown dwarf desert"
- Data: orbit table = Table 3 (deluxetable in `rv.tex`), host masses = Table 2. VizieR has only RVs.
- CSV: `wilson_2016/wilson_2016_orbits.csv`, 12 rows, CSV ready. All RV-only (SOPHIE, single Keplerian).
  - 7 BD candidates (M sin i 13.2-70 MJ): HD 283668, HD 39392, HD 77065, BD+26 1888, HD 122562,
    HD 132032, HD 134113.
  - 5 objects the paper calls "low-mass binaries" with M sin i 44-83 MJ, kept because their RV orbits
    are good ground truth: HD 6664 (80), BD+73 275 (44), HD 51813 (46), HD 56709 (72), HD 134169 (83).
    Hipparcos IAD orbits give stellar true masses for HD 6664 (0.35 Msun), BD+73 275 (0.20), HD 51813
    (0.27) and tentatively HD 56709 (0.54); these are noted per row.
  - Excluded: HD 50761 B (M sin i 145 MJ).
- Oddities: HD 122562 b has a = 4.1 AU and P = 7.6 yr. HD 77065 and HD 134113 fits include Latham+2002 RVs.
  HD 134113 host mass is from Sousa+2011 (the paper's alternate value is 0.89 Msun).
  NEA lists "HD 122562 c", not b, so check how the NEA entry relates to this object.
- Bonus: Table A.1 (`tab:appendix_table`) is a compilation of ~60 BD candidates around solar-type stars
  (m sin i, P, e, M*, Teff, references) but has no K, omega or Tp. It is useful as a checklist, not
  as an orbit source (see Leads).

### Dalal et al. 2021, A&A 651, A11 (arXiv:2105.09741), "SOPHIE XVII"
- Data: BD orbits = Table 2 (`tab:BDorbit`), host masses = Table A.2. VizieR has only RVs.
- CSV: `dalal_2021/dalal_2021_orbits.csv`, 4 rows (BD-00 4475, HD 184601, HD 205521, HD 5433;
  M sin i 25-60 MJ), CSV ready, RV-only (Yorbit).
- HD 205521 b: Hipparcos (Mc = 0.48 Msun) and Gaia DR1 excess-noise (0.20 Msun) analyses both say it is
  stellar (noted in the row).
- Not included: Table A.4 (`tab:StellarC`) has 16 SOPHIE orbits with M sin i >= 75 MJ. Four lie
  at 90-102 MJ (HD 166356 90.7, HD 187057 92.0, HD 8291 99.5, HD 153915 101.7). Add them if the
  project wants low-mass stellar companions as well. Note that the Tp column header there reads
  "BJD-2500000" but the values (e.g. 57984) are clearly BJD-2400000.
- No overlap with Kiefer+2019 or CLS.

### Díaz et al. 2012, A&A 538, A113 (arXiv:1111.1168), "SOPHIE IV: Massive companions in the planet-BD boundary"
- Data: orbit table = Table 5 (`TableParameters_v3.tex`), host masses = Table 3.
- CSV: `diaz_2012/diaz_2012_orbits.csv`, 7 rows, CSV ready, RV-only.
  HD 22781 b (13.65), HD 156279 b (9.71), HD 14651 b (47.0), HD 30246 b (55.1), HD 92320 b (59.4),
  HD 137510 b (27.3), HD 16702 b (48.7).
- HD 16702: Hipparcos astrometric orbit gives Mc = 0.35-0.45 Msun (stellar). HD 137510: Hipparcos
  limits Mc to 0.019-0.057 Msun (a true BD); its fit combines SOPHIE with McDonald, TLS and Lick RVs.
  HD 30246 has a poorly constrained orbit (K +1600/-560 m/s).
- Excluded: HD 167215 (M sin i 92 MJ, Hipparcos 0.20-0.29 Msun). It is also in Kiefer+2019 Table B1/B3.
- Overlaps: HD 22781 b and HD 156279 b are in NEA. HD 156279 is also in CLS I (status K).

### Bouchy et al. 2016, A&A (arXiv:1511.08397), "SOPHIE VIII: Follow-up of ELODIE candidates: long-period BD companions"
- A&A 585, A46 (2016); confirmed via Crossref (DOI 10.1051/0004-6361/201526347) on 2026-10-26.
- Data: MCMC orbit table = Table 4 (`\label{fit}`; P is given in years), host masses = Table 2.
- All 5 companions lie beyond 5 AU (a = 5.95-22 AU). CSV `bouchy_2016/bouchy_2016_orbits.csv` keeps
  the three borderline ones, flagged "OUTSIDE a<=5 AU": HD 14348 (48.9 MJ, 5.95 AU, P 13 yr),
  HD 72946 (60.4 MJ, 6.37 AU, P 16 yr), HD 209262 (32.3 MJ, 6.16 AU, P 15 yr).
  Excluded: HD 10844 (83.4 MJ, 10.2 AU) and HD 18757 (35.2 MJ, 22 AU).
- RV-only (ELODIE + SOPHIE). No Hipparcos detections.
- Status: CSV ready, but 0 companions strictly within the range.

### Other SOPHIE papers checked
- Bouchy et al. 2009 (arXiv:0907.3559, SOPHIE I, HD 16760 b, 13.1 MJ): in NEA and listed in the
  Wilson+2016 compilation. Not extracted.
- Demangeon+2021 (SOPHIE XVIII, arXiv:2107.05089, cold Jupiters) and Borgniet+2019 (arXiv:1809.09914,
  AF stars): not downloaded. By their titles they are mostly about planets, so they are low priority.
- Kiefer+2019 (already covered) also has orbits for HD 5470 and HD 211681 (both trend-only in Patel+2007)
  and for HD 167215.

## Keck / Lick (California Planet Search)

### Patel et al. 2007, ApJ 665, 744 (arXiv:0704.3418), "Fourteen new companions from the Keck & Lick RV survey"
- Data: `tab2.tex` (Table 2) has the only full orbit: HD 167665 b (M sin i 50.3 MJ, P 4385 d,
  a 5.47 AU). Host mass is from Table 1.
- CSV: `patel_2007/patel_2007_orbits.csv`, 1 row (borderline a), RV-only. Tp is in JD.
- The other 13 stars have trend-only constraints (Table 4 gives M sin i/P/a ranges and no orbits).
  Status: not usable as orbits. HD 5470 and HD 211681 now have orbits in Kiefer+2019.
- Overlap: HD 167665 is in CLS I (status S, 48.4 MJ, P 4433 d), so CLS already covers it.

## SDSS: MARVELS and APOGEE

### Grieves et al. 2017, MNRAS (arXiv:1702.01784), "Exploring the BD desert: new substellar companions from SDSS-III MARVELS"
- MNRAS 467, 4264 (2017); confirmed via Crossref (DOI 10.1093/mnras/stx334) on 2026-10-26.
- Data: companion table = Table 2 (two stacked tabulars: msini/P/e/K, then a/Tp/omega), host masses =
  Table 1.
- CSV: `grieves_2017/grieves_2017_orbits.csv`, 10 rows (MARVELS-8 to 13 and 16 to 19; M sin i 13.1-76 MJ), RV-only
  (RVLIN, bootstrap errors). Tp is in JD.
- Excluded: MARVELS-14b (100.2 MJ) and MARVELS-15b (97.6 MJ).
- Caveat from the paper: these are candidates. Some fits have weak phase coverage or aliases. MARVELS RV
  precision is ~50-100+ m/s.

### MARVELS discovery papers (combined into `marvels_bd/`)
CSV: `marvels_bd/marvels_bd_orbits.csv`, 7 rows, CSV ready, all RV-only.
| Object | Paper | M sin i (MJ) | Notes |
|---|---|---|---|
| TYC 2949-557-1 b | Fleming+2010, ApJ 718, 1186 (arXiv:1006.1627) | 64.3 | only Tc given (Tp blank), e~0 |
| TYC 2930-872-1 b (MARVELS-2b) | Fleming+2012, AJ 144, 72 (arXiv:1206.5514) | 68.1 | only Tc; tertiary trend; true mass likely stellar |
| TYC 2087-255-1 b (MARVELS-4b) | Ma+2013, AJ 145, 20 (arXiv:1211.6140) | 40.0 | MARVELS+TNG fit |
| HIP 67526 b (MARVELS-5b) | Jiang+2013, AJ 146, 65 (arXiv:1307.3308) | 65.0 | prior-dependent true mass 76-244 MJ |
| GSC 03546-01452 b (MARVELS-6b) | De Lee+2013, AJ 145, 155 (arXiv:1304.2815) | 31.7 | BJD_TDB |
| HD 87646 b, c (MARVELS-7b,c) | Ma+2016, AJ 152, 112 (arXiv:1608.03597) | 12.4, 57.0 | host is in a ~22 AU binary; HD 87646 b is in NEA |
- Volumes/pages for the AJ papers come from their DOIs. The Fleming 2010 reference comes from the
  arXiv journal_ref.
- Skipped: MARVELS-1b (Lee+2011; Wright+2013 showed it is a face-on SB2, arXiv:1305.0280),
  MARVELS-3b (Wisniewski+2012, ~98 MJ), and TYC 3010-1494-1 (Mack+2013, an SB2, arXiv:1306.3157).
- MARVELS-4b, 5b and 6b were re-recovered by Grieves+2017 with similar parameters, so the discovery-paper
  values are kept here.

### Troup et al. 2016, AJ 151, 85 (arXiv:1601.00688), "Companions to APOGEE stars I"
- Data: the arXiv source includes the machine-readable `goldOrbit-dr12.fits` (382 gold-sample orbits,
  the same content as VizieR J/AJ/151/85 table4). The VizieR `table4.dat` URL returned 404 (it may only
  exist gzipped), so the FITS file from arXiv was used. VizieR `ReadMe` is saved in `troup_2016/`.
- CSV: `troup_2016/troup_2016_orbits.csv`, 163 rows with 5 <= M sin i <= 85 MJ (111 in 13-80 MJ;
  46 at 5-13 MJ; 6 at 80-85 MJ). All a < 2 AU. Classes: 70 MS, 73 RG, 11 SG, 8 RC, 1 PMS.
  RV-only.
- Conversions: MSINI (Msun) was converted with 1047.57 MJ/Msun, and OMEGA from radians. TPERI is in JD
  (the time standard is not stated). Hosts are 2MASS IDs ("2MASS J..."), which SIMBAD resolves.
- Data quality: these are low-precision orbits (sigma_RV ~100-200 m/s) from 8-34 visits over at most
  ~1000 d. Aliases are common, and the fits include a global linear trend (`slope` in notes).
  SEMIAMP_ERR = 0 means the error is missing (the ReadMe says "?=0"). Many RG hosts have large jitter.
  Treat these as candidates. Many are faint (H = 8-14) and distant (up to kpc), so Gaia sensitivity
  will be low for most; the dwarfs near the bright end are the useful subset.
- Status: CSV ready.

## Giant-star surveys

### Lick giants: Quirrenbach et al. 2019, A&A 624, A18 (arXiv:1904.03557, nu Oph) and Mitchell et al. 2013, A&A (arXiv:1305.5107, tau Gem)
- CSV: `lick_giants/lick_giants_orbits.csv`, 3 rows: nu Oph b (22.2), nu Oph c (24.7, a 5.93 AU),
  tau Gem b (20.6). RV-only.
- nu Oph: the paper gives the mean anomaly M0 at JD 2451853.595. I derived Tp = epoch - M0/360*P
  (for c, this agrees with Teng+2023's independent Tp to within 2 d). Mitchell+2013: A&A
  555, A87, confirmed via Crossref (DOI 10.1051/0004-6361/201321714) on 2026-10-26.
- Quirrenbach+2019 Table `Lick_BDs` lists all Lick-sample BDs: tau Gem b, 11 Com b, nu Oph b, nu Oph c.
  That is the complete Lick-giant BD list; 11 Com b is in `teng_2023`.
- Reffert et al. 2015, A&A 574, A116 (arXiv:1412.4634): its tables contain stellar parameters for 373
  giants only, with no companion orbit table. Not extracted.
- All four objects are in NEA (< 30 MJ).

### Okayama Planet Search: Teng et al. 2023, PASJ (arXiv:2308.05343), "Revisiting planetary systems in OPSP"
- PASJ 75, 1030 (2023); confirmed via Crossref (DOI 10.1093/pasj/psad056) on 2026-10-26.
- Data: RV-only refits of 32 systems (`tab:orbpar_orbits`, omega in rad, Tp JD-2450000), host masses in
  `tab:star_pars`. Joint RV+HGCA fits for 4 systems are in the table "Joint fitting results".
- CSV: `teng_2023/teng_2023_orbits.csv`, 6 rows (RV-only parameters): 11 Com b (15.5), nu Oph b (21.5),
  nu Oph c (23.7), HD 175679 b (36.9 MJ; RV+HGCA true mass 61.8 +15/-9.4 MJ in notes),
  HD 180314 b (22.7), HD 14067 b (8.8 MJ, a 5.11 AU; RV+HGCA true mass 15.3 MJ in notes).
- Not put in the CSV: other 5-13 MJ planets in the table that are already in NEA
  (eps Tau b 7.2, 4 UMa b 7.9, HD 104985 b 5.6, gamma Lib c 5.2, 18 Del b 9.2).
- HD 175679 b (> 30 MJ) is the only new object relative to NEA.

### EAPS-Net (BOAO + OAO): Omiya et al. 2009, PASJ 61, 825 (arXiv:0906.3762)
- CSV: `okayama_eapsnet/okayama_eapsnet_orbits.csv`, 1 row: HD 119445 b (37.6 MJ, P 410 d, 3.9 Msun
  host). RV-only; the orbit table was hand-transcribed from the LaTeX. Tp is in JD.
- The proceedings paper Omiya et al. 2011 (arXiv:1101.4108) repeats the HD 119445 orbit and adds a
  1.7 MJ planet, so it was not used.
- Teng+2023 says EAPS-Net stars (HD 119445, HD 100655, HD 40956, HD 111591) will be updated in a
  "forthcoming paper". That paper is a lead.

### BOAO (Lee et al.)
- The arXiv API listing of B.-C. Lee BOAO papers (2011-2023) shows only planet-mass companions in their
  titles (alpha Ari, eps CrB, HD 66141, beta Cnc, mu Leo, beta UMi, HD 208527, HD 220074, ...).
  Lee+2023 (arXiv:2307.15897, "nine M red giants: detection of sub-stellar companions") may contain a
  BD-mass object. It was not downloaded (see Leads). No CSV.

### PTPS / TAPAS (Niedzielski et al.)
- CSV: `ptps/ptps_orbits.csv`, 5 rows, RV-only (HET/HRS, plus HARPS-N for TAPAS):
  - Niedzielski+2009, ApJ 707, 768 (arXiv:0906.1804), Table 4: BD+20 2457 b (21.4) and c (12.5).
    The host mass is only a range (1.3-4.3 Msun), so mstar is left blank. The system's dynamical
    stability is doubtful (Horner+2014, Trifonov+2014, per Quirrenbach+2019).
  - Niedzielski+2021, A&A 648, A58 (arXiv:2101.10410, TAPAS VII), Table "TableKeplerian":
    HD 4760 b (13.9), BD+02 3313 b (34.1 MJ, P 1393 d, not in NEA), and TYC 434-4538-1 b (6.1).
- T0 is given in MJD and is not defined explicitly as periastron. I converted it to JD and put it in
  Tp_bjd with a "verify" flag in each row.
- Other PTPS papers on the arXiv list (Gettel+2012, Nowak+2013, TAPAS I-IV and VI, Niedzielski+2015)
  report companions <= ~9 MJ that are in NEA. Not extracted.

### EXPRESS (Jones et al.; note that this is a SOUTHERN survey, FEROS/CHIRON)
- CSV: `express_jones/express_jones_orbits.csv`, 2 rows, RV-only:
  HIP 97233 b (20.0 MJ, a 2.55; Jones+2015, A&A 573, A3, arXiv:1409.7429) and
  HIP 67537 b (11.1 MJ, a 4.91; Jones+2017, A&A 602, A58, arXiv:1612.06252). Both are in NEA.
- Not searched further: there are probably other EXPRESS papers (e.g. Jones+2016 and the
  EXPRESS+PPPS papers by Wittenmyer et al.), but they mostly report planets.

## Summary counts (rows in CSV / of which 13-80 MJ and a <= 5 AU)
| short_id | rows | strict in-range BD (13-80 MJ, a<=5) | notes |
|---|---|---|---|
| wilson_2016 | 12 | 11 | HD 134169 at 83 MJ; several are stellar by Hipparcos |
| dalal_2021 | 4 | 4 | HD 205521 is stellar by astrometry |
| diaz_2012 | 7 | 6 | HD 156279 b is 9.7 MJ; HD 16702 is stellar |
| bouchy_2016 | 3 | 0 | all a = 5.95-6.4 AU (borderline) |
| patel_2007 | 1 | 0 | a 5.47 AU; already in CLS |
| grieves_2017 | 10 | 10 | candidates |
| marvels_bd | 7 | 6 | HD 87646 b 12.4 MJ |
| troup_2016 | 163 | ~111 | low precision, candidates |
| lick_giants | 3 | 2 | nu Oph c a 5.9 AU; all in NEA |
| teng_2023 | 6 | 4 | nu Oph c a 5.8; HD 14067 b 8.8 MJ |
| okayama_eapsnet | 1 | 1 | |
| ptps | 5 | 3 | BD+20 2457 b,c doubtful |
| express_jones | 2 | 1 | southern |
nu Oph b/c appear in both lick_giants and teng_2023 (two independent fits), so de-duplicate them.

---

## Leads
Other agents' directories already exist in `sample/raw_lit/` for feng_2022, ma_ge_2014,
sahlmann_2011, stevenson_2023, unger_2023, sb9, kiefer_2021, exoplanets_org and others. Those are not
repeated here.

1. **Wilson+2016 Table A.1** (arXiv:1602.02749, in `wilson_2016/src_1602.02749/rv.tex`,
   `tab:appendix_table`): ~60 BD candidates (10-80 MJ) with references. Use it as a checklist to find
   northern BDs missing from the pool, e.g. HD 168443 c, HD 38529 c, HD 136118 b, HD 190228 b,
   HD 191760 b, HD 189310 b, HD 52756 b, HD 89707 b, HD 91669 b, HD 174457 b, HD 30339 b, HD 162020 b,
   HD 202206 b. Many come from Sahlmann+2011, Nidever+2002, Butler+2006 or Wittenmyer+2009.
2. **Nidever et al. 2002, ApJS 141, 503** (Keck/Lick, "RVs for 889 late-type stars"): has SB1 orbits,
   including BD-desert companions (HD 174457, HD 30339, HD 65430 per Wilson's table). Check VizieR
   J/ApJS/141/503 for an orbit table; otherwise the arXiv source.
3. **Latham et al. 2002, AJ 124, 1144** (CfA, northern, "Orbits for 171 single-lined binaries"):
   machine-readable on VizieR (J/AJ/124/1144). It covers HD 114762, G9-42 (=HD 77065) and G66-65
   (=HD 134113). It gives f(m), so M sin i needs a primary mass; worth a targeted pass for small-K orbits.
4. **Wittenmyer et al. 2009, ApJS 182, 97** (HET/McDonald long-term survey): HD 91669 b (30.6 MJ) and
   HD 136118 b. Look for the arXiv source.
5. **EAPS-Net update paper** promised by Teng+2023 (HD 119445, HD 40956, HD 111591): check arXiv
   2024-2026 for Teng or Omiya.
6. **Lee et al. 2023, arXiv:2307.15897** (BOAO M giants, "sub-stellar companions"): check for BD-mass
   objects with P < 11 yr.
7. **Price-Whelan et al. 2020, ApJ 895, 2** (APOGEE DR16 binaries, The Joker): machine-readable
   catalog of ~20k SB1 orbits with posterior samples. Mostly stellar, but its low-K tail overlaps BDs and
   it supersedes many Troup+2016 candidates. A newer APOGEE DR17 catalog may also exist (not verified).
8. **Halbwachs et al. 2000, A&A 355, 581** (ELODIE/CORAVEL BD candidates with Hipparcos masses):
   northern ELODIE orbits, mostly re-observed in Bouchy+2016 and Kiefer+2019.
9. **Kane et al. 2019 / Stevenson et al. 2023**: Keck and Lick legacy refits. Stevenson+2023 already
   has a directory here.

## Data-quality concerns (summary)
- Troup+2016: low RV precision, few epochs, aliases, and large giant jitter. Many hosts are faint or
  distant.
- MARVELS (Grieves+2017 and discovery papers): ~50-200 m/s precision. Several "BDs" are stellar per
  the authors (MARVELS-2b, possibly 5b).
- Several SOPHIE "BD candidates" are stellar by Hipparcos (Wilson: HD 6664, BD+73 275, HD 51813,
  HD 56709; Díaz: HD 16702; Dalal: HD 205521). Their RV orbits are still valid ground truth, but their
  Gaia signals will be large, so do not count them as substellar.
- Time standards are mixed (JD, HJD, BJD); each row's notes say which. The PTPS T0 definition is not
  stated.
- Citations: Bouchy+2016 (A&A 585, A46), Grieves+2017 (MNRAS 467, 4264), Mitchell+2013 (A&A 555,
  A87), and Teng+2023 (PASJ 75, 1030) were confirmed via Crossref on 2026-10-26. All other
  citations were confirmed from the arXiv journal_ref or DOI.
