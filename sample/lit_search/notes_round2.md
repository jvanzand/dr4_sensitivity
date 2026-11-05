# Literature search, round 2 (2026-10-26)

Round 2 adds the sources listed as leads in `README.md`. Each source's downloads, parsing
script (`build.py`, run with `python -I build.py`) and clean table
(`<short_id>_orbits.csv`, standard 15 columns) are in `../raw_lit/<short_id>/`. Nothing here
is merged into the census yet.

"In range" means M sin i ≤ 80 M_J and (a ≤ 5 AU or P ≤ 4000 d). "BD" counts rows with
13 ≤ M sin i ≤ 80 M_J.

## Summary

| Source | Directory | Rows in range | BD | Fit type | Status |
|---|---|---:|---:|---|---|
| Ren+2026 | `ren_2026` | 1 | 1 | RV+astrometry | done (only one companion tabulated) |
| Cui+2026 | `cui_2026` | 23 | 7 | RV+astrometry | done (republishes Feng+2022) |
| Barbato+2026 (GAPS LXXV) | `barbato_2026` | 8 | 6 | RV-only | done (no Tp) |
| Jaros+2026 (PTPS/JOTA) | `jaros_2026` | 0 | 0 | RV-only | skipped: all companions ≥ 103 M_J |
| Latham+2002 | `latham_2002` | 9 | 9 | RV-only | done (M sin i from f(m)) |
| Nidever+2002 | `nidever_2002` | 6 | 6 | RV-only | done |
| SB9 with TIC/FLAME masses | `sb9_tic` | 13 | 11 | RV-only | done (3 duplicate rows, 2 sdB hosts) |
| Wittenmyer+2009 (+ HD 91669 paper) | `wittenmyer_2009` | 26 | 2 | RV-only | done |
| Lee+2023 (BOAO M giants) | `lee_2023` | 2 | 2 | RV-only | done |
| EAPS-Net update | none | 0 | 0 | | skipped: no 2024-2026 paper on arXiv |
| Halbwachs+2000 | none | 0 | 0 | | skipped: not on arXiv or VizieR |

Total: 88 rows in range, 44 of them BD. 64 are RV-only.

Duplicates within round 2:
- Wolf 272, G 42-30, BD+07 2634 and G 65-52 are in both `latham_2002` and `sb9_tic`, because
  SB9 takes those orbits from Latham+2002.
- HD 283750 appears three times in `sb9_tic`: three independent SB9 orbits of the same system.
- HD 65430 is in `nidever_2002` at 67.8 M_J. Latham fits it to 90 M_J, out of range.

## 2026 arXiv papers

### Ren et al. 2026 (`ren_2026`)
- **Citation:** Ren, L., Gan, T., Xuan, Y., Mao, S., Ida, S., Xiao, G.-Y., Feng, F., & Pallé, E.
  2026, "An Intermediate-orbit Brown Dwarf Confirmed with Astrometry+Radial Velocity Joint
  Analysis, and the Eccentricity Distribution of Brown Dwarfs Across Different Locations,"
  arXiv:2610.16662.
  - The arXiv comment says "accepted for publication in AJ". There is no journal_ref or DOI,
    so the volume and page are UNCONFIRMED.
- **Data:** LaTeX source, Table 2 (`allpriors`), in `src_2610.16662/main.tex`.
- **The 159-BD sample is not tabulated.** It appears only in a figure. It is compiled from
  Gan+2025 (transiting BDs), Barbato+2026 and Unger+2023 (RV + Gaia two-body), Feng+2022
  (RV + Hipparcos-Gaia) and Bowler+2020 (imaging). So the only companion recorded is HD 203771 B.
- **Rows:** 1 in range, 1 BD.
  - HD 203771 B (HIP 105740): P = 1140 d, true mass 69.0 M_J at i = 129.6°, so M sin i = 53.2 M_J.
- **Fit type:** RV+astrometry. The paper gives no RV-only fit.
  - I used the RV+Hip+GDR2+GDR3 column because it is the only one that reports K.
  - The RV + Gaia DR3 two-body column is summarized in the row's notes (74.2 M_J, i = 133.6°).
- **Tp:** JD 2457389 + t_p, as the table labels it, so tp_system is JD.
- **Overlaps:** a new fit of a companion already in `feng_2022` (P = 1112 d) and `stevenson_2023`.

### Cui et al. 2026 (`cui_2026`)
- **Citation:** Cui, K., Xiao, G.-Y., Feng, F., Liu, B., Nayakshin, S., Hall, C., Guo, K., Lai, D.,
  Ogihara, M., Rui, Y., Boss, A. P., Butler, R. P., & Xuan, Y. 2026, "A universal brown dwarf
  desert formed between planets and stars," PNAS, arXiv:2603.01808, DOI 10.1073/pnas.2524764123.
  - The DOI comes from arXiv. The volume and article number are UNCONFIRMED.
- **Data:** LaTeX longtable "Catalog of companions" (`tab:sample`, Supplementary Information) in
  `src_2603.01808/BDdesert.tex`.
  - The table gives true mass m, q, e, a, P (yr), I and Ω only. There is no K, ω or Tp.
- **Republication:** these are not new fits. The paper selects 55 companions from Feng+2022's
  RV + Hipparcos-Gaia catalog (ApJS 262, 21) and uses Feng+2022's MCMC posteriors.
  - Cite Feng+2022 as the original.
  - Treat this table as a cross-check of `feng_2022`. 20 of the 23 in-range hosts already appear
    there or in other tables. HD 145675, HD 150554 and HIP 22203 had no exact host-name match.
- **Rows:**
  - 55 parsed and written to `cui_2026_all.csv` with an `in_range` flag.
  - 23 in range (`cui_2026_orbits.csv`), 7 of them BD: HD 217786, HD 39213, HD 127506, HD 150554,
    HD 168443, HD 8673 and HIP 22203.
- **Fit type:** RV+astrometry.
- **Issues:**
  - M sin i is computed here as m sin I.
    - HD 39213 has a true mass of 118.9 M_J, but its M sin i is 69.3 M_J, so it passes the cut.
    - HD 217786 has a true mass of 37.4 M_J at I = 21°, so its M sin i is 13.5 M_J.
  - P is given in years to 0.01 yr (about 3.7 d) precision. Many printed P errors are 0.00 yr;
    those e_P_d cells are left blank.
  - The stellar mass is not tabulated (the paper uses TIC v8.2), so mstar is blank.

### Barbato et al. 2026, GAPS LXXV (`barbato_2026`)
- **Citation:** Barbato, D., Pinamonti, M., Sozzetti, A., Desidera, S., D'Orazi, V., Maldonado, J.,
  et al. (25 authors) 2026, "The GAPS programme at TNG LXXV. Validating and confirming Gaia
  substellar astrometric candidates with HARPS-N," arXiv:2606.29001.
  - The arXiv comment says "accepted for publication on A&A". There is no journal_ref or DOI,
    so the volume and page are UNCONFIRMED.
- **Data:** LaTeX `table:rvfit` (RV fit block) and `table:star-validated` (Gaia DR3 IDs and SED
  masses) in `src_2606.29001/aa57585-25.tex`.
- **Rows:** 8 in range, 6 BD.
  - BDs: 2MASS J20570943+1406590 (30.3 M_J), HIP 105707 (20.3), HIP 60865 (43.9),
    TYC 3456-251-1 (59.5), TYC 4998-437-1 (22.8) and StKM 1-1494 (63.4).
  - Planets: TYC 198-493-1 (8.6) and PM J13580+3141 (10.3).
  - All are Gaia DR3 NSS astrometric candidates. The other 6 targets turned out to be SB2s.
- **Fit type:** RV-only (PyORBIT).
  - Most fits use HARPS-N RVs only.
  - PM J13580+3141 (= Gaia-4 b of Stefansson+2025) adds literature HPF, NEID and FIES RVs and a
    quadratic trend. The adopted solution and the single-Keplerian alternative are both in notes.
- **Republication:** 6 of the 8 are first RV orbits.
  - PM J13580+3141 was first characterized by Stefansson+2025 (RV + Gaia).
  - StKM 1-1494 was confirmed by GRAVITY imaging (Winterhalder+2024).
- **Issues:**
  - **Tp is blank.** The fits use mean longitude λ0, and the paper never states its reference epoch.
  - **Time system:** the RV table uses BJD, but this is not needed with Tp blank.
  - **Minor inconsistency in the paper:** its text says minimum masses run "8 to 62 M_J" and
    a "0.75/0.76 to 1.42 au". The table gives 8.63–63.35 M_J and 0.75–1.47 au. The table values
    are used.
  - **Already in `exoplanet_eu`:** HIP 60865 and TYC 198-493-1 appear there with identical P and
    e, apparently taken from this paper. HIP 105707 appears with a different, earlier orbit.
  - **Disturbed astrometry:** RUWE is 1.2–3.8. These are by construction Gaia DR3 NSS orbit
    stars, so they are good DR4 validation targets.

### Jaros et al. 2026, PTPS/JOTA (`jaros_2026`)
- **Citation:** Jaros, R., Niedzielski, A., Ali, S., Srivastava, D., Sierzputowska, J., Adamów, M.,
  Lillo-Box, J., Orikhovskyi, D., Wolszczan, A., Villaver, E., Pribulla, T., & Vaňko, M. 2026,
  "Joint Observations of PTPS Targets (JOTA). Combined spectroscopic and photometric analysis of
  16 SB1 systems," arXiv:2608.30546.
  - The journal is UNCONFIRMED: arXiv gives no journal_ref, DOI or comment. The `aa.cls` file
    suggests A&A.
- **Data:** LaTeX Table "Results of Keplerian modeling" in `src_2608.30546/main.tex`.
- **Rows:** 0 in range, so no CSV was written.
  - All 16 companions have m sin i between 0.0983 and 1.38 M_sun, which is at least 103 M_J.
  - The lowest is TYC 2267-00101-1 at 0.0983 M_sun.
- **Republication:** 5 systems are updates of previously known binaries. The others are new
  detections, most around giants.
- **Status:** skipped (out of range). The source is kept. Note that the table's epoch column is
  T_c, not periastron.

## Older SB1 catalogs

### Latham et al. 2002 (`latham_2002`)
- **Citation:** Latham, D. W., Stefanik, R. P., Torres, G., Davis, R. J., Mazeh, T., Carney, B. W.,
  Laird, J. B., & Morse, J. A. 2002, AJ, 124, 1144, "A survey of proper-motion stars. XVI.
  Orbital solutions for 171 single-lined spectroscopic binaries."
  - The journal, volume, page and bibcode 2002AJ....124.1144L are confirmed from the VizieR ReadMe.
  - No arXiv ID was looked up.
- **Data:** VizieR J/AJ/124/1144.
  - `table4.dat`: 156 robust orbits.
  - `table5.dat`: 15 preliminary orbits.
  - `table1.dat`: names from Paper XII.
- **Method:**
  - Names were resolved through SIMBAD/Sesame (171 of 171).
  - M1 is the Gaia DR3 FLAME mass, else the TIC v8 mass.
  - HD 169822 and Ross 783 have neither. Their M1 falls back to the paper's own Table 6 value,
    and their notes say so.
  - M2 sin i is solved from the published f(m) with the minimum-mass approximation.
  - Tp is the catalog's t0 + 2400000. VizieR labels t0 as HJD, so tp_system is HJD.
  - All 171 orbits are in `latham_2002_all_msini.csv`.
- **Rows:** 9 in range, all 9 BD, all RV-only. G 42-30 is a preliminary (Table 5) orbit.
  - Spot check: HD 77065 (G 9-42) matches table4.dat for P, K, e, ω and T0.
- **Issues:**
  - **Primary masses differ from the paper's.** The paper's Table 6 masses for these metal-poor
    stars are 0.6–0.8 M_sun, and FLAME and TIC usually give more.
    - G 178-27, G 15-6, HD 118742 and HD 210631 rise to 82–90 M_J and drop out, although the
      paper's masses would keep them below 80 M_J.
  - **Disturbed astrometry:** many kept hosts have RUWE > 1.4 or a non-single-star flag
    (HD 217580 RUWE 18, G 42-30 8.2, G 65-52 4.9). Their FLAME masses may be biased.
  - **Evolved host:** HD 134113 is post-main-sequence per FLAME.
  - **No mass:** 10 stars have no usable mass.
  - **Velocity system:** K is on the CfA native system.

### Nidever et al. 2002 (`nidever_2002`)
- **Citation:** Nidever, D. L., Marcy, G. W., Butler, R. P., Fischer, D. A., & Vogt, S. S. 2002,
  ApJS, 141, 503, "Radial Velocities for 889 Late-type Stars," arXiv:astro-ph/0112477.
  Confirmed from the arXiv journal_ref and the VizieR ReadMe.
- **Data:** VizieR J/ApJS/141/503 `table3.dat` (15 SB1 orbits). The arXiv LaTeX is in `arxiv_src/`.
- **Method:**
  - M2 sin i, M1 and a are the paper's own Table 3 values.
  - The notes also give M2 sin i recomputed with FLAME or TIC masses. The two agree within about 15%.
- **Rows:** 6 in range, 6 BD, all RV-only.
  - In range: HD 30339 (77.8), HD 65430 (67.8), HD 140913 (43.2), HD 174457 (65.8), GJ 595 (60.0)
    and GJ 623 (42.0).
  - HD 4747 (42.3 M_J) is out because P = 6832 d and a = 6.7 AU.
- **Issues:**
  - **Time system:** T0 is given as "T0 − 2450000" without a time system, so tp_system is unknown.
  - **HD 140913:** e is fixed in the fit.
  - **Earlier orbits:** HD 140913 and GJ 623 have earlier published orbits; the notes list the
    references.
  - **Disturbed hosts:** HD 174457 is evolved (log g 3.74) with RUWE 14. HD 65430 and GJ 595 have
    RUWE > 1.4.

### SB9 orbits without a primary mass (`sb9_tic`)
- **Citation:** Pourbaix, D., et al. 2004, A&A, 424, 727 (SB9; VizieR B/sb9).
  - This was taken from round 1 and not re-verified.
  - Each row's table_ref gives the bibcode of the original orbit paper.
- **Data:** the 242 orbits in `../sb9/sb9_excluded_no_primary_mass.csv`, with full elements from
  `../sb9/tap_*.csv`. These were only read; `../sb9/` is unchanged.
- **Method:**
  - M1 is the FLAME mass, else the TIC v8 mass. The TIC ID comes from SIMBAD, else from the
    Gaia DR2 ID.
  - f(m) comes from P, K and e, and M2 sin i is solved with the minimum-mass approximation.
  - Tp is SB9's T0, labelled "Periastron time (JD)" in the ReadMe, so tp_system is JD. It is
    blank where SB9 flags T0 as not a periastron epoch.
  - All 242 orbits are in `sb9_tic_all.csv`.
- **Rows:**
  - Only 72 of the 242 got a mass, 24 of them from FLAME. TIC v8 has no mass for giants or hot
    stars, which make up most of the list.
  - 13 rows are in range, 11 of them BD, all RV-only.
- **Issues:**
  - **HD 283750 appears three times;** remove two.
  - **Duplicates of Latham+2002:** 4 rows.
  - **G 178-27 sits at the threshold.** The SB9 orbit (Latham+1992) gives 79.5 M_J; Latham+2002's
    refit gives 81.6 M_J.
  - **Two hot-subdwarf hosts:** BPS CS 22169-0001 and CPD-64 481. Their TIC mass of 0.47 M_sun is
    unreliable; drop them.
  - **OGLE-TR 113:** two rows for the transiting planet (1.0–1.3 M_J).
  - **Unreliable giant mass:** 2MASS J19413531+4011459 (RUWE 30).

## Individual papers

### Wittenmyer et al. 2009 (`wittenmyer_2009`)
- **Citations:** this directory holds two papers.
  - Wittenmyer, R. A., Endl, M., Cochran, W. D., Levison, H. F., & Henry, G. W. 2009, ApJS, 182, 97,
    "A Search for Multi-Planet Systems Using the Hobby-Eberly Telescope," arXiv:0903.0652.
    Confirmed from DOI 10.1088/0067-0049/182/1/97.
  - HD 91669 b comes from Wittenmyer, R. A., Endl, M., Cochran, W. D., Ramírez, I., Reffert, S.,
    MacQueen, P. J., & Shetrone, M. 2009, AJ, 137, 3529, arXiv:0812.4094. Confirmed from
    DOI 10.1088/0004-6256/137/3/3529.
- **Data:** ApJS LaTeX Table 3 (25 companions in 22 systems), with M* from Table 1. AJ paper
  Tables 1–2.
- **Rows:** 26, all in range. 2 are BD: HD 38529 c (13.38 M_J) and HD 91669 b (30.6 M_J). The rest
  are known planets.
- **Fit type:** RV-only.
  - The ApJS fits are new combined fits of HET, McDonald 2.7 m and published RVs.
  - HD 91669 b uses 2.7 m data only.
- **Issues:**
  - **Tp definition:** Tp is "T_0 (JD−2400000)" + 2400000, with tp_system JD. Neither paper defines
    T_0 as periastron. It is consistent with periastron via Kiefer+2021 (HD 136118) and Unger+2023
    (HD 91669).
  - **HD 128311 c:** T_0 = 2456987 lies after the data span; it is kept as printed.
  - **Overlaps:** many, including kiefer_2021, unger_2023, feng_2022, xiao_2023 and
    stevenson_2023. Most hosts are NEA planet hosts.

### Lee et al. 2023, BOAO M giants (`lee_2023`)
- **Citation:** Lee, B.-C., Do, H.-J., Park, M.-G., Lim, B., Choi, Y.-H., Koo, J.-R., Bang, T.-Y.,
  Oh, H.-I., Han, I., & Chang, H.-Y. 2023, "Long-period radial velocity variations of nine
  M red giants: The detection of sub-stellar companions around HD 6860 and HD 112300,"
  arXiv:2307.15897.
  - The journal is UNCONFIRMED: arXiv gives no journal_ref or DOI. The file name AA43725.tex
    suggests A&A.
- **Data:** LaTeX Table 3. M* comes from the text.
- **Rows:** 2 in range, 2 BD, both RV-only (BOES).
  - HD 6860 b (β And): 28.3 M_J, P = 664 d.
  - HD 112300 b (δ Vir): 15.8 M_J, P = 467 d.
  - Tp is "T_periastron (JD)", so tp_system is JD.
- **Issues:**
  - **Possible intrinsic signals:** both hosts are pulsating M giants with jitter of about
    100 m/s.
  - **HD 112300:** the fit includes a linear trend.
  - **Bright hosts:** V = 2.05 and 3.38, so their Gaia astrometry is likely unreliable.
  - **Already in `exoplanet_eu`,** without K or Tp.

### EAPS-Net update after Teng+2023
- **Search:** the arXiv API, by author (Teng_H, Teng_Huan, Omiya_M), by target (HD 119445,
  HD 40956, HD 111591) and by keyword. No 2024–2026 update paper exists.
- **Checked and skipped:**
  - Xiao, Teng et al. 2023 (arXiv:2312.01388): Jupiter-mass planets with RV + HGCA fits.
  - Jeong et al. 2018 (arXiv:1711.07173): planets of about 2.7 M_J.
- **Status:** skipped.

### Halbwachs et al. 2000, A&A 355, 581
- **Not accessible:** it is not on arXiv, and VizieR J/A+A/355/581 returns 404. aanda.org is blocked.
- **Status:** skipped. The citation is unverified.

## Leads

1. **Gan et al. 2025 transiting brown dwarf compilation** (cited by Ren+2026 for its 159-BD
   sample). Transiting BDs have RV orbits at short P. They are below Gaia's astrometric
   sensitivity but useful as a control sample. Not yet checked on arXiv.
2. **Stefansson et al. 2025:** RV follow-up of Gaia DR3 NSS substellar candidates (including
   Gaia-4 b, per Barbato+2026). Barbato+2026 quotes it for about 32% binary contamination. Its RV-only and joint orbits
   would add Gaia-validated companions. Not yet checked on arXiv.
3. **Better masses for SB9 giants:** 170 of the 242 SB9 orbits still have no M1, mostly giants.
   Asteroseismic or spectroscopic giant masses would rescue some.
4. **A Latham variant list:** with the paper's own Table 6 masses, G 178-27, G 15-6, HD 118742
   and HD 210631 fall below 80 M_J. They would enter the list if the paper's metal-poor masses
   are preferred over FLAME or TIC.
5. **A cross-check for HD 140913:** Mazeh+1996 (VizieR J/ApJ/466/415, also in SB9) gives an
   independent CfA orbit to compare with Nidever+2002.
6. **Feng+2022 RV-only refits:** Cui+2026 and Ren+2026 both reuse Feng+2022 posteriors. RV-only
   refits of the 23 in-range Cui companions remain the way to make them ground truth.
