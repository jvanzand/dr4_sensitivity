# Literature search for RV companions

Updated 2026-10-30 (RV+astrometry fits admitted). Goal: a comprehensive pool of RV
companions, especially brown dwarfs above the NASA Exoplanet Archive's (NEA's) 30 M_J
limit, to map Gaia DR4 sensitivity.

Four search passes produced the source tables:
1. Southern surveys, northern surveys, and aggregate catalogs (2026-10-24; notes in
   `notes_southern.md`, `notes_northern.md`, and `notes_aggregate.md`).
2. The remaining accessible sources (2026-10-26; notes in `notes_round2.md`).

Each source's raw downloads and transcribed orbit table are in `../raw_lit/<source>/`.
`literature.py` reads them, and `census.py` merges them with the NEA.

## Summary

- **Sources:** 41 sources and 1484 transcribed orbit rows, from arXiv LaTeX sources,
  VizieR, the Stevenson+2023 GitHub catalog, and the Open Exoplanet Catalogue's mirror of
  exoplanet.eu.
- **Audit:** every table was checked against its raw source. Periods, K, eccentricities,
  and M sin i were matched to the LaTeX, VizieR, or FITS files. All mismatches were
  correct unit conversions. A row-completeness check of Wilson+2016 found 2 omitted
  companions, now added.
- **Census:** after removing 745 NEA duplicates and applying the cuts (M sin i ≤ 80 M_J,
  SMA ≤ 5 AU), the literature adds **294 companions**. The census holds **2667**, including
  **338 brown dwarfs**. **189** are plausibly recoverable (S/N ≥ 1.5, P ≤ 5.5 yr): 94 from
  the literature and 95 from the NEA.
- **exoplanet.eu is down for everyone,** not only on our networks. On 2026-10-26 it timed
  out from this machine and refused connections from Anthropic's fetch servers, and the
  whole Observatoire de Paris domain (obspm.fr) was unreachable. The nightly mirror stands
  in for it.

## Rules applied when merging

- **RV+astrometry fits are allowed** (amended 2026-10-30; previously RV-only). Ground
  truth may come in any published form, since the recovery stage uses only DR4. This
  admits Xiao+2023, Ren+2026, Cui+2026, and the RV+astrometry rows of Unger+2023 and
  Stevenson+2023. Their `msini_mj` is the published M sin i or the true mass times sin i,
  and `fit_type` records the fit. Of the 192 usable RV+astrometry rows (Feng+2022 aside),
  only 5 enter the census: 105 duplicate NEA planets (101 of them in the census), 81
  duplicate a higher-ranked literature row (often the same paper's RV-only fit, as in
  Unger+2023), and 1 fails the cuts. Of Feng+2022's 114 rows, 54 match NEA planets.
- **Feng+2022 is excluded** (`EXCLUDED_SOURCES` in `literature.py`; `census.py` also skips
  its NEA solutions), because its fits are considered suspect.
- **Precedence**, when a companion appears in several sources (same Gaia DR3 source and
  period within 10% of any published value):
  1. NEA;
  2. survey papers with full orbits;
  3. RV+astrometry papers (Xiao+2023, Ren+2026, Cui+2026);
  4. low-precision candidates (MARVELS, APOGEE/Troup+2016);
  5. compilations, starting with those that list full elements;
  6. single-lined binary catalogs whose M sin i comes from the mass function with an
     assumed primary mass (Latham+2002, SB9).
- **Planet regime (< 13 M_J):** the NEA's vetting governs. Compilations, low-precision
  lists, and mass-function catalogs add a planet only as a stand-in for an NEA planet that
  lacks a usable NEA solution. 136 other entries are set aside in
  `../raw_lit/review_planet_regime.csv`.
- **Exclusions:** exoplanet.eu hosts that are white dwarfs, cataclysmic variables,
  pulsars, or hot subdwarfs are dropped, as are two hot-subdwarf SB9 hosts (BPS
  CS 22169-0001, CPD-64 481). Where exoplanet.eu lists an inclination, its true mass is
  converted back to M sin i.
- **Periastron times:** each source's Tp is tagged with its stated time system (BJD_TDB,
  BJD, BJD_UTC, HJD, JD, or unknown). Older tables are tagged per source in
  `literature.py` (`TP_SYSTEM_BY_SOURCE`); newer tables carry a `tp_system` column.
  `census.py` converts everything to BJD_TDB.
- **Host matching:** hosts are matched to Gaia DR3 through SIMBAD (exact identifiers, then
  the name-normalizing script interface), then through Gaia's 2MASS crossmatch. When CDS
  doesn't answer, `literature.py` uses SIMBAD's Harvard mirror (used on 2026-10-30).
- **Stellar masses:** from the source, else the NEA, else Gaia DR3 FLAME.

## Source inventory

"Rows" counts the transcribed table; "BD" counts rows with 13 ≤ M sin i ≤ 80 M_J; "Added"
counts census companions after de-duplication and cuts; "Recoverable" uses S/N ≥ 1.5 and
P ≤ 5.5 yr.

| Source | Tier | Rows | BD | Added | Recoverable |
|---|---|---:|---:|---:|---:|
| CLS I (Rosenthal+2021) | survey (VizieR) | 227 | 21 | 5 | 2 |
| Kiefer+2019 SOPHIE | survey (VizieR) | 11 | 9 | 7 | 4 |
| Sahlmann+2011 CORALIE | survey | 20 | 19 | 13 | 10 |
| Unger+2023 CORALIE/HARPS | survey | 61 (31 RV-only, 30 RV+astrometry) | 46 | 15 | 15 |
| Barbato+2023 CORALIE XIX | survey | 29 | 28 | 8 | 5 |
| Barbato+2026 GAPS LXXV | survey | 8 | 6 | 7 | 7 |
| Wilson+2016 SOPHIE | survey | 14 | 13 | 10 | 8 |
| Díaz+2012, Dalal+2021, Bouchy+2016 SOPHIE | survey | 14 | 13 | 5 | 4 |
| Nidever+2002 Keck/Lick | survey | 6 | 6 | 3 | 2 |
| Wittenmyer+2009 HET | survey | 26 | 2 | 0 | 0 |
| Parc+2026, Frensch+2023, Jenkins+2009/2015, Rickman+2019, Ségransan+2010 | survey | 20 | 12 | 5 | 2 |
| Teng+2023, Omiya+2009, PTPS, Lick giants, EXPRESS, Patel+2007, Lee+2023 | survey | 20 | 16 | 3 | 2 |
| Grieves+2017, MARVELS discovery papers | low precision | 17 | 16 | 16 | 4 |
| Troup+2016 APOGEE | low precision | 163 | 111 | 107 | 4 |
| Xiao+2023, Ren+2026, Cui+2026 | RV+astrometry | 100 | 30 | 4 | 1 |
| exoplanet.eu (mirror) | compilation | 229 | 155 | 37 | 12 |
| Stevenson+2023 | compilation | 153 (74 RV-only, 79 RV+astrometry) | 105 | 23 | 4 |
| Ma & Ge 2014, Kiefer+2021, exoplanets.org, Gaia DR3 validation | compilation | 168 | 75 | 4 | 1 |
| Latham+2002, SB9, SB9 with TIC masses | mass function | 84 | 71 | 22 | 7 |
| Feng+2022 | excluded (suspect) | 114 | 20 | 0 | 0 |
| Jaros+2026 (all companions above 100 M_J) | — | 0 | 0 | 0 | 0 |

## Reproducibility

- **Every source is scripted.** Each `../raw_lit/<source>/build.py` regenerates
  `<source>_orbits.csv` from the raw files in its directory (arXiv LaTeX sources, VizieR
  tables, the Stevenson+2023 GitHub catalog, the exoplanet.eu mirror in `../raw_lit/oec/`).
  The scripts are self-contained (standard library, numpy, pandas, astropy, requests) and
  run from any directory. Each docstring lists its inputs, their URLs, and the
  hand-written per-row text (notes, table references) the script keeps.
- **Rebuild everything** from the project root with `python -I sample/build_all.py`. It runs
  every source script, then `literature.py --cached`, `census.py --cached`, and the plots,
  without using the network. Rebuilding twice gives byte-identical tables.
- **Flags.** The first-round scripts take `--download` (fetch missing raw files) and, where
  they look anything up, `--refresh`. The second-round scripts (`barbato_2026`, `cui_2026`,
  `latham_2002`, `lee_2023`, `nidever_2002`, `ren_2026`, `sb9_tic`, `wittenmyer_2009`) take
  no arguments. They read files already in their directories, and query SIMBAD, Gaia, or
  the TIC only for names missing from their caches. `build_all.py` passes each flag only
  to the scripts that accept it.
- **Lookup caches.** Live lookups are stored in the source directories and used offline.
  - The aggregate sources (`exoplanet_eu`, `exoplanets_org`, `stevenson_2023`,
    `gaia_dr3_val`, `feng_2022`, `xiao_2023`, `sb9`) keep `simbad_cache.csv` (main ID,
    object type, identifiers), `nea_pscomppars_cache.csv` (the NEA rows they name-match),
    and `literature_companions_cache.csv` (CLS I and Kiefer+2019 rows).
  - `sb9` also keeps its Gaia DR3 FLAME and astrometry tables and the Mamajek
    spectral-type table, and reads the VizieR TAP dumps `tap_main.csv`, `tap_orbits.csv`,
    and `tap_alias.csv` (`sb9_tic` reads them too).
  - `literature.py` keeps `../raw_lit/simbad_gaia_ids.csv`. Hosts that SIMBAD can't
    resolve are cached with an empty ID, so cached runs never retry them.
  - The caches were seeded from the original 2026-10-24 query results, not re-queried, so
    they reproduce the tables exactly.
- **`--refresh` changes results.** It re-queries the NEA, Gaia, VizieR, and SIMBAD
  (through `build_all.py --refresh`). Today's NEA can change the aggregate tables' NEA-overlap
  notes, and which exoplanet.eu and exoplanets.org rows are kept. The NEA and Gaia
  downloads are dated, with checksums, in `../raw/SNAPSHOT.json`.
- **Acceptance check (2026-10-26).** Every rebuilt table matched the transcribed table it
  replaced: the same rows, numbers within a relative 1e-9, and the same text. No
  transcription errors were found. The Okayama, PTPS, and EXPRESS tables were typed by
  hand; their scripts now parse the papers' LaTeX tables and give the same values.
- **Original code.** The first-pass agents' parsing code is kept in `build_scripts/` for
  reference. It hard-codes temporary paths and isn't run.

## Remaining leads

1. **APOGEE binary catalogs** (Price-Whelan et al. 2020, ApJ 895, 2; any DR17 successor
   is unverified). This is the largest untapped pool, with mostly distant hosts.
2. **Feng+2022 companions**, excluded as suspect: 114 transcribed rows, plus 28 NEA
   planets whose only NEA solution is Feng+2022. They could return through our own fits
   or a later decision about Feng+2022.
3. **Stefansson+2025** (RV follow-up of Gaia DR3 orbit candidates) and **Gan+2025** (a
   transiting brown-dwarf compilation), flagged in `notes_round2.md`.
4. **Better primary masses for SB9 giants:** 170 SB9 orbits still lack a usable primary
   mass.
5. **exoplanet.eu directly**, once it's back up (K, Tp, and separate M sin i and true-mass
   columns).

## Known issues

- **Time systems:** Barbato+2023 Tp values were reconstructed from mean longitudes at an
  inferred epoch (BJD 2455500, consistent with 16 independent Tp values). Barbato+2026
  and Parc+2026 have no Tp, since their papers don't state the reference epoch.
- **Low precision:** Troup+2016 and the MARVELS orbits come from few epochs at 50-200 m/s
  precision and may be aliases.
- **Stellar companions:** several survey "brown dwarfs" are stellar in true mass, by
  Hipparcos (e.g., HD 6664, BD+73 275, HD 51813, HD 16702, HD 205521). Their RV orbits
  remain valid ground truth.
- **Pulsating giants:** Lee+2023's companions orbit pulsating M giants (β And, δ Vir),
  whose RV signals may be intrinsic.
- **Citations:** all first-pass citations are confirmed, the flagged ones via Crossref on
  2026-10-26. Still unconfirmed: volume and page for Ren+2026, Cui+2026, and
  Barbato+2026, and the journal for Jaros+2026 and Lee+2023 (see `notes_round2.md`).
