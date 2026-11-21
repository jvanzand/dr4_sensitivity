"""Census of RV companions for the Gaia DR4 sensitivity study.

Queries the NASA Exoplanet Archive for planets with RV data, picks one solution per
planet (RV-only or RV+astrometry; newest paper per system unless another has a 2x
better period; Feng+2022 and Gaia-only orbits excluded), adds brown dwarfs from survey catalogs (literature.py), cross-matches the
hosts to Gaia DR3, and forecasts each companion's astrometric signal-to-noise and whether
DR4 could plausibly recover it.

Usage (from the project root; run literature.py first):
    python -I sample/census.py            # query the archives, then build the census
    python -I sample/census.py --cached   # reuse the downloads in sample/raw/
"""

import argparse
import datetime
import hashlib
import io
import json
import os

import astropy.units as u
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(HERE, 'raw')
NEA_COMP_FILE = os.path.join(RAW_DIR, 'nea_pscomppars_rv.csv')
NEA_PS_FILE = os.path.join(RAW_DIR, 'nea_ps_rv.csv')
GAIA_FILE = os.path.join(RAW_DIR, 'gaia_dr3_hosts.csv')
GAIA_CONE_FILE = os.path.join(RAW_DIR, 'gaia_dr3_cone_no_id.csv')
GAIA_FLAME_FILE = os.path.join(RAW_DIR, 'gaia_dr3_flame_mass.csv')
## Error columns added 2026-11-06 (kept in separate files so the 2026-10-30 NEA snapshot
## stays as it was): eccentricity, RV semi-amplitude, and stellar-mass errors per NEA
## solution, the composite table's stellar-mass errors, and Gaia FLAME mass bounds
NEA_PS_ERR_FILE = os.path.join(RAW_DIR, 'nea_ps_errors.csv')
NEA_COMP_ERR_FILE = os.path.join(RAW_DIR, 'nea_pscomppars_st_mass_err.csv')
FLAME_BOUNDS_FILE = os.path.join(RAW_DIR, 'gaia_dr3_flame_bounds.csv')
STANDIN_FILE = os.path.join(HERE, 'st_mass_err_standin.csv')
LIT_FILE = os.path.join(HERE, 'raw_lit', 'literature_companions.csv')
CUT_FILE = os.path.join(HERE, 'census_cut_companions.csv')     # cut companions of census hosts
REVIEW_FILE = os.path.join(HERE, 'raw_lit', 'review_planet_regime.csv')

## Literature tiers (set in literature.py) whose sub-13 M_J entries need NEA vetting
PLANET_REGIME_MJ = 13.0
UNVETTED_TIERS = ('low-precision survey', 'compilation', 'mass-function M sin i')

## NEA solutions to skip (decided 2026-10-30): Feng et al. 2022 (ApJS 262, 21), whose
## RV+astrometry fits are considered suspect, and Gaia Collaboration (Arenou) et al. 2023
## (A&A 674, A34), whose Gaia DR3 orbits use no RVs
EXCLUDED_NEA_REFS = ('FENG_ET_AL__2022', 'GAIA_COLLABORATION_ET_AL__2023')

GAIA_COLUMNS = ('source_id, ra, dec, phot_g_mean_mag, parallax, parallax_error, ruwe, '
                'non_single_star')

NEA_TAP = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync'
GAIA_TAP = 'https://gea.esac.esa.int/tap-server/tap/sync'

## Sample cuts: upper bounds only. The mass bound (the hydrogen-burning limit) admits
## brown dwarfs; the NEA itself
## stops at 30 M_J, so heavier companions come from the literature catalogs.
MSINI_MAX_MJ = 80.0
SMA_MAX_AU = 5.0

## A published SMA that disagrees with Kepler's third law for its own P and M_* by more
## than this fraction is replaced by the Kepler value: the period is closer to the
## observable (decided 2026-11-04)
SMA_KEPLER_TOL = 0.2

## DR4 covers 66 months
T_BASELINE_YR = 5.5

## "Plausibly recoverable": alpha/sigma_FOV >= 1.5, the DR4 detection threshold of
## Lammers & Winn (2025, arXiv:2511.04673), and P within the DR4 time span. This is a
## rough estimate only. Lammers & Winn adopt P <= 4 yr instead: in their simulated DR4
## fits, the recovered period's mean error reaches 20% at P ~ 4 yr. That limit is kept
## for reference in the figures.
SNR_THRESHOLD = 1.5
P_MAX_YR = T_BASELINE_YR
P_ACCURATE_YR = 4.0

MJUP_MSUN = 9.545e-4

COMP_COLUMNS = ['pl_name', 'hostname', 'discoverymethod', 'sy_pnum', 'tran_flag',
                'pl_controv_flag', 'ra', 'dec', 'st_mass', 'sy_dist', 'sy_plx',
                'sy_gaiamag', 'sy_vmag', 'gaia_dr3_id', 'gaia_dr2_id']
PS_COLUMNS = ['pl_name', 'default_flag', 'pl_pubdate', 'pl_refname', 'soltype',
              'pl_orbper', 'pl_orbpererr1', 'pl_orbpererr2', 'pl_orbsmax',
              'pl_orbsmaxerr1', 'pl_orbsmaxerr2', 'pl_orbeccen', 'pl_msinij',
              'pl_msinijerr1', 'pl_msinijerr2', 'pl_bmassj', 'pl_bmassjerr1',
              'pl_bmassjerr2', 'pl_bmassprov', 'st_mass', 'pl_rvamp', 'pl_orblper',
              'pl_orbtper', 'pl_tsystemref']

PS_ERR_KEY = ['pl_name', 'pl_refname', 'pl_pubdate', 'soltype']
PS_ERR_COLUMNS = PS_ERR_KEY + ['pl_orbeccenerr1', 'pl_orbeccenerr2', 'pl_rvamperr1',
                               'pl_rvamperr2', 'st_masserr1', 'st_masserr2']

## Precision rule (decided 2026-10-27; missing errors pass, 2026-11-06): a companion is
## "precisely characterized" if its orbit (P, or SMA when P has no error) and its mass
## are measured at >= this many sigma. A missing or zero error never excludes a companion
## ("innocent until proven guilty"); such companions are flagged instead.
PRECISION_SIGMA = 3.0

## Stellar-mass bins (M_sun) for the stand-in fractional error given to stars without
## a published error
STANDIN_BINS = [0, 0.6, 0.9, 1.2, 1.6, 3.0, np.inf]
STANDIN_MIN_STARS = 10

## NEA time-reference labels for the epoch of periastron, mapped to the labels used for
## literature sources (see literature.py)
NEA_TSYSTEM = {'BJD-TDB': 'BJD_TDB', 'BJD': 'BJD', 'BJD-UTC': 'BJD_UTC', 'BJD-TCB': 'BJD_TCB',
               'HJD': 'HJD', 'JD': 'JD', 'JD-UTC': 'JD', 'TJD-TDB': 'TJD_TDB'}


def tap_query(url, adql):
    """Run a synchronous TAP query and return a DataFrame."""
    resp = requests.post(url, data={'REQUEST': 'doQuery', 'LANG': 'ADQL',
                                    'FORMAT': 'csv', 'QUERY': adql}, timeout=300)
    resp.raise_for_status()
    if resp.text.lstrip().startswith('<'):
        raise RuntimeError(f'TAP query failed:\n{resp.text[:2000]}')
    return pd.read_csv(io.StringIO(resp.text), dtype={'source_id': str})


def record_snapshot(raw_dir, sources, from_mtime=False, appended=False):
    """Write raw_dir/SNAPSHOT.json: for each downloaded file, where it came from, when it
    was retrieved (UTC), its row count, and its SHA-256. Archive queries return whatever
    the archive holds that day, so the snapshot is what makes a cached rebuild
    reproducible. With from_mtime, the retrieval time is the file's modification time.
    Files whose SHA-256 is unchanged keep their earlier record. With appended (rows added
    to a cached file), a changed file keeps its original retrieval time and gains an
    'appended_utc' time."""
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    path = os.path.join(raw_dir, 'SNAPSHOT.json')
    old = {}
    if os.path.exists(path):
        with open(path) as fh:
            old = json.load(fh)
    record = {}
    for name, origin in sources.items():
        file = os.path.join(raw_dir, name)
        with open(file, 'rb') as fh:
            data = fh.read()
        sha = hashlib.sha256(data).hexdigest()
        prev = old.get(name, {})
        if prev.get('sha256') == sha:
            record[name] = prev
            continue
        when = (datetime.datetime.fromtimestamp(os.path.getmtime(file), datetime.timezone.utc)
                .isoformat(timespec='seconds') if from_mtime else now)
        record[name] = {'source': origin, 'retrieved_utc': when,
                        'lines': data.count(b'\n'), 'sha256': sha}
        if appended and prev:
            record[name]['retrieved_utc'] = prev['retrieved_utc']
            record[name]['appended_utc'] = when
    with open(path, 'w') as fh:
        json.dump(record, fh, indent=2)


## What each cached download in raw/ is, for SNAPSHOT.json
RAW_SOURCES = {
    os.path.basename(NEA_COMP_FILE): f'{NEA_TAP}: pscomppars, rv_flag = 1',
    os.path.basename(NEA_PS_FILE): f'{NEA_TAP}: ps, rv_flag = 1',
    os.path.basename(GAIA_FILE): f'{GAIA_TAP}: gaiadr3.gaia_source, NEA and literature hosts',
    os.path.basename(GAIA_FLAME_FILE): f'{GAIA_TAP}: gaiadr3.astrophysical_parameters '
                                       '(mass_flame)',
    os.path.basename(GAIA_CONE_FILE): f'{GAIA_TAP}: gaiadr3.gaia_source cone searches '
                                      'around NEA hosts without a DR3 ID',
}


def query_gaia(ids):
    """Gaia DR3 astrometry/photometry and FLAME masses for a list of source IDs."""
    rows, flame = [], []
    for i in range(0, len(ids), 500):
        id_list = ', '.join(ids[i:i+500])
        rows.append(tap_query(GAIA_TAP, f'select {GAIA_COLUMNS} from gaiadr3.gaia_source '
                                        f'where source_id in ({id_list})'))
        flame.append(tap_query(GAIA_TAP, 'select source_id, mass_flame from '
                                         'gaiadr3.astrophysical_parameters '
                                         f'where source_id in ({id_list}) '
                                         'and mass_flame is not null'))
    return pd.concat(rows), pd.concat(flame)


ERR_SOURCES = {
    os.path.basename(NEA_PS_ERR_FILE): f'{NEA_TAP}: ps, rv_flag = 1 (e, K, M_* errors)',
    os.path.basename(NEA_COMP_ERR_FILE): f'{NEA_TAP}: pscomppars, rv_flag = 1 (M_* errors)',
    os.path.basename(FLAME_BOUNDS_FILE): f'{GAIA_TAP}: gaiadr3.astrophysical_parameters '
                                         '(mass_flame_lower, mass_flame_upper)',
}


def fetch_error_columns():
    """Download the error columns of ERR_SOURCES for any that the cache lacks. The NEA
    keys are matched to the cached solutions on planet, reference, date, and solution
    type; rows the archive added since the snapshot are ignored."""
    if not os.path.exists(NEA_PS_ERR_FILE):
        tap_query(NEA_TAP, f"select {', '.join(PS_ERR_COLUMNS)} from ps "
                           "where rv_flag = 1").to_csv(NEA_PS_ERR_FILE, index=False)
    if not os.path.exists(NEA_COMP_ERR_FILE):
        tap_query(NEA_TAP, 'select pl_name, st_masserr1, st_masserr2 from pscomppars '
                           'where rv_flag = 1').to_csv(NEA_COMP_ERR_FILE, index=False)
    if not os.path.exists(FLAME_BOUNDS_FILE):
        ids = pd.read_csv(GAIA_FLAME_FILE, dtype={'source_id': str}).source_id.tolist()
        chunks = [tap_query(GAIA_TAP, 'select source_id, mass_flame_lower, mass_flame_upper '
                                      'from gaiadr3.astrophysical_parameters where source_id '
                                      f"in ({', '.join(ids[i:i+200])})")
                  for i in range(0, len(ids), 200)]
        pd.concat(chunks).to_csv(FLAME_BOUNDS_FILE, index=False)
    record_snapshot(RAW_DIR, {**RAW_SOURCES, **ERR_SOURCES})


def add_missing_gaia():
    """For a cached rebuild: fetch Gaia rows only for literature hosts that the cached
    Gaia files lack (e.g., after adding a literature source), append them, and update
    the snapshot. Rows already cached are never re-queried."""
    if not os.path.exists(LIT_FILE):
        return
    lit_ids = pd.read_csv(LIT_FILE, dtype={'source_id': str}).source_id.dropna()
    cached = set(pd.read_csv(GAIA_FILE, dtype={'source_id': str}).source_id)
    missing = sorted(set(lit_ids) - cached)
    if not missing:
        return
    print(f'Fetching Gaia DR3 rows for {len(missing)} new literature hosts')
    rows, flame = query_gaia(missing)
    rows.to_csv(GAIA_FILE, mode='a', header=False, index=False)
    flame.to_csv(GAIA_FLAME_FILE, mode='a', header=False, index=False)
    record_snapshot(RAW_DIR, RAW_SOURCES, appended=True)


def download():
    os.makedirs(RAW_DIR, exist_ok=True)

    ## Every planet with RV data, and every published solution for those planets
    comp = tap_query(NEA_TAP, f"select {', '.join(COMP_COLUMNS)} from pscomppars "
                              "where rv_flag = 1")
    comp.to_csv(NEA_COMP_FILE, index=False)
    tap_query(NEA_TAP, f"select {', '.join(PS_COLUMNS)} from ps "
                       "where rv_flag = 1").to_csv(NEA_PS_FILE, index=False)

    ids = comp.gaia_dr3_id.dropna().str.replace('Gaia DR3 ', '', regex=False)
    if os.path.exists(LIT_FILE):
        lit_ids = pd.read_csv(LIT_FILE, dtype={'source_id': str}).source_id.dropna()
        ids = pd.concat([ids, lit_ids])
    ## Gaia astrometry and photometry, plus FLAME stellar masses, a fallback for hosts
    ## whose orbit source gives no mass
    rows, flame = query_gaia(list(ids.unique()))
    rows.to_csv(GAIA_FILE, index=False)
    flame.to_csv(GAIA_FLAME_FILE, index=False)

    ## Hosts without a DR3 ID: cone search around the archive position. Use a wide
    ## radius for nearby stars, whose proper motions can carry them far from J2000
    ## positions; a narrow one for distant stars in crowded fields.
    no_id = comp[comp.gaia_dr3_id.isna()].drop_duplicates('hostname')
    cones = []
    for _, row in no_id.iterrows():
        radius_deg = (120.0 if row.sy_dist < 100 else 5.0)/3600
        adql = (f'select {GAIA_COLUMNS} from gaiadr3.gaia_source where 1 = contains('
                f"point('ICRS', ra, dec), circle('ICRS', {row.ra}, {row.dec}, {radius_deg}))")
        cone = tap_query(GAIA_TAP, adql)
        cone['hostname'] = row.hostname
        cones.append(cone)
    pd.concat(cones).to_csv(GAIA_CONE_FILE, index=False)
    record_snapshot(RAW_DIR, RAW_SOURCES)
    for f in (NEA_PS_ERR_FILE, NEA_COMP_ERR_FILE, FLAME_BOUNDS_FILE):
        if os.path.exists(f):
            os.remove(f)
    fetch_error_columns()


def match_by_position(comp):
    """Pick the cone-search source that matches each ID-less host.

    First choice: parallax within 20% of the archive parallax (sy_plx, which is more
    reliable than sy_dist for some nearby stars) and G within 1.5 mag of the expected G
    (archive G, else V - 0.5). Fallback: a source with no parallax (a two-parameter
    solution) within 5 arcsec of the archive position and 0.75 mag of the expected G.
    Returns {hostname: (source_id, match label)}."""
    cones = pd.read_csv(GAIA_CONE_FILE, dtype={'source_id': str})
    hosts = comp.drop_duplicates('hostname').set_index('hostname')
    matches = {}
    for host, cands in cones.groupby('hostname'):
        row = hosts.loc[host]
        expected_G = row.sy_gaiamag if np.isfinite(row.sy_gaiamag) else row.sy_vmag - 0.5
        expected_plx = row.sy_plx if np.isfinite(row.sy_plx) else 1000/row.sy_dist
        dG = (cands.phot_g_mean_mag - expected_G).abs()
        ok = ((cands.parallax - expected_plx).abs() < 0.2*expected_plx) & (dG < 1.5)
        if ok.any():
            matches[host] = (cands.loc[dG[ok].idxmin(), 'source_id'], 'position')
            continue
        sep_arcsec = 3600*np.hypot((cands.ra - row.ra)*np.cos(np.radians(row.dec)),
                                   cands.dec - row.dec)
        ok = cands.parallax.isna() & (sep_arcsec < 5) & (dG < 0.75)
        if ok.any():
            matches[host] = (cands.loc[dG[ok].idxmin(), 'source_id'],
                             'position (no DR3 parallax)')
    return matches


def sigma_fov_lw(G):
    """Lammers & Winn (2025) per-FOV-transit error [uas]: flat for G < 14."""
    return 54.0*np.maximum(1.0, 10**(0.2*(G - 14.0)))


def sigma_fov_gaiamock(G):
    """Per-FOV-transit error [uas] from gaiamock (El-Badry et al. 2024): the per-CCD
    "EDR3 adjusted" curve of Holl et al. (2023, arXiv:2206.05439, Fig. 3) divided by
    sqrt(8), as in Lammers & Winn's comparison. Unlike the Lammers & Winn model, it rises
    for bright stars. Held constant beyond its G = 4-20 range."""
    G_vals = [4, 5, 6, 7, 8.2, 8.4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    sigma_ccd_mas = [0.4, 0.35, 0.15, 0.17, 0.23, 0.13, 0.13, 0.135, 0.125, 0.13, 0.15,
                     0.23, 0.36, 0.63, 1.05, 2.05, 4.1]
    return 1e3*np.interp(G, G_vals, sigma_ccd_mas)/np.sqrt(8)


def sym_err(e1, e2):
    """Symmetrized archive error: the mean absolute value of the two sides (one side
    alone if the other is missing). Zero or missing errors (fixed or derived values)
    become NaN."""
    err = pd.concat([e1.abs(), e2.abs()], axis=1).mean(axis=1)
    return err.where(err > 0)


def solution_candidates(ps, comp):
    """Every usable solution for each planet, labeled with its fit type.

    A solution that reports M sin i (pl_bmassprov = 'Msini') is RV-only; one that reports
    a true mass for a non-transiting planet used astrometry (e.g., Hipparcos-Gaia or Gaia
    DR3 with RVs). Solutions from EXCLUDED_NEA_REFS, and solutions with neither P nor SMA,
    are skipped. A transiting planet with no such solution falls back to its transit+RV
    solutions (true mass, sin i ~ 1).
    """
    transiting = set(comp.pl_name[comp.tran_flag == 1])
    ps = ps[~ps.pl_refname.str.contains('|'.join(EXCLUDED_NEA_REFS), na=False)]
    true_mass = (ps.pl_bmassprov == 'Mass') & ps.pl_bmassj.notna()
    has_orbit = ps.pl_orbper.notna() | ps.pl_orbsmax.notna()
    cand = ps.assign(fit_type=np.select(
        [(ps.pl_bmassprov == 'Msini') & ps.pl_msinij.notna() & has_orbit,
         true_mass & ~ps.pl_name.isin(transiting) & has_orbit,
         true_mass & ps.pl_name.isin(transiting) & has_orbit],
        ['RV-only', 'RV+astrometry', 'transit+RV'], ''))
    cand = cand[cand.fit_type != ''].copy()
    rv = cand.fit_type != 'transit+RV'
    cand = cand[rv | ~cand.pl_name.isin(cand.pl_name[rv])]
    cand['mass_mj'] = cand.pl_msinij.fillna(cand.pl_bmassj)
    cand['mass_type'] = np.where(cand.pl_msinij.notna(), 'M sin i', 'true mass')
    ## Fractional period error; missing or zero errors (fixed or derived periods) can't
    ## be compared
    err = (cand.pl_orbpererr1.abs() + cand.pl_orbpererr2.abs())/2
    cand['period_err_d'] = err.where(err > 0)
    cand['frac_p_err'] = cand.period_err_d/cand.pl_orbper
    ## Errors on the other quantities (symmetrized; zero or missing -> NaN)
    cand['mass_err_pub'] = np.where(cand.pl_msinij.notna(),
                                    sym_err(cand.pl_msinijerr1, cand.pl_msinijerr2),
                                    sym_err(cand.pl_bmassjerr1, cand.pl_bmassjerr2))
    cand['sma_err_pub'] = sym_err(cand.pl_orbsmaxerr1, cand.pl_orbsmaxerr2)
    cand['ecc_err'] = sym_err(cand.pl_orbeccenerr1, cand.pl_orbeccenerr2)
    cand['K_err'] = sym_err(cand.pl_rvamperr1, cand.pl_rvamperr2)
    cand['st_mass_err_sol'] = sym_err(cand.st_masserr1, cand.st_masserr2)
    cand['ref'] = cand.pl_refname.str.extract(r'refstr=(\S+)', expand=False)
    return cand.merge(comp[['pl_name', 'hostname']], on='pl_name')


def pick_solution(ps, comp):
    """One solution per planet (decided 2026-10-30; replaces the archive default):

    1. Per system, the newest paper with solutions for the most of its planets. Planets
       that paper doesn't cover take their own newest solution.
    2. If another solution (older or newer) has a period error less than half the chosen
       one's, take it instead (the newest such solution), and repeat. Solutions without
       a nonzero period error are never compared.
    3. Gaia-only solutions are excluded (EXCLUDED_NEA_REFS), so every orbit uses RVs.

    solution_rule records which step chose each solution.
    """
    cand = solution_candidates(ps, comp)
    ## Newest first; within a month, by reference name, so rebuilds are deterministic
    cand = cand.sort_values(['pl_pubdate', 'ref'], ascending=False)
    chosen = []
    for host, sys in cand.groupby('hostname', sort=True):
        papers = (sys.groupby('ref').agg(n=('pl_name', 'nunique'), date=('pl_pubdate', 'max'))
                  .reset_index().sort_values(['n', 'date', 'ref'], ascending=False))
        paper = papers.ref.iloc[0]
        for planet, sols in sys.groupby('pl_name', sort=True):
            if paper in set(sols.ref):
                sol, rule = sols[sols.ref == paper].iloc[0], 'newest paper for the system'
            else:
                sol, rule = sols.iloc[0], 'newest paper for the planet'
            step1 = sol
            while np.isfinite(sol.frac_p_err):
                better = sols[sols.frac_p_err < sol.frac_p_err/2]
                if better.empty:
                    break
                sol = better.iloc[0]
                rule = 'more precise P (>2x)'
            chosen.append(sol.to_frame().T.assign(
                solution_rule=rule, step1_ref=step1.ref, step1_period_d=step1.pl_orbper,
                step1_period_err_d=step1.period_err_d))
    chosen = pd.concat(chosen, ignore_index=True)
    chosen['mass_source'] = chosen.fit_type
    return chosen.infer_objects()


def add_literature(df, ps):
    """Append companions from the survey catalogs (literature.py), skipping any already
    present: same Gaia DR3 source and period within 10% of any published period for that
    companion. NEA entries take precedence, then catalogs in the order literature.py lists
    them. An NEA planet with no usable solution doesn't count as present, so a
    catalog's orbit stands in for it."""
    lit = pd.read_csv(LIT_FILE, dtype={'source_id': str})
    lit = lit.rename(columns={'st_mass': 'st_mass_sol'}).assign(
        discoverymethod='Radial Velocity', tran_flag=0, pl_controv_flag=0,
        mass_source=lit.catalog, pl_refname=lit.catalog, mass_type='M sin i',
        gaia_match=np.where(lit.source_id.notna(), 'SIMBAD', 'not in DR3'))

    periods = ps.dropna(subset=['pl_orbper']).groupby('pl_name').pl_orbper.apply(list)
    def published(mask):
        return [(s, p) for name, s in zip(df.pl_name[mask], df.source_id[mask])
                for p in periods.get(name, [])]
    def matches(s, p, pairs):
        return any(s == s0 and abs(p/p0 - 1) < 0.1 for s0, p0 in pairs)
    known = published(df.mass_mj.notna() & df.source_id.notna())
    nea_no_rv = published(df.mass_mj.isna() & df.source_id.notna())

    keep, planet_regime = [], []
    for s, p, m, tier in zip(lit.source_id, lit.pl_orbper, lit.mass_mj, lit.catalog_status):
        dup = matches(s, p, known)
        ## Below 13 M_J the NEA's vetting governs: compilations, low-precision candidate
        ## lists, and SB9 add a planet only as a stand-in for an NEA planet that lacks an
        ## usable solution. Their other sub-13 M_J entries are often unconfirmed or
        ## retracted planets.
        unvetted = (m < PLANET_REGIME_MJ and tier in UNVETTED_TIERS
                    and not matches(s, p, nea_no_rv))
        keep.append(not dup and not unvetted)
        planet_regime.append(not dup and unvetted)
        if not dup and not unvetted:
            known.append((s, p))
    lit[planet_regime].to_csv(REVIEW_FILE, index=False)
    print(f'Unvetted sub-{PLANET_REGIME_MJ:.0f} M_J literature entries set aside: '
          f'{sum(planet_regime)} (see {os.path.relpath(REVIEW_FILE, HERE)})')
    n_dup = len(keep) - sum(keep) - sum(planet_regime)
    lit = lit[keep].copy()
    lit['sy_pnum'] = lit.groupby('hostname').pl_name.transform('size')
    return pd.concat([df.assign(catalog='NEA'), lit], ignore_index=True), n_dup


def tp_to_bjd_tdb(tp, system, ra, dec):
    """Convert periastron times to BJD_TDB. Returns (bjd_tdb, flag) arrays.

    JD and HJD are taken as UTC; JD as geocentric. The barycentric (and heliocentric)
    light-travel times use the host's coordinates. 'BJD' with no stated time scale and
    'unknown' systems are taken as BJD_TDB; the error is at most ~70 s for a BJD_UTC
    value and ~0.006 d (the light-travel time across 1 AU) for a JD or HJD value. Values
    that look like reduced dates (< 2,000,000) with no stated offset are left blank."""
    tp = np.asarray(tp, float)
    system = np.asarray(system, object)
    bjd = np.full(len(tp), np.nan)
    flag = np.full(len(tp), 'no Tp', dtype=object)
    geocenter = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
    has = np.isfinite(tp)
    tjd = has & (system == 'TJD_TDB') & (tp < 1e5)
    bjd[tjd], flag[tjd] = tp[tjd] + 2457000, 'TESS JD (BJD_TDB - 2457000), offset added'
    reduced = has & ~tjd & (tp < 2e6)
    flag[reduced] = 'reduced date with unstated offset; not converted'
    ok = has & ~tjd & ~reduced
    as_is = {'BJD_TDB': 'BJD_TDB as published',
             'BJD': 'BJD with unstated time scale; taken as BJD_TDB',
             'unknown': 'time system not stated; taken as BJD_TDB'}
    for sys_label, note in as_is.items():
        m = ok & (system == sys_label)
        bjd[m], flag[m] = tp[m], note
    for sys_label, scale in [('BJD_UTC', 'utc'), ('BJD_TCB', 'tcb')]:
        m = ok & (system == sys_label)
        if m.any():
            bjd[m] = Time(tp[m], format='jd', scale=scale).tdb.jd
            flag[m] = f'converted from {sys_label}'
    coords_ok = np.isfinite(np.asarray(ra, float)) & np.isfinite(np.asarray(dec, float))
    for sys_label in ['JD', 'HJD']:
        m = ok & (system == sys_label)
        no_coord = m & ~coords_ok
        bjd[no_coord] = tp[no_coord]
        flag[no_coord] = f'{sys_label} not converted (no coordinates); taken as BJD_TDB'
        m &= coords_ok
        if not m.any():
            continue
        coord = SkyCoord(np.asarray(ra, float)[m]*u.deg, np.asarray(dec, float)[m]*u.deg)
        t = Time(tp[m], format='jd', scale='utc', location=geocenter)
        if sys_label == 'HJD':
            t = t - t.light_travel_time(coord, kind='heliocentric')
        bjd[m] = (t.tdb + t.light_travel_time(coord, kind='barycentric')).jd
        flag[m] = f'converted from {sys_label} (UTC assumed)'
    return bjd, flag


def add_errors(df):
    """Uncertainties on stellar mass, SMA, and mass, and the 3-sigma precision flags.

    Stellar mass: the error that goes with the mass actually used (the orbit source's, the
    NEA composite's, or Gaia FLAME's half-interval). Stars without one get a stand-in
    fractional error: the median among stars in the same mass bin that have a published
    error (STANDIN_FILE records the medians). Stand-in errors are labeled as such.

    SMA: the published error where the published SMA is kept; otherwise derived from the
    period error and the stellar-mass error (a^3 ~ M P^2, so
    sigma_a/a = sqrt((2 sigma_P/P)^2 + (sigma_M/M)^2)/3).

    Mass: the published error; else, for an M sin i with a published K error, derived from
    M sin i ~ K sqrt(1 - e^2) P^(1/3) M_*^(2/3), adding in quadrature the relative errors of
    K, P (1/3 power), M_* (2/3 power), and e (e sigma_e/(1 - e^2)), each used where known
    and treated as independent. Companions with neither get no mass error.

    Precision (PRECISION_SIGMA): orbit sigma is P/sigma_P, else the published SMA over its
    published error. A missing error counts as passing, and the companion is flagged
    precision_unverified. Returns the table and the stand-in table."""
    fb = (pd.read_csv(FLAME_BOUNDS_FILE, dtype={'source_id': str})
          .drop_duplicates('source_id'))
    fb['st_mass_err_flame'] = (fb.mass_flame_upper - fb.mass_flame_lower)/2
    df = df.merge(fb[['source_id', 'st_mass_err_flame']], on='source_id', how='left')

    ## Stellar mass
    src = df.st_mass_source
    published = pd.Series(np.select([src == 'orbit source', src == 'NEA composite',
                                     src == 'Gaia FLAME'],
                                    [df.st_mass_err_sol, df.st_mass_err_comp,
                                     df.st_mass_err_flame], np.nan), index=df.index)
    published = published.where(published > 0)
    stars = (df.assign(frac=published/df.st_mass).dropna(subset=['frac'])
             .drop_duplicates('hostname'))
    overall = stars.frac.median()
    idx = pd.cut(stars.st_mass, STANDIN_BINS, labels=False)
    standin = pd.DataFrame({
        'st_mass_lo': STANDIN_BINS[:-1], 'st_mass_hi': STANDIN_BINS[1:],
        'n_stars': [int((idx == i).sum()) for i in range(len(STANDIN_BINS) - 1)],
        'median_frac_err': [stars.frac[idx == i].median()
                            if (idx == i).sum() >= STANDIN_MIN_STARS else overall
                            for i in range(len(STANDIN_BINS) - 1)]})
    standin['overall_median_frac_err'] = overall
    row_bin = pd.cut(df.st_mass, STANDIN_BINS, labels=False)
    frac = row_bin.map(standin.median_frac_err)
    df['st_mass_err'] = published.fillna(frac*df.st_mass)
    df['st_mass_err_source'] = np.select(
        [published.notna(), df.st_mass_err.notna()],
        [np.where(src == 'orbit source', 'published (orbit source)',
                  np.where(src == 'NEA composite', 'published (NEA composite)',
                           'Gaia FLAME bounds')), 'stand-in (median fractional error)'],
        'none')

    ## Period error: published, else from the published SMA error (sigma_P/P =
    ## 1.5 sigma_a/a at fixed mass). period_err_d itself stays as published, since the
    ## 2x period-precision rule that picks the NEA solutions uses it.
    P, M, sM = df.period_d, df.st_mass, df.st_mass_err
    pub_sma = (df.sma_source == 'published') & df.sma_err_pub.notna()
    p_from_sma = 1.5*P*df.sma_err_pub/df.pl_orbsmax
    sP = df.period_err_d.where(df.period_err_d.notna(), p_from_sma.where(pub_sma))
    df['period_err_filled_d'] = sP
    df['period_err_source'] = np.select(
        [df.period_err_d.notna(), sP.notna()], ['published', 'derived from SMA error'], 'none')

    ## SMA
    derived_sma = df.sma_au*np.sqrt((2*sP/P)**2 + (sM/(M + df.mass_mj*MJUP_MSUN))**2)/3
    df['sma_err_au'] = np.where(pub_sma, df.sma_err_pub, derived_sma)
    df['sma_err_source'] = np.select(
        [pub_sma, np.isfinite(derived_sma)], ['published', 'derived from P and M_*'], 'none')

    ## Mass
    K = df.K_ms
    rel2 = ((df.K_err/K)**2 + (2/3*sM/M)**2 + ((sP/P)/3).fillna(0)**2
            + (df.pl_orbeccen*df.ecc_err/(1 - df.pl_orbeccen**2)).fillna(0)**2)
    ## An M sin i, or a transiting planet's mass (sin i ~ 1); not an RV+astrometry true
    ## mass, whose inclination uncertainty K cannot supply
    can_derive = (df.mass_err_pub.isna()
                  & ((df.mass_type == 'M sin i') | (df.fit_type == 'transit+RV'))
                  & (K > 0) & (df.K_err > 0) & sM.notna())
    derived_mass = df.mass_mj*np.sqrt(rel2)
    df['mass_err_mj'] = np.where(df.mass_err_pub.notna(), df.mass_err_pub,
                                 np.where(can_derive, derived_mass, np.nan))
    df['mass_err_source'] = np.select(
        [df.mass_err_pub.notna(), can_derive],
        ['published', 'derived from K (with P, e, M_* terms)'], 'none')
    df['msini_err_from_K'] = df.mass_err_source.str.startswith('derived')

    ## Precision
    df['p_sigma'] = P/sP
    df['orbit_sigma'] = df.p_sigma.where(df.p_sigma.notna(), df.pl_orbsmax/df.sma_err_pub)
    df['mass_sigma'] = df.mass_mj/df.mass_err_mj
    ## Errors that no route (published, derived from a related quantity, literature) could
    ## supply. These are concerning: the paper gave neither the quantity's error nor a
    ## related one's (K, M sin i, P, SMA).
    df['missing_errors'] = np.select(
        [df.orbit_sigma.isna() & df.mass_sigma.isna(), df.orbit_sigma.isna(),
         df.mass_sigma.isna()], ['orbit and mass', 'orbit', 'mass'], '')
    df['orbit_precise'] = df.orbit_sigma.isna() | (df.orbit_sigma >= PRECISION_SIGMA)
    df['mass_precise'] = df.mass_sigma.isna() | (df.mass_sigma >= PRECISION_SIGMA)
    df['precise'] = df.orbit_precise & df.mass_precise
    df['precision_unverified'] = df.precise & (df.orbit_sigma.isna() | df.mass_sigma.isna())
    return df, standin


def build_census():
    comp = pd.read_csv(NEA_COMP_FILE)
    ps = pd.read_csv(NEA_PS_FILE)
    gaia = pd.read_csv(GAIA_FILE, dtype={'source_id': str})
    ## Error columns (fetch_error_columns): matched to the cached solutions by key
    n_ps = len(ps)
    ps = ps.merge(pd.read_csv(NEA_PS_ERR_FILE).drop_duplicates(PS_ERR_KEY),
                  on=PS_ERR_KEY, how='left')
    assert len(ps) == n_ps
    cerr = pd.read_csv(NEA_COMP_ERR_FILE).drop_duplicates('pl_name')
    comp = comp.merge(cerr.assign(st_mass_err_comp=sym_err(cerr.st_masserr1, cerr.st_masserr2))
                      [['pl_name', 'st_mass_err_comp']], on='pl_name', how='left')

    sol = pick_solution(ps, comp)
    sol = sol.rename(columns={'st_mass': 'st_mass_sol', 'pl_rvamp': 'K_ms',
                              'pl_orblper': 'omega_deg', 'pl_orbtper': 'tp_raw'})
    sol['tp_system'] = np.where(sol.tp_raw.isna(), 'none',
                                sol.pl_tsystemref.map(NEA_TSYSTEM).fillna('unknown'))
    sol = sol[['pl_name', 'mass_mj', 'mass_type', 'fit_type', 'mass_source',
               'solution_rule', 'step1_ref', 'step1_period_d', 'step1_period_err_d',
               'pl_refname', 'pl_pubdate', 'pl_orbper', 'period_err_d',
               'pl_orbsmax', 'pl_orbeccen', 'K_ms', 'omega_deg', 'tp_raw', 'tp_system',
               'st_mass_sol', 'mass_err_pub', 'sma_err_pub', 'ecc_err', 'K_err',
               'st_mass_err_sol']]
    df = comp.merge(sol, on='pl_name', how='left')
    df['pl_refname'] = df.pl_refname.str.extract(r'refstr=(\S+)', expand=False)

    ## Gaia DR3: keep source IDs as strings, since they exceed float precision. Hosts
    ## without an archive ID are matched by position.
    df['source_id'] = df.gaia_dr3_id.str.replace('Gaia DR3 ', '', regex=False).str.strip()
    df['gaia_match'] = np.where(df.gaia_dr3_id.notna(), 'archive ID', 'not in DR3')
    for host, (source_id, label) in match_by_position(comp).items():
        df.loc[df.hostname == host, ['source_id', 'gaia_match']] = [source_id, label]
    df, n_lit_dup = add_literature(df, ps)
    cone = pd.read_csv(GAIA_CONE_FILE, dtype={'source_id': str})
    gaia = pd.concat([gaia, cone.drop(columns='hostname')]).drop_duplicates('source_id')
    df = df.merge(gaia.rename(columns={'ra': 'gaia_ra', 'dec': 'gaia_dec'}),
                  on='source_id', how='left')

    ## Flag matches whose Gaia photometry or parallax disagrees with the archive: likely
    ## the wrong source or a blended binary. For these, fall back to archive values.
    suspect = (df.phot_g_mean_mag - df.sy_gaiamag).abs() > 0.5
    suspect |= (df.phot_g_mean_mag - df.sy_vmag) > 0.5
    suspect |= (df.phot_g_mean_mag.notna() & ~(df.parallax > 0)
                & (df.gaia_match != 'position (no DR3 parallax)'))
    df.loc[suspect, 'gaia_match'] = df.loc[suspect, 'gaia_match'] + ' (suspect)'
    ## Without usable Gaia photometry, fall back to the archive G, then to V (a rough
    ## stand-in, mostly for very bright stars absent from DR3)
    use_gaia = df.phot_g_mean_mag.notna() & ~suspect
    df['G'] = np.where(use_gaia, df.phot_g_mean_mag, df.sy_gaiamag.fillna(df.sy_vmag))
    df['G_source'] = np.where(use_gaia, 'Gaia DR3',
                              np.where(df.sy_gaiamag.notna(), 'archive G',
                                       np.where(df.sy_vmag.notna(), 'archive V', 'none')))
    ## Distance: Gaia parallax, else the archive parallax, else the archive distance
    archive_dist = (1000/df.sy_plx).where(df.sy_plx > 0, df.sy_dist)
    df['dist_pc'] = np.where(use_gaia & (df.parallax > 0), 1000/df.parallax, archive_dist)
    df['dist_source'] = np.where(use_gaia & (df.parallax > 0), 'Gaia DR3 parallax',
                                 np.where(df.sy_plx > 0, 'archive parallax',
                                          np.where(df.sy_dist.notna(), 'archive distance',
                                                   'none')))
    df['abs_G'] = df.G - 5*np.log10(df.dist_pc/10)

    ## Stellar mass from the chosen solution, else the composite table, else Gaia DR3
    ## FLAME.
    flame = pd.read_csv(GAIA_FLAME_FILE, dtype={'source_id': str})
    df = df.merge(flame, on='source_id', how='left')
    df['st_mass_source'] = np.where(df.st_mass_sol.notna(), 'orbit source',
                                    np.where(df.st_mass.notna(), 'NEA composite',
                                             np.where(df.mass_flame.notna(), 'Gaia FLAME',
                                                      'none')))
    df['st_mass'] = df.st_mass_sol.fillna(df.st_mass).fillna(df.mass_flame)

    ## P and SMA: Kepler's third law with the total mass (M_* plus M sin i, or the true
    ## mass) fills whichever the solution lacks, and replaces a published SMA that
    ## disagrees with the period by more than SMA_KEPLER_TOL
    m_tot = df.st_mass + df.mass_mj*MJUP_MSUN
    df['period_d'] = df.pl_orbper.fillna(365.25*np.sqrt(df.pl_orbsmax**3/m_tot))
    kepler_sma = (m_tot*(df.pl_orbper/365.25)**2)**(1/3)
    df['sma_source'] = np.select(
        [df.pl_orbsmax.isna() & kepler_sma.isna(), df.pl_orbsmax.isna(),
         (df.pl_orbsmax/kepler_sma - 1).abs() > SMA_KEPLER_TOL],
        ['none', 'Kepler (no SMA given)',
         f'Kepler (published SMA off by > {SMA_KEPLER_TOL:.0%})'], 'published')
    df['sma_au'] = np.where(df.sma_source == 'published', df.pl_orbsmax, kepler_sma)
    df['p_over_T'] = df.period_d/365.25/T_BASELINE_YR
    df['period_err_d'] = df.period_err_d.where(df.period_err_d > 0)
    df, standin = add_errors(df)
    standin.to_csv(STANDIN_FILE, index=False)

    ## Periastron times in a common system, BJD_TDB
    df['Tp_bjd_tdb'], df['tp_flag'] = tp_to_bjd_tdb(
        df.tp_raw, df.tp_system, df.gaia_ra.fillna(df.ra), df.gaia_dec.fillna(df.dec))

    ## Minimum astrometric signature of the star's orbit about the barycentre [uas]
    q = df.mass_mj*MJUP_MSUN/df.st_mass
    df['alpha_min_uas'] = 1e6*q/(1 + q)*df.sma_au/df.dist_pc
    df['sigma_fov_uas'] = sigma_fov_gaiamock(df.G)
    df['snr'] = df.alpha_min_uas/df.sigma_fov_uas
    df['snr_lw'] = df.alpha_min_uas/sigma_fov_lw(df.G)
    ## Hosts absent from DR3 (mostly G < 3 giants) are assumed absent from DR4 too
    df['recoverable'] = ((df.snr >= SNR_THRESHOLD) & (df.period_d <= P_MAX_YR*365.25)
                         & (df.gaia_match != 'not in DR3'))

    ## Cut flow
    steps = [
        ('NEA rv_flag = 1, plus new literature companions', np.ones(len(df), bool)),
        ('M sin i or mass from a usable solution', df.mass_mj.notna()),
        (f'M sin i (or true mass) <= {MSINI_MAX_MJ} M_J', df.mass_mj <= MSINI_MAX_MJ),
        (f'SMA <= {SMA_MAX_AU} AU', df.sma_au <= SMA_MAX_AU),
        ('stellar mass and distance known', df.st_mass.notna() & df.dist_pc.notna()),
    ]
    keep = np.ones(len(df), bool)
    cutflow = []
    cut_reason = pd.Series('', index=df.index)
    for label, mask in steps:
        cut_reason[keep & ~np.asarray(mask)] = label
        keep &= np.asarray(mask)
        cutflow.append((label, int(keep.sum()), df[keep].hostname.nunique()))
    sample = df[keep].copy()

    out_cols = ['pl_name', 'hostname', 'catalog', 'catalog_status', 'discoverymethod',
                'sy_pnum', 'tran_flag',
                'pl_controv_flag', 'mass_mj', 'mass_type', 'fit_type', 'mass_source',
                'solution_rule', 'pl_refname', 'pl_pubdate', 'sma_au', 'sma_source',
                'period_d', 'period_err_d', 'pl_orbeccen', 'K_ms', 'omega_deg', 'tp_raw', 'tp_system',
                'Tp_bjd_tdb', 'tp_flag', 'p_over_T', 'st_mass', 'st_mass_source', 'dist_pc',
                'dist_source', 'G', 'abs_G',
                'G_source', 'gaia_match', 'source_id', 'ruwe', 'non_single_star', 'sy_gaiamag',
                'sy_vmag', 'parallax', 'sy_plx', 'alpha_min_uas', 'sigma_fov_uas', 'snr',
                'snr_lw', 'recoverable',
                ## Added 2026-11-06: errors, how each was obtained, and precision flags
                'period_err_filled_d', 'period_err_source', 'missing_errors',
                'st_mass_err', 'st_mass_err_source', 'sma_err_au', 'sma_err_source',
                'mass_err_mj', 'mass_err_source', 'msini_err_from_K', 'p_sigma',
                'orbit_sigma', 'mass_sigma', 'orbit_precise', 'mass_precise', 'precise',
                'precision_unverified']
    sample = sample[out_cols].sort_values('snr', ascending=False)
    ## Full precision, since periods and periastron times need more than a few digits
    sample.to_csv(os.path.join(HERE, 'census_table.csv'), index=False)
    ## Companions that fail a cut but orbit a census host (Judah, 2026-11-20): not part of
    ## the sample, but the assisted RV fit of a system puts them back in (no period or SMA
    ## limit), so its baseline RV model is complete. cut_reason = the first cut failed.
    cut = df[~keep & df.hostname.isin(set(sample.hostname))].assign(cut_reason=cut_reason[~keep])
    cut[out_cols + ['cut_reason']].sort_values(['hostname', 'period_d']).to_csv(
        CUT_FILE, index=False)
    return df, sample, cutflow, n_lit_dup


def tables(sample):
    g_bins = [-5, 6, 8, 10, 12, 14, 30]
    g_labels = ['<6', '6-8', '8-10', '10-12', '12-14', '>14']
    p_bins = [0, 0.1, 0.5, 1, 2, np.inf]
    p_labels = ['<0.1', '0.1-0.5', '0.5-1', '1-2', '>2']
    s_bins = [0, 0.1, 0.5, 1.5, 5, np.inf]
    s_labels = ['<0.1', '0.1-0.5', '0.5-1.5', '1.5-5', '>5']

    s = sample.copy()
    s['G bin'] = pd.cut(s.G, g_bins, labels=g_labels, right=False)
    s['P/T bin'] = pd.cut(s.p_over_T, p_bins, labels=p_labels, right=False)
    s['S/N bin'] = pd.cut(s.snr, s_bins, labels=s_labels, right=False)
    return {
        'G vs P/T': pd.crosstab(s['G bin'], s['P/T bin'], margins=True,
                                margins_name='total'),
        'S/N vs P/T': pd.crosstab(s['S/N bin'], s['P/T bin'], margins=True,
                                  margins_name='total'),
        'Recoverable (S/N >= 1.5, P <= 5.5 yr) by G': pd.crosstab(
            s['G bin'], s.recoverable, margins=True, margins_name='total'),
        'Gaia match (hosts)': s.drop_duplicates('hostname').gaia_match.value_counts(),
        'Discovery method': s.discoverymethod.value_counts(),
        'Mass source': s.mass_source.value_counts(),
        'NEA solution rule': s.solution_rule.value_counts(),
        'SMA source': s.sma_source.value_counts(),
        'Mass (M sin i, M_J) by catalog': pd.crosstab(
            s.catalog, pd.cut(s.mass_mj, [0, 13, 20, 30, 50, 80]), margins=True,
            margins_name='total'),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cached', action='store_true',
                        help='reuse the downloads in sample/raw/')
    args = parser.parse_args()

    if args.cached:
        fetch_error_columns()
        add_missing_gaia()
    else:
        download()
    df, sample, cutflow, n_lit_dup = build_census()

    print(f'Literature companions already in the NEA (dropped): {n_lit_dup}')
    print('Cut flow (companions, hosts):')
    for label, n_pl, n_host in cutflow:
        print(f'  {label:48s} {n_pl:5d} {n_host:5d}')
    for name, table in tables(sample).items():
        print(f'\n{name}\n{table}')

    ## Pathological solutions, noted but not corrected (QUESTIONS.md, "To remember")
    s = df[df.pl_name.isin(sample.pl_name)]
    cols = ['pl_name', 'pl_refname', 'period_d', 'sma_au', 'st_mass', 'snr']
    off = s[s.sma_source.str.contains('off')].assign(published_sma=s.pl_orbsmax)
    print(f'\nPublished SMA replaced by the Kepler value ({len(off)}):')
    print(off[cols + ['published_sma']].sort_values('snr', ascending=False)
          .to_string(index=False))
    sw = s[s.solution_rule == 'more precise P (>2x)']
    z = ((sw.period_d - sw.step1_period_d).abs()
         / np.hypot(sw.period_err_d, sw.step1_period_err_d))
    tension = sw[z > 3].assign(sigma=z[z > 3].round(1))
    print(f'\n2x period-precision switches whose period disagrees with the step-1 '
          f'(newest) solution at > 3 sigma ({len(tension)} of {len(sw)}):')
    print(tension[['pl_name', 'pl_refname', 'period_d', 'period_err_d', 'step1_ref',
                   'step1_period_d', 'step1_period_err_d', 'sigma', 'snr']]
          .sort_values('snr', ascending=False).to_string(index=False))
