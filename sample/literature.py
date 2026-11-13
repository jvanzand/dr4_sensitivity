"""RV companions from survey catalogs that the NASA Exoplanet Archive doesn't cover.

The NEA lists objects up to 30 M_J, so brown dwarfs come from survey papers instead. Two
RV-only orbit catalogs are read straight from VizieR:

    CLS I (Rosenthal et al. 2021, ApJS 255, 8; J/ApJS/255/8): every signal in the
        California Legacy Survey, with RV-only fits. Signals CLS classifies as known
        planets (K), candidates (C), or stellar companions (S) are kept; activity (A),
        rotation (R), and artifacts (N) are dropped. "S" includes brown dwarfs.
    Kiefer et al. 2019 (A&A 631, A125; J/A+A/631/A125), Table B5: SOPHIE orbits of 11
        companions with M sin i = 20-90 M_J.

Other sources were transcribed from papers (mostly arXiv LaTeX sources), VizieR, or
compiled catalogs into raw_lit/<short_id>/<short_id>_orbits.csv; see SOURCE_TIERS below
and the notes in sample/lit_search/.

Hosts are resolved to Gaia DR3 source IDs through SIMBAD. Writes
sample/raw_lit/literature_companions.csv.

Usage (from the project root):
    python -I -W ignore sample/literature.py            # download, then build
    python -I -W ignore sample/literature.py --cached   # reuse sample/raw_lit/
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import requests
from astropy.io import ascii

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census import GAIA_TAP, record_snapshot, tap_query

RAW_DIR = os.path.join(HERE, 'raw_lit')
OUT_FILE = os.path.join(RAW_DIR, 'literature_companions.csv')
SIMBAD_FILE = os.path.join(RAW_DIR, 'simbad_gaia_ids.csv')

## SIMBAD at CDS, else its Harvard mirror (CDS was unreachable on 2026-10-30)
SIMBAD_HOSTS = ['https://simbad.cds.unistra.fr', 'https://simbad.harvard.edu']
_simbad_host = []
CDS_FTP = 'https://cdsarc.cds.unistra.fr/ftp'

CATALOGS = {
    'cls_2021': ('J/ApJS/255/8', ['ReadMe', 'planets.dat', 'table2.dat']),
    'kiefer_2019': ('J/A+A/631/A125', ['ReadMe', 'tableb1.dat', 'tableb3.dat',
                                        'tableb5.dat']),
}

## Columns of CLS planets.dat, which is '|'-separated (the ReadMe's byte-by-byte header
## for this file is malformed, so astropy's CDS reader can't parse it)
CLS_PLANET_COLUMNS = {1: 'CPS', 2: 'Status', 3: 'Npl', 4: 'Msini', 8: 'a', 12: 'Per',
                      20: 'omega', 24: 'K', 28: 'e'}

## Time system of each source's periastron times (Tp), as the papers state them, for
## sources whose table predates the tp_system column. Values: BJD_TDB, BJD (no time scale
## given), BJD_UTC, HJD, JD, or unknown. census.py converts every Tp to BJD_TDB.
## marvels_bd mixes systems, so its rows are read from their notes.
TP_SYSTEM_BY_SOURCE = {
    'barbato_2023': 'BJD', 'bouchy_2016': 'unknown', 'dalal_2021': 'BJD',
    'diaz_2012': 'BJD', 'exoplanets_org': 'unknown', 'express_jones': 'JD',
    'feng_2022': 'BJD', 'frensch_2023': 'BJD', 'grieves_2017': 'JD', 'jenkins_2009': 'JD',
    'jenkins_2015': 'JD', 'kiefer_2021': 'unknown', 'lick_giants': 'JD',
    'okayama_eapsnet': 'JD', 'patel_2007': 'JD', 'ptps': 'JD', 'rickman_2019': 'JD',
    'sahlmann_2011': 'JD', 'sb9': 'unknown', 'teng_2023': 'JD', 'troup_2016': 'unknown',
    'unger_2023': 'BJD', 'wilson_2016': 'BJD', 'xiao_2023': 'unknown',
}


def tp_system_from_notes(note):
    """Time system named in a row's notes, for sources that mix systems."""
    for pattern, system in [(r'BJD[_ ]?TDB', 'BJD_TDB'), (r'HJD', 'HJD'), (r'BJD', 'BJD'),
                            (r'(?<![A-Z])JD', 'JD')]:
        if re.search(pattern, note):
            return system
    return 'unknown'


def download():
    for subdir, (catalog, files) in CATALOGS.items():
        os.makedirs(os.path.join(RAW_DIR, subdir), exist_ok=True)
        for f in files:
            resp = requests.get(f'{CDS_FTP}/{catalog}/{f}', timeout=300)
            resp.raise_for_status()
            with open(os.path.join(RAW_DIR, subdir, f), 'wb') as fh:
                fh.write(resp.content)
        record_snapshot(os.path.join(RAW_DIR, subdir),
                        {f: f'{CDS_FTP}/{catalog}/{f}' for f in files})


def simbad_name(raw):
    """Turn a catalog identifier (e.g., '11964a', 'hip63510', 'gl876', 'bd-103166',
    'BD+291539', 'HD217850') into a SIMBAD-resolvable name."""
    raw = raw.strip()
    m = re.fullmatch(r'(\d+)([a-zA-Z]?)', raw)
    if m:
        return f'HD {m.group(1)}{m.group(2).upper()}'
    m = re.fullmatch(r'(?i)(hd|hip|gl|gj)\s*(\d+)([a-zA-Z]?)', raw)
    if m:
        prefix = {'gl': 'GJ'}.get(m.group(1).lower(), m.group(1).upper())
        return f'{prefix} {m.group(2)}{m.group(3).upper()}'
    m = re.fullmatch(r'(?i)bd([+-]\d{2})(\d+)', raw)
    if m:
        return f'BD{m.group(1)} {int(m.group(2))}'
    return raw


def load_cls():
    d = os.path.join(RAW_DIR, 'cls_2021')
    pl = pd.read_csv(os.path.join(d, 'planets.dat'), sep='|', header=None)
    pl = pl[list(CLS_PLANET_COLUMNS)].rename(columns=CLS_PLANET_COLUMNS)
    pl['CPS'] = pl.CPS.str.strip()
    pl['Status'] = pl.Status.str.strip()
    pl = pl[pl.Status.isin(['K', 'C', 'S'])]

    stars = ascii.read(os.path.join(d, 'table2.dat'), readme=os.path.join(d, 'ReadMe'),
                       format='cds').to_pandas()
    stars['CPS'] = stars.CPS.str.strip()
    ## Take only stellar mass from CLS. Its 'Vmag' column holds what look like absolute
    ## magnitudes (HD 107148: 4.47, vs. V = 8.0), and at least one parallax is wrong
    ## (HD 18445: 6.98 arcsec), so photometry and parallax come from Gaia instead.
    pl = pl.merge(stars[['CPS', 'Mass']], on='CPS', how='left')

    status = {'K': 'known planet', 'C': 'candidate', 'S': 'stellar/substellar companion'}
    host = pl.CPS.map(simbad_name)
    return pd.DataFrame({
        'pl_name': host + ' [CLS ' + pl.Npl.astype(int).astype(str) + ']',
        'hostname': host,
        'catalog': 'CLS I (Rosenthal+2021)',
        'catalog_status': pl.Status.map(status),
        'mass_mj': pl.Msini,
        'pl_orbper': pl.Per,
        'pl_orbsmax': pl.a,
        'pl_orbeccen': pl.e,
        'K_ms': pl.K,
        'omega_deg': np.degrees(pl.omega),
        'tp_raw': np.nan,
        'tp_system': 'none',
        'st_mass': pl.Mass,
        'sy_vmag': np.nan,
        'sy_plx': np.nan,
    })


def load_kiefer():
    d = os.path.join(RAW_DIR, 'kiefer_2019')
    read = lambda f: ascii.read(os.path.join(d, f), readme=os.path.join(d, 'ReadMe'),
                                format='cds').to_pandas()
    orbits = read('tableb5.dat').merge(read('tableb1.dat')[['Name', 'Vmag', 'Plx']],
                                       on='Name').merge(read('tableb3.dat')[['Name', 'Mass']],
                                                        on='Name')
    host = orbits.Name.map(simbad_name)
    return pd.DataFrame({
        'pl_name': host + ' [Kiefer+2019]',
        'hostname': host,
        'catalog': 'Kiefer+2019',
        'catalog_status': 'substellar companion',
        'mass_mj': orbits.M2sini,
        'pl_orbper': orbits.Per,
        'pl_orbsmax': orbits.a2,
        'pl_orbeccen': orbits.e,
        'K_ms': orbits.K,
        'omega_deg': orbits.omega,
        ## Tp is listed as Tp-2400000; rows flagged '*' give a transit time instead.
        ## SOPHIE time stamps in this catalog are BJD.
        'tp_raw': np.where(orbits.n_Tp.fillna('').str.strip() == '*', np.nan,
                           orbits.Tp + 2400000),
        'tp_system': 'BJD',
        'st_mass': orbits.Mass,
        'sy_vmag': orbits.Vmag,
        'sy_plx': orbits.Plx,
    })


## Sources in the standard format, in order of precedence when a companion appears in
## several (census.add_literature keeps the first). Survey papers with full RV orbits come
## first, then few-epoch, low-precision survey candidates, then compilations of published
## orbits (those listing full elements first), then single-lined binary catalogs whose
## M sin i is derived from the mass function with an assumed primary mass. Papers with
## only RV+astrometry fits (Xiao+2023, Ren+2026, Cui+2026) rank just below the survey
## papers. Notes on each source are in sample/lit_search/.
SOURCE_TIERS = {
    'survey paper': ['sahlmann_2011', 'unger_2023', 'barbato_2023', 'segransan_2010',
                     'rickman_2019', 'jenkins_2009', 'jenkins_2015', 'frensch_2023',
                     'parc_2026', 'wilson_2016', 'dalal_2021', 'diaz_2012', 'bouchy_2016',
                     'patel_2007', 'lick_giants', 'teng_2023', 'okayama_eapsnet', 'ptps',
                     'express_jones', 'barbato_2026', 'nidever_2002', 'wittenmyer_2009',
                     'lee_2023'],
    'RV+astrometry paper': ['xiao_2023', 'ren_2026', 'cui_2026'],
    'low-precision survey': ['grieves_2017', 'marvels_bd', 'troup_2016'],
    'compilation': ['kiefer_2021', 'exoplanets_org', 'ma_ge_2014', 'stevenson_2023',
                    'gaia_dr3_val', 'exoplanet_eu'],
    'mass-function M sin i': ['latham_2002', 'sb9', 'sb9_tic'],
}

## Sources left out entirely: Feng+2022's RV+astrometry fits are considered suspect
## (decided 2026-10-30)
EXCLUDED_SOURCES = ['feng_2022']

## Rows to drop, by source and host: hot-subdwarf primaries, whose assumed masses make the
## derived M sin i unreliable
EXCLUDED_HOSTS = {'sb9_tic': ['BPS CS 22169-0001', 'CPD-64 481']}

## exoplanet.eu includes companions of white dwarfs, cataclysmic variables, pulsars, and
## hot subdwarfs; drop hosts with these SIMBAD object types
EU_EXCLUDED_OTYPES = ('CV*', 'WD*', 'WD?', 'Psr', 'HS*')


def optional_error(t, column, notes=None, pattern=None):
    """An error column of an orbit table as numbers (NaN where blank or <= 0), falling
    back to a regex on the notes column when the table has no such column or a row is
    blank there."""
    out = (pd.to_numeric(t[column], errors='coerce') if column in t
           else pd.Series(np.nan, index=t.index))
    if pattern is not None and notes is not None:
        out = out.fillna(pd.to_numeric(notes.str.extract(pattern, expand=False), errors='coerce'))
    return out.where(out > 0)


def load_standard():
    """Orbit tables transcribed from papers into the standard format, one per source:
    raw_lit/<short_id>/<short_id>_orbits.csv with columns host, companion, P_d, e_P_d,
    K_ms, e, omega_deg, Tp (with tp_system; older tables have Tp_bjd instead), msini_mj,
    a_au, mstar_msun, fit_type, table_ref, notes.
    RV-only and RV+astrometry fits with an M sin i and a period are used; for the latter,
    msini_mj is the published M sin i or the true mass times sin i."""
    order = [s for tier in SOURCE_TIERS.values() for s in tier]
    tier_of = {s: tier for tier, sources in SOURCE_TIERS.items() for s in sources}
    found = sorted(os.path.basename(os.path.dirname(p)) for p in
                   glob.glob(os.path.join(RAW_DIR, '*', '*_orbits.csv'))
                   if os.path.basename(p) == os.path.basename(os.path.dirname(p)) + '_orbits.csv')
    found = [s for s in found if s not in EXCLUDED_SOURCES]
    unknown = [s for s in found if s not in tier_of]
    if unknown:
        print('Sources not in SOURCE_TIERS (added last):', unknown)

    frames = []
    for short_id in [s for s in order if s in found] + unknown:
        t = pd.read_csv(os.path.join(RAW_DIR, short_id, f'{short_id}_orbits.csv'))
        t = t[t.fit_type.str.strip().isin(['RV-only', 'RV+astrometry']) & t.msini_mj.notna()
              & t.P_d.notna()
              & ~t.host.str.strip().isin(EXCLUDED_HOSTS.get(short_id, []))]
        notes = t.notes.fillna('')
        mass = t.msini_mj.copy()
        if short_id == 'exoplanet_eu':
            otype = notes.str.extract(r'otype ([^\s;)]+)', expand=False)
            t, notes, mass, otype = [x[~otype.isin(EU_EXCLUDED_OTYPES)]
                                     for x in (t, notes, mass, otype)]
            ## Its mass column is the true mass when an inclination is listed; convert
            ## back to M sin i so all sources share the same quantity
            incl = notes.str.extract(r'inclination ([\d.]+) deg listed',
                                     expand=False).astype(float)
            mass = mass.where(incl.isna(), mass*np.sin(np.radians(incl)))
        ## Periastron times: newer tables have Tp and tp_system columns; older ones have
        ## Tp_bjd, whose system is set per source above
        tp = t['Tp'] if 'Tp' in t else t['Tp_bjd']
        if 'tp_system' in t:
            tp_system = t.tp_system.fillna('unknown')
        elif short_id in TP_SYSTEM_BY_SOURCE:
            tp_system = pd.Series(TP_SYSTEM_BY_SOURCE[short_id], index=t.index)
        else:
            tp_system = notes.map(tp_system_from_notes)
        tp_system = tp_system.where(tp.notna(), 'none')
        host = t.host.str.strip()
        frames.append(pd.DataFrame({
            'pl_name': host + ' [' + short_id + ' ' + t.companion.astype(str).str.strip() + ']',
            'hostname': host,
            'catalog': short_id,
            'catalog_status': tier_of.get(short_id, 'unclassified'),
            'fit_type': t.fit_type.str.strip(),
            'mass_mj': mass,
            'pl_orbper': t.P_d,
            'period_err_d': pd.to_numeric(t.e_P_d, errors='coerce'),
            ## Published 1-sigma errors on the other quantities, where the source gives
            ## them (optional columns of the orbit tables; zero or blank = not published).
            ## K errors may also sit in the notes ("e_K=12.3 m/s"). They fill in the mass
            ## and SMA errors in census.add_errors.
            'K_err': optional_error(t, 'e_K_ms', notes, r'e_K\s*=\s*([0-9.]+(?:[eE][-+]?\d+)?)'),
            'mass_err_pub': optional_error(t, 'e_msini_mj'),
            'ecc_err': optional_error(t, 'e_e', notes, r'e_e\s*=\s*([0-9.]+(?:[eE][-+]?\d+)?)'),
            'sma_err_pub': optional_error(t, 'e_a_au'),
            'pl_orbsmax': t.a_au,
            'pl_orbeccen': t.e,
            'K_ms': t.K_ms,
            'omega_deg': t.omega_deg,
            'tp_raw': tp,
            'tp_system': tp_system,
            'st_mass': t.mstar_msun,
            'sy_vmag': np.nan,
            'sy_plx': np.nan,
        }))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def simbad_url(path):
    """URL of a SIMBAD service on the first host that answers (checked once per run)."""
    if not _simbad_host:
        for host in SIMBAD_HOSTS:
            try:
                requests.get(host + '/simbad/', timeout=20)
            except requests.RequestException:
                continue
            _simbad_host.append(host)
            break
        else:
            raise RuntimeError(f'No SIMBAD host answered: {SIMBAD_HOSTS}')
    return _simbad_host[0] + path


def resolve_gaia_script(names):
    """Gaia DR3 source IDs via SIMBAD's script interface, which normalizes names the
    identifier table matches only exactly (e.g., '11 Com', 'nu Oph', 'TYC 0123-00549-1')."""
    script = ('output console=off script=off\n'
              'format object "%OBJECT|%IDLIST(Gaia DR3)"\n'
              + '\n'.join(f'query id {n}' for n in names))
    resp = requests.post(simbad_url('/simbad/sim-script'), data={'script': script},
                         timeout=300)
    resp.raise_for_status()
    rows = [line.split('|') for line in resp.text.splitlines() if '|Gaia DR3 ' in line]
    return pd.DataFrame({'key': [r[0].strip() for r in rows],
                         'source_id': [r[1].replace('Gaia DR3', '').strip() for r in rows]})


def resolve_gaia_2mass(names):
    """Gaia DR3 source IDs for '2MASS Jhhmmss...' names through Gaia's own 2MASS
    crossmatch, for APOGEE stars SIMBAD doesn't link to Gaia."""
    ids = {n: n.replace('2MASS J', '').strip() for n in names if n.startswith('2MASS J')}
    if not ids:
        return pd.DataFrame(columns=['key', 'source_id'])
    quoted = ', '.join(f"'{v}'" for v in ids.values())
    res = tap_query(GAIA_TAP, 'select original_ext_source_id, source_id from '
                              'gaiadr3.tmass_psc_xsc_best_neighbour '
                              f'where original_ext_source_id in ({quoted})')
    back = {v: k for k, v in ids.items()}
    return pd.DataFrame({'key': res.original_ext_source_id.map(back),
                         'source_id': res.source_id.astype(str)})


def resolve_gaia(names):
    """Gaia DR3 source IDs for host names, via SIMBAD's identifier table."""
    chunks = []
    for i in range(0, len(names), 300):
        quoted = ', '.join("'" + n.replace("'", "''") + "'" for n in names[i:i+300])
        adql = ('select i1.id as name, i2.id as gaia from ident as i1 join ident as i2 '
                'on i1.oidref = i2.oidref '
                f"where i1.id in ({quoted}) and i2.id like 'Gaia DR3 %'")
        chunks.append(tap_query(simbad_url('/simbad/sim-tap/sync'), adql))
    res = pd.concat(chunks)
    ## SIMBAD returns its own spacing (e.g., 'HD  11964'); normalize to match the input
    res['key'] = res.name.str.split().str.join(' ')
    res['source_id'] = res.gaia.str.replace('Gaia DR3 ', '', regex=False).str.strip()
    return res[['key', 'source_id']].drop_duplicates('key')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cached', action='store_true',
                        help='reuse the downloads in sample/raw_lit/')
    args = parser.parse_args()

    if not args.cached:
        download()
    lit = pd.concat([load_cls(), load_kiefer(), load_standard()], ignore_index=True)
    ## CLS I and Kiefer+2019 (read from VizieR) are RV-only fits
    lit['fit_type'] = lit.fit_type.fillna('RV-only')
    ## Published errors for those two (raw_lit/<source>/errors.py writes <source>_errors.csv)
    for folder, name in (('cls_2021', 'cls_2021_errors.csv'), ('kiefer_2019', 'kiefer_2019_errors.csv')):
        e = pd.read_csv(os.path.join(RAW_DIR, folder, name)).drop_duplicates('pl_name_key')
        e = e.set_index('pl_name_key')
        has = lit.pl_name.isin(e.index)
        for col, src in (('period_err_d', 'e_P_d'), ('K_err', 'e_K_ms'),
                         ('mass_err_pub', 'e_msini_mj'), ('ecc_err', 'e_e'),
                         ('sma_err_pub', 'e_a_au')):
            if col not in lit:
                lit[col] = np.nan
            vals = lit.loc[has, 'pl_name'].map(e[src])
            lit.loc[has, col] = vals.where(vals > 0)

    ## Resolve hosts not already in the SIMBAD cache. Hosts that no lookup resolves are
    ## cached with an empty source_id, so --cached runs stay offline; a full run retries them.
    ids = (pd.read_csv(SIMBAD_FILE, dtype={'source_id': str})
           if os.path.exists(SIMBAD_FILE) else pd.DataFrame(columns=['key', 'source_id']))
    if not args.cached:
        ids = ids[ids.source_id.notna()]
    new = [n for n in lit.hostname.unique() if n not in set(ids.key)]
    if new:
        found = resolve_gaia(new)
        ## Retry unresolved component names (e.g., 'HD 147379A') without the suffix
        missing = [n for n in new
                   if n not in set(found.key) and re.search(r'\d[A-Z]$', n)]
        if missing:
            retry = resolve_gaia([n[:-1] for n in missing])
            retry['key'] = retry.key.map({n[:-1]: n for n in missing})
            found = pd.concat([found, retry])
        missing = [n for n in new if n not in set(found.key)]
        if missing:
            found = pd.concat([found, resolve_gaia_script(missing)])
        missing = [n for n in new if n not in set(found.key)]
        if missing:
            found = pd.concat([found, resolve_gaia_2mass(missing)])
        unresolved = pd.DataFrame({'key': [n for n in new if n not in set(found.key)]})
        ids = pd.concat([ids, found, unresolved]).drop_duplicates('key')
        ids.to_csv(SIMBAD_FILE, index=False)
    lit = lit.merge(ids, left_on='hostname', right_on='key', how='left').drop(columns='key')
    lit.to_csv(OUT_FILE, index=False)

    print(f'{len(lit)} companions from {lit.hostname.nunique()} hosts')
    print(lit.groupby(['catalog', 'catalog_status']).size())
    print('Hosts without a Gaia DR3 ID from SIMBAD:',
          sorted(lit[lit.source_id.isna()].hostname.unique()))
