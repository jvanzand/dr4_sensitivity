"""Pull public DACE RVs for every system in targets.csv into one table keyed by system.

What the Jump case taught, applied here:
- *Naming.* DACE resolves some spellings itself ("HD111232" = "HD 111232", "TOI-1339" =
  HD 191939) but not others (HIP numbers are slow or fail). Each system is queried under
  its target name, SIMBAD main ID, and the aliases in DACE_ALIAS_PREFIXES, and the
  results are combined; a row returned under several names is kept once (by `rv_id`).
  Each row keeps DACE's own name for the star (`source_record_name`) and the name that
  found it. Rows whose DACE position is far from SIMBAD's are flagged.
- *Instruments and reductions.* DACE holds several RV series per instrument: its own
  DRS reductions (possibly several DRS versions and modes), archive reductions it
  labels in its publication fields (e.g. "HR DRS-3.1" for SOPHIE, "TACOS" for ELODIE),
  and series copied from publications. Each becomes its own dataset, with its own offset and jitter. For DRS
  reductions only the latest DRS version is kept (`is_latest_drs`). Instrument upgrades
  are already separate instruments in DACE (HARPS03/HARPS15, CORALIE98/07/14, ...).
- *The same exposure twice.* A publication's RVs often repeat DACE's own reduction of
  the same spectra (HD 43691: da Silva+2007's SOPHIE RVs). Every series is kept here,
  with `within_rank` = 0 for a reduction of the spectra (DACE's DRS or an archive
  reduction) and 1 for a publication; merge_rvs.py keeps one copy of each exposure.
- *Provenance.* `rv_source` = DACE; `rv_source_detail` = the DRS version and mode, the
  archive reduction's label, or "publication"; `study` = the publication's bibcode
  (or reference) where there is one.
- *Times.* DACE gives `rjd` = BJD - 2400000. Its time scale isn't stated, so for each
  instrument with a known site the script compares rjd with the exposure midpoint
  (from the raw file name's UTC start plus half the exposure time) and reports the
  median offset from the BJD_TDB model; ~0 s means TDB, ~69 s means UTC. Times are
  converted to BJD_TCB on the BJD_TDB assumption, and publication series, which may
  use other conventions, are flagged.
- *Flags, not cuts.* `drs_qc` and DACE's public status are columns.

Raw responses are cached in rv_data/dace_raw/ (one CSV per system, plus manifest.csv with
the names queried, row counts, and retrieval times) and reused unless --refresh.

Outputs (in rv_data/): dace_rvs.csv, dace_time_check.csv, and a
printed summary. The columns of dace_rvs.csv match jump_rvs.csv.

Usage (from the project root, in the gaia_dr4_sensitivity conda environment):
    python -I rv_data/pull_dace.py            # query only uncached systems
    python -I rv_data/pull_dace.py --refresh  # re-query every system
"""

import argparse
import concurrent.futures as cf
import datetime
import logging
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
logging.disable(logging.CRITICAL)  # dace-query logs a warning per call in public mode

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(HERE, 'dace_raw')
MANIFEST = os.path.join(RAW_DIR, 'manifest.csv')

## Aliases worth trying in DACE, besides the target name and SIMBAD's main ID
## (Durchmusterung names, e.g. "BD+66 1273", time out in DACE, so they are not tried.)
DACE_ALIAS_PREFIXES = ('HD ', 'HIP ', 'HR ', 'GJ ', 'Gl ', 'TOI-', 'WASP-', 'HAT-P-', 'KELT-',
                       'K2-', 'Kepler-', 'CoRoT-', 'NGTS-', 'LHS ', 'LP ', 'Wolf ', 'Ross ',
                       'NAME ')
CALL_TIMEOUT_S = 90
POSITION_TOL_ARCSEC = 30.0
POSITION_YEARS = 30.0

## Observatory per DACE instrument (prefix match), for the time-scale check
SITES = {
    'HARPS': EarthLocation.from_geodetic(-70.7346*u.deg, -29.2584*u.deg, 2400*u.m),     # La Silla
    'CORALIE': EarthLocation.from_geodetic(-70.7346*u.deg, -29.2584*u.deg, 2400*u.m),
    'ESPRESSO': EarthLocation.from_geodetic(-70.4045*u.deg, -24.6272*u.deg, 2635*u.m),  # Paranal
    'SOPHIE': EarthLocation.from_geodetic(5.7133*u.deg, 43.9308*u.deg, 650*u.m),        # OHP
    'ELODIE': EarthLocation.from_geodetic(5.7133*u.deg, 43.9308*u.deg, 650*u.m),
    'HARPN': EarthLocation.from_geodetic(-17.8890*u.deg, 28.7540*u.deg, 2387*u.m),      # La Palma
}

KEEP = ['rv_id', 'spectrum_id', 'rjd', 'rv', 'rv_err', 'texp', 'instrument_name', 'ins_mode',
        'drs_qc', 'is_latest_drs', 'version_major', 'version_minor', 'version_patch',
        'rv_extraction_method', 'pub_reference', 'pub_bibcode', 'program_id',
        'program_code', 'status_public', 'target_catname', 'ra_deg', 'dec_deg',
        'file_rootname', 'spectro_smw', 'spectro_rhk', 'ccf_fwhm', 'ccf_bispan']


def query_names(target):
    """DACE names to try for one system, most likely first."""
    names = [target.target, re.sub(r'\s+', ' ', str(target.simbad_main_id)).strip()]
    names += [a for a in target.aliases.split('; ') if a.startswith(DACE_ALIAS_PREFIXES)]
    return list(dict.fromkeys(names))


def fetch_system(target, pool):
    """All DACE rows for one system, over every name, with a per-name status."""
    from dace_query.spectroscopy import Spectroscopy

    def call(name):
        return Spectroscopy.get_timeseries(name, sorted_by_instrument=False,
                                           output_format='pandas')

    futures = {name: pool.submit(call, name) for name in query_names(target)}
    frames, log = [], []
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    for name, fut in futures.items():
        try:
            d = fut.result(timeout=CALL_TIMEOUT_S)
            status = 'ok'
        except cf.TimeoutError:
            d, status = pd.DataFrame(), f'timeout ({CALL_TIMEOUT_S} s)'
        except Exception as err:  # dace-query raises plain Exceptions for "no data"
            d, status = pd.DataFrame(), f'error: {err!r}'[:200]
        log.append({'system_id': target.system_id, 'query_name': name, 'n_rows': len(d),
                    'status': status, 'retrieved_utc': now})
        if len(d):
            frames.append(d[[c for c in KEEP if c in d.columns]].assign(query_name=name))
    raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=KEEP)
    return raw, pd.DataFrame(log)


def time_check(rows, targets):
    """Median offset [s] of DACE's BJD from the BJD_TDB model, per instrument and series,
    from the raw file's UTC start time plus half the exposure time."""
    out = []
    t_start = pd.to_datetime(rows.file_rootname.str.extract(
        r'(\d{4}-\d\d-\d\dT\d\d[:\-]\d\d[:\-]\d\d(?:\.\d+)?)')[0].str.replace(
        r'T(\d\d)-(\d\d)-(\d\d)', r'T\1:\2:\3', regex=True), errors='coerce')
    rows = rows.assign(t_start=t_start).dropna(subset=['t_start', 'texp'])
    for (system_id, ins, detail), g in rows.groupby(['system_id', 'instrument_name',
                                                    'rv_source_detail']):
        site = next((loc for key, loc in SITES.items() if ins.upper().startswith(key)), None)
        if site is None:
            continue
        tg = targets.loc[system_id]
        t = Time(g.t_start.values, scale='utc', location=site) + (g.texp.values/2)*u.s
        tdb = (t.tdb + t.light_travel_time(SkyCoord(tg.ra*u.deg, tg.dec*u.deg))).jd
        off = (g.rjd.values + 2400000 - tdb)*86400
        out.append({'system_id': system_id, 'instrument': ins, 'series': detail,
                    'n': len(g), 'median_offset_s': round(float(np.median(off)), 1),
                    'p16_s': round(float(np.percentile(off, 16)), 1),
                    'p84_s': round(float(np.percentile(off, 84)), 1)})
    return pd.DataFrame(out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--refresh', action='store_true', help='re-query every system')
    args = parser.parse_args()

    os.makedirs(RAW_DIR, exist_ok=True)
    targets = pd.read_csv(os.path.join(HERE, 'targets.csv'), dtype={'system_id': str})
    manifest = (pd.read_csv(MANIFEST, dtype={'system_id': str}) if os.path.exists(MANIFEST)
                else pd.DataFrame(columns=['system_id']))

    ## Query DACE for systems not yet cached
    todo = targets if args.refresh else targets[~targets.system_id.isin(manifest.system_id)]
    pool = cf.ThreadPoolExecutor(max_workers=4)
    for t in todo.itertuples():
        raw, log = fetch_system(t, pool)
        raw.to_csv(os.path.join(RAW_DIR, f'{t.system_id}.csv'), index=False)
        manifest = pd.concat([manifest[manifest.system_id != t.system_id], log],
                             ignore_index=True)
        manifest.to_csv(MANIFEST, index=False)
        print(f'{t.target}: ' + ', '.join(f'{r.query_name} {r.n_rows}'
                                          + ('' if r.status == 'ok' else f' ({r.status})')
                                          for r in log.itertuples()))
    pool.shutdown(wait=False, cancel_futures=True)

    ## Combine the cached responses
    frames = []
    for t in targets.itertuples():
        path = os.path.join(RAW_DIR, f'{t.system_id}.csv')
        if os.path.exists(path) and os.path.getsize(path) > 1:
            d = pd.read_csv(path, dtype={'rv_id': str, 'spectrum_id': str})
            if len(d):
                frames.append(d.assign(system_id=t.system_id, target=t.target))
    rv = pd.concat(frames, ignore_index=True)
    ## A row found under several names is kept once
    key = rv.rv_id.where(rv.rv_id.notna(), rv.instrument_name + '|' + rv.rjd.astype(str)
                         + '|' + rv.rv.astype(str))
    rv = rv.assign(_key=key).groupby(['system_id', '_key'], as_index=False).agg(
        lambda s: s.iloc[0] if s.name != 'query_name' else '; '.join(sorted(set(s))))

    ## Series: DACE's DRS reductions (latest version only), archive reductions that DACE
    ## labels in its publication fields, or publications (a bibcode or a reference)
    for col in ['pub_reference', 'pub_bibcode']:
        rv[col] = rv[col].str.replace(r'\s+', ' ', regex=True).str.strip()
    label = rv.pub_bibcode.where(rv.pub_bibcode.notna(), rv.pub_reference)
    is_paper = (label.str.match(r'^\d{4}[A-Za-z&.]', na=False)
                | label.str.contains(r'et al|\(\d{4}\)', na=False))
    archive = label.notna() & ~is_paper
    pub = is_paper
    drs = ('DRS-' + rv.version_major.astype('Int64').astype(str) + '.'
           + rv.version_minor.astype('Int64').astype(str) + '.'
           + rv.version_patch.astype('Int64').astype(str))
    rv['rv_source_detail'] = np.where(
        pub, 'publication', np.where(archive, 'archive reduction: ' + label.astype(str),
        drs + np.where(rv.rv_extraction_method.notna(), '-' + rv.rv_extraction_method.astype(str), '')
        + np.where(rv.ins_mode.notna(), ' ' + rv.ins_mode.astype(str), '')))
    old_drs = ~pub & ~archive & (rv.is_latest_drs.astype(str).str.lower() == 'false')
    print(f'Older DRS versions left out: {old_drs.sum()}')
    rv = rv[~old_drs].copy()
    pub, archive = pub[~old_drs], archive[~old_drs]
    rv['study'] = np.where(pub, rv.pub_bibcode.where(rv.pub_bibcode.notna(), rv.pub_reference),
                           None)
    series = np.where(pub, rv.study.astype(str),
                      rv.rv_source_detail.str.replace('archive reduction: ', '', regex=False))
    rv['dataset'] = ('dace:' + rv.instrument_name.str.replace(r'\s+', '_', regex=True) + ':'
                     + pd.Series(series, index=rv.index).str.replace(' ', '_'))

    ## Positions against SIMBAD
    tg = targets.set_index('system_id')
    t_ra, t_dec = tg.ra.reindex(rv.system_id).values, tg.dec.reindex(rv.system_id).values
    pm = np.hypot(tg.pmra.fillna(0), tg.pmdec.fillna(0)).reindex(rv.system_id).values
    sep = 3600*np.hypot(((rv.ra_deg - t_ra + 180) % 360 - 180)*np.cos(np.radians(t_dec)),
                        rv.dec_deg - t_dec)
    rv['far'] = (sep > POSITION_TOL_ARCSEC + pm/1000*POSITION_YEARS).values

    ## Reductions of the spectra rank above publications of the same exposures
    rv['within_rank'] = np.where(pub, 1, 0)

    ## Times
    checks = time_check(rv, tg)
    checks.to_csv(os.path.join(HERE, 'dace_time_check.csv'), index=False)
    rv['bjd_tdb'] = rv.rjd + 2400000
    rv['bjd_tcb'] = Time(rv.bjd_tdb.values, format='jd', scale='tdb').tcb.jd
    rv['time_flag'] = np.where(rv.rv_source_detail.str.startswith('DRS'),
                               'DACE rjd; time scale taken as BJD_TDB',
                               np.where(rv.rv_source_detail == 'publication',
                               'published time stamps; time scale taken as BJD_TDB',
                               'archive reduction; time scale taken as BJD_TDB'))

    ## Common columns (as jump_rvs.csv)
    out = pd.DataFrame({
        'system_id': rv.system_id, 'target': rv.target, 'rv_source': 'DACE',
        'rv_source_detail': rv.rv_source_detail, 'study': rv.study,
        'source_record_name': rv.target_catname, 'obs_id': rv.rv_id,
        'dataset': rv.dataset, 'instrument': rv.instrument_name,
        'bjd_tdb': rv.bjd_tdb, 'bjd_tcb': rv.bjd_tcb, 'time_flag': rv.time_flag,
        'rv_ms': rv.rv, 'rv_err_ms': rv.rv_err,
        'isjunk': rv.drs_qc.astype(str).str.lower() == 'false',
        'quality_flag': ('drs_qc=' + rv.drs_qc.astype(str)
                         + np.where(rv.far, '; position far from SIMBAD', '')),
        's_index': rv.spectro_smw,
        'program': rv.program_code.where(rv.program_code.notna(), rv.program_id),
        'access': 'status_public=' + rv.status_public.astype(str),
        'exptime_s': rv.texp, 'query_name': rv.query_name, 'within_rank': rv.within_rank,
        ## Exposure-level identity for duplicate checks (DUPLICATE_RULES.md): DACE's DRS
        ## rows name their raw file; uploaded and published series don't
        'file_id': np.where(rv.rv_source_detail.str.startswith('DRS'), rv.file_rootname, None),
        'time_quality': np.where(rv.rv_source_detail.str.startswith('DRS'), 'exposure',
                                 'published'),
    }).sort_values(['target', 'dataset', 'bjd_tcb'])
    out.to_csv(os.path.join(HERE, 'dace_rvs.csv'), index=False)

    print(f'Rows far from SIMBAD: {int(rv.far.sum())}')
    print(f'{len(out)} RVs for {out.target.nunique()} systems (copies of one exposure in '
          'several series are resolved by merge_rvs.py)')
    summary = out.groupby(['target', 'dataset']).agg(
        n=('rv_ms', 'size'), qc_fail=('isjunk', 'sum'),
        first=('bjd_tcb', 'min'), last=('bjd_tcb', 'max'),
        dace_names=('source_record_name', lambda s: ', '.join(sorted(set(s.astype(str))))))
    summary['baseline_yr'] = ((summary['last'] - summary['first'])/365.25).round(1)
    print(summary.drop(columns=['first', 'last']).to_string())
    print('\nTime check (offset of DACE BJD from the BJD_TDB model; ~0 s = TDB, ~69 s = UTC):')
    print(checks.to_string(index=False) if len(checks) else '  no rows with raw file times')
    no_rv = sorted(set(targets.target) - set(out.target))
    if no_rv:
        print(f'No DACE RVs: {no_rv}')
    sys.stdout.flush()
    os._exit(0)  # don't wait for a hung DACE call


if __name__ == '__main__':
    main()
