"""Turn a Jump export (query B of jump_queries.sql) into one RV table keyed by system.

Every Jump star record matched to a target (query A) belongs to that target's system,
so duplicate Jump records of one star (e.g., HD 191939 also filed as T001339 and
TIC269701147) are merged here, and each row keeps its Jump star name
(`source_record_name`). A Jump record
matched to two targets is an error. When query A's export has Jump's coordinates, each
matched record is checked against the SIMBAD position in targets.csv.

Rules (decided 2026-11-04):
- KPF RVs are left out (the instrument isn't yet stable over the long term).
- Times become BJD_TCB, the study's standard. CPS-pipeline HIRES and APF times are
  BJD_TDB (checked); imported datasets are converted only when their time scale is
  BJD_TDB, and flagged otherwise.
- Each dataset gets its own offset and jitter in the fits: CPS data by instrument
  (hires_k before the 2004 upgrade, hires_j after, apf); imports by Jump dataset, since
  different reductions of one instrument can differ in zero point.
- Every dataset is kept here, with `within_rank` = 0 for the CPS pipeline, 1 for the
  imported dataset Jump marks preferred, and 2 for other imports; merge_rvs.py keeps
  one copy of each exposure.
- Junk and quality flags are kept as columns, not applied.

Outputs (in rv_data/): jump_rvs.csv and a printed summary.

Usage (from the project root, in the gaia_dr4_sensitivity conda environment):
    python -I rv_data/ingest_jump.py --rvs rv_data/jump_raw/jump_rvs_2026-11-04.csv \
        --matches rv_data/jump_raw/jump_matches_2026-11-04.csv
"""

import argparse
import os

import numpy as np
import pandas as pd
from astropy.time import Time

HERE = os.path.dirname(os.path.abspath(__file__))
## Matched Jump records farther than this from the SIMBAD position (plus proper motion
## over POSITION_YEARS) are flagged
POSITION_TOL_ARCSEC = 30.0
POSITION_YEARS = 30.0


def check_positions(matches, targets):
    """Separation [arcsec] of each matched Jump star record from the target's SIMBAD
    position. Jump's RA may be in hours or degrees, so both are tried."""
    m = matches.dropna(subset=['star_id', 'ra', 'dec']).merge(targets, on='target',
                                                             suffixes=('_jump', ''))
    if m.empty:
        return m
    ## Small-angle separations, which is all a sanity check needs
    cos_dec = np.cos(np.radians(m.dec))
    sep = [3600*np.hypot(((ra_j - m.ra + 180) % 360 - 180)*cos_dec, m.dec_jump - m.dec)
           for ra_j in (m.ra_jump, 15*m.ra_jump)]
    m['sep_arcsec'] = np.minimum(*sep)
    pm = np.hypot(m.pmra.fillna(0), m.pmdec.fillna(0))/1000*POSITION_YEARS
    m['position_ok'] = m.sep_arcsec <= POSITION_TOL_ARCSEC + pm
    return m


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rvs', required=True, help="query B's export (CSV)")
    parser.add_argument('--matches', required=True, help="query A's export (CSV)")
    args = parser.parse_args()

    targets = pd.read_csv(os.path.join(HERE, 'targets.csv'), dtype={'system_id': str})
    rv = pd.read_csv(args.rvs, na_values=['null'])
    matches = pd.read_csv(args.matches, na_values=['null'])

    ## Identity: every Jump record belongs to exactly one target
    found = matches.dropna(subset=['star_id'])
    clash = found.groupby('star_id').target.nunique()
    if (clash > 1).any():
        raise SystemExit(f'Jump records matched to more than one target: '
                         f'{found[found.star_id.isin(clash.index[clash > 1])].to_string()}')
    unknown = set(rv.target) - set(targets.target)
    if unknown:
        raise SystemExit(f'targets missing from targets.csv: {sorted(unknown)}')
    if {'ra', 'dec'} <= set(matches.columns):
        pos = check_positions(matches, targets)
        bad = pos[~pos.position_ok]
        print(f'Position check: {len(pos) - len(bad)} of {len(pos)} Jump records within '
              f'{POSITION_TOL_ARCSEC:.0f}" (+ proper motion) of SIMBAD')
        if len(bad):
            print(bad[['target', 'jump_star_name', 'sep_arcsec']].to_string(index=False))
    else:
        print("Position check skipped: query A's export has no Jump coordinates "
              '(rerun query A from the current jump_queries.sql)')

    ## KPF out
    kpf = rv.instrument.str.lower().str.startswith('kpf')
    print(f'KPF RVs left out: {kpf.sum()}')
    rv = rv[~kpf].copy()

    ## Provenance tiers and fitting datasets
    cps = rv.database == 'Jump: CPS pipeline'
    rv['rv_source'] = 'Jump'
    rv['rv_source_detail'] = rv.database.str.replace(r'^Jump( import)?: ', '', regex=True).str.strip()
    rv['study'] = rv.publication
    rv['instrument'] = rv.instrument.str.strip()
    rv['dataset'] = np.where(cps, rv.instrument, rv.obs_id.str.split(':').str[0])
    rv['within_rank'] = np.select(
        [cps, rv.quality_flag.str.contains('preferred=true', na=False)], [0, 1], 2)
    rv['isjunk'] = rv.quality_flag.str.contains('isjunk=true', na=False)
    rv['source_record_name'] = rv.jump_star_name
    ## Exposure-level identity for duplicate checks (DUPLICATE_RULES.md): Jump's own
    ## observation IDs, and the raw file names of imported archive datasets
    rv['file_id'] = np.where(cps, rv.obs_id, rv.obs_id.str.split(':', n=1).str[1])
    rv['time_quality'] = np.where(
        cps | rv.file_id.str.contains(r'\d{4}-\d\d-\d\dT\d\d', na=False), 'exposure', 'published')

    ## Times: BJD_TDB to BJD_TCB
    tdb = cps | rv.time_scale.str.lower().isin(['bjd_tdb', 'tdb'])
    rv['bjd_tdb'] = np.where(tdb, rv.time_bjd, np.nan)
    rv['bjd_tcb'] = np.nan
    rv.loc[tdb, 'bjd_tcb'] = Time(rv.bjd_tdb[tdb].values, format='jd', scale='tdb').tcb.jd
    rv['time_flag'] = np.where(tdb, '', 'time scale ' + rv.time_scale.astype(str)
                               + '; not converted')

    out = rv.merge(targets[['target', 'system_id']], on='target')
    cols = ['system_id', 'target', 'rv_source', 'rv_source_detail', 'study',
            'source_record_name', 'obs_id', 'dataset', 'instrument', 'bjd_tdb', 'bjd_tcb',
            'time_flag', 'rv_ms', 'rv_err_ms', 'isjunk', 'quality_flag', 's_index',
            'program', 'access', 'exptime_s', 'within_rank', 'file_id', 'time_quality']
    out = out[cols].sort_values(['target', 'dataset', 'bjd_tcb'])
    out.to_csv(os.path.join(HERE, 'jump_rvs.csv'), index=False)

    print(f'{len(out)} RVs for {out.target.nunique()} systems (copies of one exposure in '
          'several datasets are resolved by merge_rvs.py)')
    summary = out.groupby(['target', 'dataset']).agg(
        n=('rv_ms', 'size'), junk=('isjunk', 'sum'),
        first=('bjd_tcb', 'min'), last=('bjd_tcb', 'max'),
        jump_records=('source_record_name', lambda s: ', '.join(sorted(set(s.astype(str))))),
        study=('study', lambda s: ', '.join(sorted(set(s.dropna()))) or '-'))
    summary['baseline_yr'] = ((summary['last'] - summary['first'])/365.25).round(1)
    print(summary.drop(columns=['first', 'last']).to_string())
    no_rv = sorted(set(targets.target) - set(out.target))
    if no_rv:
        print(f'No Jump RVs: {no_rv}')


if __name__ == '__main__':
    main()
