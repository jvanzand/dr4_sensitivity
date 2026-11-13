"""Merge the per-source RV tables (jump_rvs.csv, dace_rvs.csv) into one table per system,
keeping each exposure once.

The rules for duplicates are in DUPLICATE_RULES.md; in short, within an instrument
family, datasets are ranked (source, then reduction before publication, then size), and
an RV of a lower-ranked dataset is dropped as a copy if it shares a raw-frame ID with a
higher-ranked RV (E1) or is paired with one in time (E2: one to one, nearest first,
within 2 min for exposure-level times or 15 min when a published series is involved).
Pairs are then checked dataset by dataset: true copies differ by a constant offset, so
pairs far from it are rejected (both RVs kept) and flagged.

Each kept RV keeps its source columns: rv_source (database), rv_source_detail
(reduction), study, source_record_name (the source's own name for the star), and dataset
(own offset and jitter in the fits).

Outputs (in rv_data/):
    rvs_merged.csv                     every kept RV, all systems
    systems/<system_id>_<target>.csv   the same, one file per system
    rv_coverage.csv                    per system and dataset: RVs, time span, sources,
                                       studies, and programs that may need permission
    rv_duplicates.csv                  each dropped RV, with its evidence and partner
    rv_duplicate_pairs.csv             the checks, per pair of datasets

Usage (from the project root, in the gaia_dr4_sensitivity conda environment):
    python -I rv_data/merge_rvs.py
"""

import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = ['jump_rvs.csv', 'dace_rvs.csv']
WINDOW_EXPOSURE_S = 120.0   # E2 window when both datasets have exposure-level times
WINDOW_PUBLISHED_S = 900.0  # E2 window when either is a published or uploaded series
REJECT_SIGMA = 5.0
INCONSISTENT_RATIO = 3.0
MIN_PAIRS_TO_CHECK = 3
SMALL_REMAINDER = 5
## An RV with an error this large (m/s) is a placeholder, not a measurement (HD 43691 has a
## SOPHIE row with RV 999001.8 and error 999000 from DACE, not flagged by its QC); it is
## marked isjunk. A general outlier rule is still to be decided (QUESTIONS.md, 2026-11-06).
PLACEHOLDER_ERR_MS = 1e5
## Outlier rule (decided 2026-11-09): per system, center every dataset on its own median,
## take the pooled MAD of those residuals as a robust sigma (1.4826 x median |deviation|;
## equal to the standard deviation for Gaussian noise, but one wild value cannot inflate it),
## and mark an RV isjunk if it lies more than OUTLIER_SIGMA of these from the median.
OUTLIER_SIGMA = 10.0
OUTLIER_MIN_N = 8       # smaller series have no meaningful MAD
## Jump's CPS programs that are public; any other CPS program may need permission
PUBLIC_PROGRAMS = {'JUMP_LEGACY_PUBLIC', 'APF_PUBLIC'}

## Instrument families: every source's name for one spectrograph, all upgrade eras
FAMILIES = {
    'hires': ['hires_j', 'hires_k', 'HIRES', 'HIRES-PRE04', 'HIRES-POST04'],
    'apf': ['apf', 'APF'],
    'harps': ['HARPS03', 'HARPS15', 'harps03', 'harps15', 'HARPS'],
    'harpn': ['HARPN', 'harpn', 'HARPS-N'],
    'coralie': ['CORALIE', 'CORALIE98', 'CORALIE07', 'CORALIE14', 'coralie98', 'coralie07',
                'coralie14'],
    'coravel': ['CORAVEL'],
    'sophie': ['SOPHIE', 'sophie', 'sophie06', 'sophie11'],
    'elodie': ['ELODIE', 'elodie'],
    'espresso': ['ESPRESSO', 'ESPRESSO18', 'ESPRESSO19', 'espresso18', 'espresso19'],
    'mike': ['MIKE ECHELLE SPECTROGRAPH', 'mike'],
    'neid': ['neid', 'neid_era4', 'neid_era5', 'neid_era6', 'neid_era7'],
    'hamilton': ['hamilton', 'HAMILTON'],
    'ucles': ['UCLES'],
    'ces': ['CES'],
    'carmenes': ['carmenes_vis', 'CARMENES'],
}
FAMILY_OF = {name: fam for fam, names in FAMILIES.items() for name in names}

TIMESTAMP = re.compile(r'(\d{4}-\d\d-\d\d)T(\d\d)[:\-](\d\d)[:\-](\d\d)')


def family(instrument):
    name = str(instrument).strip()
    return FAMILY_OF.get(name, FAMILY_OF.get(name.split()[0] if name else name))


def exposure_key(file_id):
    """Raw-frame timestamp from a file name, to the second (None if there is none)."""
    m = TIMESTAMP.search(str(file_id))
    return f'{m.group(1)}T{m.group(2)}:{m.group(3)}:{m.group(4)}' if m else None


def rank_datasets(rv):
    """Rank the datasets of each system and family (0 = kept first)."""
    rv['tier'] = np.select([rv.rv_source.eq('Jump') & rv.rv_source_detail.eq('CPS pipeline'),
                            rv.rv_source.eq('DACE')], [0, 1], 2)
    rv['within_rank'] = rv.within_rank.fillna(0)
    rv['n_dataset'] = rv.groupby(['system_id', 'dataset']).rv_ms.transform('size')
    order = (rv.drop_duplicates(['system_id', 'family', 'dataset'])
             .sort_values(['system_id', 'family', 'tier', 'within_rank', 'n_dataset', 'dataset'],
                          ascending=[True, True, True, True, False, True]))
    order['rank'] = order.groupby(['system_id', 'family']).cumcount()
    return rv.merge(order[['system_id', 'family', 'dataset', 'rank']],
                    on=['system_id', 'family', 'dataset'])


def match_family(g):
    """Pairs (copy index, kept index, evidence, candidates) within one system and family."""
    pairs = []
    dropped = set()
    for rank in sorted(g['rank'].unique())[1:]:
        low = g[(g['rank'] == rank) & ~g.index.isin(dropped)]
        high = g[(g['rank'] < rank) & ~g.index.isin(dropped)]
        if low.empty or high.empty:
            continue
        used_low, used_high = set(), set()
        ## E1: the same raw frame
        keys = high.dropna(subset=['exposure_key']).groupby('exposure_key').head(1)
        key_to_high = dict(zip(keys.exposure_key, keys.index))
        for i, key in low.exposure_key.dropna().items():
            j = key_to_high.get(key)
            if j is not None and j not in used_high:
                pairs.append((i, j, 'E1 exposure ID', 1))
                used_low.add(i)
                used_high.add(j)
        ## E2: time, one to one, nearest first
        cand = []
        for i, row in low[~low.index.isin(used_low)].iterrows():
            published_low = row.time_quality != 'exposure'
            dt = (high.bjd_tdb - row.bjd_tdb).abs()*86400
            window = np.where(published_low | (high.time_quality != 'exposure'),
                              WINDOW_PUBLISHED_S, WINDOW_EXPOSURE_S)
            ok = dt[dt.values < window]
            for j, d in ok.items():
                cand.append((d, i, j, len(ok)))
        for d, i, j, n in sorted(cand):
            if i in used_low or j in used_high:
                continue
            pairs.append((i, j, 'E2 time', n))
            used_low.add(i)
            used_high.add(j)
        dropped |= used_low
    return pairs


def resolve_duplicates(rv):
    """Apply DUPLICATE_RULES.md to a combined RV table. Returns the kept RVs, the dropped
    copies (with their evidence), the per-dataset-pair checks, the pairs rejected for
    review, and the number of exact repeats dropped."""
    rv = rv.copy()
    rv['family'] = rv.instrument.map(family)
    unknown = rv.family.isna()
    if unknown.any():
        print('WARNING: instruments without a family (add them to FAMILIES): '
              f'{sorted(set(rv.instrument[unknown]))}')
        rv.loc[unknown, 'family'] = rv.instrument[unknown].str.lower()
    rv['exposure_key'] = rv.file_id.map(exposure_key) if 'file_id' in rv else None
    rv['time_quality'] = rv.time_quality.fillna('published')

    ## Exact repeats within one dataset
    repeat = rv.assign(t=rv.bjd_tdb.round(5)).duplicated(['system_id', 'dataset', 't', 'rv_ms'])
    rv = rank_datasets(rv[~repeat].copy())

    ## Copies across datasets
    records = []
    for _, g in rv.groupby(['system_id', 'family']):
        for i, j, evidence, n in match_family(g):
            records.append({'copy': i, 'kept': j, 'evidence': evidence, 'candidates': n})
    p = pd.DataFrame(records, columns=['copy', 'kept', 'evidence', 'candidates'])
    a, b = rv.loc[p['copy']].reset_index(drop=True), rv.loc[p['kept']].reset_index(drop=True)
    p = p.assign(target=a.target.values, dataset=a.dataset.values, kept_dataset=b.dataset.values,
                 dt_s=(a.bjd_tdb.values - b.bjd_tdb.values)*86400,
                 drv_ms=a.rv_ms.values - b.rv_ms.values,
                 err_ms=np.hypot(a.rv_err_ms.values, b.rv_err_ms.values))

    ## Check each pair of datasets: copies differ by a constant offset
    p['flag'] = np.where(p.candidates > 1, 'ambiguous', '')
    checks = []
    for (target, ds, kept_ds), q in p.groupby(['target', 'dataset', 'kept_dataset']):
        offset = q.drv_ms.median()
        scatter = 1.4826*(q.drv_ms - offset).abs().median()
        expected = q.err_ms.median()
        flags = []
        if len(q) >= MIN_PAIRS_TO_CHECK:
            bad = (q.drv_ms - offset).abs() > REJECT_SIGMA*max(scatter, expected)
            p.loc[q.index[bad], 'flag'] = 'rejected: RV difference far from the offset'
            if scatter > INCONSISTENT_RATIO*expected:
                flags.append('scatter > 3x expected: review')
        else:
            bad = pd.Series(False, index=q.index)
            e2 = q.evidence.str.startswith('E2')
            p.loc[q.index[e2], 'flag'] = (p.loc[q.index[e2], 'flag'] + '; unvalidated').str.strip('; ')
            flags.append(f'fewer than {MIN_PAIRS_TO_CHECK} pairs: unvalidated')
        n_low = int(((rv.target == target) & (rv.dataset == ds)).sum())
        checks.append({'target': target, 'dataset': ds, 'kept_dataset': kept_ds,
                       'n_pairs': len(q), 'n_rejected': int(bad.sum()),
                       'fraction_of_dataset': round(len(q)/n_low, 2),
                       'evidence': '; '.join(sorted(set(q.evidence))),
                       'offset_ms': round(offset, 2), 'scatter_ms': round(scatter, 2),
                       'expected_ms': round(expected, 2),
                       'max_abs_dt_s': round(q.dt_s.abs().max(), 1),
                       'max_candidates': int(q.candidates.max()), 'flags': '; '.join(flags)})
    checks = pd.DataFrame(checks)
    accepted = p[~p.flag.str.startswith('rejected')]

    ## Remainders: a lower-ranked dataset's RVs that copy nothing (counted over all the
    ## datasets it copies)
    if len(checks):
        n_copied = accepted.groupby(['target', 'dataset']).size()
        for k, row in checks.iterrows():
            n_low = int(((rv.target == row.target) & (rv.dataset == row.dataset)).sum())
            remainder = n_low - int(n_copied.get((row.target, row.dataset), 0))
            if 0 < remainder < SMALL_REMAINDER:
                checks.at[k, 'flags'] = '; '.join(filter(None, [
                    row['flags'], f'small remainder ({remainder} RVs; offset to '
                               f'{row.kept_dataset} {row.offset_ms:.1f} m/s)']))

    ## Outputs
    drop = rv.loc[accepted['copy']].assign(
        evidence=accepted.evidence.values, copy_of=rv.loc[accepted['kept'], 'obs_id'].values,
        copy_dataset=accepted.kept_dataset.values, copy_dt_s=accepted.dt_s.values,
        copy_drv_ms=accepted.drv_ms.values, copy_candidates=accepted.candidates.values,
        copy_flag=accepted.flag.values)
    review = p[p.flag.str.startswith('rejected')]
    out = rv.drop(index=accepted['copy']).drop(
        columns=['tier', 'n_dataset', 'rank', 'exposure_key'], errors='ignore')
    out['dup_review'] = np.where(out.index.isin(review['copy']) | out.index.isin(review['kept']),
                                 'rejected copy pair: review', '')
    out = out.sort_values(['target', 'dataset', 'bjd_tcb'])
    return out, drop, checks, review, int(repeat.sum())


def flag_outliers(rv, nsigma=OUTLIER_SIGMA):
    """MAD outlier rule (see OUTLIER_SIGMA). Marks rows isjunk in place, appends the reason to
    quality_flag, and returns (per-system statistics, the flagged rows)."""
    stats, flagged = [], []
    ok = ~rv.isjunk.astype(bool) & rv.rv_ms.notna() & rv.rv_err_ms.notna()
    for sid, g in rv[ok].groupby('system_id'):
        dev = g.rv_ms - g.groupby('dataset').rv_ms.transform('median')
        center = np.median(dev)
        sigma = 1.4826*np.median(np.abs(dev - center))
        if len(g) < OUTLIER_MIN_N or not sigma > 0:
            stats.append({'system_id': sid, 'target': g.target.iloc[0], 'n': len(g),
                          'mad_sigma_ms': sigma, 'max_z': np.nan, 'n_flagged': 0})
            continue
        z = np.abs(dev - center)/sigma
        bad = z.index[z > nsigma]
        for i in bad:
            rv.loc[i, 'isjunk'] = True
            rv.loc[i, 'quality_flag'] = (str(rv.loc[i, 'quality_flag']).replace('nan', '')
                                         + f'; outlier: {z[i]:.1f} MAD-sigma from the dataset '
                                           f'median (limit {nsigma:g})')
        flagged.append(rv.loc[bad, ['system_id', 'target', 'dataset', 'bjd_tcb', 'rv_ms',
                                    'rv_err_ms']].assign(z_mad=z[bad].values))
        stats.append({'system_id': sid, 'target': g.target.iloc[0], 'n': len(g),
                      'mad_sigma_ms': sigma, 'max_z': float(z.max()), 'n_flagged': len(bad)})
    return (pd.DataFrame(stats),
            pd.concat(flagged) if flagged else pd.DataFrame(
                columns=['system_id', 'target', 'dataset', 'bjd_tcb', 'rv_ms', 'rv_err_ms', 'z_mad']))


def selftest():
    """A clean Gaussian series with two offset datasets and one wild point."""
    rng = np.random.default_rng(0)
    n = 60
    rv = pd.DataFrame({'system_id': 'x', 'target': 'x', 'dataset': ['a', 'b']*(n//2),
                       'bjd_tcb': np.arange(n, dtype=float),
                       'rv_ms': rng.normal(0, 5, n) + np.where(np.arange(n) % 2, 3e4, 0.0),
                       'rv_err_ms': 2.0, 'isjunk': False, 'quality_flag': ''})
    rv.loc[7, 'rv_ms'] += 400.0
    stats, flagged = flag_outliers(rv)
    assert list(flagged.index) == [7] and bool(rv.loc[7, 'isjunk']) and rv.isjunk.sum() == 1
    print('selftest passed:', stats.iloc[0].to_dict())


def main():
    rv = pd.concat([pd.read_csv(os.path.join(HERE, f), dtype={'system_id': str})
                    for f in SOURCES if os.path.exists(os.path.join(HERE, f))],
                   ignore_index=True)
    placeholder = rv.rv_err_ms.abs() >= PLACEHOLDER_ERR_MS
    rv.loc[placeholder, 'isjunk'] = True
    rv.loc[placeholder, 'quality_flag'] = (
        rv.loc[placeholder, 'quality_flag'].fillna('').astype(str)
        + '; invalid: RV error >= 1e5 m/s (placeholder, not a measurement)')
    print(f'{int(placeholder.sum())} placeholder RVs (error >= {PLACEHOLDER_ERR_MS:.0e} m/s) '
          f'marked isjunk')
    out_stats, out_rows = flag_outliers(rv)
    out_stats.to_csv(os.path.join(HERE, 'rv_outlier_stats.csv'), index=False)
    out_rows.to_csv(os.path.join(HERE, 'rv_outliers.csv'), index=False)
    print(f'{len(out_rows)} outliers (> {OUTLIER_SIGMA:g} MAD-sigma) marked isjunk; per-system '
          f'statistics in rv_outlier_stats.csv')
    print(out_stats.to_string(index=False))
    out, drop, checks, review, n_repeat = resolve_duplicates(rv)
    accepted = drop
    drop.to_csv(os.path.join(HERE, 'rv_duplicates.csv'), index=False)
    checks.to_csv(os.path.join(HERE, 'rv_duplicate_pairs.csv'), index=False)
    out.to_csv(os.path.join(HERE, 'rvs_merged.csv'), index=False)
    out.to_csv(os.path.join(HERE, 'rvs_merged.csv'), index=False)
    os.makedirs(os.path.join(HERE, 'systems'), exist_ok=True)
    for (sid, target), g in out.groupby(['system_id', 'target']):
        g.to_csv(os.path.join(HERE, 'systems', f"{sid}_{re.sub(r'[^A-Za-z0-9+-]', '', target)}.csv"),
                 index=False)

    ## Coverage, with programs that may need permission
    restricted = ((out.rv_source.eq('Jump') & out.rv_source_detail.eq('CPS pipeline')
                   & ~out.program.isin(PUBLIC_PROGRAMS))
                  | out.access.astype(str).str.contains('status_public=False'))
    cov = out.assign(restricted=np.where(restricted, out.program.astype(str), None)).groupby(
        ['system_id', 'target', 'dataset']).agg(
        n=('rv_ms', 'size'), flagged=('isjunk', 'sum'),
        first_bjd_tcb=('bjd_tcb', 'min'), last_bjd_tcb=('bjd_tcb', 'max'),
        source=('rv_source', 'first'), detail=('rv_source_detail', 'first'),
        study=('study', lambda s: '; '.join(sorted(set(s.dropna().astype(str))))),
        source_record_names=('source_record_name', lambda s: '; '.join(sorted(set(s.astype(str))))),
        restricted_programs=('restricted', lambda s: '; '.join(sorted(set(s.dropna())))))
    cov['baseline_yr'] = ((cov.last_bjd_tcb - cov.first_bjd_tcb)/365.25).round(2)
    cov.reset_index().to_csv(os.path.join(HERE, 'rv_coverage.csv'), index=False)

    print(f'{len(rv)} RVs in; {n_repeat} exact repeats and '
          f'{len(accepted)} copies dropped ({(accepted.evidence.str.startswith("E1")).sum()} by '
          f'exposure ID, {(accepted.evidence.str.startswith("E2")).sum()} by time); '
          f'{len(review)} pairs rejected for review')
    print(f'{len(out)} RVs kept for {out.system_id.nunique()} systems')
    print(cov.drop(columns=['first_bjd_tcb', 'last_bjd_tcb', 'source_record_names'])
          .reset_index(level=0, drop=True).to_string())
    print('\nCopy checks, per pair of datasets (rv_duplicate_pairs.csv):')
    print(checks.drop(columns=['target']).to_string(index=False) if len(checks) else '  none')


if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main()
