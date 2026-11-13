"""How well does a (blind) RV search recover the literature companions?

For each system folder in fits_rv/<results-dir>/ this compares the companions of the final
model of `companion_search.py` with the companions the literature reports for the same star
(sample/stage1_catalog.csv), one to one. A literature companion counts as RECOVERED if some
model companion has a period within `RECOVERY_TOL` (10%) of it; each model companion can match
only one literature companion (nearest in period first). Literature companions the search
could not have found are labeled rather than counted as misses: those with a period below the
blind period floor (2 x mean RV spacing) or above the period cap (2 x baseline).

Model companions that match nothing are "extras": real companions the literature missed,
spurious signals (aliases, activity, noise), or curvature stand-ins; the table cannot tell
which.

Usage (project root, gaia_dr4_sensitivity environment):
    python -I fits_rv/recovery.py                         # fits_rv/results, writes recovery.csv
    python -I fits_rv/recovery.py --results-dir results_WithLitPlanets --out /tmp/r.csv
    python -I fits_rv/recovery.py --selftest
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import companion_search as cs

RECOVERY_TOL = 0.10        # |P_model / P_lit - 1| for a match


def match(p_lit, p_model, tol=RECOVERY_TOL):
    """One-to-one matching in period: all (literature, model) pairs within `tol`, nearest in
    ln P first, each companion used once. Returns {literature index: model index}."""
    pairs = sorted((abs(np.log(pm/pl)), i, j) for i, pl in enumerate(p_lit)
                   for j, pm in enumerate(p_model) if abs(pm/pl - 1) <= tol)
    used_i, used_j, out = set(), set(), {}
    for _, i, j in pairs:
        if i not in used_i and j not in used_j:
            out[i] = j
            used_i.add(i)
            used_j.add(j)
    return out


def model_companions(res, data):
    """The final model's companions in parameter order, each with the Delta-BIC of the step
    that added it (the first adopted step with more companions than its index)."""
    final = [s for s in res['steps'] if s['adopted']][-1]
    pl = cs.unpack(np.array(final['theta']), final['n_planets'], final['trend_order'],
                   data.nds)[0]
    out = []
    for j, p in enumerate(pl):
        d = next((s['dBIC'] for s in res['steps'] if s['adopted'] and s['n_planets'] > j),
                 np.nan)
        out.append({'P': float(np.exp(p[0])), 'K': float(np.exp(p[1])),
                    'e': float(min(p[2]**2 + p[3]**2, cs.E_MAX)), 'dBIC': float(d)})
    return out, final['trend_order']


def table(results_dir):
    rows, extras = [], []
    for f in sorted(glob.glob(os.path.join(HERE, results_dir, '*', 'companion_search.json'))):
        res = json.load(open(f))
        name = os.path.basename(os.path.dirname(f))
        rv = res['rv_table'] if os.path.isabs(res['rv_table']) else os.path.join(ROOT, res['rv_table'])
        data = cs.load_system(rv)
        lit = sorted(cs.literature_priors(res['system_id']), key=lambda L: L['P0'])
        model, trend = model_companions(res, data)
        m = match([L['P0'] for L in lit], [c['P'] for c in model])
        for i, L in enumerate(lit):
            if L['P0'] < data.p_min_blind:
                where = 'below period floor'
            elif L['P0'] > data.pmax:
                where = 'above period cap'
            else:
                where = 'searchable'
            c = model[m[i]] if i in m else None
            rows.append({'system': name, 'companion': L['name'], 'P_lit_d': L['P0'],
                         'K_lit_ms': L['K0'], 'e_lit': L['e0'], 'range': where,
                         'recovered': c is not None,
                         'P_found_d': c['P'] if c else np.nan,
                         'K_found_ms': c['K'] if c else np.nan,
                         'e_found': c['e'] if c else np.nan,
                         'dBIC_at_adoption': c['dBIC'] if c else np.nan})
        for j, c in enumerate(model):
            if j not in m.values():
                extras.append({'system': name, 'P_d': c['P'], 'K_ms': c['K'], 'e': c['e'],
                               'dBIC_at_adoption': c['dBIC']})
        rows.append({'system': name, 'companion': f'(trend order {trend}; floor '
                     f'{data.p_min_blind:.0f} d, cap {data.pmax:.0f} d)'})
    return pd.DataFrame(rows), pd.DataFrame(extras)


def selftest():
    m = match([10.0, 100.0, 105.0], [102.0, 9.2, 50.0])
    assert m == {0: 1, 1: 0}, m                 # 105 loses 102 to 100 (nearer); 50 matches nothing
    assert match([10.0], [11.5]) == {}, 'a 15% offset is not a match'
    assert match([10.0], [10.9]) == {0: 0}
    print('recovery self-test passed')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--results-dir', default='results')
    ap.add_argument('--out', help='CSV (default fits_rv/<results-dir>/recovery.csv)')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest()
        sys.exit()
    rows, extras = table(a.results_dir)
    out = a.out or os.path.join(HERE, a.results_dir, 'recovery.csv')
    rows.to_csv(out, index=False)
    extras.to_csv(out.replace('.csv', '_extras.csv'), index=False)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_rows', 200)
    pd.set_option('display.max_colwidth', 48)
    print(rows.round(3).to_string(index=False))
    lit = rows[rows.P_lit_d.notna()]
    ok = lit[lit.range == 'searchable']
    print(f'\nliterature companions: {len(lit)}; in the searchable period range: {len(ok)}; '
          f'recovered (P within {RECOVERY_TOL:.0%}): {int(ok.recovered.sum())} of {len(ok)}, '
          f'and {int(lit[lit.range != "searchable"].recovered.sum())} outside that range')
    print(f'extra model companions matching no literature companion: {len(extras)}')
    print('wrote', out)
