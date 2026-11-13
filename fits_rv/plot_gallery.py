"""Gallery of the companion-search (Delta-BIC) runs: one row per system.

Columns:
  1. BIC at each step of the search (lower is better). Squares are the literature base model
     (companions under Gaussian priors) and trend moves, which are not Delta-BIC gated; filled
     circles are blind additions adopted at Delta-BIC >= 20; the open cross is the first
     rejected step. Labels give each step's action and Delta-BIC against the model before it.
  2. The final model (the last adopted step) on the RVs, with each dataset's offset removed;
     the red dashed line is the trend component, if any.
  3. Periods: literature periods (triangles) against the companions in the final model, marked
     as RV-supported or held only by their prior, or added blindly; the RV baseline and the
     period cap (2x baseline; 3x in the older runs) and the blind period floor are marked. For a
     blind run (no literature companions in the fit) the literature periods are tallied as
     recovered (a model companion within 10%), missed, or outside the searchable range.

The RV search is not a test of recovering every literature companion (transiting planets
with small K cannot be found from RVs alone); what matters for the project is whether Gaia
recovers what the RVs show.

Usage (project root, gaia_dr4_sensitivity environment):
    python -I fits_rv/plot_gallery.py                  # every folder in fits_rv/results
    python -I fits_rv/plot_gallery.py --systems TOI1694 --out /tmp/test.png   # local test
Writes fits_rv/results/search_gallery.png by default.
"""

import argparse
import glob
import json
import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import companion_search as cs
import names
import plot_search as ps
import recovery as rec


def short(action):
    """A compact step label for the BIC panel."""
    a = action.replace('forced companion', 'forced').replace('add companion and drop trend '
                                                              'order', 'add + drop trend')
    return a.replace('raise trend order', 'trend').replace('drop trend term', 'drop trend')


def literature_periods(sid, cat):
    g = cat[cat.system_id == sid]
    return sorted(g.period_d.dropna().tolist())


def bic_panel(ax, steps, title):
    from matplotlib.ticker import NullLocator
    x = [s['step'] for s in steps]
    y = [s['BIC'] for s in steps]
    ax.plot(x, y, '-', color='0.6', lw=1, zorder=1)
    kinds = [('forced' if s['forced'] else ('adopted' if s['adopted'] else 'rejected'))
             for s in steps]
    style = {'forced': dict(marker='s', mfc='C0', mec='C0', label='literature base model, or a trend move'),
             'adopted': dict(marker='o', mfc='C2', mec='C2', label=f'adopted (Delta-BIC >= {cs.DELTA_BIC:g})'),
             'rejected': dict(marker='X', mfc='none', mec='C3', label=f'rejected (Delta-BIC < {cs.DELTA_BIC:g})')}
    seen = set()
    for s, k, xi, yi in zip(steps, kinds, x, y):
        if s['step'] == 0:
            k = 'adopted'
            lab = 'baseline: offsets and jitter only'
            ax.plot([xi], [yi], marker='D', ms=9, mfc='k', mec='k', ls='', label=lab, zorder=3)
            continue
        ax.plot([xi], [yi], ms=10, ls='', zorder=3,
                **{**style[k], **({} if k not in seen else {'label': None})})
        seen.add(k)
        ax.annotate(f'{short(s["action"])}\ndBIC {s["dBIC"]:.1f}', (xi, yi),
                    textcoords='offset points', xytext=(6, 8), fontsize=12, rotation=20)
    ax.set_yscale('log')
    ax.set_xlabel('search step')
    ax.set_ylabel('BIC (lower is better)')
    ax.set_xticks(x)
    ax.xaxis.set_minor_locator(NullLocator())     # discrete steps: no minor ticks
    ax.set_xlim(-0.4, max(x) + 1.3)
    ax.set_title(title, fontsize=14)


def period_panel(ax, steps, final, res, lit, data):
    """Literature periods (triangles) above the periods in the final model (circles, with the
    Delta-BIC of the step that added each). Literature-prior runs: the first N companions of the
    model are the literature ones, drawn filled if the RVs support them (Delta-BIC for dropping
    >= the threshold) and hollow if only the prior holds them; the rest were added blindly.
    Blind runs: every model companion was found blindly; each literature period is a filled
    triangle if a model companion lies within the recovery tolerance, a red open one if the
    search should have been able to find it and did not, and a grey open one if its period lies
    below the blind floor or above the cap (not searchable)."""
    span = data.span
    cap = res.get('period_cap_baselines', 3.0)      # runs before 2026-11-09 used 3x
    pmax = min(cs.P_MAX_YR*365.25, cap*span)
    pl = cs.unpack(np.array(final['theta']), final['n_planets'], final['trend_order'],
                   data.nds)[0]
    blind = not res.get('literature', res['n_known'] > 0)
    allp = [float(np.exp(p[0])) for p in pl]
    added = lambda j: next(x['dBIC'] for x in steps if x['adopted'] and x['n_planets'] > j)
    if blind:
        m = rec.match(lit, allp)
        matched = set(m.values())
        k_lit = {'rec': dict(marker='^', mfc='k', mec='k',
                             label=f'literature period, recovered (model companion within '
                                   f'{rec.RECOVERY_TOL:.0%})'),
                 'miss': dict(marker='^', mfc='none', mec='C3', mew=2,
                              label='literature period, not recovered'),
                 'out': dict(marker='^', mfc='none', mec='0.5', mew=2,
                             label='literature period outside the searchable range')}
        seen_l = set()
        for i, P in enumerate(lit):
            k = ('rec' if i in m else 'out' if (P < data.p_min_blind or P > pmax) else 'miss')
            ax.plot([P], [1], ms=12, ls='', **{**k_lit[k], **({'label': None} if k in seen_l else {})})
            seen_l.add(k)
        kinds = {'match': dict(marker='o', mfc='C2', mec='C2',
                               label='found blindly, matches a literature period'),
                 'extra': dict(marker='o', mfc='none', mec='C2', mew=2,
                               label='found blindly, no literature match')}
        seen = set()
        for j, P in enumerate(allp):
            k = 'match' if j in matched else 'extra'
            ax.plot([P], [0], ms=11, ls='', **{**kinds[k], **({'label': None} if k in seen else {})})
            seen.add(k)
            ax.annotate(f'dBIC {added(j):.0f}', (P, 0), textcoords='offset points', xytext=(0, -26),
                        fontsize=10, ha='center')
    else:
        sup = res['final'].get('literature_fit', [])
        kinds = {'lit_ok': dict(marker='s', mfc='C0', mec='C0',
                                label=f'literature companion, RV-supported (dBIC for dropping '
                                      f'>= {cs.DELTA_BIC:g})'),
                 'lit_weak': dict(marker='s', mfc='none', mec='C3', mew=2,
                                  label=f'literature companion held by its prior (dBIC for '
                                        f'dropping < {cs.DELTA_BIC:g})'),
                 'added': dict(marker='o', mfc='C2', mec='C2',
                               label=f'added blindly by the search (dBIC >= {cs.DELTA_BIC:g})')}
        ax.plot(lit, [1]*len(lit), '^', ms=12, color='k',
                label='literature periods (catalog)' if lit else None)
        seen = set()
        for j, P in enumerate(allp):
            if j < len(sup):
                k, d = ('lit_ok' if sup[j]['supported'] else 'lit_weak'), sup[j]['dBIC_drop']
            else:
                k, d = 'added', added(j)
            ax.plot([P], [0], ms=11, ls='', **{**kinds[k], **({'label': None} if k in seen else {})})
            seen.add(k)
            ax.annotate(f'dBIC {d:.0f}', (P, 0), textcoords='offset points', xytext=(0, -26),
                        fontsize=10, ha='center')
    ax.axvline(span, color='C1', ls='--', label=f'RV baseline ({span:.0f} d)')
    ax.axvline(pmax, color='C3', ls=':', label=f'period cap, {cap:g}x baseline ({pmax:.0f} d)')
    ax.axvline(data.p_min_blind, color='C4', ls='-.',
               label=f'blind period floor, 2x mean RV spacing ({data.p_min_blind:.0f} d)')
    ax.set_xscale('log')
    ax.set_ylim(-0.7, 1.7)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['final model', 'literature'])
    ax.set_xlabel('orbital period [d]')
    ax.set_title('periods: literature vs. final model', fontsize=14)
    ax.set_xlim(0.5*min(allp + lit + [data.p_min_blind]), 3*max(pmax, max(allp + lit)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--systems', help='comma-separated display names (default: all results)')
    ap.add_argument('--results-dir', default='results', help='folder under fits_rv/')
    ap.add_argument('--out', help='default fits_rv/<results-dir>/search_gallery.png')
    a = ap.parse_args()
    a.out = a.out or os.path.join(HERE, a.results_dir, 'search_gallery.png')

    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt

    folders = (sorted(glob.glob(os.path.join(HERE, a.results_dir, '*', 'companion_search.json')))
               if not a.systems else
               [os.path.join(HERE, a.results_dir, s, 'companion_search.json')
                for s in a.systems.split(',')])
    cat = pd.read_csv(os.path.join(ROOT, 'sample', 'stage1_catalog.csv'),
                      dtype={'system_id': str}, low_memory=False)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fig, axes = plt.subplots(len(folders), 3, figsize=(32, 7*len(folders)), squeeze=False,
                             gridspec_kw={'width_ratios': [1.1, 1.5, 0.9]})
    for row, f in enumerate(folders):
        res = json.load(open(f))
        sys_name = os.path.basename(os.path.dirname(f))
        sid = names.system_id(sys_name)
        steps = res['steps']
        data = cs.load_system(res['rv_table'] if os.path.isabs(res['rv_table'])
                              else os.path.join(ROOT, res['rv_table']))
        final = [s for s in steps if s['adopted']][-1]
        lit = literature_periods(sid, cat)
        lit_used = res.get('literature', res['n_known'] > 0)
        nk = res['n_known'] if lit_used else len(lit)
        bic_panel(axes[row, 0], steps,
                  f'{sys_name}: BIC by search step (N_RV = {data.n}, {nk} literature '
                  f'companion{"s" if nk != 1 else ""}' + ('' if lit_used else ', not used') + ')')
        t = np.linspace(data.t.min() - 0.03*data.span, data.t.max() + 0.03*data.span, 6000)
        grid = SimpleNamespace(t=t, tref=data.tref, span=data.span, nds=data.nds)
        pl = cs.unpack(np.array(final['theta']), final['n_planets'], final['trend_order'],
                       data.nds)[0]
        per = ', '.join(f'{np.exp(p[0]):.4g}' for p in pl)
        ## the residual axis is dropped here; plot_search.draw wants a pair of axes
        axr = axes[row, 1].inset_axes([0, -0.0, 1, 0.0001])
        axr.set_visible(False)
        ps.draw(data, grid, colors, axes[row, 1], axr, np.array(final['theta']),
                final['n_planets'], final['trend_order'],
                f'final model (step {final["step"]}): periods [d] {per}; trend order '
                f'{final["trend_order"]}', False)
        axes[row, 1].set_xlabel('time [MJD]')
        for ln in axes[row, 1].lines:       # the model curve is dense for short periods
            if ln.get_linestyle() == '-' and ln.get_color() == 'k':
                ln.set_linewidth(0.6)
        axes[row, 1].title.set_fontsize(14)
        period_panel(axes[row, 2], steps, final, res, lit, data)
    fig.suptitle('Companion-search gallery: maximum-likelihood fits, BIC = -2 ln L + k ln N '
                 '(RV-only N)', fontsize=18, y=1.0)
    fig.tight_layout(h_pad=9)
    for row in range(len(folders)):          # legends last: they need the final layout
        for col, below in ((0, -0.17), (1, -0.2), (2, -0.17)):
            ps.place_legend(axes[row, col], fontsize=12, below=below)
    fig.savefig(a.out, dpi=300, bbox_inches='tight')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
