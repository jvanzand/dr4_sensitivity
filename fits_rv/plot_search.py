"""Plot maximum-likelihood fits from companion_search.py for one system.

Each column is one step of the search (a model): the RVs with every dataset's fitted offset
removed, the model curve (Keplerians plus trend), the trend component alone, and the
residuals below. Titles give the step's action, its BIC, and its Delta-BIC against the
model before it (positive = better).

Usage (project root, gaia_dr4_sensitivity environment):
    python -I fits_rv/plot_search.py HD111232 [--steps 2,3,4 | --compare-from 1]
Reads fits_rv/<results-dir>/<display name>/companion_search.json (default results), writes search_fits.png or
trend_vs_planet.png there.
"""

import argparse
import json
import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import companion_search as cs
import names


def place_legend(ax, anchor_ax=None, fontsize=11, markerscale=0.6, below=-0.2, ncol_below=2):
    """Legend at the corner or edge of `ax` where it covers the fewest data points; if every
    spot covers data, put it below `anchor_ax` (default `ax`) instead. Call after the layout
    is final (fig.tight_layout()), since the test uses the rendered sizes. Data points are
    markers (counted fully), annotations (their corners), and the points of plotted curves
    (counted lightly)."""
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return None
    fig = ax.figure
    rend = fig.canvas.get_renderer()
    inv = ax.transAxes.inverted()
    marks, curves = [], []
    for ln in ax.lines:
        xy = np.asarray(ln.get_xydata(), float)
        if len(xy) < 3 and ln.get_marker() in ('None', None, ''):
            if len(xy) == 2 and xy[0, 0] == xy[1, 0]:     # axvline: x in data, y in axes
                xa = inv.transform(ax.transData.transform([[xy[0, 0], ax.get_ylim()[0]]]))[0, 0]
                curves.append(np.column_stack([np.full(50, xa), np.linspace(0, 1, 50)]))
            continue
        pts = inv.transform(ax.transData.transform(xy[np.isfinite(xy).all(axis=1)]))
        (marks if ln.get_marker() not in ('None', None, '') else curves).append(pts)
    for t in ax.texts:
        bb = t.get_window_extent(rend)
        marks.append(inv.transform(np.array([[bb.x0, bb.y0], [bb.x1, bb.y0], [bb.x0, bb.y1],
                                             [bb.x1, bb.y1], [(bb.x0 + bb.x1)/2, (bb.y0 + bb.y1)/2]])))
    marks = np.vstack(marks) if marks else np.zeros((0, 2))
    curves = np.vstack(curves) if curves else np.zeros((0, 2))

    def covered(box, pts):
        (x0, y0), (x1, y1) = box
        pad = 0.01
        return int(np.sum((pts[:, 0] > x0 - pad) & (pts[:, 0] < x1 + pad)
                          & (pts[:, 1] > y0 - pad) & (pts[:, 1] < y1 + pad)))
    best = None
    for loc in ('upper right', 'upper left', 'lower right', 'lower left', 'upper center',
                'lower center', 'center right', 'center left', 'center'):
        leg = ax.legend(handles, labels, loc=loc, fontsize=fontsize, markerscale=markerscale,
                        framealpha=0.9)
        box = inv.transform(leg.get_window_extent(rend).get_points())
        score = (covered(box, marks), covered(box, curves)/max(1, len(curves)))
        leg.remove()
        if best is None or score < best[0]:
            best = (score, loc)
    if best[0][0] == 0 and best[0][1] <= 0.03:   # no data and little of any curve covered
        return ax.legend(handles, labels, loc=best[1], fontsize=fontsize,
                         markerscale=markerscale, framealpha=0.9)
    host = anchor_ax or ax
    return host.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, below),
                       ncol=min(ncol_below, len(handles)), fontsize=fontsize,
                       markerscale=markerscale)


def draw(data, grid, colors, ax, axr, th, npl, tord, title, legend):
    """One model on a pair of axes: data minus offsets with the model curve, and residuals."""
    import numpy as np
    _, off, resid = cs.neg_lnL(th, data, npl, tord, True)
    model = cs.model_rv(th, grid, npl, tord)
    coefs = cs.unpack(th, npl, tord, data.nds)[1]
    for d, name in enumerate(data.names):
        idx = data.members[d]
        ax.errorbar(data.t[idx] - 2400000.5, data.rv[idx] - off[d], data.err[idx], fmt='o', ms=4,
                    color=colors[d % len(colors)], alpha=0.8, label=f'{name} (n={len(idx)})')
        axr.errorbar(data.t[idx] - 2400000.5, resid[idx], data.err[idx], fmt='o', ms=4,
                     color=colors[d % len(colors)], alpha=0.8)
    ax.plot(grid.t - 2400000.5, model, 'k-', lw=1.5, label='maximum-likelihood model')
    if tord:
        ax.plot(grid.t - 2400000.5, cs.trend_model(coefs, grid), 'r--', lw=1.5,
                label=f'trend component (order {tord})')
    ax.set_title(title, fontsize=13)
    ax.set_ylabel('RV minus dataset offset [m/s]')
    axr.axhline(0, c='0.5', lw=0.8)
    axr.set_ylabel('residual [m/s]')
    axr.set_xlabel('time [MJD]')
    if legend:      # an in-axes legend; callers wanting no overlap use place_legend afterwards
        ax.legend(fontsize=9, loc='best')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('system', help='display name (HD111232) or Gaia DR3 source_id')
    ap.add_argument('--steps', help='comma-separated step numbers (default: the last forced '
                                    'step and the next two)')
    ap.add_argument('--compare-from', type=int,
                    help='instead, refit from this step and compare: as is, + linear trend, '
                         '+ quadratic trend, + one more companion (the trend-or-planet question)')
    ap.add_argument('--results-dir', default='results', help='folder under fits_rv/')
    ap.add_argument('--out')
    a = ap.parse_args()

    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt

    sid = names.system_id(a.system)
    sys_name = names.display_name(sid)
    folder = os.path.join(HERE, a.results_dir, sys_name)
    with open(os.path.join(folder, 'companion_search.json')) as fh:
        res = json.load(fh)
    steps = res['steps']
    lit_pri = res.get('literature_priors', [])    # priors the search used (JSON lists)
    data = cs.load_system(os.path.join(ROOT, res['rv_table']) if not os.path.isabs(res['rv_table'])
                          else res['rv_table'])
    t_dense = np.linspace(data.t.min() - 0.03*data.span, data.t.max() + 0.03*data.span, 6000)
    grid = SimpleNamespace(t=t_dense, tref=data.tref, span=data.span, nds=data.nds)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    panels = []     # (theta, npl, tord, title)
    if a.compare_from is not None:
        s = steps[a.compare_from]
        th0, npl, t0 = np.array(s['theta']), s['n_planets'], s['trend_order']
        assert t0 == 0, 'compare-from expects a step without a trend'
        floors = [cs.P_MIN_KNOWN_D]*npl
        base_bic = cs.bic(s['lnL'], npl, 0, data)
        models = [('as is', th0, s['lnL'], npl, 0)]
        th_prev = th0
        for order, label in ((1, '+ linear trend'), (2, '+ quadratic trend')):
            ## each order starts from the previous order's solution (the lower order is
            ## nested in the higher one, so the fit can only improve)
            th, ll = cs.fit(cs.seed_trend(th_prev, npl, order - 1, order, data), data, npl,
                            order, floors, lit_pri)
            th_prev = th
            models.append((label, th, ll, npl, order))
        th, ll = cs.add_planet(th0, npl, 0, data, data.p_min_blind, lit_pri)
        models.append((f'+ one more companion (period cap {cs.PERIOD_CAP_BASELINES:g}x baseline)', th, ll, npl + 1, 0))
        for label, th, ll, n, o in models:
            b = cs.bic(ll, n, o, data)
            per = ', '.join(f'{np.exp(p[0]):.4g} d' for p in cs.unpack(th, n, o, data.nds)[0])
            panels.append((th, n, o, f'{label}\nBIC = {b:.0f}, dBIC vs. "as is" = {base_bic - b:.1f}'
                           f'\ncompanion periods: {per or "none"}; trend order {o}'))
    else:
        if a.steps:
            chosen = [int(x) for x in a.steps.split(',')]
        else:
            lf = max([s['step'] for s in steps if s['forced']], default=0)
            chosen = [s['step'] for s in steps if lf <= s['step'] <= lf + 2]
        for k in chosen:
            s = steps[k]
            th = np.array(s['theta'])
            pl = cs.unpack(th, s['n_planets'], s['trend_order'], data.nds)[0]
            per = ', '.join(f'{np.exp(p[0]):.4g} d' for p in pl)
            verdict = 'forced' if s['forced'] else ('adopted' if s['adopted'] else 'rejected')
            panels.append((th, s['n_planets'], s['trend_order'],
                           f'step {k}: {s["action"]} ({verdict})\nBIC = {s["BIC"]:.0f}, dBIC vs. '
                           f'previous model = {s["dBIC"]:.1f}\ncompanion periods: '
                           f'{per or "none"}; trend order {s["trend_order"]}'))

    fig, axes = plt.subplots(2, len(panels), figsize=(7*len(panels), 10), squeeze=False,
                             gridspec_kw={'height_ratios': [3, 1]}, sharex='col')
    for col, (th, n, o, title) in enumerate(panels):
        draw(data, grid, colors, axes[0, col], axes[1, col], th, n, o, title, False)
    fig.suptitle(f'{sys_name} (Gaia DR3 {sid}): maximum-likelihood fits, N_RV = {data.n}, '
                 f'baseline {data.span:.0f} d (period cap {data.pmax:.0f} d)', fontsize=16)
    fig.tight_layout()
    place_legend(axes[0, 0], anchor_ax=axes[1, 0], below=-0.5)
    out = a.out or os.path.join(folder, 'trend_vs_planet.png' if a.compare_from is not None
                                else 'search_fits.png')
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print('wrote', out)


if __name__ == '__main__':
    main()
