"""Presentation figure of one benchmark dataset: a slide-ready summary, made
from saved products only (no fitting or searching).

Panels: (a) the RVs and the reference RV fit that sets the injected planets;
(b) the blind search, step by step: each step's Delta-BIC periodogram, the
simulated threshold, the decision; (c) each planet's orbit on the sky,
injected vs recovered; (d) the star's total motion over DR4 from all
companions, injected vs the search's final model; (e) known vs recovered
(star's astrometric semi-major axis vs period); (f) a table of injected vs
recovered elements; (g, h) the along-scan data the search fits, its final
model and the residuals. The diagnostic figures (synthetic_dr4.process_figure,
benchmark.summary_figure, the step figures) carry the details.

Run:        python -I validation/presentation.py TOI1736_seed1 [out_dir]
            (default out_dir: tmp_figures/ at the repository root, not in git)
Test mode:  python -I validation/presentation.py --test
"""

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import benchmark as bm          # noqa: E402  (sets the paths of core, injection, fits_rv)
import companion_search as cs   # noqa: E402
import figtools                 # noqa: E402
import orbit                    # noqa: E402
import search                   # noqa: E402
import synthetic_dr4 as sd      # noqa: E402

OUT_DIR = os.path.join(ROOT, 'tmp_figures')
YR = lambda t_mjd: 2000 + (np.asarray(t_mjd) - 51544.5)/365.25      # noqa: E731


def _minus_line(t_fit, xy_fit, t, xy):
    """xy(t) minus the straight line fitted to xy_fit(t_fit) (least squares,
    each coordinate): removes motion that looks like proper motion."""
    A = np.column_stack([np.ones_like(t_fit), t_fit - t_fit.mean()])
    out = []
    for vf, v in zip(xy_fit, xy):
        c = np.linalg.lstsq(A, vf, rcond=None)[0]
        out.append(v - c[0] - c[1]*(t - t_fit.mean()))
    return out


def recovered_model(stack, s):
    """The search's final model refit to the data: per recovered companion
    its (P, e, M0) and Thiele-Innes constants; the trend coefficients (sky
    terms in dt^2, dt^3 about the channel mean epoch, dt in years)."""
    st = stack.with_trend(s.get('trend', 0))
    known = [(f['P_d'], f['e'], math.radians(f['M0_deg'])) for f in s['found']]
    coef = st.fit_free_jitter(st.design(known))[0]
    names = [n.split(':')[1] for n in st.baseline_names]
    trend = {n: float(coef[i]) for i, n in enumerate(names) if n.startswith(('acc', 'jerk'))}
    comps = [(P, e, M0, coef[st.k_base + 4*i:st.k_base + 4*i + 4]) for i, (P, e, M0) in enumerate(known)]
    return st, comps, trend


equal_with_headroom = figtools.equal_with_headroom     # shared helper (blind_search/figtools.py)


def sky_of(t, P, e, M0, ti, tref):
    A, B, F, G = ti
    X, Y = orbit.xy(orbit.mean_anomaly(t, P, M0, tref), e)
    return B*X + G*Y, A*X + F*Y


def presentation_figure(stem, out_dir=OUT_DIR):
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    s = json.load(open(os.path.join(bm.RESULTS, stem, 'summary.json')))
    stack, truth = bm.load(stem)
    comps = truth['companions']
    targets = [c for c in comps if c.get('role') != 'trend']
    tcs = [c for c in comps if c.get('role') == 'trend']
    rows = {r['name']: r for r in s['companions']}
    colors = {c['name']: f'C{[0, 3, 4, 5, 6][k % 5]}' for k, c in enumerate(targets)}
    fs = 17
    legs = []

    fig = plt.figure(figsize=(26, 38))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1.15, 1, 0.6], hspace=0.3, wspace=0.22)

    # (a) RVs and the reference fit
    fit = sd.fit_from_truth(truth)
    rv, th, n, tord = fit['rv'], fit['theta'], fit['npl'], fit['tord']
    _, off, _ = cs.neg_lnL(th, rv, n, tord, True)
    jit = cs.unpack(th, n, tord, rv.nds)[2]
    ax = fig.add_subplot(gs[0, 0])
    for d, nm in enumerate(rv.names):
        m = rv.ds == d
        ax.errorbar(YR(rv.t[m] - 2400000.5), rv.rv[m] - off[d], np.hypot(rv.err[m], jit[d]), fmt='o', ms=5,
                    lw=1, alpha=0.75, color=f'C{d % 10}', label=f"{nm.split(':')[1] if ':' in nm else nm} ({m.sum()})")
    tt = np.linspace(rv.t.min() - 100, rv.t.max() + 100, 20000)
    ax.plot(YR(tt - 2400000.5), cs.model_rv(th, sd.grid_for(tt, rv), n, tord), 'k-', lw=1.3, zorder=5,
            label='reference fit (sets the injected planets)')
    ax.set_xlabel('year')
    ax.set_ylabel('RV (m/s)')
    trend_txt = ['', ', linear trend', ', quadratic trend'][tord]
    ax.set_title(f"(a) RVs and the reference fit: {n} planet{'' if n == 1 else 's'}{trend_txt}", fontsize=fs + 3)
    legs.append((ax, dict(fontsize=fs - 3, ncol=2)))

    # (b) the search, step by step
    k = 1
    steps = []
    while os.path.exists(os.path.join(bm.RESULTS, stem, f'step{k}.npz')):
        steps.append(search.load_step(os.path.join(bm.RESULTS, stem, f'step{k}.npz')))
        k += 1
    sub = gs[0, 1].subgridspec(len(steps), 1, hspace=0.12)
    for i, (dec, meta, a, sim) in enumerate(steps):
        ax = fig.add_subplot(sub[i])
        per, dbic = a['periods'], a['dbic']
        ax.plot(per, dbic, color='C0', lw=1.8, label=r'$\Delta$BIC of adding one companion')
        ax.axhline(dec['thr'], color='C1', ls='--', lw=2.2, label='detection threshold (FAP 0.1%)')
        j = int(np.argmax(dbic))
        hit = dec['detected']
        ax.plot(per[j], dbic[j], 'o', ms=13, color='C3' if hit else '0.4', zorder=5,
                label='highest peak' if i == 0 else None)
        for p in meta.get('truth', []):        # short ticks at the bottom: legends stay free above
            ax.plot([p, p], [0, 0.12], color='C2', lw=4, transform=ax.get_xaxis_transform(), zorder=6,
                    label='injected period' if i == 0 and p == meta['truth'][0] else None)
        ax.set_xscale('log')
        ax.set_xlim(meta['pmin'], meta['pmax'])
        ax.set_ylim(max(dbic.min(), -25) - 3, max(dbic.max(), dec['thr'])*1.35 + 3)
        tr = dec.get('trend_before') or 0
        nk = len(a['known_before'])
        model_txt = (f"{nk} companion{'' if nk == 1 else 's'}"
                     f"{['', ' + acceleration', ' + acceleration + jerk'][tr]}")
        P_new = a['known_after'][-1][0] if hit and a['known_after'] else dec['P']   # after the refit
        verdict = (f"peak {dbic[j]:.0f} > threshold {dec['thr']:.0f}: companion added at "
                   f"{P_new:.0f} d" if hit else
                   f"peak {dbic[j]:.0f} < threshold {dec['thr']:.0f}: stop")
        ax.text(0.01, 0.95, f"step {dec['step']} (model: {model_txt})\n{verdict}", transform=ax.transAxes,
                va='top', fontsize=fs - 1, bbox=dict(fc='w', ec='0.7', alpha=0.9))
        ax.set_ylabel(r'$\Delta$BIC')
        if i == 0:
            ax.set_title('(b) Blind search of the DR4 astrometry, step by step\n'
                         f"before step 1: {bm.trend_summary(dict(trend_log=[x for x in s.get('trend_log', []) if 'adopted' in x], trend=steps[0][0].get('trend_before') or 0)).rsplit('; final', 1)[0] or 'no trend test'}",
                         fontsize=fs + 3)
        if i < len(steps) - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel('trial period (d)')
        if i == 0:                             # one legend for all steps
            legs.append((ax, dict(fontsize=fs - 4, grow=8)))

    # (c) each planet's orbit on the sky: injected vs recovered
    st, rec, trend = recovered_model(stack, s)
    ax = fig.add_subplot(gs[1, 0])
    found_by_P = {round(f['P_d'], 6): f for f in s['found']}
    for c in targets:
        r = rows.get(c['name'], {})
        P = c['rv_based']['P_d']
        tf = np.linspace(0, P, 3000) + stack.t.min()
        x, y = sd.companion_sky(tf, c)
        a_s = c['derived']['a_star_mas']
        size = f'{a_s:.3g} mas' if a_s >= 0.1 else f'{a_s*1e3:.2g} $\\mu$as'
        ax.plot(x, y, '--', color=colors[c['name']], lw=2.5,
                label=f"{sd._pretty(c['name'])}: injected (P = {P:.4g} d, a = {size})")
        f = r.get('found')
        if f:
            P2, e2, M0 = f['P_d'], f['e'], math.radians(f['M0_deg'])
            ti = next(t_ for (p_, e_, m_, t_) in rec if abs(p_ - P2) < 1e-6*P2)
            tf = np.linspace(0, P2, 3000) + stack.t.min()
            ax.plot(*sky_of(tf, P2, e2, M0, ti, stack.tref), '-', color=colors[c['name']], lw=2,
                    label=f"{sd._pretty(c['name'])}: recovered (P = {P2:.4g} d)")
    ax.plot(0, 0, '+', color='k', ms=18, mew=2.5, label='barycentre')
    equal_with_headroom(ax)
    ax.set_xlabel(r'$\Delta\alpha*$ (mas)')
    ax.set_ylabel(r'$\Delta\delta$ (mas)')
    ax.set_title("(c) Each planet's orbit on the sky (star's reflex motion)", fontsize=fs + 3)
    legs.append((ax, dict(fontsize=fs - 3)))

    # (d) the star's total motion over DR4: injected vs the search's model
    ax = fig.add_subplot(gs[1, 1])
    t_tr = np.sort(stack.t)
    tg = np.linspace(t_tr.min(), t_tr.max(), 4000)
    xi, yi = np.zeros_like(tg), np.zeros_like(tg)
    for c in comps:
        x, y = sd.companion_sky(tg, c, full=True)
        if c.get('role') == 'trend':          # its linear motion looks like proper motion
            x, y = _minus_line(t_tr, sd.companion_sky(t_tr, c, full=True), tg, (x, y))
        xi, yi = xi + x, yi + y
    xr, yr = np.zeros_like(tg), np.zeros_like(tg)
    for P, e, M0, ti in rec:
        x, y = sky_of(tg, P, e, M0, ti, stack.tref)
        xr, yr = xr + x, yr + y
    if trend:
        dt = lambda t: (t - stack.t.mean())/365.25                        # noqa: E731
        poly = lambda t, kk: sum(trend.get(f'{nm}_{kk}', 0.0)*dt(t)**p    # noqa: E731
                                 for nm, p in (('acc', 2), ('jerk', 3)))
        xr_t, yr_t = _minus_line(t_tr, (poly(t_tr, 'ra'), poly(t_tr, 'dec')), tg, (poly(tg, 'ra'), poly(tg, 'dec')))
        xr, yr = xr + xr_t, yr + yr_t
    ax.plot(xi, yi, '-', color='0.6', lw=6, label='injected: all companions' +
            (' (trend companion: minus its linear motion)' if tcs else ''))
    ax.plot(xr, yr, '-', color='k', lw=1.6, label='recovered: the search\'s final model' +
            (' (planets + acceleration' + (' + jerk)' if 'jerk_ra' in trend else ')') if trend else ''))
    ax.plot(0, 0, '+', color='k', ms=18, mew=2.5, label='barycentre')
    equal_with_headroom(ax)
    ax.set_xlabel(r'$\Delta\alpha*$ (mas)')
    ax.set_ylabel(r'$\Delta\delta$ (mas)')
    ax.set_title("(d) The star's motion over DR4 from all companions", fontsize=fs + 3)
    legs.append((ax, dict(fontsize=fs - 3)))

    # (e) known vs recovered
    ax = fig.add_subplot(gs[2, 0])
    pmin, pmax = s['settings']['pmin_d'], 2*stack.span
    sig = float(np.median(stack.err))
    ax.axvspan(1e-3, pmin, color='0.88', zorder=0, label='periods not searched', gid='guide')
    ax.axvspan(pmax, 1e7, color='0.88', zorder=0, gid='guide')
    ax.axhline(sig, color='C2', ls='--', lw=2, label=f'per-transit error ({sig*1e3:.0f} $\\mu$as)')
    allP, alla = [pmin, pmax], [sig]
    for c in targets:
        r = rows[c['name']]
        f = r['found']
        P, a_s = r['P_true'], r['a_star_true_mas']
        ax.plot(P, a_s, 'o', ms=22, mfc='none', mec=colors[c['name']], mew=3, zorder=4,
                label=f"{sd._pretty(c['name'])}: injected")
        if f:
            ax.plot(f['P_d'], f['a_star_mas'], '*', ms=26, color=colors[c['name']], zorder=5,
                    label=f"{sd._pretty(c['name'])}: recovered")
            allP.append(f['P_d'])
            alla.append(f['a_star_mas'])
        allP.append(P)
        alla.append(a_s)
    for f in s['spurious']:
        ax.plot(f['P_d'], f['a_star_mas'], 'X', ms=20, color='k', zorder=5, label='spurious detection')
        allP.append(f['P_d'])
        alla.append(f['a_star_mas'])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(min(allP)/2.5, max(allP)*2.5)
    ax.set_ylim(min(alla)/10, max(alla)*10)
    ax.set_xlabel('period (d)')
    ax.set_ylabel("star's astrometric semi-major axis (mas)")
    nf = sum(rows[c['name']]['found'] is not None for c in targets)
    ax.set_title(f"(e) Known vs recovered: {nf} of {len(targets)} planets, {len(s['spurious'])} spurious",
                 fontsize=fs + 3)
    legs.append((ax, dict(fontsize=fs - 3, markerscale=0.5)))

    # (f) table
    ax = fig.add_subplot(gs[2, 1])
    ax.axis('off')
    head = ['', 'P (d)', 'e', "star's a (mas)", 'i (deg)']
    cells, cols = [], []
    for c in targets:
        r, rb, dv = rows[c['name']], c['rv_based'], c['derived']
        cells.append([f"{sd._pretty(c['name'])} injected", f"{rb['P_d']:.1f}", f"{rb['e']:.2f}",
                      f"{dv['a_star_mas']:.3g}", f"{c['random']['inc_deg']:.0f}"])
        cols.append(['w']*5)
        f = r['found']
        if f:
            cells.append(['   recovered', f"{f['P_d']:.1f} ({f['P_d']/rb['P_d'] - 1:+.1%})", f"{f['e']:.2f}",
                          f"{f['a_star_mas']:.3g} ({f['a_star_mas']/dv['a_star_mas'] - 1:+.0%})",
                          f"{f['inc_deg']:.0f}"])
            cols.append(['#dbe9f6']*5)
        else:
            cells.append(['   not recovered', '', '', '', ''])
            cols.append(['#f8dcdc']*5)
    tab = ax.table(cellText=cells, colLabels=head, cellColours=cols, loc='upper center', cellLoc='center',
                   colWidths=[0.3, 0.2, 0.1, 0.24, 0.1])
    tab.auto_set_font_size(False)
    tab.set_fontsize(fs)
    tab.scale(1, 3.0)
    note = ''
    if tcs:
        tc = tcs[0]
        note = (f"Trend: the RVs' {['', 'linear', 'quadratic'][tord]} trend was injected as a distant companion "
                f"(P = {tc['rv_based']['P_d']/365.25:.0f} yr, {tc['derived']['mass_mj']:.0f} M$_J$; not a target).\n"
                f"Search: {bm.trend_summary(s)}.")
    ax.set_title('(f) Injected vs recovered' + (f'\n{note}' if note else ''), fontsize=fs + 1)

    # (g, h) the data the search fits (along scan, one point per transit), its
    # final model and the residuals (benchmark.al_panels)
    bm.al_panels(fig, gs[3, 0], gs[3, 1], stack, s, legs)
    for a_, kw in legs:
        figtools.legend_clear(a_, **kw)
    figtools.check_legends(fig)
    fig.suptitle(f"{s['label']}: blind search of synthetic Gaia DR4 astrometry (seed {s['seed']})",
                 fontsize=fs + 12, y=0.925)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f'{stem}_presentation.png')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return out


if __name__ == '__main__':
    if '--test' in sys.argv:
        import tempfile
        stem = bm.datasets(bm.config())[0][3]
        with tempfile.TemporaryDirectory() as d:
            out = presentation_figure(stem, d)
            assert os.path.getsize(out) > 0
        print(f'presentation.py test passed ({stem})')
    else:
        args = [a for a in sys.argv[1:] if not a.startswith('-')]
        print('wrote', presentation_figure(args[0], *(args[1:2] or [OUT_DIR])))
