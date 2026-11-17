"""Diagnostic figure of one search step: `plot_step_sim` draws the real
periodogram with example noise-sampled periodograms, the simulated threshold,
how the threshold is built, and the residuals. `search.run_search(out_dir=...)`
writes one per step.

Test mode (writes a step figure to /tmp):  python -I blind_search/plots.py
--test
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import figtools  # noqa: E402


def _plt():
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    return plt


def _residual_panels(fig, cell, stack, residual_models):
    """Normalized residuals of each model in `residual_models` (formal errors
    only, so leftover signal shows), one panel per model, inside the grid cell
    `cell`."""
    sub = cell.subgridspec(len(residual_models), 1, hspace=0.55)
    plt = _plt()
    z0 = np.zeros(len(stack.channels))
    for q, pl in enumerate(residual_models):
        ax3 = fig.add_subplot(sub[q])
        coef, r, chi2 = stack.fit(pl, z0)
        z = r/stack.sigma(z0)
        dof = stack.n - stack.n_params(len(pl))
        for i, ch in enumerate(stack.channels):
            m = stack.chan == i
            ax3.plot(stack.t[m] - stack.tref, z[m], '.', ms=9, color=f'C{i % 10}',
                     label=f'{ch.name} ({m.sum()} points)' if (q == 0 and i < 8) else None)
        ax3.axhline(0, color='0.5', lw=1)
        ax3.set_ylabel('residual / error')
        ax3.set_title(f'residuals to the {len(pl)}-companion model: reduced $\\chi^2$ = $\\chi^2$/dof = '
                      f'{chi2/dof:.2f}  ($\\chi^2$ = {chi2:.0f}, dof = N $-$ {stack.n_params(len(pl))} '
                      f'= {dof}; formal errors)', fontsize=12)
        if q == len(residual_models) - 1:
            ax3.set_xlabel('time - mid-baseline (d, MJD)')
        if q == 0:
            figtools.legend_clear(ax3, fontsize=10, frameon=False)


def redraw_step(stack, path, out):
    """Redraw one search step's figure from its saved products
    (`search.save_step`), without searching again. `stack` is the same data
    the search ran on."""
    import search
    dec, meta, a, sim = search.load_step(path)
    if sim is None:
        raise ValueError(f'{path}: no simulated threshold saved (method {dec["method"]!r})')
    known = a['known_before']
    models = [known] + ([a['known_after']] if dec['detected'] else [])
    tau = dec.get('trend_before') or 0
    stack = stack.with_trend(tau)            # the step's trend order (residuals of both models drawn at it)
    return plot_step_sim(stack, a['periods'], a['dbic'], sim, known, meta['pmin'], meta['pmax'], out,
                         search.step_title(meta['label'], dec['step'], len(known), tau,
                                           'rv' if stack.has_rv and not stack.has_astro else 'astro'),
                         truth=meta['truth'], truth_label=meta['truth_label'],
                         residual_models=models, fap=meta['fap'])


def plot_step_sim(stack, periods, dbic, sim, planets, pmin, pmax, out, title, truth=(), truth_label=None,
                  residual_models=None, fap=1e-3, nshow=20, stat=None):
    """One search step with the simulated threshold
    (simthreshold.simulated_threshold): top, the real data's Delta-BIC
    periodogram, `nshow` of the simulated no-further-companion periodograms
    (grey) with the highest peak of each marked, and the threshold; bottom left,
    how the threshold is built (histogram of the highest peak of every simulated
    dataset, the Gumbel fit and the thresholds, with the real data's highest
    peak); bottom right, the residuals of the models in `residual_models`. `sim`
    is the dict returned by simulated_threshold(..., keep=nshow); `planets` are
    the companions already in the model. `stat` is the statistic compared with
    the threshold when it is not the periodogram's highest peak (rv_search.py:
    the larger of that peak and the trend gain)."""
    from scipy.stats import gumbel_r
    plt = _plt()
    residual_models = [planets] if residual_models is None else residual_models
    maxima, thr = sim['maxima'], sim['thr'][fap]
    M, ex = len(maxima), sim['examples'][:nshow]
    k = int(np.argmax(dbic))
    peak = float(dbic[k])
    stat = peak if stat is None else float(stat)
    unit = 'm/s' if stack.has_rv and not stack.has_astro else 'mas'
    jit = float(np.sqrt(sim['jitter2'][0]))
    fig = plt.figure(figsize=(20, 15.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1], hspace=0.62)
    ax = fig.add_subplot(gs[0, :])
    for i, d in enumerate(ex):
        ax.plot(periods, d, color='0.72', lw=0.9, alpha=0.9, zorder=1,
                label=f'$\\Delta$BIC periodogram of {len(ex)} of the {M} noise-sampled datasets (no further companion)' if i == 0 else None)
        j = int(np.argmax(d))
        ax.plot(periods[j], d[j], 'o', color='0.35', ms=5.5, zorder=3,
                label='highest peak of each of those noise-sampled periodograms (the number recorded)' if i == 0 else None)
    ax.plot(periods, dbic, color='C0', lw=2, zorder=4, label='$\\Delta$BIC of adding one circular companion: real data')
    ax.plot(periods[k], peak, 'o', color='C3', ms=12, zorder=5,
            label=f'highest peak of the real data: P = {periods[k]:.4g} d, $\\Delta$BIC = {peak:.1f}')
    ax.axhline(thr, color='C1', ls='--', lw=2.2, zorder=2,
               label=f'simulated threshold (FAP = {fap:g}): {thr:.1f}; only the highest peak is compared with it')
    ax.axhline(0, color='0.5', lw=1)
    for i, p in enumerate(planets):
        ax.axvline(p[0], color='C4', ls='-.', lw=1.5,
                   label='period of a companion already in the model' if i == 0 else None)
    for i, p in enumerate(truth):
        if isinstance(p, (tuple, list)):                     # (period, label): one coloured line per known signal
            p, lab, col = p[0], p[1], ['C2', 'C9', 'C8'][i % 3]
        else:
            lab, col = ((truth_label or 'known (injected or literature) period') if i == 0 else None), 'C2'
        ax.axvline(p, color=col, ls='-', lw=3, alpha=0.4, label=lab)
    lo = min(dbic.min(), min(d.min() for d in ex)) - 3 if ex else dbic.min() - 3
    if peak > 3*max(thr, 10):
        ax.set_yscale('symlog', linthresh=max(30.0, 2*thr))
        ax.set_ylabel('$\\Delta$BIC (linear below %.0f, log above)' % max(30.0, 2*thr), fontsize=18)
    else:
        ax.set_ylabel('$\\Delta$BIC')
    ax.set_ylim(lo, max(peak, thr)*1.25 if peak > 0 else max(thr, 5)*1.4)
    ax.set_xscale('log')
    ax.set_xlim(pmin, pmax)
    ax.set_xlabel('trial period (d)')
    nge = int(np.sum(maxima >= stat))
    verdict = ('DETECTED' if stat > thr else 'NOT DETECTED')
    ax.set_title(f'{title}\n{verdict}: {"highest peak" if stat == peak else "trend gain"} $\\Delta$BIC = {stat:.1f} {">" if stat > thr else "<"} threshold {thr:.1f}; '
                 f'{nge} of {M} noise-sampled datasets have a highest peak at or above the real one', fontsize=15)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=12, frameon=True, framealpha=0.95)

    ## how the threshold is built
    ax2 = fig.add_subplot(gs[1, 0])
    top = max(maxima.max() + 3, thr*1.3)
    bins = np.linspace(maxima.min() - 1, maxima.max() + 1, 22)
    ax2.hist(maxima, bins=bins, color='0.75', edgecolor='k', lw=1,
             label=f'highest peak of each noise-sampled dataset ({M})')
    x = np.linspace(bins[0], top, 300)
    ax2.plot(x, M*(bins[1] - bins[0])*gumbel_r.pdf(x, sim['loc'], sim['scale']), color='C5', lw=2.5,
             label='Gumbel (extreme-value) fit to those peaks')
    ytop = ax2.get_ylim()[1]*1.45
    ax2.set_ylim(0, ytop)
    for f_, ls in ((1e-2, ':'), (1e-3, '--')):
        if f_ in sim['thr']:
            ax2.axvline(sim['thr'][f_], color='C1', ls=ls, lw=2.2)
            ax2.annotate(f"threshold, FAP = {f_:g}: {sim['thr'][f_]:.1f}", xy=(sim['thr'][f_], 0.5*ytop), xytext=(-7, 0),
                         textcoords='offset points', rotation=90, color='C1', ha='right', va='center', fontsize=12,
                         fontweight='bold')
    if stat <= top:
        ax2.axvline(stat, color='C3', ls='-', lw=2.8)
        ax2.annotate(f'highest peak of the real data: {stat:.1f}', xy=(stat, 0.68*ytop), xytext=(7, 0), textcoords='offset points',
                     rotation=90, color='C3', ha='left', va='center', fontsize=12, fontweight='bold')
    else:
        ax2.annotate(f'real data:\n{stat:.0f}, off scale', xy=(top, 0.5*ytop),
                     xytext=(bins[0] + 0.665*(top - bins[0]), 0.80*ytop), fontsize=12, color='C3', fontweight='bold', ha='left',
                     arrowprops=dict(arrowstyle='-|>', color='C3', lw=2.5))
    ax2.set_xlim(bins[0], top)
    ax2.set_xlabel('highest $\\Delta$BIC of one periodogram')
    ax2.set_ylabel('number of noise-sampled datasets')
    ax2.set_title('how the threshold is built: each noise-sampled dataset contributes its highest peak;\nFAP = fraction of noise-sampled datasets whose highest peak exceeds the threshold', fontsize=12.5)

    _residual_panels(fig, gs[1, 1], stack, residual_models)
    figtools.legend_clear(ax2, fontsize=10.5, frameon=False)   # after all panels exist
    sub = (f'Noise-sampled datasets: at every real transit, the predicted value of the fitted {len(planets)}-companion model plus Gaussian '
           f'noise of standard deviation $\\sqrt{{\\sigma_i^2+s^2}}$ (formal error $\\sigma_i$, fitted jitter $s$ = {jit:.3f} {unit}); each is refit '
           f'and searched like the real data; only its highest peak is kept (grey dots).')
    fig.text(0.5, 0.012, sub, ha='center', fontsize=12.5)
    figtools.check_legends(fig)          # no legend may cover data
    fig.savefig(out, dpi=130, bbox_inches='tight')
    plt.close(fig)
    return out


if __name__ == '__main__':
    if '--test' in sys.argv:
        import data
        import model
        import periodogram as pg
        import simthreshold as st
        sys.path.insert(0, os.path.join(ROOT, 'tests'))
        import test_blind_search as T
        g = data.load_gaia_dr4(T.GAIA4)
        s = model.Stack([g])
        T.fill(s, [(300.0, 0.0, 0.0)], [[0.25, -0.15, 0.2, 0.3]], 11)
        per = pg.period_grid(s.span, 20.0, 2*s.span)
        d = pg.delta_bic(s, [], None, per)
        sim = st.simulated_threshold(s, [], M=10, seed=1, periods=per, keep=5)
        plot_step_sim(s, per, d, sim, [], 20.0, 2*s.span,
                      '/tmp/blind_step_test.png', 'plots.py test',
                      truth=[300.0])
        print('wrote /tmp/blind_step_test.png')
