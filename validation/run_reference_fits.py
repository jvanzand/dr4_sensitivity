"""Run the reference (assisted) RV fit of one or more systems, in parallel.

The reference RV fit (`synthetic_dr4.reference_rv_fit`) is the baseline RV model
of a system: its literature companions (plus the context companions that failed a
census cut) under priors, a trend or long-period Keplerian if the RVs need one,
on the nightly-binned RVs. This script runs it per system and saves, in
`--out` (default validation/results/reference_fits/):
  <system>.log          everything the fit printed
  <system>_fit.pkl      the full fit (load with `load_fit`)
  <system>_summary.json companions (P, K, e, context or sample, Delta-BIC for
                        dropping each), trend order, extra Keplerian, jitters
  <system>_rv.png       RVs, model, long-period part, residuals (with --figures)
  summary.csv           one row per system: status, N_RV, companions, trend

Usage: set SYSTEMS, JOBS, FIGURES and OUT in the `if __name__ == '__main__':` block
at the bottom, then run the module from the project root:
  python -I validation/run_reference_fits.py
  python -I validation/run_reference_fits.py --dry-run    # list what would run (local test mode)
A progress bar counts the finished systems.

Interactive use:
  import sys; sys.path.insert(0, 'validation')
  import run_reference_fits as rf
  fit = rf.fit_system('HD111232')            # one fit, in this process
  fit = rf.load_fit('HD111232')              # a saved fit
  rf.plot_fit(fit, 'HD111232.png')
"""

import contextlib
import json
import os
import pickle
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for d in ('validation', 'fits_rv', 'blind_search'):
    sys.path.insert(0, os.path.join(ROOT, d))
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np  # noqa: E402

OUT_DIR = os.path.join(HERE, 'results', 'reference_fits')
RV_DIR = os.path.join(ROOT, 'rv_data', 'systems')


def all_systems():
    """System keys of every RV table in rv_data/systems/ (e.g. 'TOI1736')."""
    return sorted(f.split('_', 1)[1][:-4].replace('-', '') for f in os.listdir(RV_DIR) if f.endswith('.csv'))


def fit_system(system):
    """The reference RV fit of one system (prints its progress)."""
    import synthetic_dr4 as sd
    return sd.reference_rv_fit(system)


def summarize(fit):
    """A JSON-ready summary of a fit."""
    import rv_model as cs
    th, P = fit['theta'], cs.PARAMS_PER_PLANET
    support = {s['index']: s['dBIC_drop'] for s in fit['support']}
    comps = []
    for j, L in enumerate(fit['lit'] + [{'name': 'long-period Keplerian (extra)'}]*fit['n_extra']):
        k = P*j
        comps.append(dict(name=L['name'], role=('extra' if j >= len(fit['lit']) else
                                                 'context' if L.get('context') else 'sample'),
                          P_d=float(np.exp(th[k])), P_lit_d=L.get('P0'), K_ms=float(np.exp(th[k + 1])),
                          e=float(th[k + 2]**2 + th[k + 3]**2), dBIC_drop=support.get(j)))
    jit = np.exp(th[-fit['rv'].nds:])
    return dict(system=fit['system'], kind=fit['kind'], n_rv=int(fit['rv'].n), trend_order=int(fit['tord']),
                n_extra=int(fit['n_extra']), lnL=float(fit['lnL']), companions=comps,
                jitter_ms=dict(zip(map(str, fit['rv'].names), map(float, jit))),
                jitter_at_bound=[str(n) for n, j in zip(fit['rv'].names, jit) if j < 0.11 or j > 29.9])


def load_fit(system, out=OUT_DIR):
    with open(os.path.join(out, f'{system}_fit.pkl'), 'rb') as fh:
        return pickle.load(fh)


def plot_fit(fit, path):
    """RVs and model (a), RVs minus the companions with P < 2000 d against the
    longer-period ones (b), residuals (c)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import figtools
    import rv_model as cs
    import synthetic_dr4 as sd
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    rv, th, n, tord = fit['rv'], fit['theta'], fit['npl'], fit['tord']
    P = cs.PARAMS_PER_PLANET
    names = [L['name'] for L in fit['lit']] + ['long-period Keplerian (extra)']*fit['n_extra']
    _, off, res = cs.neg_lnL(th, rv, n, tord, True)
    comp = lambda t, j: cs.keplerian(t, rv.tref, *th[P*j:P*(j + 1)])     # noqa: E731
    y = rv.rv - off[rv.ds]
    long_j = [j for j in range(n) if np.exp(th[P*j]) > 2000]
    short_j = [j for j in range(n) if j not in long_j]
    k = P*n
    trend = lambda t: cs.trend_model(th[k:k + tord], sd.grid_for(t, rv)) if tord else 0*t   # noqa: E731
    y_long = y - sum((comp(rv.t, j) for j in short_j), 0*rv.t)
    tm = np.linspace(rv.t.min() - 200, rv.t.max() + 200, 20000)
    yr = lambda t: 2000 + (t - 2451544.5)/365.25                          # noqa: E731
    cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, ax = plt.subplots(3, 1, figsize=(18, 22), gridspec_kw=dict(height_ratios=[1.2, 1.2, 0.8]))
    for d, idx in enumerate(rv.members):
        lab = str(rv.names[d]).split(':')[1] if ':' in str(rv.names[d]) else str(rv.names[d])
        jit = np.exp(th[-rv.nds + d])
        e = np.hypot(rv.err[idx], jit)
        kw = dict(fmt='o', ms=4, color=cols[d % len(cols)], elinewidth=1)
        ax[0].errorbar(yr(rv.t[idx]), y[idx], e, label=f'{lab} (N={len(idx)}, jitter {jit:.1f} m/s)', **kw)
        ax[1].errorbar(yr(rv.t[idx]), y_long[idx], e, label=lab, **kw)
        ax[2].errorbar(yr(rv.t[idx]), res[idx], e, label=lab, **kw)
    ax[0].plot(yr(tm), cs.model_rv(th, sd.grid_for(tm, rv), n, tord), 'k-', lw=0.6,
               label='full model', zorder=5)
    ax[0].set(ylabel='RV - offset (m/s)', title=f"(a) {fit['system']}: RVs (nightly bins) and the reference RV fit")
    ax[1].plot(yr(tm), sum((comp(tm, j) for j in long_j), 0*tm) + trend(tm), 'k-', lw=2,
               label='long-period companions' + (f' + trend (order {tord})' if tord else ''), zorder=5)
    for i, j in enumerate(long_j):
        role = 'extra' if j >= len(fit['lit']) else ('context' if fit['lit'][j].get('context') else 'sample')
        ax[1].plot(yr(tm), comp(tm, j), ['--', ':'][i % 2], color='0.35', lw=2,
                   label=f'{names[j]} ({role}): P {np.exp(th[P*j]):.0f} d, K {np.exp(th[P*j + 1]):.1f} m/s, '
                         f'e {th[P*j + 2]**2 + th[P*j + 3]**2:.2f}')
    ax[1].set(ylabel='RV - offset - companions with P < 2000 d (m/s)',
              title='(b) The long-period signal (companions with P < 2000 d removed)')
    ax[2].axhline(0, color='0.5', lw=1, gid='guide')
    ax[2].set(xlabel='year', ylabel='residual (m/s)', title='(c) Residuals (error bars include jitter)')
    for a, v in ((ax[1], y_long), (ax[2], res)):
        lo, hi = np.percentile(v, [0.5, 99.5])
        a.set_ylim(lo - 0.15*(hi - lo), hi + 0.15*(hi - lo))
    for a in ax:
        a.set_xlim(yr(tm[0]), yr(tm[-1]))
        figtools.legend_clear(a, fontsize=11)
    fig.tight_layout()
    figtools.check_legends(fig)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def _run_one(system, out, figures):
    """Worker: fit one system, its output going to <system>.log."""
    t0 = time.time()
    with open(os.path.join(out, f'{system}.log'), 'w') as log, contextlib.redirect_stdout(log), \
            contextlib.redirect_stderr(log):
        try:
            fit = fit_system(system)
            with open(os.path.join(out, f'{system}_fit.pkl'), 'wb') as fh:
                pickle.dump(fit, fh)
            s = summarize(fit)
            json.dump(s, open(os.path.join(out, f'{system}_summary.json'), 'w'), indent=1)
            if figures:
                plot_fit(fit, os.path.join(out, f'{system}_rv.png'))
            status = 'ok'
        except Exception as ex:                                   # recorded, the other systems go on
            traceback.print_exc()
            s, status = None, f'FAILED: {type(ex).__name__}: {ex}'
    return dict(system=system, status=status, minutes=round((time.time() - t0)/60, 1),
                n_rv=s and s['n_rv'], kind=s and s['kind'], trend_order=s and s['trend_order'],
                n_extra=s and s['n_extra'],
                companions=s and '; '.join(f"{c['name']} ({c['role']}) P {c['P_d']:.1f} K {c['K_ms']:.1f}"
                                           for c in s['companions']))


def progress(done, total, last='', width=30):
    """One-line progress bar, redrawn in place: [#####-----] done/total."""
    n = int(width*done/total)
    sys.stdout.write(f"\r[{'#'*n}{'-'*(width - n)}] {done}/{total} finished  {last:<50.50s}")
    sys.stdout.flush()


def run(systems, jobs, out, figures):
    """Fit `systems` in parallel (`jobs` at a time) into `out`, with a progress
    bar; writes summary.csv and returns its rows."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import pandas as pd
    os.makedirs(out, exist_ok=True)
    rows = []
    progress(0, len(systems))
    with ProcessPoolExecutor(min(jobs, len(systems))) as ex:
        futs = [ex.submit(_run_one, s, out, figures) for s in systems]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            progress(len(rows), len(systems), f"{r['system']}: {r['status'][:30]} ({r['minutes']} min)")
    print()
    rows.sort(key=lambda r: r['system'])
    for r in rows:
        print(f"{r['system']:10s} {r['status'][:70]:70s} {r['minutes']:5.1f} min")
    pd.DataFrame(rows).to_csv(os.path.join(out, 'summary.csv'), index=False)
    print(f"wrote {os.path.join(out, 'summary.csv')}")
    return rows


PROJECT_PYTHON = os.path.expanduser('~/opt/anaconda3/envs/gaia_dr4_sensitivity/bin/python')


if __name__ == '__main__':
    ## The project environment (Python 3.14): run from any other Python (e.g. the base
    ## Anaconda 3.7, whose scipy fails in the worker processes), the module restarts itself
    ## with the project interpreter
    if sys.version_info < (3, 12):
        if os.path.exists(PROJECT_PYTHON):
            print(f'restarting with {PROJECT_PYTHON} (this is Python {sys.version.split()[0]})', flush=True)
            os.execv(PROJECT_PYTHON, [PROJECT_PYTHON, '-I', os.path.abspath(__file__), *sys.argv[1:]])
        sys.exit(f'needs the gaia_dr4_sensitivity environment (Python 3.14); this is {sys.executable}')

    ## Settings
    SYSTEMS = all_systems()              # or e.g. ['HD111232', 'HD160691']
    JOBS = max(1, (os.cpu_count() or 2) - 1)   # fits running at once
    FIGURES = True                       # draw <system>_rv.png per system
    OUT = OUT_DIR                        # validation/results/reference_fits/
    DRY_RUN = '--dry-run' in sys.argv    # list what would run, fit nothing

    unknown = sorted(set(SYSTEMS) - set(all_systems()))
    if unknown:
        sys.exit(f'no RV table in rv_data/systems/ for {unknown}')
    print(f'{len(SYSTEMS)} systems, {min(JOBS, len(SYSTEMS))} at a time -> {OUT}: {" ".join(SYSTEMS)}')
    if not DRY_RUN:
        run(SYSTEMS, JOBS, OUT, FIGURES)
