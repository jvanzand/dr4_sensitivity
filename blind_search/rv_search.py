"""Iterative blind search of radial velocities with the simulated threshold.

The loop of `search.run_search` on RV channels: companions are refit with free
eccentricity (e up to `fit.E_MAX`), although the periodogram template is
circular, and the model may carry a linear or quadratic trend, handled as in
RVSearch (trend.py; Judah, 2026-11-16):
  1. before the first companion, a slope is added if the data warrant it, then
     a curvature if they warrant that too;
  2. at each step the new companion may replace trend terms (in the
     periodogram, and at its eccentric refit, where the lowest-BIC trend
     order is kept), so a signal first fitted as a trend can come back as a
     companion;
  3. after a companion is added, a trend term the data no longer warrant is
     dropped.
"Warranted" = likelihood-ratio test at the search's FAP per decision. The
threshold for a companion comes from M noise-sampled datasets drawn from the
current model (companions, trend, offsets and fitted jitter of each
instrument, plus white noise), each refit and searched like the real data.

No RV noise model: the noise-sampled datasets carry white noise only.

Test mode:  python -I blind_search/rv_search.py --test   (a short synthetic
check). The HD 191939 run is validation/rv_hd191939.py.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'tests'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np  # noqa: E402

import data  # noqa: E402
import model  # noqa: E402
import search  # noqa: E402


def amplitudes(stack, known):
    """RV semi-amplitude K (m/s) of each companion in the joint linear fit."""
    coef = stack.fit_free_jitter(stack.design(known))[0]
    k0, kc = stack.k_base, stack.k_comp_rv
    return [float(np.hypot(*coef[k0 + kc*i:k0 + kc*(i + 1)]))
            for i in range(len(known))]


def run_rv_search(channels, executor, M=50, fap=1e-3, max_steps=12, pmin=2.0,
                  max_trend=2, seed=0, label='', out_dir=None, truth=(),
                  verbose=True):
    """Run the search (`search.run_search` with trends up to `max_trend`).
    Returns dict(known, trend, trend_log, steps, amplitudes, periods): `steps`
    has one dict per step (peak, period, threshold, detected, trend order
    before and after, trend moves; the arrays are dropped). With `out_dir` a
    step figure is written for each step (`truth` periods are drawn on it)."""
    res = search.run_search(model.Stack(channels), executor, M=M, fap=fap,
                            max_companions=max_steps, seed=seed, truth=truth,
                            label=label, out_dir=out_dir, pmin=pmin,
                            truth_label='known planet period', verbose=verbose,
                            max_trend=max_trend)
    steps = [{k: v for k, v in st_.items() if k not in ('dbic', 'maxima')}
             for st_ in res['steps']]
    return dict(known=res['known'], trend=res['trend'],
                trend_log=res['trend_log'], steps=steps,
                amplitudes=amplitudes(res['stack'], res['known']),
                periods=res['periods'])


def selftest():
    """Synthetic RV system: two circular planets plus a linear trend. Checks
    that the loop adopts the trend, finds both planets and stops."""
    from concurrent.futures import ProcessPoolExecutor
    import test_blind_search as T
    with ProcessPoolExecutor(3) as ex:
        s = model.Stack(T.rv_channels(5, n=160, span=3000.0), rv_trend=1)
        base = np.zeros(s.k_base)
        base[-1] = 12.0                                   # trend term, m/s
        T.fill(s, [(40.0, 0.0, 0.0), (210.0, 0.0, 0.0)],
               [[9.0, 4.0], [6.0, -5.0]], 3, base=base)
        ch = [data.Channel('rv', c.name, c.t, s.y[s.chan == i], c.err)
              for i, c in enumerate(s.channels)]
        r = run_rv_search(ch, ex, M=10, pmin=5.0, max_steps=6, seed=1)
        Ps = sorted(p[0] for p in r['known'])
        ok = (len(Ps) == 2 and abs(Ps[0]/40 - 1) < 0.03
              and abs(Ps[1]/210 - 1) < 0.05 and r['trend'] >= 1)
        assert ok, (Ps, r['trend'])
        print('rv_search selftest passed: planets', [round(p, 1) for p in Ps],
              'trend order', r['trend'])


if __name__ == '__main__':
    if '--test' in sys.argv:
        selftest()
    else:
        print(__doc__)
