"""Simulated per-step detection threshold for the blind search (replaces eFAP;
QUESTIONS.md, 2026-11-11).

At search step k the current model is the (k - 1)-planet model: the reference
motion plus the k - 1 companions already found, with their refitted (P, e, M0)
and the jitter, all fitted to the real data. M noise-sampled datasets are
generated from it: at every real transit (real epoch, scan angle, parallax
factor) take the model's predicted along-scan position and add an independent
Gaussian random number of standard deviation sqrt(err^2 + s^2), s = the fitted
jitter; each noise-sampled dataset uses new random numbers. Each noise-sampled
dataset goes through the same pipeline as the real one: the known companions are
refitted (Nelder-Mead, `fit._polish`), then the Delta-BIC periodogram is
computed (`periodogram.delta_bic`, jitter free). Only its highest peak is kept.
The M highest peaks show how large a false alarm gets when there is nothing more
to find; the detection threshold for a false-alarm probability FAP is their (1 -
FAP) quantile. Because M is small, the quantile is read from a Gumbel
(extreme-value) distribution fitted to the M maxima (two parameters: location,
scale); the direct 0.95 and 0.99 quantiles are returned as well.

The real data's own periodogram body is never used, so a real unknown companion
cannot raise the threshold through its leakage (the failure of eFAP). The noise
is white by construction.

Test mode:  python -I blind_search/simthreshold.py --test
"""

import copy
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit  # noqa: E402
import periodogram as pg  # noqa: E402

PMIN = 20.0       # d; lower end of the trial-period grid (same as the Gaia-4 tests)


def null_model(stack, known):
    """Prediction (native units) and jitter^2 of the (k - 1)-planet model with
    only the `known` companions, fitted to the stack's data with free jitter."""
    X = stack.design(known)
    coef, _, _, s2 = stack.fit_free_jitter(X)
    return X@coef, s2


def one_sim(args):
    """One noise-sampled dataset through the pipeline; returns (its highest
    Delta-BIC peak, its whole periodogram if `keep_full` else None)."""
    stack, known, mean, s2, seed, periods, refit, noise, window, keep_full = args
    rng = np.random.default_rng(seed)
    s = copy.copy(stack)
    sd = np.sqrt(stack.err**2 + s2[stack.chan]) if noise == 'fitted' else stack.err
    s.y = mean + rng.normal(0.0, sd)
    kn = list(known)
    if kn and refit:
        kn, _ = fit._polish(s, [], kn, window)
    d = pg.delta_bic(s, kn, None, periods)
    return float(d.max()), (d if keep_full else None)


def gumbel_fit(maxima):
    from scipy.stats import gumbel_r
    loc, scale = gumbel_r.fit(maxima)
    return float(loc), float(scale)


def gumbel_quantile(loc, scale, fap):
    """Value exceeded with probability `fap` by the highest peak: loc - scale
    ln(-ln(1 - fap))."""
    return loc - scale*np.log(-np.log1p(-fap))


def simulated_threshold(stack, known, M=50, seed=0, executor=None, refit=True, noise='fitted', window=0.1,
                        periods=None, faps=(1e-2, 1e-3), boot=200, keep=0):
    """Simulated threshold for the next companion given the `known` ones.

    Returns a dict: maxima (the M highest peaks), loc/scale (Gumbel fit), thr
    (FAP -> Gumbel threshold), thr_sd (FAP -> bootstrap standard deviation of
    that threshold), q95/q99 (direct quantiles of the maxima), noise ('fitted' =
    white Gaussian noise at the formal errors plus the fitted jitter; 'formal' =
    formal errors only), refit (known companions refitted in every simulation),
    periods, and `examples` = the whole Delta-BIC periodograms of the first
    `keep` noise-sampled datasets (for plotting)."""
    periods = pg.period_grid(stack.span, PMIN, 2*stack.span) if periods is None else periods
    mean, s2 = null_model(stack, known)
    args = [(stack, known, mean, s2, 1000003*seed + i, periods, refit, noise, window, i < keep) for i in range(M)]
    if isinstance(executor, str) and executor == 'serial':           # run in this process (for callers that parallelise elsewhere)
        res = [one_sim(a) for a in args]
    elif executor is None:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(M, max(1, (os.cpu_count() or 2) - 1))) as ex:
            res = list(ex.map(one_sim, args, chunksize=2))
    else:
        res = list(executor.map(one_sim, args, chunksize=2))
    maxima = np.array([r[0] for r in res])
    examples = [r[1] for r in res if r[1] is not None]
    loc, scale = gumbel_fit(maxima)
    rng = np.random.default_rng(seed + 1)
    bfits = [gumbel_fit(rng.choice(maxima, len(maxima))) for _ in range(boot)]
    return dict(maxima=maxima, loc=loc, scale=scale, jitter2=s2, noise=noise, refit=refit, periods=periods, examples=examples,
                thr={f: float(gumbel_quantile(loc, scale, f)) for f in faps},
                thr_sd={f: float(np.std([gumbel_quantile(l, sc, f) for l, sc in bfits])) for f in faps},
                q95=float(np.quantile(maxima, 0.95)), q99=float(np.quantile(maxima, 0.99)))


if __name__ == '__main__':
    if '--test' in sys.argv:
        import data
        import model
        sys.path.insert(0, os.path.join(data.ROOT, 'tests'))
        import test_blind_search as T
        g = data.load_gaia_dr4(T.GAIA4)
        s = T.fill(model.Stack([g]), [], [], 3)
        r = simulated_threshold(s, [], M=12, seed=1)
        assert np.all(np.isfinite(r['maxima'])) and np.isfinite(r['thr'][1e-3]) and r['thr'][1e-3] > r['thr'][1e-2]
        print('simthreshold.py test passed: median max %.1f, Gumbel thresholds %s' % (np.median(r['maxima']), r['thr']))
