"""Linear Delta-BIC periodogram for the blind search.

For every trial period P a circular companion (e = 0, M0 = 0; the phase is
absorbed into the linear amplitudes) is added to the current model: baseline
plus the companions already found (their P, e, M0 held fixed, their amplitudes
refit). The gain in chi2 comes from a small linear problem solved after
projecting out the current model's columns (partial regression), so each period
costs only a few small matrix products. Jitter is refitted at every trial period
(default) or held at the value passed in.

    Delta-BIC(P) = Delta chi2(P) - k_new ln N

where k_new = 1 (the period) + the new linear amplitudes (4 astrometric, 2 RV;
both when both data types are selected).

Eccentricity and phase enter only in the refit of the best peak (search.py), as
in RVSearch. The period grid is rvsearch's: uniform in frequency with spacing
1/(2 pi T), T the baseline.

Test mode:  python -I blind_search/periodogram.py --test
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def period_grid(span, pmin, pmax, oversampling=1.0):
    """Periods uniform in frequency from 1/pmax to 1/pmin, spacing 1/(2 pi
    span)/oversampling."""
    fmin, fmax = 1/pmax, 1/pmin
    nf = int(((fmax - fmin)*2*np.pi*span + 1)*oversampling)
    return 1/np.linspace(fmax, fmin, nf)


def _delta_bic_refit_known(stack, planets, periods, window):
    """Delta-BIC with the known companions' (P, e, M0) re-optimised at every
    trial period, as RVSearch does (its trial fits keep every known planet's
    per, tc, k, secosw, sesinw free). The null model is likewise re-optimised
    (RVSearch's null is a full maximum-likelihood fit). Each optimisation is
    Nelder-Mead on the known companions' 3 nonlinear parameters each, started
    from the values in `planets`, with the same bounds as fit.py (P within
    +-window of its start, 0 <= e <= E_MAX); linear parameters and jitter are
    profiled out exactly as in the fixed-known path."""
    import fit
    from scipy.optimize import minimize
    n = len(planets)
    th0 = np.array([[np.log(P), e, M0] for P, e, M0 in planets]).ravel()
    starts = [p[0] for p in planets]
    pen = (1 + stack.k_comp_astro + stack.k_comp_rv)*np.log(stack.n)
    s20 = stack.fit_free_jitter(stack.design(planets))[3]

    def best_lnL(extra):
        def nl(th):
            pl = fit._unpack(th, n)
            for i in range(n):
                if not (abs(pl[i][0]/starts[i] - 1) <= window and 0 <= pl[i][1] <= fit.E_MAX):
                    return 1e30
            X = stack.design(pl) if extra is None else np.hstack([stack.design(pl), extra])
            return -stack.fit_free_jitter(X, s20)[2]
        r = minimize(nl, th0, method='Nelder-Mead',
                     options={'xatol': 1e-5, 'fatol': 1e-4, 'maxiter': 400*n, 'adaptive': True})
        return -min(r.fun, nl(th0))

    lnL0 = best_lnL(None)
    return np.array([2*(best_lnL(stack.companion_columns(P, 0.0, 0.0)) - lnL0) - pen for P in periods])


def delta_bic(stack, planets, jit, periods, refit_known=False, window=0.1):
    """Delta-BIC of adding one circular companion at each of `periods`, to the
    model with the baseline and `planets` (a list of (P, e, M0)).

    jit=None (the default of the search; RVSearch's formulation): every
    channel's jitter is refitted at every trial period, in the model with and
    without the new companion. The null model then absorbs a strong signal into
    its jitter, so the likelihood gain of a wrong period is bounded (about the
    number of fitted amplitudes) instead of growing with the signal power; this
    is what keeps the eFAP tail clean. With jit an array (jitter held fixed) the
    gain at wrong periods grows as (amplitudes/N) x signal power and swamps eFAP
    for strong signals (QUESTIONS.md, 2026-11-10).

    refit_known=False (default): the known companions' (P, e, M0) are held
    fixed, only their linear amplitudes are refit. refit_known=True: their (P,
    e, M0) are re-optimised at every trial period (RVSearch's convention),
    within +-`window` in P; slower (see `_delta_bic_refit_known`)."""
    k_new = 1 + stack.k_comp_astro + stack.k_comp_rv
    pen = k_new*np.log(stack.n)
    out = np.empty(len(periods))
    if refit_known and planets:
        # RVSearch convention: known companions' (P, e, M0) free in every fit (jitter free, as above).
        if jit is not None:
            raise ValueError('refit_known needs jit=None (free jitter)')
        return _delta_bic_refit_known(stack, planets, periods, window)
    if jit is None:
        X0 = stack.design(planets)
        _, _, lnL0, s20 = stack.fit_free_jitter(X0)
        for i, P in enumerate(periods):
            X = np.hstack([X0, stack.companion_columns(P, 0.0, 0.0)])
            out[i] = 2*(stack.fit_free_jitter(X, s20)[2] - lnL0) - pen
        return out
    sig = stack.sigma(jit)
    w = 1/sig
    X0 = stack.design(planets)*w[:, None]
    U, s, _ = np.linalg.svd(X0, full_matrices=False)
    Q = U[:, s > s[0]*1e-10]
    yw = stack.y*w
    r = yw - Q@(Q.T@yw)
    for i, P in enumerate(periods):
        C = stack.companion_columns(P, 0.0, 0.0)*w[:, None]
        C = C - Q@(Q.T@C)
        G, g = C.T@C, C.T@r
        out[i] = g@np.linalg.lstsq(G, g, rcond=1e-10)[0] - pen
    return out


if __name__ == '__main__':
    if '--test' in sys.argv:
        t = period_grid(1000.0, 10.0, 2000.0)
        assert t[0] > 9.99 and t[-1] < 2000.1 and np.all(np.diff(t) > 0)
        print('periodogram.py test passed (full tests: tests/test_blind_search.py)')
