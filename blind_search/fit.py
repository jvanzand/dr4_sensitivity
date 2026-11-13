"""Nonlinear refit of companions' (P, e, M0) by variable projection with the
jitter free.

For fixed (P, e, M0) of every companion all linear parameters (reference motion,
offsets, trend, Thiele-Innes constants, K cos w, K sin w) and the per-channel
jitter are fitted exactly/jointly by `Stack.fit_free_jitter`; only (ln P, e, M0)
per companion are searched (Nelder-Mead). P is kept within +-`window` of its
start and e in [0, E_MAX].

`add_and_refit(stack, planets, P0)`: add a companion near the periodogram peak
P0 (several starting eccentricities and phases, others held fixed), then polish
all companions together.

Test mode:  python -I blind_search/fit.py --test
"""

import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

E_MAX = 0.7        # blind eccentricity cap (as in the RV search)


def _unpack(th, n):
    return [(float(np.exp(th[3*i])), float(th[3*i + 1]), float(np.mod(th[3*i + 2], 2*np.pi))) for i in range(n)]


def _neg_lnL(th, stack, fixed, starts, window):
    n = len(th)//3
    pl = _unpack(th, n)
    for i in range(n):
        if not (abs(pl[i][0]/starts[i] - 1) <= window and 0 <= pl[i][1] <= E_MAX):
            return 1e30
    return -stack.fit_free_jitter(stack.design(list(fixed) + pl))[2]


def _polish(stack, fixed, planets, window):
    th0 = np.array([[np.log(P), e, M0] for P, e, M0 in planets]).ravel()
    starts = [p[0] for p in planets]
    r = minimize(_neg_lnL, th0, args=(stack, fixed, starts, window), method='Nelder-Mead',
                 options={'xatol': 1e-6, 'fatol': 1e-7, 'maxiter': 1500*len(planets), 'adaptive': True})
    return _unpack(r.x, len(planets)), -r.fun


def add_and_refit(stack, planets, P0, window=0.1):
    """Add a companion near P0 and refit. Returns (planets, lnL)."""
    best = (None, -np.inf)
    for e0, M00 in ((0.0, 0.0), (0.3, 0.0), (0.3, 2.1), (0.3, 4.2), (0.6, 1.0), (0.6, 3.5)):
        pl, l = _polish(stack, list(planets), [(P0, e0, M00)], window)
        if l > best[1]:
            best = (pl, l)
    new = best[0][0]
    allp, l = _polish(stack, [], list(planets) + [new], window) if planets else (best[0], best[1])
    return allp, l


def refit_circular(stack, planets, P0, window=0.5, n_scan=400):
    """Refit the period of ONE new CIRCULAR companion (e = 0; the phase lives in
    the free linear amplitudes) with the companions in `planets` held fixed: a
    scan of `n_scan` log-spaced periods over [P0/(1 + window), P0 (1 + window)]
    (profile likelihood: the linear parameters and the jitter are refitted at
    every period) followed by a bounded 1-D polish around the best point. The
    window must be wide enough to reach the neighbouring periodogram grid
    periods (the grid is uniform in frequency, so at long periods adjacent
    periods differ by up to a factor 1.45 for the Gaia-4 injection map). Returns
    (period, lnL). TO REVISIT: eccentric orbits."""
    from scipy.optimize import minimize_scalar
    lo, hi = np.log(P0/(1 + window)), np.log(P0*(1 + window))

    def nl(lnP):
        return -stack.fit_free_jitter(stack.design(list(planets) + [(float(np.exp(lnP)), 0.0, 0.0)]))[2]
    grid = np.linspace(lo, hi, n_scan)
    vals = np.array([nl(x) for x in grid])
    j = int(np.argmin(vals))
    r = minimize_scalar(nl, bounds=(grid[max(j - 1, 0)], grid[min(j + 1, n_scan - 1)]), method='bounded', options={'xatol': 1e-8})
    best, val = (r.x, r.fun) if r.fun <= vals[j] else (grid[j], vals[j])
    return float(np.exp(best)), float(-val)


if __name__ == '__main__':
    if '--test' in sys.argv:
        import data
        import model
        sys.path.insert(0, os.path.join(data.ROOT, 'tests'))
        import test_blind_search as T
        g = data.load_gaia_dr4(T.GAIA4)
        s = model.Stack([g])
        T.fill(s, [(300.0, 0.3, 1.0)], [[0.25, -0.15, 0.2, 0.3]], 1)
        pl, l = add_and_refit(s, [], 295.0)
        assert abs(pl[0][0]/300 - 1) < 0.02 and abs(pl[0][1] - 0.3) < 0.1, pl
        print('fit.py test passed: refit', np.round(pl[0], 3))
