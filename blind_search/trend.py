"""Trend handling of the blind searches (RVSearch's logic; Judah, 2026-11-16).

A trend is a polynomial in time added to the model: for RVs a slope, then a
curvature (`Stack.rv_trend`); for astrometry an acceleration, then its rate of
change (`Stack.astro_trend`; Gaia's 7- and 9-parameter solutions).

  1. Before the first companion (`raise_trend`): the first trend term is added
     if the data warrant it, then the second if they warrant that too.
  2. At every companion step the periodogram lets the new companion replace
     trend terms (`periodogram.delta_bic`), and when a companion is detected
     the trend order is chosen at the eccentric refit
     (`search.add_with_trend_choice`: the lowest BIC over orders 0 .. current),
     so a signal first fitted as a trend can come back as a companion.
  3. After a companion is added (`settle_trend`): the highest trend term is
     dropped if the data no longer warrant it, and the check repeats.

"Warranted": a likelihood-ratio test of the term (jitter free, companions'
P, e, M0 held fixed) at the search's false-alarm probability per decision:
the gain 2 Delta lnL is compared with a chi-square distribution with k degrees
of freedom, k = the columns the term adds (1 for RVs, 2 per astrometric
channel); the term is warranted if P(chi2_k > gain) < fap.

Test mode:  python -I blind_search/trend.py --test
"""

import os
import sys

import numpy as np
from scipy.stats import chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def term_test(stack, planets):
    """Likelihood-ratio test of the highest trend term of `stack` (order >= 1)
    against the model one order lower, both with `planets` (P, e, M0 fixed).
    Returns (gain = 2 Delta lnL, k = columns of the term, false-alarm
    probability P(chi2_k > gain))."""
    lo = stack.with_trend(stack.trend - 1)
    l1 = stack.fit_free_jitter(stack.design(planets))[2]
    l0 = lo.fit_free_jitter(lo.design(planets))[2]
    gain, k = 2*(l1 - l0), stack.k_base - lo.k_base
    return float(gain), int(k), float(chi2.sf(max(gain, 0.0), k))


def raise_trend(stack, planets, fap, max_trend=2):
    """Step 1: raise the trend order one term at a time while each new term is
    warranted (up to `max_trend`). Returns (stack at the adopted order, log:
    one dict per test)."""
    log = []
    while stack.trend < max_trend:
        hi = stack.with_trend(stack.trend + 1)
        gain, k, p = term_test(hi, planets)
        ok = p < fap
        log.append(dict(action=f'raise trend order {stack.trend} -> {hi.trend}', gain_2dlnL=gain,
                        dof=k, fap=p, adopted=ok, n_companions=len(planets)))
        if not ok:
            break
        stack = hi
    return stack, log


def settle_trend(stack, planets, fap, window=0.1):
    """Step 3: drop the highest trend term while it is no longer warranted;
    after a drop the companions' (P, e, M0) are re-polished at the lower
    order. Returns (stack, planets, log)."""
    import fit
    log = []
    while stack.trend > 0:
        gain, k, p = term_test(stack, planets)
        drop = p >= fap
        log.append(dict(action=f'test trend term {stack.trend} (keep or drop)', gain_2dlnL=gain,
                        dof=k, fap=p, dropped=drop, n_companions=len(planets)))
        if not drop:
            break
        stack = stack.with_trend(stack.trend - 1)
        if planets:
            planets = fit._polish(stack, [], list(planets), window)[0]
    return stack, list(planets), log


if __name__ == '__main__':
    if '--test' in sys.argv:
        import data
        import model
        sys.path.insert(0, os.path.join(data.ROOT, 'tests'))
        import test_blind_search as T
        g = data.load_gaia_dr4(T.GAIA4)
        s = model.Stack([g], astro_trend=1)
        base = np.zeros(s.k_base)
        base[5:7] = [0.08, -0.05]                    # acceleration, mas/yr^2
        T.fill(s, [], [], 2, base=base)
        a, log = raise_trend(s.with_trend(0), [], 1e-3)
        flat = s.with_trend(0)
        T.fill(flat, [], [], 2)
        b, _ = raise_trend(flat, [], 1e-3)
        assert a.trend >= 1 and b.trend == 0, (a.trend, b.trend, log)
        print(f'trend.py test passed: acceleration adopted (order {a.trend}); noise alone: order {b.trend}')
