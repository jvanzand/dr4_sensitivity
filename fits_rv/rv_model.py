"""RV model library: Keplerian RVs, likelihood, maximum-likelihood / -a-posteriori fits, BIC
and literature priors (Judah, 2026-11-20: the blind RV search that lived here as
companion_search.py was removed; the blind RV search is blind_search/rv_search.py).

Used by the reference RV fit of the benchmarks (validation/synthetic_dr4.py) and by
validation/presentation.py.

Parameters per companion: ln P, ln K, sqrt(e) cos w, sqrt(e) sin w, M0 (at the RV
reference time `RVData.tref`); then the trend coefficients (order 0-2, in m/s per
baseline^k) and one ln jitter per dataset. Each dataset's offset is profiled out
exactly in `neg_lnL`. The prior on e >= E_MAX (0.99) is zero (Judah, 2026-11-20):
a steep wall in `neg_lnL` keeps every fit below it.

Literature companions enter under Gaussian priors on ln P, ln K and e around the
quoted values (`literature_priors`; widths: the published error, never narrower than
1% in P, 10% in K and 0.1 in e; 20% in P and K where no error is known). A companion
added beyond them (the reference fit's long-period Keplerian) is capped at
e <= E_BLIND_MAX by `fit` (callers may raise the cap).

BIC = -2 lnL + k ln N_RV, with lnL the data log-likelihood (prior terms excluded) and
k counting every parameter (5 per Keplerian, 1 per trend term, offset and jitter per
dataset).

Test mode:  python -I fits_rv/rv_model.py --test
"""

import math
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.signal import argrelmax

## Keplerian curves come from kepler.c (compiled on first use, ~2-3x faster than numpy),
## else from the numpy solver below
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kepler_c

DELTA_BIC = 20.0     # raised from 10 on 2026-11-06 (Judah), after the first pilot searches
JITTER_BOUNDS = (0.1, 30.0)          # m/s, as in the MCMC prior
E_MAX = 0.99                         # the prior on e >= E_MAX is zero (Judah, 2026-11-20)
E_WALL = 1e6                         # steepness of the wall in `neg_lnL` that enforces it
PERIOD_CAP_BASELINES = 2.0           # period cap in RV baselines (2026-11-09; was 3)
E_BLIND_MAX = 0.7                    # e cap of companions beyond the literature ones
P_MIN_KNOWN_D = 0.05                 # floor for the N forced companions
P_MAX_YR = 200.0
PARAMS_PER_PLANET = 5
MAX_TREND_ORDER = 2                  # linear, then quadratic
ECC_STARTS = [(0.05, 0.05), (0.55, 0.3), (-0.4, 0.5)]   # (sqrt(e) cos w, sqrt(e) sin w)
N_JOINT = 4                          # best single-companion seeds refit jointly
## Prior widths for the literature companions (accepted 2026-11-09)
LIT_FLOOR_P, LIT_FLOOR_K, LIT_SIGMA_E, LIT_DEFAULT_REL = 0.01, 0.10, 0.10, 0.20
CATALOG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'sample', 'stage1_catalog.csv')
CUT_CATALOG = os.path.join(os.path.dirname(CATALOG), 'census_cut_companions.csv')
OVERSAMPLE = 5


class RVData:
    """RVs in days and m/s, with an integer dataset index per point."""

    def __init__(self, t, rv, err, dataset):
        order = np.argsort(t)
        self.t = np.asarray(t, float)[order]
        self.rv = np.asarray(rv, float)[order]
        self.err = np.asarray(err, float)[order]
        names, idx = np.unique(np.asarray(dataset)[order], return_inverse=True)
        self.names = list(names)
        self.ds = idx
        self.n = len(self.t)
        self.nds = len(self.names)
        self.tref = 0.5*(self.t.min() + self.t.max())
        self.span = self.t.max() - self.t.min()
        self.pmax = min(P_MAX_YR*365.25, PERIOD_CAP_BASELINES*self.span)   # period cap
        self.members = [np.where(self.ds == d)[0] for d in range(self.nds)]


def kepler_E(M, e):
    M = np.mod(M, 2*np.pi)
    E = np.where(e < 0.8, M, np.pi)
    for _ in range(60):
        dE = -(E - e*np.sin(E) - M)/(1 - e*np.cos(E))
        E = E + dE
        if np.max(np.abs(dE)) < 1e-12:
            break
    return E


def keplerian(t, tref, lnP, lnK, s1, s2, M0):
    e = min(s1*s1 + s2*s2, E_MAX)
    w = np.arctan2(s2, s1)
    P = np.exp(lnP)
    E = kepler_E(2*np.pi*(t - tref)/P + M0, e)
    nu = 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
    return np.exp(lnK)*(np.cos(nu + w) + e*np.cos(w))


def add_keplerian(m, t, tref, lnP, lnK, s1, s2, M0):
    """m + one Keplerian (accumulated in place when kepler.c is available)."""
    if kepler_c.AVAILABLE and m.flags['C_CONTIGUOUS'] and t.flags['C_CONTIGUOUS']:
        e = min(s1*s1 + s2*s2, E_MAX)
        kepler_c.keplerian_rv_add(m, t, tref, math.exp(lnP), math.exp(lnK), e,
                                  math.atan2(s2, s1), M0)
        return m
    return m + keplerian(t, tref, lnP, lnK, s1, s2, M0)


def trend_model(coefs, data):
    """Polynomial trend in x = (t - tref)/span: coefs[0]*x + coefs[1]*x^2 (m/s)."""
    x = (data.t - data.tref)/data.span
    m = np.zeros_like(x)
    for k, c in enumerate(coefs):
        m = m + c*x**(k + 1)
    return m


def unpack(theta, npl, tord, nds):
    pl = np.asarray(theta[:PARAMS_PER_PLANET*npl]).reshape(npl, PARAMS_PER_PLANET)
    k = PARAMS_PER_PLANET*npl
    coefs = np.asarray(theta[k:k + tord])
    k += tord
    return pl, coefs, np.exp(theta[k:k + nds])


def model_rv(theta, data, npl, tord):
    pl, coefs, _ = unpack(theta, npl, tord, data.nds)
    m = trend_model(coefs, data)
    for p in pl:
        m = add_keplerian(m, data.t, data.tref, *p)
    return m


def neg_lnL(theta, data, npl, tord, return_parts=False):
    """Negative log likelihood with each dataset's offset profiled out exactly."""
    pl, coefs, jit = unpack(theta, npl, tord, data.nds)
    m = trend_model(coefs, data)
    for p in pl:
        m = add_keplerian(m, data.t, data.tref, *p)
    r = data.rv - m
    var = data.err**2 + jit[data.ds]**2
    w = 1/var
    off = np.zeros(data.nds)
    for d, idx in enumerate(data.members):
        off[d] = np.sum(w[idx]*r[idx])/np.sum(w[idx])
    res = r - off[data.ds]
    nll = 0.5*np.sum(res**2*w + np.log(2*np.pi*var)) + E_WALL*e_excess(theta, npl)
    return (nll, off, res) if return_parts else nll


def e_excess(theta, npl):
    """Sum over companions of (e - E_MAX)^2 where e > E_MAX. The Keplerian clips e
    at E_MAX, so without a wall the likelihood is flat beyond it and a poorly
    constrained companion drifts there (HD 191939 g reached e = 1.002); the wall
    (E_WALL x this, in `neg_lnL`) makes e > E_MAX impossible in every fit, and
    `cap_e` moves a fit's final values strictly below E_MAX. Zero for any valid
    fit, so data lnL and BIC are unchanged."""
    tot = 0.0
    for j in range(npl):
        k = PARAMS_PER_PLANET*j
        tot += max(theta[k + 2]**2 + theta[k + 3]**2 - E_MAX, 0.0)**2
    return tot


def bounds_for(npl, tord, nds, floors, pmax, p_hi=None):
    """Parameter bounds. A companion's period is capped at `pmax` (2x the RV baseline,
    at most 200 yr): with less than a third of an orbit covered, a Keplerian is not
    distinguishable from curvature."""
    b = []
    for k in range(npl):
        top = pmax if p_hi is None or k >= len(p_hi) else p_hi[k]
        b += [(np.log(floors[k]), np.log(top)), (np.log(0.01), np.log(1e5)),
              (-1, 1), (-1, 1), (None, None)]
    b += [(-1e5, 1e5)]*tord
    b += [(np.log(JITTER_BOUNDS[0]), np.log(JITTER_BOUNDS[1]))]*nds
    return b


def prior_penalty(theta, pri):
    """Negative log of the Gaussian priors on ln P, ln K, and e of the first len(pri)
    companions (additive constants dropped)."""
    tot = 0.0
    for j, L in enumerate(pri):
        b = theta[PARAMS_PER_PLANET*j:PARAMS_PER_PLANET*(j + 1)]
        tot += 0.5*((b[0] - L['lnP'][0])/L['lnP'][1])**2
        if L['lnK'] is not None:
            tot += 0.5*((b[1] - L['lnK'][0])/L['lnK'][1])**2
        e = min(b[2]**2 + b[3]**2, E_MAX)
        tot += 0.5*((e - L['e'][0])/L['e'][1])**2
    return tot


def ecc_penalty(theta, n_lit, npl):
    """Soft wall keeping the blind companions (those after the first `n_lit`) at
    e <= E_BLIND_MAX; `fit` projects onto the cap afterwards, so the wall only has to steer."""
    tot = 0.0
    for j in range(n_lit, npl):
        b = theta[PARAMS_PER_PLANET*j:PARAMS_PER_PLANET*(j + 1)]
        tot += 1e5*max(b[2]**2 + b[3]**2 - E_BLIND_MAX, 0.0)**2
    return tot


def cap_e(theta, n_lit, npl):
    """Rescale (sqrt(e) cos w, sqrt(e) sin w) onto the cap of any companion above it:
    E_BLIND_MAX for blind companions (those after the first `n_lit`), just below
    E_MAX for the others (the prior on e >= E_MAX is zero)."""
    theta = theta.copy()
    for j in range(npl):
        k = PARAMS_PER_PLANET*j
        cap = E_BLIND_MAX if j >= n_lit else E_MAX*(1 - 1e-6)
        e = theta[k + 2]**2 + theta[k + 3]**2
        if e > cap:
            theta[k + 2:k + 4] *= math.sqrt(cap/e)
    return theta


def fit(theta0, data, npl, tord, floors, pri=None):
    """Maximum-a-posteriori fit (maximum likelihood where there are no priors). `pri` holds
    the priors of the first companions (the literature ones); any companion after them is
    blind and capped at e <= E_BLIND_MAX. Returns (theta, data lnL at theta): the prior terms
    are not in the returned lnL, so BIC compares data fits only."""
    pri = list(pri or [])[:npl]
    p_hi = [max(data.pmax, math.exp(L['lnP'][0] + 4*L['lnP'][1])) for L in pri]
    res = minimize(lambda th: neg_lnL(th, data, npl, tord) + prior_penalty(th, pri)
                   + ecc_penalty(th, len(pri), npl), theta0,
                   method='L-BFGS-B',
                   bounds=bounds_for(npl, tord, data.nds, floors, data.pmax, p_hi),
                   options={'maxiter': 400, 'maxfun': 20000})
    th = cap_e(res.x, len(pri), npl)
    return th, -neg_lnL(th, data, npl, tord)


def n_params(npl, tord, nds):
    return PARAMS_PER_PLANET*npl + tord + 2*nds


def bic(lnL, npl, tord, data):
    return -2*lnL + n_params(npl, tord, data.nds)*np.log(data.n)


def set_trend_order(theta, npl, old, new):
    """Rearrange a parameter vector for a different trend order (new terms start at 0)."""
    k = PARAMS_PER_PLANET*npl
    coefs = list(theta[k:k + old])[:new] + [0.0]*max(0, new - old)
    return np.concatenate([theta[:k], coefs, theta[k + old:]])


def seed_trend(theta, npl, old, new, data):
    """`theta` with the trend order changed from `old` to `new`; a new top term starts at the
    weighted least-squares coefficient that best removes the current residuals (starting it
    at 0 can leave the optimizer at a point from which the trend never gets going)."""
    th = set_trend_order(theta, npl, old, new)
    if new <= old:
        return th
    _, _, res = neg_lnL(theta, data, npl, old, True)
    jit = unpack(theta, npl, old, data.nds)[2]
    w = 1/(data.err**2 + jit[data.ds]**2)
    x = ((data.t - data.tref)/data.span)**new
    for idx in data.members:
        x[idx] -= np.sum(w[idx]*x[idx])/np.sum(w[idx])
    th[PARAMS_PER_PLANET*npl + old] = float(np.clip(np.sum(w*x*res)/np.sum(w*x*x), -1e5, 1e5))
    return th


def periodogram(data, resid, jitter, periods):
    """Generalized Lomb-Scargle Delta chi^2 of the best circular sinusoid at each trial
    period, with every dataset's offset projected out (weights 1/(err^2 + jitter^2))."""
    w = 1/(data.err**2 + jitter[data.ds]**2)
    y = resid.copy()
    for idx in data.members:
        y[idx] -= np.sum(w[idx]*y[idx])/np.sum(w[idx])
    out = np.zeros(len(periods))
    omega = 2*np.pi/periods
    step = max(1, int(2e6//data.n))
    for a in range(0, len(periods), step):
        arg = np.outer(omega[a:a + step], data.t - data.tref)
        S, C = np.sin(arg), np.cos(arg)
        for idx in data.members:
            wi = w[idx]/np.sum(w[idx])
            S[:, idx] -= (S[:, idx]@wi)[:, None]
            C[:, idx] -= (C[:, idx]@wi)[:, None]
        Sw, Cw = S*w, C*w
        ss, cc, sc = np.sum(Sw*S, 1), np.sum(Cw*C, 1), np.sum(Sw*C, 1)
        sy, cy = Sw@y, Cw@y
        det = ss*cc - sc**2
        det = np.where(det > 1e-12*np.maximum(ss*cc, 1e-30), det, np.inf)
        out[a:a + step] = (cc*sy**2 - 2*sc*sy*cy + ss*cy**2)/det
    return out


def sinusoid_phase(data, resid, jitter, P):
    """Amplitude and phase phi of the best circular sinusoid A cos(2 pi (t - tref)/P - phi)
    at period P, with every dataset's offset projected out."""
    w = 1/(data.err**2 + jitter[data.ds]**2)
    y = resid.copy()
    x = 2*np.pi*(data.t - data.tref)/P
    S, C = np.sin(x), np.cos(x)
    for idx in data.members:
        wi = w[idx]/np.sum(w[idx])
        y[idx] -= np.sum(wi*y[idx])
        S[idx] -= np.sum(wi*S[idx])
        C[idx] -= np.sum(wi*C[idx])
    A = np.array([[np.sum(w*S*S), np.sum(w*S*C)], [np.sum(w*S*C), np.sum(w*C*C)]])
    b = np.array([np.sum(w*S*y), np.sum(w*C*y)])
    a_s, a_c = np.linalg.solve(A + 1e-12*np.eye(2), b)
    return np.hypot(a_s, a_c), np.arctan2(a_s, a_c)


def make_prior(P0, K0=None, e0=None, omega0=None, sigP=None, sigK_rel=None, name='',
               notes=None):
    """One literature companion's priors: Gaussian in ln P (sd = max(published error / P,
    1%), or 20% if none), in ln K (sd = max(relative K or mass error, 10%), or 20% if none;
    omitted if K is unknown), and in e (sd 0.1 around the value; 0 if unreported, i.e.
    assumed circular). Phase and omega are free (omega0 only seeds the start)."""
    sdP = max(sigP/P0, LIT_FLOOR_P) if sigP else LIT_DEFAULT_REL
    sdK = max(sigK_rel, LIT_FLOOR_K) if sigK_rel else LIT_DEFAULT_REL
    return {'name': name, 'P0': float(P0), 'K0': None if K0 is None else float(K0),
            'e0': None if e0 is None else float(e0), 'omega0': omega0,
            'lnP': (math.log(P0), float(sdP)),
            'lnK': None if K0 is None else (math.log(K0), float(sdK)),
            'e': (0.0 if e0 is None else float(e0), LIT_SIGMA_E),
            'notes': notes or {}}


def _prior_from_row(r):
    notes, P0 = {}, float(r.period_d)
    sigP = r.period_err_filled_d if pd.notna(r.period_err_filled_d) else None
    notes['P'] = (f'error {r.period_err_source}' if sigP else 'no error known: 20%')
    e0 = float(r.pl_orbeccen) if pd.notna(r.pl_orbeccen) else None
    notes['e'] = 'published' if e0 is not None else 'unreported: circular assumed'
    K0 = float(r.K_ms) if pd.notna(r.K_ms) else None
    notes['K'] = 'published' if K0 is not None else None
    if K0 is None and pd.notna(r.mass_mj) and pd.notna(r.st_mass):
        K0 = (28.4329*r.mass_mj*r.st_mass**(-2/3)*(P0/365.25)**(-1/3)
              / math.sqrt(1 - (e0 or 0.0)**2))
        notes['K'] = f'derived from {r.mass_type} (sin i = 1)'
    sigK = (float(r.mass_err_mj/r.mass_mj) if pd.notna(r.mass_err_mj) and r.mass_mj > 0
            else None)
    notes['K_width'] = 'mass-error proxy' if sigK else 'no error known: 20%'
    return make_prior(P0, K0, e0, math.radians(r.omega_deg) if pd.notna(r.omega_deg) else None,
                      sigP, sigK, r.pl_name, notes)


def literature_priors(system_id, catalog=CATALOG, context=True):
    """Priors for the system's literature companions (sample/stage1_catalog.csv), strongest
    K first. K comes from the catalog, else from M sin i, P, e, and M_* (flagged). The relative
    K width uses the mass error (the catalog carries no K errors). Notes record the source of
    each value.

    With `context` (default; Judah, 2026-11-20), the companions of the same host that failed
    a census cut (sample/census_cut_companions.csv: SMA > 5 AU, mass > 80 M_J, no mass, ...)
    are added too, marked 'context': True with their 'cut_reason', so the assisted RV fit has
    every known companion and no period limit. They are not part of the sample. A cut
    companion with neither a period nor a K or mass is skipped (noted in 'skipped')."""
    cat = pd.read_csv(catalog, dtype={'system_id': str}, low_memory=False)
    rows = cat[cat.system_id == str(system_id)]
    out = [_prior_from_row(r) for _, r in rows.iterrows()]
    if context and len(rows) and os.path.exists(CUT_CATALOG):
        cut = pd.read_csv(CUT_CATALOG, low_memory=False)
        for _, r in cut[cut.hostname.isin(set(rows.hostname))].iterrows():
            if pd.isna(r.period_d) or (pd.isna(r.K_ms) and pd.isna(r.mass_mj)):
                continue
            L = _prior_from_row(r)
            if L['K0'] is None:
                continue
            L.update(context=True, cut_reason=r.cut_reason)
            out.append(L)
    out.sort(key=lambda L: -(L['K0'] or 0.0))
    return out


def seed_literature_planet(theta, j, data, lit, tord=0, n_have=None):
    """Insert the j-th literature companion into `theta` (which holds n_have planets, by
    default the first j). Start periods: the quoted P and the three highest periodogram peaks
    of the current residuals within +-3 sigma of it (a short-period companion's likelihood has
    many narrow peaks, and the quoted P is not always on the best one). Each start has the
    phase from the residuals' sinusoid and several (e, omega) starts (the quoted ones if known);
    the best by posterior is returned. Other parameters stay fixed here; a joint fit follows."""
    n_have = j if n_have is None else n_have
    _, _, res = neg_lnL(theta, data, n_have, tord, True)
    _, _, jit = unpack(theta, n_have, tord, data.nds)
    L = lit[j]
    mu, sd = L['lnP']
    lo, hi = math.exp(mu - 3*sd), math.exp(mu + 3*sd)
    periods = 1/np.arange(1/hi, 1/lo + 1e-12, 1/(OVERSAMPLE*data.span))[::-1]
    p_starts = [L['P0']]
    if len(periods) >= 3:
        power = periodogram(data, res, jit, periods)
        pk = argrelmax(power)[0]
        p_starts += [float(periods[i]) for i in pk[np.argsort(power[pk])[::-1]][:3]]
    e_start = max(L['e0'] or 0.0, 0.01)
    oms = [L['omega0']] if L['omega0'] is not None else [0.0, np.pi/2, np.pi, 1.5*np.pi]
    ecc = [(math.sqrt(e_start)*math.cos(om), math.sqrt(e_start)*math.sin(om)) for om in oms]
    ecc += [(0.05, 0.05)]
    k0 = PARAMS_PER_PLANET*j
    pri = lit[:n_have + 1]
    bnds = bounds_for(n_have + 1, tord, data.nds, [P_MIN_KNOWN_D]*(n_have + 1), data.pmax,
                      [max(data.pmax, math.exp(l['lnP'][0] + 4*l['lnP'][1])) for l in pri]
                      )[k0:k0 + PARAMS_PER_PLANET]
    best = None
    for P0 in p_starts:
        A0, phi = sinusoid_phase(data, res, jit, P0)
        K_start = L['K0'] or max(A0, 0.05)
        for s1, s2 in ecc:
            w = math.atan2(s2, s1)
            new = np.array([math.log(P0), math.log(K_start), s1, s2, -phi - w])
            full = np.concatenate([theta[:k0], new, theta[k0:]])

            def f(x, full=full):
                th = full.copy()
                th[k0:k0 + PARAMS_PER_PLANET] = x
                return neg_lnL(th, data, n_have + 1, tord) + prior_penalty(th, pri)
            r = minimize(f, new, method='L-BFGS-B', bounds=bnds, options={'maxiter': 200})
            if best is None or r.fun < best[0]:
                th = full.copy()
                th[k0:k0 + PARAMS_PER_PLANET] = r.x
                best = (r.fun, th)
    return best[1]


def polish_literature(theta, npl, tord, data, lit, rounds=2):
    """Coordinate polish: re-seed each literature companion in turn (see
    seed_literature_planet) with everything else fixed, then refit jointly; keep the result
    only if the posterior improves. Returns (theta, data lnL)."""
    floors = [P_MIN_KNOWN_D]*npl
    theta, lnL = fit(theta, data, npl, tord, floors, lit)
    post = lnL - prior_penalty(theta, lit)
    for _ in range(rounds):
        improved = False
        for j in range(len(lit)):
            k = PARAMS_PER_PLANET*j
            th = seed_literature_planet(np.delete(theta, np.s_[k:k + PARAMS_PER_PLANET]), j, data,
                                        lit, tord, npl - 1)
            th, ll = fit(th, data, npl, tord, floors, lit)
            if ll - prior_penalty(th, lit) > post + 1e-6:
                theta, lnL, post, improved = th, ll, ll - prior_penalty(th, lit), True
        if not improved:
            break
    return theta, lnL


def rebuild_with_trend(data, lit, theta, old, new):
    """The literature companions refit from scratch with a trend of order `new` already in the
    model (the trend term seeded by least squares, then each companion re-seeded at its quoted
    period). `theta` holds exactly the literature companions at trend order `old`. An
    alternative start for a trend raise: a refit from the old solution can stay in a poor
    optimum, since the periods of short-period companions have many narrow peaks."""
    n = len(lit)
    th = seed_trend(theta, n, old, new, data)
    coefs_jit = th[PARAMS_PER_PLANET*n:]
    for j in range(n):
        coefs_jit = seed_literature_planet(coefs_jit, j, data, lit, new)
    return polish_literature(coefs_jit, n, new, data, lit)


def initial_theta(npl, tord, data):
    r = data.rv.copy()
    for idx in data.members:
        r[idx] -= np.mean(r[idx])
    jit0 = np.log(max(np.std(r)*0.3, JITTER_BOUNDS[0]*2))
    return np.concatenate([np.zeros(PARAMS_PER_PLANET*npl), [0.0]*tord,
                           np.full(data.nds, min(jit0, np.log(JITTER_BOUNDS[1])*0.9))])


def synthetic(seed, planets, noise=3.0, n=120, nds=2, span=2500.0, trend=0.0, curv=0.0):
    """Synthetic RVs: Keplerian planets (P, K, e, omega, M0), per-dataset offsets and
    jitter, heteroscedastic errors, a random observing pattern, and an optional linear
    (m/s per day) and quadratic (m/s per day^2) trend."""
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, span, n))
    ds = rng.integers(0, nds, n)
    tref = 0.5*(t.min() + t.max())
    rv = trend*(t - tref) + curv*(t - tref)**2
    for P, K, e, om, M0 in planets:
        s1, s2 = np.sqrt(e)*np.cos(om), np.sqrt(e)*np.sin(om)
        rv = rv + keplerian(t, tref, np.log(P), np.log(K), s1, s2, M0)
    err = rng.uniform(1.0, 3.0, n)
    rv = rv + rng.normal(0, np.sqrt(err**2 + noise**2)) + 50*ds*rng.choice([-1, 1])
    return RVData(t, rv, err, np.array([f'inst{d}' for d in ds]))



def selftest():
    """kepler.c against the numpy solver; two companions recovered under literature
    priors that are slightly off, with every fitted e below E_MAX."""
    t0 = time.time()
    if kepler_c.AVAILABLE:
        rng = np.random.default_rng(0)
        t = np.sort(rng.uniform(0, 6000, 410))
        for P, K, s1, s2, M0 in [(30, 10, .2, .1, 1), (0.5, 5, .9, .3, 2), (800, 50, 0, 0, .3),
                                 (12.3, 7, .99, .1, 5), (3.1, 2, .05, -.02, 0)]:
            c_out = add_keplerian(np.zeros_like(t), t, 3000.0, np.log(P), np.log(K), s1, s2, M0)
            ref = keplerian(t, 3000.0, np.log(P), np.log(K), s1, s2, M0)
            assert np.max(np.abs(c_out - ref)) < 1e-8*K, (P, K, s1, s2)
        print('kepler.c agrees with the numpy solver')
    else:
        print('kepler.c unavailable (no C compiler): numpy solver only')
    d = synthetic(1, [(11.7, 25.0, 0.1, 1.0, 0.5), (430.0, 15.0, 0.3, 2.0, 2.0)])
    lit = sorted([make_prior(11.7*1.002, 27.0, 0.1), make_prior(430.0*0.99, 14.0, 0.3)],
                 key=lambda L: -L['K0'])
    theta = initial_theta(0, 0, d)
    for j in range(len(lit)):
        theta = seed_literature_planet(theta, j, d, lit)
    theta, lnL = polish_literature(theta, len(lit), 0, d, lit)
    pl = unpack(theta, len(lit), 0, d.nds)[0]
    P = sorted(float(np.exp(p[0])) for p in pl)
    assert abs(P[0]/11.7 - 1) < 0.01 and abs(P[1]/430 - 1) < 0.03, P
    assert all(p[2]**2 + p[3]**2 < E_MAX for p in pl), pl
    print(f'rv_model.py test passed in {time.time() - t0:.0f} s: periods {P[0]:.2f}, {P[1]:.1f} d')


if __name__ == '__main__':
    if '--test' in sys.argv:
        selftest()
