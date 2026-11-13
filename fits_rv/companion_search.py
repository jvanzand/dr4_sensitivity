"""RV-only companion search: periodogram, maximum likelihood, and a Delta-BIC test.

For one system's RVs (rv_data/systems/<system_id>_<name>.csv) this answers: how many
companions do the RVs support beyond the literature's? The scheme was decided on
2026-11-05/06 and revised on 2026-11-09:

0. By default the search is fully blind: no literature companions (the toggle is
   `--literature`; decided 2026-11-09, so the result is comparable to what Gaia sees, and
   the literature can then be compared with what the RVs alone recover). With the toggle on:

1. The N companions the literature reports (the census count) are in the model from the
   start, with Gaussian priors on P, K, and e around the quoted values (widths: the published
   error, never narrower than 1% in P, 10% in K, and 0.1 in e; 20% in P and K where no error
   is known; see `literature_priors`). Phase and omega are free. They are refit under these
   priors (maximum a posteriori), not demanded to be found blindly: a transiting planet with
   a small K could never be found from the RVs alone. After the fit, each one's Delta-BIC for
   dropping it is recorded, so a companion held only by its prior is visible.
2. Then test one more companion, found blindly, and adopt it if Delta-BIC >= 20. Repeat until
   the next one is rejected.
3. A polynomial trend (order 0, 1 = linear, or 2 = quadratic) competes at every step. First
   the trend is simplified if its highest term no longer earns its keep (dropping it costs
   Delta-BIC < 20). Then the candidates are: add a companion; add a companion and lower the
   trend order at once (a companion that replaces a trend); raise the trend order. The best
   one with Delta-BIC >= 20 is adopted. So a signal that could be either a planet or a
   trend is not settled by test order.
   A blind companion needs period <= 2x the RV baseline (`PERIOD_CAP_BASELINES`, set from 3x on
   2026-11-09); a longer-period signal can only be fit with the trend terms. (A literature
   companion's period may reach its prior's upper tail beyond that.)
   A blind companion also needs P >= 2 x the mean RV spacing (span/(N-1); decided 2026-11-09):
   with sparser sampling a shorter period is not constrained (aliased), and the fit tends to
   put a spurious, often highly eccentric, signal there. Literature companions are exempt.
   Blind companions are also capped at e <= 0.7 (`E_BLIND_MAX`, decided 2026-11-09): without
   it, long-period "companions" at e ~ 0.99 with periastron in a data gap mimic curvature
   (HD 111232, TOI-1736). Literature companions keep e <= 0.99 and their e prior.
4. Only maximum-likelihood or maximum-a-posteriori fits are used; no MCMC. In the BIC, lnL is
   the data log-likelihood at the fitted point (the prior terms are not in it), k counts every
   parameter, and N counts RVs only.

Delta-BIC = 2 (lnL_new - lnL_old) - (extra parameters) ln N_RV, with 5 parameters per
Keplerian (P, K, sqrt(e) cos w, sqrt(e) sin w, M0) and 1 per trend term. Each dataset has
its own offset (profiled out exactly) and white-noise jitter (fitted, bounded to
JITTER_BOUNDS, the same range as the MCMC prior). The threshold is fixed, not calibrated
against the data's noise, so sparse series don't get a noise-dependent threshold; the
price is that a peak chosen as the best of thousands of trial periods is less
significant than a single-trial Delta-BIC of 10 suggests (the look-elsewhere effect).
`--selftest` measures that false-add rate on noise-only data.

The periodogram is a generalized Lomb-Scargle: for each trial period the best circular
sinusoid is fit with every dataset's offset projected out, which is linear and fast. It
ranks candidate periods; the Delta-BIC test uses full eccentric fits seeded at the
highest peaks.

Usage (from the project root, in the gaia_dr4_sensitivity environment):
    python -I fits_rv/companion_search.py rv_data/systems/959971546241608576_HD43691.csv
    (blind; --literature adds the literature companions from sample/stage1_catalog.csv as
    priors; --results-dir names the output folder under fits_rv/, default "results")
    python -I fits_rv/companion_search.py --selftest
"""

import argparse
import json
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
import names

DELTA_BIC = 20.0     # raised from 10 on 2026-11-06 (Judah), after the first pilot searches
JITTER_BOUNDS = (0.1, 30.0)          # m/s, as in the MCMC prior
E_MAX = 0.99
P_MIN_SEARCH_D = 1.0                 # extra-companion search floor (decided 2026-11-06)
PERIOD_CAP_BASELINES = 2.0            # blind period cap in RV baselines (2026-11-09; was 3)
E_BLIND_MAX = 0.7                    # eccentricity cap for blind companions (2026-11-09)
NYQUIST_FACTOR = 2.0                # blind companions: P >= this x the mean RV spacing (2026-11-09)
P_MIN_KNOWN_D = 0.05                 # floor for the N forced companions
P_MAX_YR = 200.0
PARAMS_PER_PLANET = 5
MAX_TREND_ORDER = 2                  # linear, then quadratic
N_PEAKS = 8                          # periodogram peaks seeded as full fits per step
ECC_STARTS = [(0.05, 0.05), (0.55, 0.3), (-0.4, 0.5)]   # (sqrt(e) cos w, sqrt(e) sin w)
N_JOINT = 4                          # best single-companion seeds refit jointly
## Prior widths for the literature companions (accepted 2026-11-09)
LIT_FLOOR_P, LIT_FLOOR_K, LIT_SIGMA_E, LIT_DEFAULT_REL = 0.01, 0.10, 0.10, 0.20
CATALOG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'sample', 'stage1_catalog.csv')
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
        ## Blind period floor: twice the mean RV spacing (span/(N-1)), and at least
        ## P_MIN_SEARCH_D. Sampling slower than this cannot pin down a shorter period.
        self.mean_gap = self.span/max(self.n - 1, 1)
        self.p_min_blind = max(P_MIN_SEARCH_D, NYQUIST_FACTOR*self.mean_gap)
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
    nll = 0.5*np.sum(res**2*w + np.log(2*np.pi*var))
    return (nll, off, res) if return_parts else nll


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


def cap_blind_e(theta, n_lit, npl):
    """Rescale (sqrt(e) cos w, sqrt(e) sin w) of any blind companion above the cap onto it."""
    theta = theta.copy()
    for j in range(n_lit, npl):
        k = PARAMS_PER_PLANET*j
        e = theta[k + 2]**2 + theta[k + 3]**2
        if e > E_BLIND_MAX:
            theta[k + 2:k + 4] *= math.sqrt(E_BLIND_MAX/e)
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
    th = cap_blind_e(res.x, len(pri), npl)
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


def trial_periods(data, p_min):
    p_max = data.pmax
    df = 1/(OVERSAMPLE*data.span)
    freq = np.arange(1/p_max, 1/p_min, df)
    return 1/freq[::-1]


def seeds_from_peaks(data, resid, jitter, p_min, n_peaks=N_PEAKS):
    periods = trial_periods(data, p_min)
    power = periodogram(data, resid, jitter, periods)
    peaks = argrelmax(power)[0]
    peaks = peaks[np.argsort(power[peaks])[::-1]][:n_peaks]
    w = 1/(data.err**2 + jitter[data.ds]**2)
    amp = np.sqrt(max(np.average((resid - np.average(resid, weights=w))**2, weights=w), 0.01))
    return [(periods[k], power[k]) for k in peaks], amp


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


def add_planet(theta, npl, tord, data, p_min, pri=None):
    """Best (npl + 1)-companion maximum-likelihood fit. Seeds sit at the highest
    periodogram peaks of the current residuals (period, amplitude, and phase from the
    sinusoid fit), with several eccentricity starts each. Each seed is first fit
    alone (the new companion's five parameters, everything else fixed), and the best few
    are then refit jointly (the literature companions under their priors, `pri`). Returns
    (theta, lnL)."""
    _, off, res = neg_lnL(theta, data, npl, tord, True)
    _, _, jit = unpack(theta, npl, tord, data.nds)
    peaks, amp = seeds_from_peaks(data, res, jit, p_min)
    k0 = PARAMS_PER_PLANET*npl
    floors = [P_MIN_KNOWN_D]*npl + [p_min]
    seeds = []
    for P0, _ in peaks:
        A0, phi = sinusoid_phase(data, res, jit, P0)
        for (s1, s2) in ECC_STARTS:
            w = np.arctan2(s2, s1)
            ## circular limit: rv = K cos(x + M0 + w) = A cos(x - phi)  =>  M0 = -phi - w
            new = np.array([np.log(P0), np.log(max(A0, 0.05)), s1, s2, -phi - w])
            full = np.concatenate([theta[:k0], new, theta[k0:]])
            # fit only the new companion's parameters
            def f(x, full=full):
                th = full.copy()
                th[k0:k0 + PARAMS_PER_PLANET] = x
                return neg_lnL(th, data, npl + 1, tord) + 1e5*max(x[2]**2 + x[3]**2
                                                                  - E_BLIND_MAX, 0.0)**2
            bnds = bounds_for(npl + 1, tord, data.nds, floors, data.pmax)[k0:k0 + PARAMS_PER_PLANET]
            r = minimize(f, new, method='L-BFGS-B', bounds=bnds, options={'maxiter': 200})
            th = full.copy()
            th[k0:k0 + PARAMS_PER_PLANET] = r.x
            seeds.append((-r.fun, th))
    seeds.sort(key=lambda s: -s[0])
    best = (None, -np.inf)
    for _, th in seeds[:N_JOINT]:
        th, ll = fit(th, data, npl + 1, tord, floors, pri)
        if ll > best[1]:
            best = (th, ll)
    return best


def summarize(theta, npl, tord, data):
    pl, coefs, jit = unpack(theta, npl, tord, data.nds)
    rows = []
    for p in pl:
        e = min(p[2]**2 + p[3]**2, E_MAX)
        rows.append({'P_d': float(np.exp(p[0])), 'K_ms': float(np.exp(p[1])),
                     'e': float(e), 'omega_deg': float(np.degrees(np.arctan2(p[3], p[2])) % 360),
                     'baseline_over_P': float(data.span/np.exp(p[0]))})
    rows.sort(key=lambda r: r['P_d'])
    lo, hi = JITTER_BOUNDS
    return {'planets': rows, 'trend_order': int(tord),
            ## per day and per day^2, in physical units
            'slope_ms_per_d': float(coefs[0]/data.span) if tord >= 1 else 0.0,
            'curvature_ms_per_d2': float(coefs[1]/data.span**2) if tord >= 2 else 0.0,
            'jitter_ms': dict(zip(data.names, map(float, jit))),
            'jitter_at_bound': [n for n, j in zip(data.names, jit)
                                if j < lo*1.01 or j > hi*0.99]}


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


def literature_priors(system_id, catalog=CATALOG):
    """Priors for the system's literature companions (sample/stage1_catalog.csv), strongest
    K first. K comes from the catalog, else from M sin i, P, e, and M_* (flagged). The relative
    K width uses the mass error (the catalog carries no K errors). Notes record the source of
    each value."""
    cat = pd.read_csv(catalog, dtype={'system_id': str}, low_memory=False)
    rows = cat[cat.system_id == str(system_id)]
    out = []
    for _, r in rows.iterrows():
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
        out.append(make_prior(P0, K0, e0,
                              math.radians(r.omega_deg) if pd.notna(r.omega_deg) else None,
                              sigP, sigK, r.pl_name, notes))
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


def fit_literature(data, lit, theta):
    """The base model: every literature companion under its priors, no trend. Returns
    (theta, data lnL)."""
    for j in range(len(lit)):
        theta = seed_literature_planet(theta, j, data, lit)
    return polish_literature(theta, len(lit), 0, data, lit)


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


def floors_for(npl, n_known, data):
    """Period floors of the model's companions: the first `n_known` are literature companions
    (P_MIN_KNOWN_D), the rest were added blindly (the sampling floor `data.p_min_blind`)."""
    k = min(npl, n_known)
    return [P_MIN_KNOWN_D]*k + [data.p_min_blind]*(npl - k)


def support_check(theta, npl, tord, data, lit, lnL):
    """For each literature companion: Delta-BIC of dropping it (refit without it, the others
    still under their priors). Small values mean the RVs do not require it and the prior is
    what holds it."""
    pl = unpack(theta, npl, tord, data.nds)[0]
    out = []
    for j, L in enumerate(lit):
        k = PARAMS_PER_PLANET*j
        th = np.delete(theta, np.s_[k:k + PARAMS_PER_PLANET])
        th, ll = fit(th, data, npl - 1, tord, floors_for(npl - 1, len(lit) - 1, data),
                     lit[:j] + lit[j + 1:])
        d = bic(ll, npl - 1, tord, data) - bic(lnL, npl, tord, data)
        e = min(pl[j][2]**2 + pl[j][3]**2, E_MAX)
        out.append({'name': L['name'], 'P_lit': L['P0'], 'P_fit': float(np.exp(pl[j][0])),
                    'K_lit': L['K0'], 'K_fit': float(np.exp(pl[j][1])), 'e_lit': L['e0'],
                    'e_fit': float(e), 'dBIC_drop': float(d), 'supported': bool(d >= DELTA_BIC)})
    return out


def initial_theta(npl, tord, data):
    r = data.rv.copy()
    for idx in data.members:
        r[idx] -= np.mean(r[idx])
    jit0 = np.log(max(np.std(r)*0.3, JITTER_BOUNDS[0]*2))
    return np.concatenate([np.zeros(PARAMS_PER_PLANET*npl), [0.0]*tord,
                           np.full(data.nds, min(jit0, np.log(JITTER_BOUNDS[1])*0.9))])


def search(data, lit, max_extra=8, verbose=True):
    """Run the scheme in the module docstring; `lit` is the list of literature priors (see
    `literature_priors`; empty for a fully blind search). Returns the list of steps (each
    with the fitted parameter vector, for plotting) and the final model summary, which
    includes the literature companions' support check. dBIC in the steps is BIC(before) -
    BIC(after): positive when the new model is better."""
    log = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
    N = data.n
    lit = sorted(lit, key=lambda L: -(L['K0'] or 0.0))      # fixed order = order in theta
    n_known = len(lit)
    npl, tord = 0, 0
    theta = initial_theta(0, 0, data)
    theta, lnL = fit(theta, data, 0, 0, [])
    steps = []

    def record(action, th, ll, n, tr, adopted, forced, d_bic):
        s = {'step': len(steps), 'action': action, 'n_planets': n, 'trend_order': tr,
             'lnL': ll, 'BIC': bic(ll, n, tr, data), 'dBIC': d_bic, 'adopted': adopted,
             'forced': forced, 'theta': [float(x) for x in th], 'tref': data.tref,
             'span': data.span, 'n_rv': N, 'datasets': data.names}
        s.update(summarize(th, n, tr, data))
        steps.append(s)
        log(f'  step {s["step"]}: {action}: dBIC = {d_bic:.1f}  '
            f'{"ADOPTED" if adopted else "rejected"}' + (' (forced)' if forced else '')
            + f'  trend order {tr}  periods ' + ', '.join(f'{p["P_d"]:.2f}' for p in s['planets']))

    log(f'N_RV = {N}, datasets = {data.names}, baseline span = {data.span:.0f} d, '
        f'period cap = {data.pmax:.0f} d, blind-period floor = {data.p_min_blind:.1f} d '
        f'({NYQUIST_FACTOR:g} x mean RV spacing of {data.mean_gap:.1f} d), '
        f'{n_known} literature companions under priors')
    record('baseline (offsets and jitter only)', theta, lnL, 0, 0, True, False, np.nan)

    def settle_trend(theta, lnL, npl, tord):
        """Trend moves until stable: drop the top term if it is worth < DELTA_BIC, raise the
        order if that gains >= DELTA_BIC. Without them, curvature from a long-period signal
        would be fit as a spiky high-eccentricity companion (seen in TOI-1736)."""
        for _ in range(6):
            cur = bic(lnL, npl, tord, data)
            if tord > 0:
                th, ll = fit(set_trend_order(theta, npl, tord, tord - 1), data, npl, tord - 1,
                             floors_for(npl, n_known, data), lit)
                d = cur - bic(ll, npl, tord - 1, data)
                if d > -DELTA_BIC:
                    theta, lnL, tord = th, ll, tord - 1
                    record(f'drop trend term (order {tord + 1} -> {tord}); it was worth only '
                           f'dBIC = {-d:.1f}', theta, lnL, npl, tord, True, False, d)
                    continue
            if tord < MAX_TREND_ORDER:
                th, ll = fit(seed_trend(theta, npl, tord, tord + 1, data), data, npl, tord + 1,
                             floors_for(npl, n_known, data), lit)
                if npl == n_known:      # only the literature companions so far: second start
                    th2, ll2 = rebuild_with_trend(data, lit, theta, tord, tord + 1)
                    if ll2 - prior_penalty(th2, lit) > ll - prior_penalty(th, lit):
                        th, ll = th2, ll2
                d = cur - bic(ll, npl, tord + 1, data)
                if d >= DELTA_BIC:
                    theta, lnL, tord = th, ll, tord + 1
                    record(f'raise trend order {tord - 1} -> {tord}', theta, lnL, npl, tord,
                           True, False, d)
                    continue
            break
        return theta, lnL, tord

    ## 1. The literature companions, under their priors, then the trend moves
    if n_known:
        before = bic(lnL, 0, 0, data)
        theta, lnL = fit_literature(data, lit, theta)
        npl = n_known
        record(f'{n_known} literature companions under Gaussian priors', theta, lnL, npl, 0,
               True, True, before - bic(lnL, npl, 0, data))
        theta, lnL, tord = settle_trend(theta, lnL, npl, tord)

    ## 2-3. Extra (blind) companions and the trend
    for _ in range(3*max_extra):
        cur = bic(lnL, npl, tord, data)
        ## First simplify: drop the highest trend term if it no longer pays for itself
        if tord > 0:
            th0 = set_trend_order(theta, npl, tord, tord - 1)
            th, ll = fit(th0, data, npl, tord - 1, floors_for(npl, n_known, data), lit)
            d = cur - bic(ll, npl, tord - 1, data)
            if d > -DELTA_BIC:
                theta, lnL, tord = th, ll, tord - 1
                record(f'drop trend term (order {tord + 1} -> {tord}); it was worth only '
                       f'dBIC = {-d:.1f}', theta, lnL, npl, tord, True, False, d)
                continue
        ## Then the additive candidates
        cands = []
        blind_ok = data.p_min_blind < 0.5*data.pmax     # room between the floor and the cap
        if blind_ok:
            th, ll = add_planet(theta, npl, tord, data, data.p_min_blind, lit)
            cands.append(('add companion', th, ll, npl + 1, tord))
        if tord > 0 and blind_ok:        # a companion that replaces a trend term
            th0 = set_trend_order(theta, npl, tord, tord - 1)
            th, ll = add_planet(th0, npl, tord - 1, data, data.p_min_blind, lit)
            cands.append((f'add companion and drop trend order {tord} -> {tord - 1}',
                          th, ll, npl + 1, tord - 1))
        if tord < MAX_TREND_ORDER:
            th0 = seed_trend(theta, npl, tord, tord + 1, data)
            th, ll = fit(th0, data, npl, tord + 1, floors_for(npl, n_known, data), lit)
            cands.append((f'raise trend order {tord} -> {tord + 1}', th, ll, npl, tord + 1))
        if not cands:
            break
        scored = [(cur - bic(c[2], c[3], c[4], data),) + c for c in cands]
        best = max(scored, key=lambda c: c[0])
        for c in scored:
            if c is not best:
                log(f'  (also tried: {c[1]}, dBIC = {c[0]:.1f})')
        d, action, th, ll, n_new, t_new = best
        if d < DELTA_BIC:
            record(action, th, ll, n_new, t_new, False, False, d)
            break
        theta, lnL, npl, tord = th, ll, n_new, t_new
        record(action, theta, lnL, npl, tord, True, False, d)

    final = summarize(theta, npl, tord, data)
    if n_known:
        final['literature_fit'] = support_check(theta, npl, tord, data, lit, lnL)
        for r in final['literature_fit']:
            log(f'  literature {r["name"]}: P {r["P_lit"]:.4g} -> {r["P_fit"]:.4g} d, '
                f'dBIC for dropping it = {r["dBIC_drop"]:.1f}'
                + ('' if r['supported'] else '  (held by its prior)'))
    return steps, final


def load_system(path):
    t = pd.read_csv(path, low_memory=False)
    ok = ~t.isjunk.astype(bool) & t.rv_ms.notna() & t.rv_err_ms.notna()
    t = t[ok]
    return RVData(t.bjd_tcb.to_numpy(), t.rv_ms.to_numpy(), t.rv_err_ms.to_numpy(),
                  t.dataset.to_numpy())


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


def selftest(n_noise_trials=10):
    t0 = time.time()
    print('0. kepler.c against the numpy solver' if kepler_c.AVAILABLE
          else '0. kepler.c unavailable (no C compiler): using the numpy solver')
    if kepler_c.AVAILABLE:
        rng = np.random.default_rng(0)
        t = np.sort(rng.uniform(0, 6000, 410))
        for P, K, s1, s2, M0 in [(30, 10, .2, .1, 1), (0.5, 5, .9, .3, 2), (800, 50, 0, 0, .3),
                                 (12.3, 7, .99, .1, 5), (3.1, 2, .05, -.02, 0)]:
            c_out = add_keplerian(np.zeros_like(t), t, 3000.0, np.log(P), np.log(K), s1, s2, M0)
            ref = keplerian(t, 3000.0, np.log(P), np.log(K), s1, s2, M0)
            assert np.max(np.abs(c_out - ref)) < 1e-8*K, (P, K, s1, s2)
    print()
    print('1. Two planets with literature priors slightly off, then expect no extra companion')
    d = synthetic(1, [(11.7, 25.0, 0.1, 1.0, 0.5), (430.0, 15.0, 0.3, 2.0, 2.0)])
    steps, final = search(d, [make_prior(11.7*1.002, 27.0, 0.1), make_prior(430.0*0.99, 14.0, 0.3)])
    found = sorted(p['P_d'] for p in final['planets'])
    assert len(found) == 2, final
    assert abs(found[0]/11.7 - 1) < 0.01 and abs(found[1]/430 - 1) < 0.03, found

    print('\n2. One literature planet, a second signal that must be added')
    d = synthetic(2, [(23.4, 20.0, 0.05, 0.5, 1.0), (260.0, 18.0, 0.1, 1.0, 3.0)])
    steps, final = search(d, [make_prior(23.4*0.999, 22.0, 0.05)])
    found = sorted(p['P_d'] for p in final['planets'])
    assert len(found) == 2, final
    assert abs(found[0]/23.4 - 1) < 0.01 and abs(found[1]/260 - 1) < 0.03, found

    print('\n3. Signal plus a linear trend: the linear trend should be adopted')
    d = synthetic(3, [(40.0, 15.0, 0.1, 0.0, 0.0)], trend=0.02)
    steps, final = search(d, [make_prior(40.0*1.001, 14.0, 0.1)])
    assert final['trend_order'] == 1 and abs(final['slope_ms_per_d'] - 0.02) < 0.005, final

    print('\n4. Signal plus a quadratic trend: order 2 should be adopted, with no long-period '
          'companion (the period cap, 2x the baseline, rules one out)')
    d = synthetic(4, [(40.0, 15.0, 0.1, 0.0, 0.0)], trend=0.01, curv=1.2e-5)
    steps, final = search(d, [make_prior(40.0*1.001, 14.0, 0.1)])
    assert final['trend_order'] == 2, final
    assert all(p['baseline_over_P'] >= 1/PERIOD_CAP_BASELINES - 1e-9 for p in final['planets']), final

    print(f'\n5. Noise only, N = 0: false-add rate over {n_noise_trials} seeds')
    adds = 0
    for s in range(n_noise_trials):
        d = synthetic(100 + s, [], noise=3.0)
        _, final = search(d, [], verbose=False)
        adds += int(len(final['planets']) > 0 or final['trend_order'] > 0)
    print(f'   {adds} of {n_noise_trials} noise-only datasets got a spurious companion or trend')
    assert adds <= max(2, n_noise_trials//5), adds
    print('\n6. A literature companion that is not in the data: it must stay near its prior and '
          'be reported as held by the prior')
    d = synthetic(300, [], noise=3.0)
    _, final = search(d, [make_prior(17.0, 5.0, 0.0, sigP=0.05)], verbose=False)
    lf = final['literature_fit'][0]
    assert not lf['supported'] and abs(lf['P_fit']/17.0 - 1) < 0.1, lf

    print('\n7. Sparse sampling: a 5.3 d signal in 30 RVs over 2500 d (mean spacing ~86 d) lies '
          'below the blind floor, so no companion shorter than the floor may be added')
    d = synthetic(7, [(5.3, 15.0, 0.1, 0.5, 1.0)], noise=3.0, n=30)
    assert abs(d.p_min_blind - 2*d.span/(d.n - 1)) < 1e-9, d.p_min_blind
    _, final = search(d, [], verbose=False)
    assert all(p['P_d'] >= d.p_min_blind*0.999 for p in final['planets']), final

    print('\n8. Eccentricity cap: a blind search on an e = 0.9 signal may not return e > 0.7')
    d = synthetic(8, [(60.0, 25.0, 0.9, 0.5, 1.0)], noise=3.0)
    _, final = search(d, [], verbose=False)
    assert final['planets'] and all(p['e'] <= E_BLIND_MAX + 1e-9 for p in final['planets']), final

    print(f'\nselftest passed in {time.time() - t0:.0f} s')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('rv_csv', nargs='?')
    ap.add_argument('--literature', action='store_true',
                    help='put the literature companions in the model under Gaussian priors '
                         '(default: off, a fully blind search)')
    ap.add_argument('--results-dir', default='results',
                    help='output folder under fits_rv/ (default "results")')
    ap.add_argument('--max-extra', type=int, default=8)
    ap.add_argument('--out', help='JSON of the steps (default fits_rv/results/<display name>/'
                                  'companion_search.json)')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest()
        sys.exit()
    if not a.rv_csv:
        ap.error('give an RV table or --selftest')
    t_start = time.time()
    sid_arg = os.path.basename(a.rv_csv).split('_')[0]
    lit = literature_priors(sid_arg) if a.literature else []
    steps, final = search(load_system(a.rv_csv), lit, a.max_extra)
    print(json.dumps(final, indent=1))
    out = a.out
    if out is None:
        sid = os.path.basename(a.rv_csv).split('_')[0]
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), a.results_dir,
                           names.display_name(sid), 'companion_search.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as fh:
        json.dump({'rv_table': a.rv_csv, 'system_id': os.path.basename(a.rv_csv).split('_')[0],
                   'display_name': names.display_name(os.path.basename(a.rv_csv).split('_')[0]),
                   'n_known': len(lit), 'literature': bool(a.literature),
                   'period_cap_baselines': PERIOD_CAP_BASELINES, 'literature_priors': lit,
                   'final': final,
                   'seconds': time.time() - t_start, 'steps': steps}, fh, indent=1,
                  default=float)
    print(f'steps written to {out}')
