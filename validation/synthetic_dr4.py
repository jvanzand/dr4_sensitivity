"""Synthetic Gaia DR4 epoch astrometry of a known RV system, written in the
archive CSV layout that `data.load_gaia_dr4` reads.

Steps (`make_system`):
  1. Reference RV fit (`reference_rv_fit`): the system's literature companions
     fitted to its RVs (rv_data/systems/) by maximum likelihood (or, where
     that has no proper maximum, with the literature priors), with a
     polynomial trend, per-dataset offsets (profiled) and jitters, using the
     likelihood of fits_rv/rv_model.py. The trend order (none by
     default, linear or quadratic) is chosen step by step (`select_trend`):
     tested before any companion is in the model, and re-checked after each
     companion is added (a term that stops paying Delta-BIC >= 20 is
     dropped). The literature-prior fit seeds the final fit, which has no
     priors. No uncertainties are estimated (Judah, 2026-11-21: they come
     from the Octofitter fits later).
  2. Ground truth (`draw_truth`): the best fit gives the RV-based
     parameters (P, K, e, omega, M0); the unknown ones are drawn at
     random: cos i uniform in [-1, 1] (isotropic, i in [0, 180] deg) and Omega
     uniform in [0, 360) deg. The true mass follows from K, P, e, i and the
     host mass (dark companion: the photocentre is the star). If the fit has
     a trend, a trend companion is drawn to reproduce it
     (`draw_trend_companion`) and its acceleration is injected too: it stands
     in for the long-term motion that real DR4 data will carry, and is not a
     target of the search.
  3. DR4 epochs (`scan_epochs`): gaiamock's scanning law
     (`get_gost_one_position`, El-Badry et al. 2024) at the star's position,
     DR4 window; a random 10% of transits removed, as gaiamock does.
  4. Astrometry: the star's along-scan position = position and proper motion
     (Gaia DR3, reference epoch J2016.0) + parallax x along-scan parallax
     factor + the companions' orbits, plus noise: per CCD, white noise at
     gaiamock's per-CCD error sigma_ccd(G) (Holl et al. 2023 "EDR3 adjusted"
     curve); per transit, gaiamock's extra unmodelled noise for G < 13 (one
     sigma per dataset, uniform in [0, 0.04] mas, not in the formal errors).
     Each transit has 9 CCD rows (8 on CCD row 4), matching the pre-release
     DR4 data of Gaia-4 and HD 114762 (9 used CCDs per transit; 8 in about
     1/7 of transits).
The zero-planet dataset uses the same epochs and noise draw as the dataset
with the companions, so the two differ only by the injected orbits.

Conventions are those of orbit.py and injection.py: omega is the star's (the
RV omega), Thiele-Innes constants of the star's orbit. The sign link between
the RV direction and the sky orientation is not verified (injection.py), which
does not matter while i and Omega are drawn at random.

Outputs (validation/benchmarks/, in git: the frozen benchmark inputs):
<name>_seed<k>.csv (with companions), <name>_seed<k>_noplanet.csv,
<name>_seed<k>_truth.json. An existing dataset is never overwritten (`force`).
Figure (`process_figure`): validation/results/figures/synthetic_<name>_seed<k>.png
(not in git); legends
are placed by figtools.legend_clear so they cover no data. `--figure-only`
redraws it from the saved files (the truth JSON keeps the fit's theta).

Run:        python -I validation/synthetic_dr4.py --system HD39091 --seed 1
            (usually through validation/benchmark.py generate)
Test mode:  python -I validation/synthetic_dr4.py --test
"""

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _d in ('blind_search', 'injection', 'fits_rv'):   # core, injection, RV likelihood
    sys.path.insert(0, os.path.join(ROOT, _d))
sys.path.insert(0, HERE)
import rv_model as cs   # noqa: E402
import data as bsdata          # noqa: E402
import figtools                # noqa: E402
import injection               # noqa: E402
import orbit                   # noqa: E402

BJD_MJD = 2400000.5
J2016_MJD = 57388.5            # Gaia DR3 reference epoch J2016.0 (TCB)
GAIA_CACHE = os.path.join(ROOT, 'data', 'gaia_dr3_astrometry.csv')
OUT_DIR = os.path.join(HERE, 'benchmarks')           # benchmark inputs (in git)
FIG_DIR = os.path.join(HERE, 'results', 'figures')   # not in git
REJECT_FRAC = 0.10             # gaiamock: a random 10% of transits removed
EXTRA_NOISE_MAX = 0.04         # mas; gaiamock's unmodelled noise, G < 13
N_CCD, N_CCD_ROW4 = 9, 8


# ---- inputs --------------------------------------------------------------------------

def host_astrometry(source_id):
    """Gaia DR3 position, parallax, proper motion and G of one source, cached
    in data/gaia_dr3_astrometry.csv (fetched from the Gaia archive once)."""
    sid = int(source_id)
    if os.path.exists(GAIA_CACHE):
        c = pd.read_csv(GAIA_CACHE)
        if sid in set(c.source_id):
            return c[c.source_id == sid].iloc[0].to_dict()
    from astroquery.gaia import Gaia
    r = Gaia.launch_job(
        'select source_id, ra, dec, parallax, pmra, pmdec, ref_epoch, phot_g_mean_mag '
        f'from gaiadr3.gaia_source where source_id = {sid}').get_results().to_pandas()
    r.columns = [x.lower() for x in r.columns]
    os.makedirs(os.path.dirname(GAIA_CACHE), exist_ok=True)
    old = pd.read_csv(GAIA_CACHE) if os.path.exists(GAIA_CACHE) else None
    pd.concat([old, r]).to_csv(GAIA_CACHE, index=False)
    return r.iloc[0].to_dict()




def load_rvs(path, system, nightly=True):
    """The system's RVs as the RV model's data (`rv_model.RVData`), read by the
    shared loader `data.read_rv_table` in blind_search/data.py (junk flags, the
    RVs listed in rv_data/excluded_rvs.csv, nightly binning unless `nightly`
    is False). Returns (RVData, number excluded)."""
    t, y, e, ds, n_ex = bsdata.read_rv_table(path, system, nightly)
    return cs.RVData(t, y, e, ds), n_ex


def rv_table(system):
    """Path of the system's RV table (rv_data/systems/<source_id>_<name>.csv)."""
    d = os.path.join(ROOT, 'rv_data', 'systems')
    key = system.replace('-', '')                       # TOI1736 matches ..._TOI-1736.csv
    hits = [f for f in os.listdir(d) if f.replace('-', '').endswith(f'_{key}.csv')]
    if len(hits) != 1:
        raise FileNotFoundError(f'{system}: {hits}')
    return os.path.join(d, hits[0])


# ---- 1. reference RV fit ---------------------------------------------------------------

def axis_steps(f, x, hold=()):
    """Per parameter, the step along which f rises by about 0.5 (one sigma),
    by central differences. A parameter along which f stays flat (rise < 0.05
    even at a step of 1e3), or one listed in `hold`, is listed in `flat`.
    Returns (steps, flat)."""
    n, f0 = len(x), f(x)
    h = np.full(n, 1e-3)
    flat = list(hold)
    for i in range(n):
        if i in flat:
            continue
        e = np.zeros(n)
        for _ in range(40):
            e[i] = h[i]
            d = 0.5*(f(x + e) + f(x - e)) - f0
            if d < 0.05 and h[i] < 1e3:
                h[i] *= 3
            elif d > 2:
                h[i] /= 3
            else:
                break
        if d < 0.05:
            flat.append(i)
    return h, flat


def polish_scaled(f, x, hold=(), rounds=3):
    """Minimize f in coordinates scaled by `axis_steps` (each parameter in
    units of its own one-sigma step), so that parameters whose scales differ
    by orders of magnitude (ln P of a short-period planet vs a jitter)
    converge together. Flat and `hold` parameters stay fixed."""
    from scipy.optimize import minimize
    for _ in range(rounds):
        h, flat = axis_steps(f, x, hold)
        free = np.array([i for i in range(len(x)) if i not in flat])
        x0 = x.copy()

        def g(z):
            y = x0.copy()
            y[free] += z*h[free]
            return f(y)
        r = minimize(g, np.zeros(len(free)), method='BFGS', options={'gtol': 1e-5, 'maxiter': 5000})
        x = x0.copy()
        x[free] += r.x*h[free]
    return x


def to_lambda(theta, npl):
    """Replace each companion's M0 by the mean longitude lambda0 = M0 + omega
    (well defined as e -> 0, where omega and M0 are degenerate)."""
    phi = np.array(theta, float)
    for j in range(npl):
        k = cs.PARAMS_PER_PLANET*j
        phi[k + 4] = theta[k + 4] + math.atan2(theta[k + 3], theta[k + 2])
    return phi


def from_lambda(phi, npl):
    theta = np.array(phi, float)
    for j in range(npl):
        k = cs.PARAMS_PER_PLANET*j
        theta[k + 4] = phi[k + 4] - math.atan2(phi[k + 3], phi[k + 2])
    return theta


def literature_sorted(source_id):
    """The system's literature companions (rv_model priors), largest K
    first: the order in which the reference fit adds them, and the order of
    the companions in its theta."""
    return sorted(cs.literature_priors(source_id), key=lambda L: -(L['K0'] or 0.0))


def select_trend(rv, lit):
    """Trend order of the reference fit, chosen step by step (Judah,
    2026-11-16): no trend by default; with no companions yet, a linear and
    then a quadratic term are each adopted if they lower BIC by >=
    cs.DELTA_BIC (20). The literature companions are then added one at a time
    (largest K first), and after each the trend is re-checked: its top term
    is dropped if it is worth < DELTA_BIC, and the order raised if that gains
    >= DELTA_BIC. Fits are rv_model.fit (literature priors on the
    companions, BIC from the data lnL). Returns (theta, order, steps)."""
    steps = []
    floors = lambda n: [cs.P_MIN_KNOWN_D]*n        # noqa: E731

    def settle(theta, lnL, n, tord):
        for _ in range(2*cs.MAX_TREND_ORDER + 2):
            cur = cs.bic(lnL, n, tord, rv)
            if tord > 0:
                th, ll = cs.fit(cs.set_trend_order(theta, n, tord, tord - 1), rv, n, tord - 1,
                                floors(n), lit[:n])
                d = cur - cs.bic(ll, n, tord - 1, rv)
                if d > -cs.DELTA_BIC:
                    steps.append(dict(n_companions=n, action=f'drop trend term {tord} -> {tord - 1}',
                                      worth_dBIC=float(-d)))
                    theta, lnL, tord = th, ll, tord - 1
                    continue
            if tord < cs.MAX_TREND_ORDER:
                th, ll = cs.fit(cs.seed_trend(theta, n, tord, tord + 1, rv), rv, n, tord + 1,
                                floors(n), lit[:n])
                if n:                                  # second start: companions refit with the trend
                    th2, ll2 = cs.rebuild_with_trend(rv, lit[:n], theta, tord, tord + 1)
                    if ll2 - cs.prior_penalty(th2, lit[:n]) > ll - cs.prior_penalty(th, lit[:n]):
                        th, ll = th2, ll2
                d = cur - cs.bic(ll, n, tord + 1, rv)
                if d >= cs.DELTA_BIC:
                    steps.append(dict(n_companions=n, action=f'raise trend order {tord} -> {tord + 1}',
                                      dBIC=float(d)))
                    theta, lnL, tord = th, ll, tord + 1
                    continue
                steps.append(dict(n_companions=n, action=f'keep trend order {tord}',
                                  dBIC_of_raising=float(d)))
            else:
                steps.append(dict(n_companions=n, action=f'keep trend order {tord}',
                                  top_term_worth_dBIC=float(-d)))
            break
        return theta, lnL, tord

    theta, lnL = cs.fit(cs.initial_theta(0, 0, rv), rv, 0, 0, [])
    theta, lnL, tord = settle(theta, lnL, 0, 0)
    for j in range(len(lit)):
        theta = cs.seed_literature_planet(theta, j, rv, lit, tord)
        theta, lnL = cs.polish_literature(theta, j + 1, tord, rv, lit[:j + 1])
        steps.append(dict(n_companions=j + 1, action=f"add {lit[j]['name']}"))
        theta, lnL, tord = settle(theta, lnL, j + 1, tord)
    return theta, tord, steps


TREND_REPLACES_LONG = False     # option (Judah, 2026-11-20): may a trend replace a known companion with P > RV span?


def replace_long_by_trend(rv, lit, theta, tord):
    """Judah, 2026-11-20 (option, `TREND_REPLACES_LONG`): a known companion
    whose literature period exceeds the RV span may be replaced by a polynomial
    trend. For each such companion (longest period first), the model without
    it is rebuilt by `select_trend` (trend order at least 1: the trend stands
    for the companion) and replaces the current one if its BIC (data lnL) is
    lower. The replaced companion then enters the truth only through the trend
    (drawn as a trend companion), not as a target. Returns (lit, theta, tord,
    steps, tests)."""
    floors = lambda n: [cs.P_MIN_KNOWN_D]*n        # noqa: E731
    tests, steps = [], None
    for L in sorted([L for L in lit if (L['P0'] or 0) > rv.span], key=lambda L: -L['P0']):
        n = len(lit)
        bic_full = float(cs.bic(-cs.neg_lnL(theta, rv, n, tord), n, tord, rv))
        lit_r = [M for M in lit if M is not L]
        th_r, t_r, st_r = select_trend(rv, lit_r)
        if t_r == 0:
            th_r, _ = cs.fit(cs.seed_trend(th_r, n - 1, 0, 1, rv), rv, n - 1, 1, floors(n - 1), lit_r)
            t_r = 1
        bic_r = float(cs.bic(-cs.neg_lnL(th_r, rv, n - 1, t_r), n - 1, t_r, rv))
        adopt = bic_r < bic_full
        tests.append(dict(name=L['name'], P_lit_d=float(L['P0']), span_d=float(rv.span),
                          bic_with_companion=bic_full, bic_trend=bic_r, trend_order=int(t_r),
                          adopted=bool(adopt)))
        print(f"{L['name']} (P {L['P0']:.0f} d > RV span {rv.span:.0f} d): trend of order {t_r} instead, "
              f"BIC {bic_r:.1f} vs {bic_full:.1f} with it: {'replaced' if adopt else 'kept'}")
        if adopt:
            lit, theta, tord, steps = lit_r, th_r, t_r, st_r
    return lit, theta, tord, steps, tests


def profile_jitters(f, th, k, nds, lo, hi, rounds=3):
    """Each dataset's ln jitter set to its maximum of -f with everything else
    held (bounded 1-D search in [lo, hi]), `rounds` times over all datasets.
    L-BFGS-B can leave a jitter at its upper bound far from its maximum (e.g.
    HD 168443 HIRES: 30 m/s with 6 m/s residuals; found 2026-11-16)."""
    from scipy.optimize import minimize_scalar
    th = th.copy()
    for _ in range(rounds):
        for d in range(nds):
            def g(v, d=d):
                x = th.copy()
                x[k + d] = v
                return f(x)
            r = minimize_scalar(g, bounds=(lo, hi), method='bounded', options={'xatol': 1e-4})
            if r.fun < f(th):
                th[k + d] = r.x
    return th


P_EXTRA_MAX_D = cs.P_MAX_YR*365.25    # long-period Keplerian: longest period tried (200 yr)
EXTRA_P_SEEDS = 8                       # its starting periods, log-spaced from the RV baseline
EXTRA_E_MAX = 0.95                      # its eccentricity cap (Judah: an eccentric companion is fine)


def long_period_keplerian(rv, lit, theta, tord):
    """Judah, 2026-11-17: is the RVs' polynomial trend better described by one
    long-period Keplerian? A blind companion (no prior; P from the RV baseline
    to P_EXTRA_MAX_D, e <= EXTRA_E_MAX) is fitted with the literature
    companions (under their priors) at each lower trend order 0 .. tord - 1;
    starts at EXTRA_P_SEEDS periods x 4 phases x 2 eccentricities, the best 4
    refit jointly (as cs.add_planet). Adopted if its BIC is lower than the
    trend model's: the trend already shows the signal is there, so this is a
    choice between two descriptions of it, not a detection. Returns a dict:
    adopt, theta, tord, and the BICs."""
    import copy
    from scipy.optimize import minimize
    n, k0 = len(lit), cs.PARAMS_PER_PLANET*len(lit)
    rvk = copy.copy(rv)
    rvk.pmax = P_EXTRA_MAX_D                      # cs.fit caps blind periods at data.pmax
    bic_trend = float(cs.bic(-cs.neg_lnL(theta, rv, n, tord), n, tord, rv))
    floors = [cs.P_MIN_KNOWN_D]*n + [rv.span]
    k_amp = float(np.ptp(cs.trend_model(theta[k0:k0 + tord], rv)))/2 or 1.0
    best = dict(adopt=False, bic=np.inf)
    e_cap, cs.E_BLIND_MAX = cs.E_BLIND_MAX, EXTRA_E_MAX      # cs.fit's cap on blind companions
    try:
        best = _extra_search(rv, rvk, lit, theta, tord, n, k0, floors, k_amp, best)
    finally:
        cs.E_BLIND_MAX = e_cap
    P = float(np.exp(best['theta'][k0]))
    return dict(adopt=best['bic'] < bic_trend, theta=best['theta'], tord=best['tord'],
                bic_keplerian=best['bic'], bic_trend=bic_trend, trend_order_tested=int(tord),
                P_d=P, K_ms=float(np.exp(best['theta'][k0 + 1])),
                e=float(best['theta'][k0 + 2]**2 + best['theta'][k0 + 3]**2))


def _extra_search(rv, rvk, lit, theta, tord, n, k0, floors, k_amp, best):
    from scipy.optimize import minimize
    for t2 in range(tord):
        base = cs.set_trend_order(theta, n, tord, t2)
        bnd = cs.bounds_for(n + 1, t2, rv.nds, floors, P_EXTRA_MAX_D)[k0:k0 + cs.PARAMS_PER_PLANET]
        seeds = []
        for P0 in np.geomspace(rv.span, P_EXTRA_MAX_D, EXTRA_P_SEEDS):
            for ph in np.arange(4)*np.pi/2:
                for s1, s2 in cs.ECC_STARTS[:2]:
                    new = np.array([math.log(P0), math.log(k_amp), s1, s2, ph])
                    full = np.concatenate([base[:k0], new, base[k0:]])

                    def f(x, full=full):
                        th = full.copy()
                        th[k0:k0 + cs.PARAMS_PER_PLANET] = x
                        return (cs.neg_lnL(th, rv, n + 1, t2)
                                + 1e5*max(x[2]**2 + x[3]**2 - EXTRA_E_MAX, 0.0)**2)
                    r = minimize(f, new, method='L-BFGS-B', bounds=bnd, options={'maxiter': 200})
                    th = full.copy()
                    th[k0:k0 + cs.PARAMS_PER_PLANET] = r.x
                    seeds.append((r.fun, th))
        seeds.sort(key=lambda x: x[0])
        for _, th in seeds[:cs.N_JOINT]:
            th, ll = cs.fit(th, rvk, n + 1, t2, floors, lit)
            b = float(cs.bic(ll, n + 1, t2, rv))
            if b < best['bic']:
                best = dict(theta=th, tord=t2, bic=b)
    return best


def planet_support(f, phi, n, n_lit, tord, rv, bnds):
    """Judah, 2026-11-17: Delta-BIC for removing each literature companion
    (refit without it, same objective `f` in theta; data BIC). Below cs.DELTA_BIC
    (20) the RVs do not detect it (reported; the companion stays in the
    model)."""
    from scipy.optimize import minimize
    theta = from_lambda(phi, n)
    lnl = -cs.neg_lnL(theta, rv, n, tord)
    out = []
    for j in range(n_lit):
        sl = np.s_[cs.PARAMS_PER_PLANET*j:cs.PARAMS_PER_PLANET*(j + 1)]
        th = np.delete(theta, sl)
        b = [x for i, x in enumerate(bnds) if not sl.start <= i < sl.stop]
        th = minimize(lambda x: f(x, n - 1, j), th, method='L-BFGS-B', bounds=b,
                      options={'maxiter': 2000, 'maxfun': 50000}).x
        d = float(cs.bic(-cs.neg_lnL(th, rv, n - 1, tord), n - 1, tord, rv) - cs.bic(lnl, n, tord, rv))
        out.append(dict(index=j, dBIC_drop=d, detected=bool(d >= cs.DELTA_BIC)))
    return out


def reference_rv_fit(system, trend_replaces_long=None):
    """Maximum-likelihood fit of the literature companions (+ a trend chosen by
    `select_trend`, or a long-period Keplerian) to the system's RVs. The
    literature priors only guide the steps before (`select_trend`); the final
    fit has none. No uncertainties: this fit sets the model and best-fit values;
    errors come later from the Octofitter fits (Judah, 2026-11-21; the
    inverse-Hessian estimate used before failed on numerical artifacts, e.g.
    HD 191939). Jitters stay within their bounds and are profiled to their
    maxima (`profile_jitters`); a jitter that has no effect (below 10% of its
    dataset's median error) or sits at the lower bound is held fixed during the
    polish (`held_fixed`). `support`: Delta-BIC for dropping each literature
    companion (detected if >= 20). Returns a dict with the data, the literature
    priors (names, order), theta, data lnL, cov (None), held_fixed, `kind`."""
    path = rv_table(system)
    sid = os.path.basename(path).split('_')[0]
    rv, n_ex = load_rvs(path, system)
    lit = literature_sorted(sid)
    th0, trend_order, trend_steps = select_trend(rv, lit)
    n_extra, extra, replaced = 0, None, []
    if TREND_REPLACES_LONG if trend_replaces_long is None else trend_replaces_long:
        lit, th0, trend_order, st, replaced = replace_long_by_trend(rv, lit, th0, trend_order)
        trend_steps = st or trend_steps
    # is the trend one long-period Keplerian? (Judah, 2026-11-17). Not asked when a known
    # companion was just replaced by the trend: that test already chose between the two.
    if trend_order >= 1 and not any(r['adopted'] for r in replaced):
        extra = long_period_keplerian(rv, lit, th0, trend_order)
        print(f"{system}: long-period Keplerian (P {extra['P_d']/365.25:.1f} yr, K {extra['K_ms']:.3g} m/s, "
              f"e {extra['e']:.2f}, trend order {extra['tord']}) BIC {extra['bic_keplerian']:.1f} vs trend "
              f"(order {trend_order}) {extra['bic_trend']:.1f}: {'adopted' if extra['adopt'] else 'not adopted'}")
        if extra['adopt']:
            th0, trend_order, n_extra = extra['theta'], extra['tord'], 1
        extra = {k: v for k, v in extra.items() if k != 'theta'}
    n = len(lit) + n_extra
    from scipy.optimize import minimize
    bnds = cs.bounds_for(n, trend_order, rv.nds, [cs.P_MIN_KNOWN_D]*len(lit) + [rv.span]*n_extra, rv.pmax,
                         [max(rv.pmax, 2*L['P0']) for L in lit] + [P_EXTRA_MAX_D]*n_extra)
    k = cs.PARAMS_PER_PLANET*n + trend_order
    lo, hi = np.log(cs.JITTER_BOUNDS)
    ks_extra = range(cs.PARAMS_PER_PLANET*len(lit), cs.PARAMS_PER_PLANET*n, cs.PARAMS_PER_PLANET)

    def data_nll(x):         # the long-period Keplerian keeps e <= EXTRA_E_MAX, as in its test
        wall = sum(1e5*max(x[k + 2]**2 + x[k + 3]**2 - EXTRA_E_MAX, 0.0)**2 for k in ks_extra)
        return cs.neg_lnL(x, rv, n, trend_order) + wall
    th = th0.copy()
    for _ in range(3):
        th = minimize(data_nll, th, method='L-BFGS-B', bounds=bnds,
                      options={'maxiter': 2000, 'maxfun': 50000}).x
    jit = np.exp(th[k:])
    # During the polish only jitters without effect or at the LOWER bound are
    # held; one at the upper bound may only reflect a poor starting model
    # (HD 168443 HIRES: 30 m/s from a start at -lnL 2003; 6 m/s at the
    # maximum), so it stays free, and the jitters are profiled between polishes.
    hold = [k + d for d in range(rv.nds)
            if jit[d] < 0.1*np.median(rv.err[rv.ds == d]) or th[k + d] <= lo + 0.01]

    def fl(p):
        x = from_lambda(p, n)
        x[k:] = np.clip(x[k:], lo, hi)                 # jitters stay within bounds
        return data_nll(x)
    phi = polish_scaled(fl, to_lambda(th, n), hold)
    for _ in range(2):
        x = from_lambda(phi, n)
        x[k:] = np.clip(x[k:], lo, hi)
        x = profile_jitters(data_nll, x, k, rv.nds, lo, hi)
        phi = polish_scaled(fl, to_lambda(x, n), hold)
    phi[k:] = np.clip(phi[k:], lo, hi)
    x = cs.cap_e(from_lambda(phi, n), n, n)                     # e < E_MAX (zero prior beyond)
    for k_ in ks_extra:                                         # and the extra one <= EXTRA_E_MAX
        e_ = x[k_ + 2]**2 + x[k_ + 3]**2
        if e_ > EXTRA_E_MAX:
            x[k_ + 2:k_ + 4] *= math.sqrt(EXTRA_E_MAX/e_)
    phi = to_lambda(x, n)
    support = planet_support(lambda x, m, j: cs.neg_lnL(x, rv, m, trend_order),
                             phi, n, len(lit), trend_order, rv, bnds)
    for sp, L in zip(support, lit):
        sp['name'] = L['name']
    weak = [sp['name'] for sp in support if not sp['detected']]
    if weak:
        print(f'{system}: not detected by the RVs (Delta-BIC < {cs.DELTA_BIC:g} for dropping): {weak}')
    theta = from_lambda(phi, n)
    return dict(system=system, source_id=sid, rv=rv, lit=lit, npl=n, tord=trend_order,
                theta=theta, phi=phi, lnL=-data_nll(theta), cov=None, held_fixed=sorted(hold),
                kind='maximum likelihood', trend_steps=trend_steps, n_rv_excluded=n_ex, n_extra=n_extra,
                extra_test=extra, support=support, replaced_by_trend=replaced)


def planet_rows(theta, cov, fit):
    """Per companion (literature order): P, K, e, omega (star), M0 at the RV
    tref, Tp, and 1-sigma errors from 4000 draws of the inverse-Hessian Gaussian
    (`cov`, in the mean-longitude parameters of `to_lambda`)."""
    rng = np.random.default_rng(0)
    draws = None
    if cov is not None:
        draws = np.array([from_lambda(p, fit['npl']) for p in
                          rng.multivariate_normal(to_lambda(theta, fit['npl']), cov, 4000)])
    out = []
    names = [L['name'] for L in fit['lit']] + [extra_name(fit)]*fit.get('n_extra', 0)
    for j, name in enumerate(names):
        k = cs.PARAMS_PER_PLANET*j

        def conv(p):
            lnP, lnK, s1, s2, M0 = p[..., k], p[..., k + 1], p[..., k + 2], p[..., k + 3], p[..., k + 4]
            P = np.exp(lnP)
            return dict(P=P, K=np.exp(lnK), e=np.minimum(s1**2 + s2**2, cs.E_MAX),
                        omega=np.mod(np.arctan2(s2, s1), 2*np.pi), M0=np.mod(M0, 2*np.pi),
                        Tp=fit['rv'].tref - np.mod(M0, 2*np.pi)*P/(2*np.pi))
        v = {kk: float(x) for kk, x in conv(theta).items()}
        if draws is not None:
            d = conv(draws)
            for kk in list(v):
                x = d[kk]
                if kk in ('omega', 'M0'):          # circular: spread about the value
                    x = np.angle(np.exp(1j*(x - v[kk])))
                    v[f'{kk}_err'] = float(np.std(x))
                else:
                    v[f'{kk}_err'] = float(np.std(x))
        v['name'] = name
        out.append(v)
    return out


def extra_name(fit):
    """Name of the long-period Keplerian that replaced the RV trend."""
    import re
    host = re.sub(r'^.*\[\S+ (.+)\]$', r'\1', fit['lit'][0]['name']).rsplit(' ', 1)[0]
    return f'{host} long-period companion (RV)'


# ---- 2. ground truth -------------------------------------------------------------------

def true_mass_mj(P, e, K, inc, mstar):
    """Companion mass (M_J) that gives K at inclination inc (dark companion)."""
    from scipy.optimize import brentq
    return brentq(lambda m: injection.k_from_mass(P, e, inc, m, mstar) - K, 1e-8, 1e6, xtol=1e-12)


TREND_A_MAX_AU = 100.0         # trend companion: widest separation drawn
TREND_M_MAX_MSUN = 0.3         # trend companion: heaviest mass kept (dark-companion model)
TREND_DLNL_MAX = 1.0           # trend companion: largest lnL loss vs the polynomial trend
TREND_MAX_DRAWS = 20000


def draw_trend_companion(fit, th, host, mstar, rng, name):
    """A companion that reproduces the RV trend of the drawn reference model
    `th` (decided 2026-11-16): the trend's astrometric counterpart, to test
    whether long-term motion confuses the search; its parameters are not
    meant to be recovered. Circular orbit; period log-uniform from the
    longest period the RV search fits as a companion (2 x the RV baseline)
    to a separation of TREND_A_MAX_AU; cos i uniform; Omega uniform. For a
    quadratic trend, K and the phase follow from the trend's slope and
    curvature at the RV reference time; for a linear trend the phase is
    uniform and K follows from the slope (the curvature is then limited by
    the lnL test only; fixing it at zero pins the draws to the mass cap, e.g.
    HD 111232). A draw is kept if its mass is <= TREND_M_MAX_MSUN and the
    RVs fit it nearly as well as the polynomial (data lnL lower by <
    TREND_DLNL_MAX, offsets reprofiled). Returns the companion entry (role
    'trend') or None when there is no trend."""
    import copy
    rv, n, tord = fit['rv'], fit['npl'], fit['tord']
    if tord == 0:
        return None
    k = cs.PARAMS_PER_PLANET*n
    coefs = th[k:k + tord]
    slope = coefs[0]/rv.span                                    # m/s/d at rv.tref
    curv = 2*coefs[1]/rv.span**2 if tord == 2 else 0.0          # m/s/d^2
    nll_poly = cs.neg_lnL(th, rv, n, tord)
    th0 = np.delete(th, np.s_[k:k + tord])                     # same planets and jitters, no trend
    p_lo = cs.PERIOD_CAP_BASELINES*rv.span
    p_hi = 365.25*math.sqrt(TREND_A_MAX_AU**3/mstar)
    m_cap = TREND_M_MAX_MSUN*injection.GM_SUN/injection.GM_JUP
    for i in range(TREND_MAX_DRAWS):
        P = math.exp(rng.uniform(math.log(p_lo), math.log(p_hi)))
        inc = math.acos(rng.uniform(-1, 1))
        Omega = rng.uniform(0, 2*np.pi)
        nu = 2*np.pi/P                       # v = K cos(nu (t - tref) + phase): v' = -K nu sin(phase)
        if tord == 2:                        # v'(tref) = slope, v''(tref) = -K nu^2 cos(phase) = curv
            K, phase = math.hypot(slope, curv/nu)/nu, math.atan2(-slope, -curv/nu)
        else:
            phase = rng.uniform(0, 2*np.pi)
            if math.sin(phase)*slope > 0:    # sign of the slope set by the phase
                phase = -phase
            K = abs(slope/(nu*math.sin(phase)))
        try:
            m = true_mass_mj(P, 0.0, K, inc, mstar)
        except ValueError:
            continue
        if m > m_cap:
            continue
        rv2 = copy.copy(rv)
        rv2.rv = rv.rv - K*np.cos(nu*(rv.t - rv.tref) + phase)
        dlnl = nll_poly - cs.neg_lnL(th0, rv2, n, 0)
        if dlnl < -TREND_DLNL_MAX:
            continue
        m_sun = m*injection.GM_JUP/injection.GM_SUN
        a_au = ((mstar + m_sun)*(P/365.25)**2)**(1/3)
        a_star = a_au*host['parallax']*m_sun/(mstar + m_sun)
        M0 = phase % (2*np.pi)
        return dict(
            name=name, role='trend',
            rv_based=dict(P_d=P, K_ms=K, e=0.0, omega_star_deg=0.0, M0_deg=math.degrees(M0),
                          Tp_bjd_tcb=rv.tref - M0*P/(2*np.pi)),
            random=dict(inc_deg=math.degrees(inc), Omega_deg=math.degrees(Omega)),
            derived=dict(mass_mj=m, msini_mj=true_mass_mj(P, 0.0, K, np.pi/2, mstar), a_au=a_au,
                         a_star_mas=a_star,
                         thiele_innes_ABFG_mas=list(orbit.ti_from_elements(a_star, inc, 0.0, Omega))),
            trend=dict(order=int(tord), slope_ms_per_yr=float(slope*365.25),
                       curvature_ms_per_yr2=float(curv*365.25**2), period_range_d=[p_lo, p_hi],
                       dlnL_vs_polynomial=float(dlnl), n_draws=i + 1,
                       injected='orbit minus its offset and velocity at J2016.0 (see companion_sky)'),
            tref_mjd=rv.tref - BJD_MJD)
    tc = _eccentric_trend_companion(fit, th, host, mstar, rng, name, th0, nll_poly, m_cap, p_hi)
    if tc is None:
        raise RuntimeError(f'no trend companion consistent with the RVs in {TREND_MAX_DRAWS} circular and '
                           f'{TREND_MAX_DRAWS} eccentric draws')
    return tc


TREND_E_MAX = 0.9              # fallback trend companion: e uniform in [0, TREND_E_MAX]


def _eccentric_trend_companion(fit, th, host, mstar, rng, name, th0, nll_poly, m_cap, p_hi):
    """Fallback when no circular trend companion fits (Judah, 2026-11-17: an
    eccentric companion is fine): P log-uniform from 1x the RV baseline to
    p_hi, e uniform in [0, TREND_E_MAX], omega, M0, Omega uniform, cos i
    uniform; K is the weighted least-squares match of the orbit's RV curve to
    the polynomial trend at the RV epochs. Kept by the same mass and lnL tests."""
    import copy
    rv, n, tord = fit['rv'], fit['npl'], fit['tord']
    k = cs.PARAMS_PER_PLANET*n
    poly = cs.trend_model(th[k:k + tord], rv)
    w = 1/rv.err**2
    p_lo = rv.span
    for i in range(TREND_MAX_DRAWS):
        P = math.exp(rng.uniform(math.log(p_lo), math.log(p_hi)))
        e = rng.uniform(0, TREND_E_MAX)
        om, M0 = rng.uniform(0, 2*np.pi, 2)
        inc = math.acos(rng.uniform(-1, 1))
        Omega = rng.uniform(0, 2*np.pi)
        s1, s2 = math.sqrt(e)*math.cos(om), math.sqrt(e)*math.sin(om)
        shape = cs.keplerian(rv.t, rv.tref, math.log(P), 0.0, s1, s2, M0)      # K = 1
        shape = shape - np.sum(w*shape)/np.sum(w)
        pc = poly - np.sum(w*poly)/np.sum(w)
        K = float(np.sum(w*shape*pc)/np.sum(w*shape**2))
        if K < 0:                                 # -K at omega is K at omega + pi
            K, om = -K, om + np.pi
            s1, s2 = -s1, -s2
        if not K > 0:
            continue
        try:
            m = true_mass_mj(P, e, K, inc, mstar)
        except ValueError:
            continue
        if m > m_cap:
            continue
        rv2 = copy.copy(rv)
        rv2.rv = rv.rv - cs.keplerian(rv.t, rv.tref, math.log(P), math.log(K), s1, s2, M0)
        dlnl = nll_poly - cs.neg_lnL(th0, rv2, n, 0)
        if dlnl < -TREND_DLNL_MAX:
            continue
        m_sun = m*injection.GM_JUP/injection.GM_SUN
        a_au = ((mstar + m_sun)*(P/365.25)**2)**(1/3)
        a_star = a_au*host['parallax']*m_sun/(mstar + m_sun)
        om = om % (2*np.pi)
        M0 = M0 % (2*np.pi)
        return dict(
            name=name, role='trend',
            rv_based=dict(P_d=P, K_ms=K, e=e, omega_star_deg=math.degrees(om), M0_deg=math.degrees(M0),
                          Tp_bjd_tcb=rv.tref - M0*P/(2*np.pi)),
            random=dict(inc_deg=math.degrees(inc), Omega_deg=math.degrees(Omega)),
            derived=dict(mass_mj=m, msini_mj=true_mass_mj(P, e, K, np.pi/2, mstar), a_au=a_au,
                         a_star_mas=a_star,
                         thiele_innes_ABFG_mas=list(orbit.ti_from_elements(a_star, inc, om, Omega))),
            trend=dict(kind='eccentric', order=int(tord), period_range_d=[p_lo, p_hi], e_max=TREND_E_MAX,
                       dlnL_vs_polynomial=float(dlnl), n_draws=i + 1,
                       injected='orbit minus its offset and velocity at J2016.0 (see companion_sky)'),
            tref_mjd=rv.tref - BJD_MJD)
    return None


def nea_orbits_for_undetected(fit, th):
    """Judah, 2026-11-21: an NEA companion the RVs do not detect (Delta-BIC for
    dropping it < cs.DELTA_BIC) enters the truth with its NEA orbit, not the
    fitted one, which drifts without priors (HD 191939 g: fit P 422 d, e 0.80
    vs NEA 288.6 d). P and K are the NEA values, e the NEA value or 0 where
    unreported, omega the NEA value where given; otherwise the fitted omega and
    the fitted mean anomaly at the RV tref are kept (the RVs constrain neither).
    Returns (theta, {companion index: description of the change})."""
    th = np.array(th, float)
    P = cs.PARAMS_PER_PLANET
    out = {}
    for s in fit.get('support', []):
        j = s['index']
        if s['dBIC_drop'] >= cs.DELTA_BIC or j >= len(fit['lit']):
            continue
        L, k = fit['lit'][j], P*j
        if not L.get('P0') or not L.get('K0'):
            continue
        e = L['e0'] or 0.0
        om = L['omega0'] if L.get('omega0') is not None else math.atan2(th[k + 3], th[k + 2])
        old = (math.exp(th[k]), math.exp(th[k + 1]), th[k + 2]**2 + th[k + 3]**2)
        th[k], th[k + 1] = math.log(L['P0']), math.log(L['K0'])
        th[k + 2], th[k + 3] = math.sqrt(e)*math.cos(om), math.sqrt(e)*math.sin(om)
        out[j] = (f"NEA orbit (not detected by the RVs: Delta-BIC {s['dBIC_drop']:.1f}); fitted "
                  f"P {old[0]:.5g} d, K {old[1]:.4g} m/s, e {old[2]:.3f} replaced; "
                  f"e {'NEA' if L['e0'] is not None else '0 (unreported)'}, omega "
                  f"{'NEA' if L.get('omega0') is not None else 'fitted'}, phase fitted")
    return th, out


def draw_truth(fit, host, mstar, rng):
    """One ground-truth system: the reference fit's RV-based parameters (its
    best fit, Judah 2026-11-21; a covariance saved with an older fit is
    ignored), with the NEA orbit for NEA companions the RVs do not detect
    (`nea_orbits_for_undetected`); inclination and Omega drawn at random,
    derived masses and astrometric sizes, and the fixed host values."""
    th = np.array(fit['theta'], float)
    th, from_nea = nea_orbits_for_undetected(fit, th)
    rows = planet_rows(th, None, fit)
    tref_mjd = fit['rv'].tref - BJD_MJD
    comps = []
    n_lit = len(fit['lit'])
    for j, r in enumerate(rows):
        inc = math.acos(rng.uniform(-1, 1))
        Omega = rng.uniform(0, 2*np.pi)
        m = true_mass_mj(r['P'], r['e'], r['K'], inc, mstar)
        msini = true_mass_mj(r['P'], r['e'], r['K'], np.pi/2, mstar)
        m_sun = m*injection.GM_JUP/injection.GM_SUN
        a_au = ((mstar + m_sun)*(r['P']/365.25)**2)**(1/3)
        a_star = a_au*host['parallax']*m_sun/(mstar + m_sun)
        comps.append(dict(
            name=r['name'],
            rv_based=dict(P_d=r['P'], K_ms=r['K'], e=r['e'], omega_star_deg=math.degrees(r['omega']),
                          M0_deg=math.degrees(r['M0']), Tp_bjd_tcb=r['Tp'] + 0.0),
            random=dict(inc_deg=math.degrees(inc), Omega_deg=math.degrees(Omega)),
            derived=dict(mass_mj=m, msini_mj=msini, a_au=a_au, a_star_mas=a_star,
                         thiele_innes_ABFG_mas=list(orbit.ti_from_elements(a_star, inc, r['omega'], Omega))),
            tref_mjd=tref_mjd))
        if j in from_nea:
            comps[-1]['orbit_source'] = from_nea[j]
        if j >= n_lit:                   # the long-period Keplerian that replaced the RV trend
            comps[-1].update(role='trend', trend=dict(
                kind='keplerian', order=0, test=fit.get('extra_test'),
                injected='orbit minus its offset and velocity at J2016.0 (see companion_sky)'))
        elif fit['lit'][j].get('context'):  # known, but cut from the sample (Judah, 2026-11-20):
            comps[-1].update(role='context',  # injected and matchable, not scored as a target
                             cut_reason=fit['lit'][j].get('cut_reason'))
    import re
    host_name = re.sub(r'^.*\[\S+ (.+)\]$', r'\1', rows[0]['name']).rsplit(' ', 1)[0]
    tc = draw_trend_companion(fit, th, host, mstar, rng, f'{host_name} trend companion')
    if tc is not None:
        comps.append(tc)
    k = cs.PARAMS_PER_PLANET*fit['npl']
    trend = [float(c) for c in th[k:k + fit['tord']]]
    return dict(companions=comps, rv_trend_coefs=trend,
                fixed=dict(mstar_msun=mstar, ra_deg=host['ra'], dec_deg=host['dec'],
                           parallax_mas=host['parallax'], pmra_masyr=host['pmra'],
                           pmdec_masyr=host['pmdec'], G=host['phot_g_mean_mag'],
                           ref_epoch_mjd=J2016_MJD))


# ---- 3-4. epochs and astrometry --------------------------------------------------------

def scan_epochs(ra, dec, rng):
    """DR4 transits from gaiamock's scanning law: MJD (TCB), scan angle (rad),
    along-scan parallax factor, CCD row. A random REJECT_FRAC removed."""
    from gaiamock import gaiamock as gm
    tab = gm.get_gost_one_position(ra, dec, data_release='dr4')
    tab = tab[rng.uniform(0, 1, len(tab)) > REJECT_FRAC]
    psi, fal, jd, row = gm.fetch_table_element(
        ['scanAngle[rad]', 'parallaxFactorAlongScan',
         'ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]', 'CcdRow[1-7]'], tab)
    o = np.argsort(jd)
    return (np.asarray(jd)[o] - BJD_MJD, np.asarray(psi)[o], np.asarray(fal)[o],
            np.asarray(row)[o].astype(int))


def sigma_ccd_mas(G):
    from gaiamock import gaiamock as gm
    return float(gm.al_uncertainty_per_ccd_interp(G))


def baseline_al(t, psi, fal, fixed):
    """Along-scan position (mas) of the barycentre: proper motion from the
    J2016.0 position, plus parallax."""
    dt = (t - fixed['ref_epoch_mjd'])/365.25
    return (fixed['pmra_masyr']*dt*np.sin(psi) + fixed['pmdec_masyr']*dt*np.cos(psi)
            + fixed['parallax_mas']*fal)


def companion_sky(t, comp, full=False):
    """The star's offset (dRA*, dDec) in mas due to one companion. For the
    trend companion (role 'trend'), its offset and velocity at J2016.0 are
    subtracted unless `full`: the Gaia DR3 position and proper motion used
    for the barycentre are the star's apparent ones, which already contain
    them, so what is injected is the change in its motion (acceleration and
    higher orders). A position + proper-motion fit absorbs the subtracted
    part, so the search sees the same data either way."""
    def sky(t):
        rb = comp['rv_based']
        M = orbit.mean_anomaly(t, rb['P_d'], math.radians(rb['M0_deg']), comp['tref_mjd'])
        X, Y = orbit.xy(M, rb['e'])
        A, B, F, G = comp['derived']['thiele_innes_ABFG_mas']
        return np.array([B*X + G*Y, A*X + F*Y])
    t = np.asarray(t, float)
    out = sky(t)
    if comp.get('role') == 'trend' and not full:
        h = 1.0                                             # d
        x0 = sky(np.array([J2016_MJD]))
        v0 = (sky(np.array([J2016_MJD + h])) - sky(np.array([J2016_MJD - h])))/(2*h)
        out = out - x0 - v0*(t - J2016_MJD)
    return out[0], out[1]


def grid_for(t, rv):
    """An RVData on times `t` with the reference time and span of `rv`, to
    evaluate the RV model on a grid."""
    g = cs.RVData(t, np.zeros_like(t), np.ones_like(t), np.zeros(len(t)))
    g.tref, g.span = rv.tref, rv.span
    return g


def companions_al(t, psi, comps):
    out = np.zeros_like(t)
    for c in comps:
        x, y = companion_sky(t, c)
        out += x*np.sin(psi) + y*np.cos(psi)
    return out


def write_csv(path, t, psi, fal, row, al_ccd, sig):
    """Archive layout (one row per CCD), as data_conversion/ writes it."""
    n_ccd = np.where(row == 4, N_CCD_ROW4, N_CCD)
    tid = np.repeat(np.arange(len(t)) + 1, n_ccd)
    df = pd.DataFrame(dict(transit_id=tid, epoch=np.repeat(t, n_ccd),
                           scan_pos_angle=np.repeat(np.degrees(psi), n_ccd),
                           centroid_pos_al=al_ccd, centroid_pos_error_al=sig,
                           ipd_error_al=sig, parallax_factor_al=np.repeat(fal, n_ccd),
                           used_by_agis_al=True, outlier_flag=0))
    df.to_csv(path)


def make_astrometry(truth, rng):
    """Epochs, the per-CCD measurements with and without the companions (same
    noise), and the noise settings."""
    fx = truth['fixed']
    t, psi, fal, row = scan_epochs(fx['ra_deg'], fx['dec_deg'], rng)
    sig = sigma_ccd_mas(fx['G'])
    extra = rng.uniform(0, EXTRA_NOISE_MAX) if fx['G'] < 13 else 0.0
    n_ccd = np.where(row == 4, N_CCD_ROW4, N_CCD)
    noise = (np.repeat(extra*rng.standard_normal(len(t)), n_ccd)
             + sig*rng.standard_normal(n_ccd.sum()))
    base = np.repeat(baseline_al(t, psi, fal, fx), n_ccd)
    orb = np.repeat(companions_al(t, psi, truth['companions']), n_ccd)
    return dict(t=t, psi=psi, fal=fal, row=row, sig=np.full(n_ccd.sum(), sig),
                noplanet=base + noise, planets=base + orb + noise,
                noise=dict(sigma_ccd_mas=sig, extra_transit_sigma_mas=extra,
                           n_transits=int(len(t)), n_ccd_rows=int(n_ccd.sum()),
                           reject_frac=REJECT_FRAC))


def rng_for(system, seed):
    """Random generator of one benchmark dataset, seeded by the seed AND the
    system name, so different systems with the same seed get independent
    draws (from 2026-11-16; HD39091_seed1 was made with the seed alone)."""
    import zlib
    return np.random.default_rng([seed, zlib.crc32(system.encode())])


_FIT_CACHE = {}


def make_system(system, seed, fig=True, force=False):
    stem = os.path.join(OUT_DIR, f'{system}_seed{seed}')
    if os.path.exists(f'{stem}_truth.json') and not force:
        raise FileExistsError(f'{stem}_truth.json exists: benchmark inputs are frozen (force=True to replace)')
    rng = rng_for(system, seed)
    if system not in _FIT_CACHE:                       # one reference fit per system and process
        _FIT_CACHE[system] = reference_rv_fit(system)
    fit = _FIT_CACHE[system]
    host = host_astrometry(fit['source_id'])
    cat = pd.read_csv(cs.CATALOG, dtype={'system_id': str}, low_memory=False)
    mstar = float(cat[cat.system_id == fit['source_id']].st_mass.iloc[0])
    truth = draw_truth(fit, host, mstar, rng)
    ast = make_astrometry(truth, rng)
    os.makedirs(OUT_DIR, exist_ok=True)
    write_csv(f'{stem}.csv', ast['t'], ast['psi'], ast['fal'], ast['row'], ast['planets'], ast['sig'])
    write_csv(f'{stem}_noplanet.csv', ast['t'], ast['psi'], ast['fal'], ast['row'], ast['noplanet'],
              ast['sig'])
    rows = planet_rows(fit['theta'], fit['cov'], fit)
    _, coefs, jit = cs.unpack(fit['theta'], fit['npl'], fit['tord'], fit['rv'].nds)
    truth.update(system=system, seed=seed, noise=ast['noise'],
                 reference_fit=dict(planets=rows, lnL=fit['lnL'], trend_order=fit['tord'],
                                    trend_coefs=[float(c) for c in coefs],
                                    trend_steps=fit['trend_steps'],
                                    n_rv_excluded=fit['n_rv_excluded'], rv_binned_nightly=True,
                                    n_rv_points=int(fit['rv'].n),
                                    jitter_ms=dict(zip(fit['rv'].names, map(float, jit))),
                                    rv_tref_bjd_tcb=fit['rv'].tref, rv_span_d=fit['rv'].span,
                                    hessian_held_fixed=[int(i) for i in fit['held_fixed']],   # since 2026-11-21: held in the polish
                                    n_extra=fit['n_extra'], long_period_keplerian_test=fit['extra_test'],
                                    rv_support=fit['support'],
                                    theta=[float(x) for x in fit['theta']],
                                    kind=fit['kind'],
                                    method=(f"{fit['kind']}; errors from the inverse Hessian" if fit.get('cov') is not None
                                            else f"{fit['kind']}; no errors (from the Octofitter fits later)")))
    json.dump(truth, open(f'{stem}_truth.json', 'w'), indent=1)
    print(f'wrote {stem}.csv, _noplanet.csv, _truth.json')
    if fig:
        process_figure(fit, truth, f'{stem}.csv', f'{stem}_noplanet.csv',
                       os.path.join(FIG_DIR, f'synthetic_{system}_seed{seed}.png'))
    return fit, truth


def fit_from_truth(truth):
    """The reference fit as saved in a _truth.json (theta, no covariance),
    with its RV data: enough to redraw the figure without refitting."""
    ref = truth['reference_fit']
    path = rv_table(truth['system'])
    rv = load_rvs(path, truth['system'], nightly=ref.get('rv_binned_nightly', False))[0]   # as fitted
    lit = literature_sorted(os.path.basename(path).split('_')[0])
    return dict(rv=rv, lit=lit, npl=len(lit) + ref.get('n_extra', 0), n_extra=ref.get('n_extra', 0),
                tord=ref['trend_order'],
                theta=np.array(ref['theta']), held_fixed=ref.get('hessian_held_fixed', ref.get('laplace_held_fixed')))   # older files: laplace_*


def redraw(system, seed):
    stem = os.path.join(OUT_DIR, f'{system}_seed{seed}')
    truth = json.load(open(f'{stem}_truth.json'))
    process_figure(fit_from_truth(truth), truth, f'{stem}.csv', f'{stem}_noplanet.csv',
                   os.path.join(FIG_DIR, f'synthetic_{system}_seed{seed}.png'))


# ---- figure ----------------------------------------------------------------------------

def _pretty(name):
    """Figure label of a companion: 'HD 89707 [sahlmann_2011 HD 89707 b]' (a
    catalog name with its source tag) -> 'HD 89707 b'; pi Men in Greek."""
    import re
    name = re.sub(r'^.*\[\S+ (.+)\]$', r'\1', name)
    return name.replace('HD 39091', r'$\pi$ Men').replace('pi Men', r'$\pi$ Men')


def process_figure(fit, truth, csv_planets, csv_noplanet, out):
    """Four rows: the RVs with the reference fit; the ground-truth values (RV-
    based vs random vs derived vs fixed); the zero-planet model with the DR4
    epochs and errors; the synthetic data with the companions injected."""
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    rv, th, n, tord = fit['rv'], fit['theta'], fit['npl'], fit['tord']
    rows = truth['reference_fit']['planets']
    legs = []                      # (axes, legend kwargs), placed by figtools at the end
    fx, comps = truth['fixed'], truth['companions']
    name = _pretty(comps[0]['name']).rsplit(' ', 1)[0]      # host name

    fig = plt.figure(figsize=(20, 44))
    gs = fig.add_gridspec(8, 2, height_ratios=[1.25, 0.5, 0.95, 1.0, 0.45, 1.0, 0.3, 0.85], hspace=0.42,
                          wspace=0.22)                     # row 6: room for the legend of 4b

    # Row 1: RVs and the reference fit, and the residuals (below)
    sub1 = gs[0, :].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.2)
    ax = fig.add_subplot(sub1[0])
    nll, off, _ = cs.neg_lnL(th, rv, n, tord, True)
    jit = cs.unpack(th, n, tord, rv.nds)[2]
    yr = lambda t: 2000 + (t - 2451544.5)/365.25   # noqa: E731
    for d, nm in enumerate(rv.names):
        m = rv.ds == d
        lab = nm.split(':')[1] if ':' in nm else nm
        k = cs.PARAMS_PER_PLANET*n + tord + d
        jl = (f'jitter {jit[d]:.1f} m/s (held fixed: errors {np.median(rv.err[m]):.0f} m/s)'
              if k in fit['held_fixed'] else f'jitter {jit[d]:.1f} m/s')
        ax.errorbar(yr(rv.t[m]), rv.rv[m] - off[d], np.hypot(rv.err[m], jit[d]), fmt='o', ms=4,
                    color=f'C{d % 10}', alpha=0.7, lw=1, label=f'{lab}: {m.sum()} RVs, {jl}')
    tt = np.linspace(rv.t.min() - 200, rv.t.max() + 200, 30000)
    grid = cs.RVData(tt, np.zeros_like(tt), np.ones_like(tt), np.zeros(len(tt)))
    grid.tref, grid.span = rv.tref, rv.span
    kind = truth['reference_fit'].get('kind', 'maximum likelihood')
    ax.plot(yr(tt), cs.model_rv(th, grid, n, tord), 'k-', lw=1.2, zorder=5,
            label=f'reference fit ({kind}):')
    for r in rows:
        ax.plot([], [], ' ', label=(f"{_pretty(r['name'])}: P = {r['P']:.4g} $\\pm$ {r['P_err']:.2g} d, "
                                    f"K = {r['K']:.3g} $\\pm$ {r['K_err']:.2g} m/s, "
                                    f"e = {r['e']:.3f} $\\pm$ {r['e_err']:.3f}, "
                                    f"$\\omega_\\star$ = {math.degrees(r['omega']):.1f} $\\pm$ "
                                    f"{math.degrees(r['omega_err']):.2g}$^\\circ$, "
                                    f"T$_p$ = BJD {r['Tp']:.0f} $\\pm$ {r['Tp_err']:.2g}"))
    tcs = [c for c in comps if c.get('role') == 'trend']
    if tord:
        cf = truth['reference_fit']['trend_coefs']
        lab = (f"trend ({'linear' if tord == 1 else 'quadratic'}; chosen step by step, "
               f"$\\Delta$BIC $\\geq$ {cs.DELTA_BIC:g}): {cf[0]/rv.span*365.25:.3g} m/s/yr")
        if tord == 2:
            lab += f", curvature {2*cf[1]/rv.span**2*365.25**2:.3g} m/s/yr$^2$"
        tt2 = np.linspace(rv.t.min(), rv.t.max(), 3000)
        ax.plot(yr(tt2), cs.trend_model(np.array(cf), grid_for(tt2, rv)), '--', color='0.4', lw=2,
                zorder=4, label=lab + ' (dashed)')
    else:
        ax.plot([], [], ' ', label=f'no trend (raising the order gains $\\Delta$BIC < {cs.DELTA_BIC:g})')
    for c in tcs:
        if c['trend'].get('kind') == 'keplerian':          # part of the reference fit (black curve)
            continue
        rb = c['rv_based']
        tt2 = np.linspace(rv.t.min(), rv.t.max(), 3000)
        w, e = math.radians(rb['omega_star_deg']), rb['e']
        kep = lambda t: cs.keplerian(t, rv.tref, math.log(rb['P_d']), math.log(rb['K_ms']),   # noqa: E731
                                     math.sqrt(e)*math.cos(w), math.sqrt(e)*math.sin(w),
                                     math.radians(rb['M0_deg']))
        ax.plot(yr(tt2), kep(tt2) - kep(np.array([rv.tref]))[0], ':', color='C3', lw=2.5,
                zorder=4, label="trend companion's RV, drawn to match the trend (dotted; injected "
                                "into the astrometry)")
    ax.set_ylabel('RV $-$ offset (m/s)')
    ax.tick_params(labelbottom=False)
    ax.set_title(f'1. {name}: RVs and the reference fit (offsets and jitters per dataset; errors include jitter)',
                 fontsize=22)
    hl = ax.get_legend_handles_labels()
    xl = ax.get_xlim()
    ax = fig.add_subplot(sub1[1], sharex=ax)
    resid = rv.rv - off[rv.ds] - cs.model_rv(th, rv, n, tord)
    for d in range(rv.nds):
        m = rv.ds == d
        ax.errorbar(yr(rv.t[m]), resid[m], np.hypot(rv.err[m], jit[d]), fmt='o', ms=4, color=f'C{d % 10}',
                    alpha=0.7, lw=1)
    ax.axhline(0, color='0.4', lw=1, gid='guide')
    ax.set_xlim(xl)
    ax.set_xlabel('year')
    ax.set_ylabel('residual (m/s)')
    ax.set_title('Residuals: data minus the reference fit (colours as above; rms '
                 + ', '.join(f"{(rv.names[d].split(':')[1] if ':' in rv.names[d] else rv.names[d])} "
                             f"{np.std(resid[rv.ds == d]):.3g}" for d in range(rv.nds)) + ' m/s)', fontsize=15)
    ax = fig.add_subplot(gs[1, :])
    ax.axis('off')
    ax.legend(*hl, fontsize=13, ncol=2, loc='upper left', frameon=False)

    # Row 2: ground truth list
    ax = fig.add_subplot(gs[2, :])
    ax.axis('off')
    ax.set_title('2. Ground truth for the synthetic DR4 data (one draw)', fontsize=22, loc='left')
    targets = [c for c in comps if c.get('role') != 'trend']
    lines = [('RV-based: the reference fit' + (' (one draw from its inverse-Hessian Gaussian)'
              if 'inverse Hessian' in truth['reference_fit'].get('method', '') else ' (best fit)'), 'C0'),]
    for c in targets:
        rb = c['rv_based']
        lines.append((f"   {_pretty(c['name'])}:  P = {rb['P_d']:.5g} d,  K = {rb['K_ms']:.4g} m/s,  "
                      f"e = {rb['e']:.3f},  $\\omega_\\star$ = {rb['omega_star_deg']:.1f}$^\\circ$,  "
                      f"T$_p$ = BJD {rb['Tp_bjd_tcb']:.1f}", 'C0'))
    lines.append(('Random: isotropic inclination (cos i uniform), $\\Omega$ uniform in [0, 360)$^\\circ$', 'C3'))
    for c in targets:
        rd = c['random']
        lines.append((f"   {_pretty(c['name'])}:  i = {rd['inc_deg']:.1f}$^\\circ$,  "
                      f"$\\Omega$ = {rd['Omega_deg']:.1f}$^\\circ$", 'C3'))
    lines.append(('Derived (from the above and the host mass)', '0.3'))
    for c in targets:
        dv = c['derived']
        lines.append((f"   {_pretty(c['name'])}:  M = {dv['mass_mj']:.4g} M$_J$ (M sin i = "
                      f"{dv['msini_mj']:.4g}),  a = {dv['a_au']:.3f} AU,  star's astrometric "
                      f"semi-major axis = {dv['a_star_mas']*1e3:.4g} $\\mu$as", '0.3'))
    for c in tcs:
        rb, rd, dv, tr = c['rv_based'], c['random'], c['derived'], c['trend']
        kind = tr.get('kind', 'circular')
        if kind == 'keplerian':
            t = tr['test']
            lines.append(('Long-period companion (replaces the RV trend; injected, not a target): one blind '
                          f"Keplerian fitted to the RVs, BIC {t['bic_keplerian']:.1f} vs {t['bic_trend']:.1f} for "
                          f"the order-{t['trend_order_tested']} trend", 'C2'))
            lines.append(("   P, K, e, omega, M0 from the reference fit; cos i, $\\Omega$ uniform; its J2016.0 offset "
                          "and velocity are in the DR3 position and proper motion, so only its acceleration is "
                          "injected", 'C2'))
        else:
            orb = ('circular, P log-uniform from 2x the RV baseline' if kind == 'circular' else
                   f"eccentric (no circular orbit fit), e uniform in [0, {tr['e_max']:g}], P log-uniform "
                   "from 1x the RV baseline")
            lines.append(('Trend companion (reproduces the RV trend; injected, not a target): '
                          f'{orb} to {TREND_A_MAX_AU:g} AU', 'C2'))
            lines.append((f"   cos i, $\\Omega${', phase' if tr['order'] == 1 or kind != 'circular' else ''} uniform; K"
                          f"{'' if tr['order'] == 1 or kind != 'circular' else ' and phase'} set by the trend; kept if "
                          f"M $\\leq$ {TREND_M_MAX_MSUN:g} M$_\\odot$ and the RVs fit it as well as the polynomial "
                          f"($\\Delta$lnL > $-${TREND_DLNL_MAX:g}); its J2016.0 offset and velocity are in the "
                          f"DR3 position and proper motion, so only its acceleration is injected", 'C2'))
        tail = (f";  $\\Delta$lnL = {tr['dlnL_vs_polynomial']:.2f} (draw {tr['n_draws']})"
                if kind != 'keplerian' else '')
        lines.append((f"   P = {rb['P_d']/365.25:.4g} yr,  K = {rb['K_ms']:.4g} m/s,  e = {rb['e']:.2f},  "
                      f"i = {rd['inc_deg']:.1f}$^\\circ$,  $\\Omega$ = {rd['Omega_deg']:.1f}$^\\circ$,  "
                      f"M = {dv['mass_mj']:.4g} M$_J$,  a = {dv['a_au']:.3g} AU,  star's astrometric "
                      f"semi-major axis = {dv['a_star_mas']:.4g} mas{tail}", 'C2'))
    nz = truth['noise']
    lines.append((f"Fixed: M$_\\star$ = {fx['mstar_msun']} M$_\\odot$ (catalog);  Gaia DR3 (J2016.0): "
                  f"$\\varpi$ = {fx['parallax_mas']:.3f} mas, $\\mu_{{\\alpha*}}$ = {fx['pmra_masyr']:.2f}, "
                  f"$\\mu_\\delta$ = {fx['pmdec_masyr']:.2f} mas/yr, G = {fx['G']:.2f}", 'k'))
    lines.append((f"Noise: $\\sigma_{{CCD}}$ = {nz['sigma_ccd_mas']:.3f} mas (gaiamock), "
                  f"extra per-transit noise {nz['extra_transit_sigma_mas']*1e3:.1f} $\\mu$as (not in the "
                  f"formal errors);  {nz['n_transits']} transits after removing "
                  f"{nz['reject_frac']:.0%}", 'k'))
    for i, (s, col) in enumerate(lines):
        ax.text(0.0, 1 - i/len(lines), s, color=col, fontsize=16, transform=ax.transAxes, va='top',
                fontweight='bold' if not s.startswith('   ') else 'normal')

    # Rows 3-5: astrometry (per transit, as the search sees it)
    g0 = bsdata.load_gaia_dr4(csv_noplanet)
    g1 = bsdata.load_gaia_dr4(csv_planets)
    psi = np.arctan2(g0.sinpsi, g0.cospsi)
    base = baseline_al(g0.t, psi, g0.fal, fx)
    t_yr = lambda t: 2000 + (t - 51544.5)/365.25   # noqa: E731
    tg = np.linspace(g0.t.min() - 30, g0.t.max() + 30, 4000)
    dtg = (tg - fx['ref_epoch_mjd'])/365.25
    # the parallax ellipse needs both parallax factors: the sky track is drawn by
    # interpolating the transits' (dRA*, dDec) parallax displacement, so only the
    # proper motion part is drawn continuously
    ax = fig.add_subplot(gs[3, 0])
    ax.plot(fx['pmra_masyr']*dtg, fx['pmdec_masyr']*dtg, '-', color='0.6', lw=1.5,
            label='proper motion alone')
    dt0 = (g0.t - fx['ref_epoch_mjd'])/365.25
    # position at each transit: proper motion + along-scan parallax shift (the
    # only component Gaia measures), drawn along the scan direction
    x0, y0 = fx['pmra_masyr']*dt0, fx['pmdec_masyr']*dt0
    s = fx['parallax_mas']*g0.fal
    ax.plot(x0 + s*g0.sinpsi, y0 + s*g0.cospsi, 'o', color='C0', ms=6,
            label=f'DR4 transits ({len(g0)}): proper motion + along-scan parallax shift')
    ax.set_xlabel(r'$\Delta\alpha*$ from J2016.0 (mas)')
    ax.set_ylabel(r'$\Delta\delta$ from J2016.0 (mas)')
    ax.set_title('3a. Zero-planet model: sky motion over DR4', fontsize=20)
    legs.append((ax, dict(fontsize=13)))

    ax = fig.add_subplot(gs[3, 1])
    pm_part = base - fx['parallax_mas']*g0.fal
    ax.errorbar(t_yr(g0.t), g0.y - pm_part, g0.err, fmt='o', color='C0', ms=5, lw=1,
                label=f'zero-planet synthetic data, minus proper motion (median error '
                      f'{np.median(g0.err)*1e3:.0f} $\\mu$as)')
    ax.plot(t_yr(g0.t), fx['parallax_mas']*g0.fal, 'x', color='k', ms=7,
            label=r'model: $\varpi$ $\times$ along-scan parallax factor')
    ax.set_xlabel('year')
    ax.set_ylabel('along-scan position (mas)')
    ax.set_title('3b. Zero-planet model: along-scan parallax signal', fontsize=20)
    legs.append((ax, dict(fontsize=13)))

    ax = fig.add_subplot(gs[4, :])
    ax.errorbar(t_yr(g0.t), (g0.y - base)*1e3, g0.err*1e3, fmt='o', color='C0', ms=5, lw=1,
                label='zero-planet synthetic data minus the true zero-planet model (noise only)')
    ax.axhline(0, color='k', lw=1)
    ax.set_xlabel('year')
    ax.set_ylabel(r'residual ($\mu$as)')
    ax.set_title('3c. Zero-planet synthetic data: realistic DR4 epochs and errors', fontsize=20)
    legs.append((ax, dict(fontsize=13)))

    ax = fig.add_subplot(gs[5, 0])
    sig_al = companions_al(g1.t, psi, comps)
    ax.errorbar(t_yr(g1.t), g1.y - base, g1.err, fmt='o', color='C1', ms=5, lw=1,
                label='synthetic data with companions, minus the zero-planet model')
    ax.plot(t_yr(g1.t), sig_al, 'x', color='k', ms=7, label='injected companion signal (along scan)')
    ax.axhline(0, color='0.5', lw=1)
    ax.set_xlabel('year')
    ax.set_ylabel('along-scan offset (mas)')
    ax.set_title('4a. Synthetic data with the companions injected', fontsize=20)
    legs.append((ax, dict(fontsize=13)))

    ax = fig.add_subplot(gs[5, 1])
    # Each companion's offset over the DR4 window, and their sum. Planets: their
    # reflex orbits about the barycentre as they are (an orbit adds no net
    # motion, and removing a fitted line from one whose period is close to the
    # window shears it open: pi Men b). Trend companion: minus its linear
    # motion (best-fit position and proper motion over the DR4 transits, least
    # squares), which cannot be told apart from the star's proper motion
    # (Judah, 2026-11-16: the linear term is proper motion and is subtracted;
    # the quadratic and cubic terms are astrophysical and kept). The injected
    # signal holds no parallax.
    A = np.column_stack([np.ones_like(g1.t), g1.t - g1.t.mean()])

    def minus_line(c):
        xs, ys = companion_sky(g1.t, c, full=True)
        x, y = companion_sky(tg, c, full=True)
        if c.get('role') != 'trend':
            return x, y, xs, ys
        lin = [np.linalg.lstsq(A, v, rcond=None)[0] for v in (xs, ys)]
        at = lambda t, k: lin[k][0] + lin[k][1]*(t - g1.t.mean())     # noqa: E731
        return x - at(tg, 0), y - at(tg, 1), xs - at(g1.t, 0), ys - at(g1.t, 1)
    parts = [minus_line(c) for c in comps]
    for k, (c, (x, y, _, _)) in enumerate(zip(comps, parts)):
        span = max(np.ptp(x), np.ptp(y))
        size = f'{span:.3g} mas' if span >= 0.1 else f'{span*1e3:.3g} $\\mu$as'
        if c.get('role') == 'trend':
            lab = f"trend companion ({c['rv_based']['P_d']/365.25:.0f} yr; {size})"
            ax.plot(x, y, '--', color='C2', lw=2.5, zorder=3, label=lab)
        else:
            lab = f"{_pretty(c['name'])} ({c['rv_based']['P_d']:.4g} d; {size})"
            ax.plot(x, y, '-', color=f'C{[0, 4, 5, 6][k % 4]}', lw=1.5, zorder=3, label=lab)
    x, y = sum(p[0] for p in parts), sum(p[1] for p in parts)
    xs, ys = sum(p[2] for p in parts), sum(p[3] for p in parts)
    ax.plot(x, y, '-', color='0.55', lw=5, zorder=1, label='sum')
    ax.plot(xs, ys, 'o', color='C1', ms=5, zorder=4, label='sum at the DR4 transits')
    ax.plot(0, 0, '+', color='k', ms=14, mew=2, zorder=5, label='barycentre')
    ax.set_aspect('equal', 'datalim')
    ax.set_xlabel(r'$\Delta\alpha*$ (mas)')
    ax.set_ylabel(r'$\Delta\delta$ (mas)')
    ax.set_title("4b. The star's motion due to its companions over DR4 (trend companion: minus its\n"
                 "linear motion, i.e. proper motion; each companion alone: period, largest extent)", fontsize=16)
    legs.append((ax, dict(fontsize=11)))

    # 4c-d: the same curves as 4b, each coordinate against time. Gaia measures
    # only the along-scan projection (4a), so there are no RA or Dec data
    # points: the dots are the model at the transit epochs.
    for col, (q, nm, sym) in enumerate(((0, 'RA', r'$\Delta\alpha*$'), (1, 'Dec', r'$\Delta\delta$'))):
        ax = fig.add_subplot(gs[7, col])
        for k, (c, pt) in enumerate(zip(comps, parts)):
            if c.get('role') == 'trend':
                ax.plot(t_yr(tg), pt[q], '--', color='C2', lw=2.5, zorder=3, label='trend companion alone')
            else:
                ax.plot(t_yr(tg), pt[q], '-', color=f'C{[0, 4, 5, 6][k % 4]}', lw=1.5, zorder=3,
                        label=f"{_pretty(c['name'])} alone")
        ax.plot(t_yr(tg), (x, y)[q], '-', color='0.55', lw=5, zorder=1, label='sum')
        ax.plot(t_yr(g1.t), (xs, ys)[q], 'o', color='C1', ms=5, zorder=4, label='sum at the DR4 transits')
        ax.axhline(0, color='0.5', lw=1)
        ax.set_xlabel('year')
        ax.set_ylabel(f'{sym} (mas)')
        ax.set_title(f'4{"cd"[col]}. {nm} offset vs time (as 4b; model only: Gaia measures along scan)',
                     fontsize=16)
        legs.append((ax, dict(fontsize=12)))

    for a_, kw in legs:            # legends placed where they cover no data
        figtools.legend_clear(a_, **kw)
    figtools.check_legends(fig)
    fig.suptitle(f'Synthetic Gaia DR4 astrometry of {name} (seed {truth["seed"]})', fontsize=28, y=0.91)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out}')


# ---- test ------------------------------------------------------------------------------

def selftest():
    """A Gaia-4 b-like companion at pi Men's sky position (no RV fit): the CSV
    loads, and a least-squares fit at the true (P, e, M0) recovers the
    injected Thiele-Innes constants and parallax within 5 sigma."""
    import tempfile
    rng = np.random.default_rng(3)
    fx = dict(mstar_msun=1.0, ra_deg=84.2995, dec_deg=-80.4645, parallax_mas=20.0, pmra_masyr=100.0,
              pmdec_masyr=-50.0, G=9.0, ref_epoch_mjd=J2016_MJD)
    ti = orbit.ti_from_elements(0.5, 1.0, 0.3, 2.0)
    comp = dict(name='test b', rv_based=dict(P_d=571.3, e=0.3, M0_deg=40.0), tref_mjd=57800.0,
                derived=dict(thiele_innes_ABFG_mas=list(ti)))
    truth = dict(fixed=fx, companions=[comp])
    ast = make_astrometry(truth, rng)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 't.csv')
        write_csv(p, ast['t'], ast['psi'], ast['fal'], ast['row'], ast['planets'], ast['sig'])
        g = bsdata.load_gaia_dr4(p)
    assert len(g) == ast['noise']['n_transits'], (len(g), ast['noise'])
    dt = (g.t - fx['ref_epoch_mjd'])/365.25
    M = orbit.mean_anomaly(g.t, 571.3, math.radians(40.0), 57800.0)
    X = np.column_stack([g.sinpsi, dt*g.sinpsi, g.cospsi, dt*g.cospsi, g.fal,
                         orbit.astro_columns(M, 0.3, g.sinpsi, g.cospsi)])
    w = 1/g.err
    coef, *_ = np.linalg.lstsq(X*w[:, None], g.y*w, rcond=None)
    cov = np.linalg.inv((X*w[:, None]).T @ (X*w[:, None]))
    want = np.array([0, 100.0, 0, -50.0, 20.0, *ti])
    z = (coef - want)/np.sqrt(np.diag(cov))
    assert np.all(np.abs(z) < 5), z
    print(f'synthetic_dr4.py test passed: {len(g)} transits, max |z| = {np.max(np.abs(z)):.1f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--system', default='HD39091')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--test', action='store_true')
    ap.add_argument('--figure-only', action='store_true', help='redraw from the saved files')
    a = ap.parse_args()
    if a.test:
        selftest()
    elif a.figure_only:
        redraw(a.system, a.seed)
    else:
        make_system(a.system, a.seed)
