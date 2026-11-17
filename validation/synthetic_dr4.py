"""Synthetic Gaia DR4 epoch astrometry of a known RV system, written in the
archive CSV layout that `data.load_gaia_dr4` reads.

Steps (`make_system`):
  1. Reference RV fit (`reference_rv_fit`): the system's literature companions
     fitted to its RVs (rv_data/systems/) by maximum likelihood (or, where
     that has no proper maximum, with the literature priors), with a
     polynomial trend, per-dataset offsets (profiled) and jitters, using the
     likelihood of fits_rv/companion_search.py. The trend order (none by
     default, linear or quadratic) is chosen step by step (`select_trend`):
     tested before any companion is in the model, and re-checked after each
     companion is added (a term that stops paying Delta-BIC >= 20 is
     dropped). The literature-prior fit seeds the final fit, which has no
     priors. Its uncertainties come from the Laplace approximation: a
     Gaussian with covariance = inverse Hessian of -lnL at the best fit, in
     parameters (ln P, ln K, sqrt(e) cos w, sqrt(e) sin w, mean longitude
     M0 + w) per companion, trend coefficients and ln jitters.
  2. Ground truth (`draw_truth`): ONE draw from that Gaussian gives the
     RV-based parameters (P, K, e, omega, M0); the unknown ones are drawn at
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
import companion_search as cs   # noqa: E402
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


def _hessian_in_basis(f, x, B):
    """Second derivatives of f at x with respect to c, where x = x0 + B c (the
    columns of B are the step vectors), by central differences with unit
    steps in c."""
    m, f0 = B.shape[1], f(x)
    H = np.zeros((m, m))
    for a in range(m):
        ba = B[:, a]
        H[a, a] = f(x + ba) + f(x - ba) - 2*f0
        for b in range(a):
            bb = B[:, b]
            H[a, b] = H[b, a] = (f(x + ba + bb) - f(x + ba - bb) - f(x - ba + bb) + f(x - ba - bb))/4
    return H


def _scale_direction(f, x, d, f0):
    """Step length along direction d at which f rises by about 0.5 (None if
    f stays flat up to a length of 1e4)."""
    s = 1.0
    for _ in range(60):
        r = 0.5*(f(x + s*d) + f(x - s*d)) - f0
        if r < 0.05 and s < 1e4:
            s *= 3
        elif r > 2:
            s /= 3
        else:
            return s if r >= 0.05 else None
    return s


def laplace_cov(f, x, hold=()):
    """Inverse Hessian of f at x by central differences. First pass: steps
    along each parameter axis (`axis_steps`). Strongly correlated parameters
    (e.g. K, e, omega and phase of a highly eccentric orbit) make that matrix
    ill-conditioned, so a second pass measures the curvature along the
    first-pass eigen-directions, each with its own step (rise of about 0.5).
    A flat parameter (e.g. the jitter of a dataset whose errors dwarf it) and
    any parameter in `hold` are held fixed: zero variance. Raises if a
    direction stays flat or the Hessian is not positive definite. Returns
    (cov, indices held fixed)."""
    n, f0 = len(x), f(x)
    h, flat = axis_steps(f, x, hold)
    keep = [i for i in range(n) if i not in flat]
    B = np.zeros((n, len(keep)))
    for a, i in enumerate(keep):
        B[i, a] = h[i]
    w, V = np.linalg.eigh(_hessian_in_basis(f, x, B))
    D = B @ V                                          # eigen-directions in parameter space
    for a in range(D.shape[1]):
        s = (1/np.sqrt(w[a])) if w[a] > 1e-6 else _scale_direction(f, x, D[:, a], f0)
        if s is None:
            raise RuntimeError('Hessian not positive definite: a direction is flat')
        D[:, a] *= s
    Hc = _hessian_in_basis(f, x, D)
    if np.min(np.linalg.eigvalsh(Hc)) <= 0:
        raise RuntimeError('Hessian not positive definite: not at a maximum')
    cov = np.zeros((n, n))
    cov[np.ix_(keep, keep)] = (D @ np.linalg.inv(Hc) @ D.T)[np.ix_(keep, keep)]
    return cov, flat


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
    """The system's literature companions (companion_search priors), largest K
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
    >= DELTA_BIC. Fits are companion_search's (literature priors on the
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


def reference_rv_fit(system):
    """Fit of the literature companions (+ a trend chosen by `select_trend`) to
    the system's RVs, and its Laplace covariance. Maximum likelihood (no priors)
    by default. If that has no proper maximum (Hessian not positive definite,
    e.g. HD 89707 b: with periastron barely sampled, e and K trade off along a
    ridge), the literature-prior fit is used instead: the maximum of the
    likelihood times companion_search's Gaussian priors on ln P, ln K and e
    around the published values (they enter only this RV reference model).
    Jitters stay within their bounds; a jitter that has no effect (below 10%
    of its dataset's median error) or sits at a bound is held fixed.
    Returns a dict with the data, the literature priors (names, order),
    theta, data lnL, cov, held_fixed and `kind` (which fit was used)."""
    path = rv_table(system)
    sid = os.path.basename(path).split('_')[0]
    rv = cs.load_system(path)
    lit = literature_sorted(sid)
    th0, trend_order, trend_steps = select_trend(rv, lit)
    n = len(lit)
    from scipy.optimize import minimize
    bnds = cs.bounds_for(n, trend_order, rv.nds, [cs.P_MIN_KNOWN_D]*n, rv.pmax,
                         [max(rv.pmax, 2*L['P0']) for L in lit])
    k = cs.PARAMS_PER_PLANET*n + trend_order
    lo, hi = np.log(cs.JITTER_BOUNDS)
    data_nll = lambda x: cs.neg_lnL(x, rv, n, trend_order)      # noqa: E731
    for kind, f in (('maximum likelihood', data_nll),
                    ('literature-prior fit (maximum likelihood had no proper maximum)',
                     lambda x: data_nll(x) + cs.prior_penalty(x, lit))):
        th = th0.copy()
        for _ in range(3):
            th = minimize(f, th, method='L-BFGS-B', bounds=bnds,
                          options={'maxiter': 2000, 'maxfun': 50000}).x
        jit = np.exp(th[k:])
        hold = [k + d for d in range(rv.nds)
                if jit[d] < 0.1*np.median(rv.err[rv.ds == d]) or not lo + 0.01 < th[k + d] < hi - 0.01]

        def fl(p, f=f):
            x = from_lambda(p, n)
            x[k:] = np.clip(x[k:], lo, hi)                 # jitters stay within bounds
            return f(x)
        phi = polish_scaled(fl, to_lambda(th, n), hold)
        phi[k:] = np.clip(phi[k:], lo, hi)
        hold = sorted(set(hold) | {k + d for d in range(rv.nds) if not lo + 0.01 < phi[k + d] < hi - 0.01})
        try:
            cov, flat = laplace_cov(fl, phi, hold)
            break
        except RuntimeError:
            if f is not data_nll:
                raise
            print(f'{system}: maximum-likelihood Hessian not positive definite; using the literature-prior fit')
    theta = from_lambda(phi, n)
    return dict(system=system, source_id=sid, rv=rv, lit=lit, npl=n, tord=trend_order,
                theta=theta, phi=phi, lnL=-data_nll(theta), cov=cov, held_fixed=flat, kind=kind,
                trend_steps=trend_steps)


def planet_rows(theta, cov, fit):
    """Per companion (literature order): P, K, e, omega (star), M0 at the RV
    tref, Tp, and 1-sigma errors from 4000 draws of the Laplace Gaussian
    (`cov`, in the mean-longitude parameters of `to_lambda`)."""
    rng = np.random.default_rng(0)
    draws = None
    if cov is not None:
        draws = np.array([from_lambda(p, fit['npl']) for p in
                          rng.multivariate_normal(to_lambda(theta, fit['npl']), cov, 4000)])
    out = []
    for j, L in enumerate(fit['lit']):
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
        v['name'] = L['name']
        out.append(v)
    return out


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
    raise RuntimeError(f'no trend companion consistent with the RVs in {TREND_MAX_DRAWS} draws')


def draw_truth(fit, host, mstar, rng):
    """One ground-truth system: RV-based parameters drawn from the reference
    fit's Laplace Gaussian, inclination and Omega drawn at random, derived
    masses and astrometric sizes, and the fixed host values."""
    while True:
        th = from_lambda(rng.multivariate_normal(fit['phi'], fit['cov']), fit['npl'])
        pl = cs.unpack(th, fit['npl'], fit['tord'], fit['rv'].nds)[0]
        if all(p[2]**2 + p[3]**2 < cs.E_MAX for p in pl):
            break
    rows = planet_rows(th, None, fit)
    tref_mjd = fit['rv'].tref - BJD_MJD
    comps = []
    for r in rows:
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


def make_system(system, seed, fig=True, force=False):
    stem = os.path.join(OUT_DIR, f'{system}_seed{seed}')
    if os.path.exists(f'{stem}_truth.json') and not force:
        raise FileExistsError(f'{stem}_truth.json exists: benchmark inputs are frozen (force=True to replace)')
    rng = rng_for(system, seed)
    fit = reference_rv_fit(system)
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
                                    jitter_ms=dict(zip(fit['rv'].names, map(float, jit))),
                                    rv_tref_bjd_tcb=fit['rv'].tref, rv_span_d=fit['rv'].span,
                                    laplace_held_fixed=[int(i) for i in fit['held_fixed']],
                                    theta=[float(x) for x in fit['theta']],
                                    kind=fit['kind'],
                                    method=f"{fit['kind']}; errors from the Laplace approximation"))
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
    rv = cs.load_system(path)
    lit = literature_sorted(os.path.basename(path).split('_')[0])
    return dict(rv=rv, lit=lit, npl=len(lit), tord=ref['trend_order'],
                theta=np.array(ref['theta']), held_fixed=ref['laplace_held_fixed'])


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

    fig = plt.figure(figsize=(20, 34))
    gs = fig.add_gridspec(6, 2, height_ratios=[1.25, 0.5, 0.95, 1.0, 0.45, 1.0], hspace=0.42, wspace=0.22)

    # Row 1: RVs and the reference fit
    ax = fig.add_subplot(gs[0, :])
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
        rb = c['rv_based']
        tt2 = np.linspace(rv.t.min(), rv.t.max(), 3000)
        v = rb['K_ms']*np.cos(2*np.pi*(tt2 - rv.tref)/rb['P_d'] + math.radians(rb['M0_deg']))
        ax.plot(yr(tt2), v - rb['K_ms']*math.cos(math.radians(rb['M0_deg'])), ':', color='C3', lw=2.5,
                zorder=4, label="trend companion's RV, drawn to match the trend (dotted; injected "
                                "into the astrometry)")
    ax.set_xlabel('year')
    ax.set_ylabel('RV $-$ offset (m/s)')
    ax.set_title(f'1. {name}: RVs and the reference fit (offsets and jitters per dataset; errors include jitter)',
                 fontsize=22)
    hl = ax.get_legend_handles_labels()
    ax = fig.add_subplot(gs[1, :])
    ax.axis('off')
    ax.legend(*hl, fontsize=13, ncol=2, loc='upper left', frameon=False)

    # Row 2: ground truth list
    ax = fig.add_subplot(gs[2, :])
    ax.axis('off')
    ax.set_title('2. Ground truth for the synthetic DR4 data (one draw)', fontsize=22, loc='left')
    targets = [c for c in comps if c.get('role') != 'trend']
    lines = [('RV-based: one draw from the reference fit (Laplace Gaussian)', 'C0'),]
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
        lines.append(('Trend companion (reproduces the RV trend; injected, not a target): circular, '
                      f'P log-uniform from 2x the RV baseline to {TREND_A_MAX_AU:g} AU', 'C2'))
        lines.append((f"   cos i, $\\Omega${', phase' if tr['order'] == 1 else ''} uniform; K"
                      f"{'' if tr['order'] == 1 else ' and phase'} set by the trend; kept if M $\\leq$ "
                      f"{TREND_M_MAX_MSUN:g} M$_\\odot$ and the RVs fit it as well as the polynomial "
                      f"($\\Delta$lnL > $-${TREND_DLNL_MAX:g}); its J2016.0 offset and velocity are in the "
                      f"DR3 position and proper motion, so only its acceleration is injected", 'C2'))
        lines.append((f"   P = {rb['P_d']/365.25:.4g} yr,  K = {rb['K_ms']:.4g} m/s,  i = {rd['inc_deg']:.1f}"
                      f"$^\\circ$,  $\\Omega$ = {rd['Omega_deg']:.1f}$^\\circ$,  M = {dv['mass_mj']:.4g} M$_J$,  "
                      f"a = {dv['a_au']:.3g} AU,  star's astrometric semi-major axis = "
                      f"{dv['a_star_mas']:.4g} mas;  $\\Delta$lnL = {tr['dlnL_vs_polynomial']:.2f} "
                      f"(draw {tr['n_draws']})", 'C2'))
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
    # Each companion's offset over the DR4 window, and their sum, minus the
    # best-fit position and proper motion over the DR4 transits (a straight
    # line in each coordinate, least squares): the 5-parameter fit absorbs that
    # part, so what is left is what the search must model (orbits, and a trend
    # companion's curvature). The injected signal holds no parallax. Removing
    # a line is linear, so the sum of the parts is the drawn total.
    A = np.column_stack([np.ones_like(g1.t), g1.t - g1.t.mean()])

    def minus_line(c):
        xs, ys = companion_sky(g1.t, c)
        lin = [np.linalg.lstsq(A, v, rcond=None)[0] for v in (xs, ys)]
        at = lambda t, k: lin[k][0] + lin[k][1]*(t - g1.t.mean())     # noqa: E731
        x, y = companion_sky(tg, c)
        return x - at(tg, 0), y - at(tg, 1), xs - at(g1.t, 0), ys - at(g1.t, 1)
    parts = [minus_line(c) for c in comps]
    for k, (c, (x, y, _, _)) in enumerate(zip(comps, parts)):
        span = max(np.ptp(x), np.ptp(y))
        size = f'{span:.3g} mas' if span >= 0.1 else f'{span*1e3:.3g} $\\mu$as'
        if c.get('role') == 'trend':
            lab = f"trend companion alone (P = {c['rv_based']['P_d']/365.25:.0f} yr; spans {size})"
            ax.plot(x, y, '--', color='C2', lw=2.5, zorder=3, label=lab)
        else:
            lab = f"{_pretty(c['name'])} alone (P = {c['rv_based']['P_d']:.4g} d; spans {size})"
            ax.plot(x, y, '-', color=f'C{[0, 4, 5, 6][k % 4]}', lw=1.5, zorder=3, label=lab)
    x, y = sum(p[0] for p in parts), sum(p[1] for p in parts)
    xs, ys = sum(p[2] for p in parts), sum(p[3] for p in parts)
    ax.plot(x, y, '-', color='0.55', lw=5, zorder=1, label='sum: the star over the DR4 window')
    ax.plot(xs, ys, 'o', color='C1', ms=5, zorder=4, label='sum at the DR4 transits')
    ax.set_aspect('equal', 'datalim')
    ax.set_xlabel(r'$\Delta\alpha*$ (mas)')
    ax.set_ylabel(r'$\Delta\delta$ (mas)')
    ax.set_title("4b. The star's injected motion over DR4, minus best-fit position and proper motion",
                 fontsize=17)
    legs.append((ax, dict(fontsize=12)))

    for a_, kw in legs:            # legends placed where they cover no data
        figtools.legend_clear(a_, **kw)
    figtools.check_legends(fig)
    fig.suptitle(f'Synthetic Gaia DR4 astrometry of {name} (seed {truth["seed"]})', fontsize=28, y=0.91)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=80, bbox_inches='tight')
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
