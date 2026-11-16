"""Tests of the blind-search pieces (blind_search/): orbit geometry, data
channels, the linear model, the Delta-BIC periodogram and the eFAP threshold.
Octofitter-free; the V1 comparison with Octofitter's likelihood is a separate
Julia test.

Run from the project root, in the gaia_dr4_sensitivity conda environment:
    python -I tests/test_blind_search.py            # all tests (about a minute)
    python -I tests/test_blind_search.py --quick    # skip the noise-calibration
        test
"""

import os
import sys
import time

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'blind_search'))
sys.path.insert(0, os.path.join(ROOT, 'injection'))
import data  # noqa: E402
import efap as ef  # noqa: E402
import model  # noqa: E402
import injection  # noqa: E402
import orbit  # noqa: E402
import periodogram as pg  # noqa: E402

failures = []


def check(name, ok, detail=''):
    print(('PASS  ' if ok else 'FAIL  ') + name + (f'   [{detail}]' if detail else ''))
    if not ok:
        failures.append(name)


GAIA4 = os.path.join(ROOT, 'data_conversion/GaiaDR4_data_Gaia4.csv')


def rv_channels(seed=0, n=150, span=3000.0, err=(1.0, 3.0)):
    """Two synthetic RV instruments over `span` days (MJD 58000...), y = 0 until
    filled."""
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(58000, 58000 + span, n))
    lab = rng.integers(0, 2, n)
    e = rng.uniform(*err, n)
    return [data.Channel('rv', f'rv:inst{k}', t[lab == k], np.zeros((lab == k).sum()), e[lab == k])
            for k in (0, 1)]


def fill(stack, planets, amps, noise_seed, jit=0.0, base=None):
    """Replace stack.y by the model with the given companion amplitudes plus
    Gaussian noise (formal errors and jitter). `amps` is one amplitude vector
    per planet, in the order of Stack.companion_columns."""
    rng = np.random.default_rng(noise_seed)
    y = np.zeros(stack.n) if base is None else stack.baseline@base
    for p, a in zip(planets, amps):
        y = y + stack.companion_columns(*p)@np.asarray(a, float)
    sig = stack.sigma(np.full(len(stack.channels), jit))
    stack.y = y + rng.normal(0, sig)
    return stack


## ---- orbit geometry ----
def test_orbit():
    rng = np.random.default_rng(1)
    M = rng.uniform(-10, 10, 400)
    worst = 0.0
    for e in (0.0, 0.1, 0.5, 0.9, 0.99):
        E = orbit.kepler_E(M, e)
        worst = max(worst, np.max(np.abs(np.mod(E - e*np.sin(E) - M + np.pi, 2*np.pi) - np.pi)))
    check('Kepler equation solved to 1e-10 for e up to 0.99', worst < 1e-10, f'{worst:.1e}')

    Mg = rng.uniform(-10, 10, 20000)
    same = max(np.max(np.abs(np.mod(orbit.kepler_E(Mg, e) - orbit.kepler_E_numpy(Mg, e) + np.pi, 2*np.pi) - np.pi))
               for e in (0.01, 0.3, 0.7, 0.95, 0.99))
    t0 = time.time(); orbit.kepler_E(Mg, 0.6); tc = time.time() - t0
    t0 = time.time(); orbit.kepler_E_numpy(Mg, 0.6); tn = time.time() - t0
    check(f'C solver in use ({orbit.C_SOLVER}) and agrees with the numpy solver to 1e-10', orbit.C_SOLVER
          and same < 1e-10, f'max diff {same:.1e}; 20000 solves: C {tc*1e3:.1f} ms, numpy {tn*1e3:.1f} ms')

    a, inc, om, Om = 1.7, 1.1, 2.3, 4.0
    A, B, F, G = orbit.ti_from_elements(a, inc, om, Om)
    a2, i2, o2, O2 = orbit.elements_from_ti(A, B, F, G)
    same = np.allclose(orbit.ti_from_elements(a2, i2, o2, O2), (A, B, F, G), atol=1e-12)
    check('Thiele-Innes round trip (elements -> A,B,F,G -> elements -> A,B,F,G)', same
          and abs(a2 - a) < 1e-12 and abs(i2 - inc) < 1e-12)

    # Astrometric columns against direct geometry: star offsets from the orbit in the sky plane.
    e, P, M0, tref = 0.35, 120.0, 0.7, 500.0
    t = rng.uniform(0, 400, 60)
    psi = rng.uniform(0, 2*np.pi, 60)
    Mt = orbit.mean_anomaly(t, P, M0, tref)
    cols = orbit.astro_columns(Mt, e, np.sin(psi), np.cos(psi))
    X, Y = orbit.xy(Mt, e)
    north = A*X + F*Y
    east = B*X + G*Y
    direct = east*np.sin(psi) + north*np.cos(psi)
    check('astrometric columns reproduce the along-scan projection of the Thiele-Innes orbit',
          np.allclose(cols@np.array([A, B, F, G]), direct, atol=1e-12))

    # Radial velocity columns against the textbook expression with nu from atan2.
    K, w = 23.0, 0.9
    cn, sn = orbit.true_anomaly_terms(Mt, e)
    nu = np.arctan2(sn, cn)
    textbook = K*(np.cos(nu + w) + e*np.cos(w))
    lin = orbit.rv_columns(Mt, e)@np.array([K*np.cos(w), K*np.sin(w)])
    check('RV columns reproduce K[cos(nu+w) + e cos w]', np.allclose(lin, textbook, atol=1e-12))

    c1 = orbit.rv_columns(Mt, 1e-9)
    c0 = orbit.rv_columns(Mt, 0.0)
    check('columns are continuous as e -> 0 (no singularity in M0)', np.allclose(c1, c0, atol=1e-7))


## ---- physical injection ----
def test_injection():
    c = injection.gaia4b(57500.0)
    check('Gaia-4 b: a star orbit of 273 uas (the catalog lists 272 uas) and K = 368 m/s',
          abs(c['a_star_mas']*1e3 - 272.5) < 3 and abs(c['K'] - 368.0) < 0.5,
          f"{c['a_star_mas']*1e3:.1f} uas, K = {c['K']:.1f}")
    # a_star sin i = K sqrt(1 - e^2) P / (2 pi) ties the astrometric and RV amplitudes together
    p = injection.GAIA4B
    a_sini_au = c['K']*np.sqrt(1 - p['e']**2)*p['P']*injection.DAY_S/(2*np.pi)/injection.AU_M
    a_star_au = c['a_star_mas']/p['plx']
    check('astrometric and RV amplitudes agree (a_star sin i = K sqrt(1 - e^2) P / 2 pi)',
          abs(a_star_au*np.sin(c['inc'])/a_sini_au - 1) < 2e-3, f"{a_star_au*np.sin(c['inc'])/a_sini_au:.5f}")
    a, i, o, O = orbit.elements_from_ti(*c['astro'])
    check('Thiele-Innes constants of the injection give back the injected semi-major axis and inclination',
          abs(a - c['a_star_mas']) < 1e-9 and abs(i - c['inc']) < 1e-9)


## ---- data channels ----
def test_data():
    g = data.load_gaia_dr4(GAIA4)
    check('Gaia-4 pre-release: 93 transits after the transit-level collapse',
          len(g) == 93 and np.all(np.diff(g.t) >= 0))
    gd = data.load_gaia_dr4(GAIA4, reduction='demo')
    check('median transit error: tutorial reduction 0.038 mas, demo reduction 0.082 mas (PLAN.md)',
          abs(np.median(g.err) - 0.038) < 0.002 and abs(np.median(gd.err) - 0.082) < 0.002,
          f'{np.median(g.err):.4f}, {np.median(gd.err):.4f}')
    check('scan-angle columns are unit vectors', np.allclose(g.sinpsi**2 + g.cospsi**2, 1))
    rv = data.load_rv_simple(os.path.join(ROOT, 'data_conversion/RV_data_HD114762.csv'))
    allc = [g] + rv
    check('select: gaia_dr4 / astro / rv / all',
          [c.name for c in data.select(allc, 'gaia_dr4')] == ['gaia_dr4']
          and len(data.select(allc, 'rv')) == len(rv)
          and len(data.select(allc, 'astro')) == 1 and len(data.select(allc, 'all')) == len(allc))
    try:
        data.select(allc, 'nonsense')
        check('select: unknown name raises', False)
    except ValueError:
        check('select: unknown name raises', True)


## ---- linear model ----
def test_model():
    g = data.load_gaia_dr4(GAIA4)
    rv = rv_channels(0)
    s = model.Stack([g] + rv, rv_trend=1)
    check('stack sizes: 5 astro + 2 offsets + 1 trend baseline columns; 6 amplitudes per companion',
          s.k_base == 8 and s.companion_columns(100., 0.1, 0.).shape == (s.n, 6))
    check('astrometric amplitudes vanish on RV rows and vice versa',
          np.all(s.companion_columns(100., .1, 0.)[s.is_rv][:, :4] == 0)
          and np.all(s.companion_columns(100., .1, 0.)[s.is_astro][:, 4:] == 0))
    pl = [(250.0, 0.3, 1.0)]
    amps = [[0.4, -0.2, 0.3, 0.5, 15.0, -8.0]]
    base = np.random.default_rng(3).normal(0, 1, s.k_base)
    fill(s, pl, amps, 0, base=base)
    s.y = s.y - 0   # noisy
    # noiseless version: exact recovery of every linear parameter
    s2 = model.Stack([g] + rv_channels(0), rv_trend=1)
    s2.y = s2.baseline@base + s2.companion_columns(*pl[0])@np.array(amps[0])
    coef, r, chi2 = s2.fit(pl, np.zeros(len(s2.channels)))
    check('noiseless fit recovers all linear parameters exactly',
          np.allclose(coef, np.concatenate([base, amps[0]]), atol=1e-7) and chi2 < 1e-9,
          f'max diff {np.max(np.abs(coef - np.concatenate([base, amps[0]]))):.1e}')
    # jitter: noise of 5 m/s on top of formal errors is recovered
    s3 = model.Stack(rv_channels(5, n=400), rv_trend=0)
    fill(s3, [], [], 7, jit=5.0)
    coef, r, chi2 = s3.fit([], np.zeros(2))
    j = s3.fit_jitter(r)
    check('maximum-likelihood jitter recovers 5 m/s within 15%', np.all(np.abs(j/5.0 - 1) < 0.15),
          f'{j.round(2)}')


## ---- eFAP ----
def test_efap():
    rng = np.random.default_rng(0)
    x = rng.exponential(2.0, 5000)
    r = ef.efap(x)
    check('eFAP recovers the decay rate of exponential noise (A = 0.5)', abs(r['A'] - 0.5) < 0.1,
          f"A = {r['A']:.3f}")
    peak = np.append(x, r['threshold'] + 5)
    r2 = ef.efap(peak)
    check('detection (max > threshold) is the same as fap_min < FAP',
          (peak.max() > r2['threshold']) == (r2['fap_min'] < r2['fap']))
    # coverage on pure exponential noise: fraction of realizations exceeding the threshold
    hits = sum(np.max(v := rng.exponential(2.0, 1000)) > ef.efap(v)['threshold'] for _ in range(2000))
    check('false-alarm rate on exponential noise is at most the nominal 1e-3', hits/2000 <= 2e-3,
          f'{hits} of 2000')
    r3 = ef.efap(np.full(500, 3.0))
    check('degenerate periodogram: infinite threshold, no detection', not r3['ok']
          and np.isinf(r3['threshold']))
    # a periodogram with a few real peaks: the threshold must sit above the noise tail
    v = rng.exponential(2.0, 3000)
    v[:5] = [40, 38, 35, 30, 28]
    check('a handful of strong peaks do not drag the threshold into the noise',
          ef.efap(v)['threshold'] < 40, f"{ef.efap(v)['threshold']:.1f}")


## ---- periodogram ----
def detect(stack, jit, pmin, pmax):
    """Periodogram and eFAP; jit=None refits the jitter at every trial period
    (the search's default)."""
    periods = pg.period_grid(stack.span, pmin, pmax)
    d = pg.delta_bic(stack, [], jit, periods)
    r = ef.efap(d)
    return periods, d, r


def test_periodogram():
    g = data.load_gaia_dr4(GAIA4)
    # astrometry only, real Gaia-4 scan geometry; a circular companion at P = 300 d
    s = model.Stack([g])
    jit = None
    fill(s, [(300.0, 0.0, 0.0)], [[0.25, -0.15, 0.2, 0.3]], 1)
    t0 = time.time()
    periods, d, r = detect(s, jit, 20.0, 2*s.span)
    dt = time.time() - t0
    P = periods[np.argmax(d)]
    check('astrometry only: injected P = 300 d recovered, above the eFAP threshold',
          abs(P/300 - 1) < 0.03 and d.max() > r['threshold'], f"peak {P:.0f} d, dBIC {d.max():.0f}, "
          f"eFAP threshold {r['threshold']:.0f}; {len(periods)} periods in {dt:.2f} s")
    # eccentric injection (e = 0.4): the circular periodogram still peaks at P (or its first harmonic)
    fill(s, [(300.0, 0.4, 1.0)], [[0.25, -0.15, 0.2, 0.3]], 2)
    periods, d, r = detect(s, jit, 20.0, 2*s.span)
    P = periods[np.argmax(d)]
    check('astrometry only: e = 0.4 companion found at P or P/2, above the threshold',
          min(abs(P/300 - 1), abs(P/150 - 1)) < 0.05 and d.max() > r['threshold'],
          f"peak {P:.0f} d, dBIC {d.max():.0f}, eFAP threshold {r['threshold']:.0f}")
    # RV only
    s = model.Stack(rv_channels(2, n=200))
    fill(s, [(40.0, 0.0, 0.0)], [[10.0, 5.0]], 3)
    t0 = time.time()
    periods, d, r = detect(s, None, 2.0, 2*s.span)
    dt = time.time() - t0
    P = periods[np.argmax(d)]
    check('RV only: injected P = 40 d recovered, above the eFAP threshold',
          abs(P/40 - 1) < 0.01 and d.max() > r['threshold'], f'peak {P:.2f} d; {len(periods)} periods '
          f'in {dt:.1f} s')
    # second step: with the first companion in the model, the second shows up
    s = model.Stack(rv_channels(4, n=200))
    pl = [(40.0, 0.0, 0.0), (210.0, 0.0, 0.0)]
    fill(s, pl, [[10.0, 5.0], [6.0, -6.0]], 5)
    periods = pg.period_grid(s.span, 2.0, 2*s.span)
    d = pg.delta_bic(s, [pl[0]], None, periods)
    P = periods[np.argmax(d)]
    check('second step: with companion 1 in the model the periodogram peaks at companion 2',
          abs(P/210 - 1) < 0.03 and d.max() > ef.efap(d)['threshold'], f'peak {P:.0f} d')
    # joint: both data types selected
    s = model.Stack([g] + rv_channels(6, n=120, span=1800.0, err=(1.5, 2.5)))
    fill(s, [(300.0, 0.0, 0.0)], [[0.25, -0.15, 0.2, 0.3, 6.0, 3.0]], 7)
    periods, d, r = detect(s, None, 20.0, 2*s.span)
    P = periods[np.argmax(d)]
    check('astrometry + RVs jointly: P = 300 d recovered, above the threshold',
          abs(P/300 - 1) < 0.03 and d.max() > r['threshold'],
          f"peak {P:.0f} d, dBIC {d.max():.0f}, eFAP threshold {r['threshold']:.0f}")


def test_gaia4b_injection():
    """Gaia-4 b injected into simulated noise at the real epochs, with a
    realistic reference motion (nonzero baseline) that the search must fit
    itself."""
    g = data.load_gaia_dr4(GAIA4)
    s0 = model.Stack([g])
    c = injection.gaia4b(s0.tref)
    base = np.linalg.lstsq(s0.baseline/g.err[:, None], g.y/g.err, rcond=None)[0]
    s0.y = None
    s = model.Stack([g])
    fill(s, [c['planet']], [c['astro']], 11, base=base)
    model_y = s.companion_columns(*c['planet'])@c['astro']
    z = (s.y - s.baseline@base - model_y)/s.err
    check('injected data minus reference motion minus model is pure noise (chi2 about N)',
          abs(np.sum(z**2) - s.n) < 4*np.sqrt(2*s.n), f'chi2 {np.sum(z**2):.0f} for {s.n}')
    periods = pg.period_grid(s.span, 20.0, 2*s.span)
    d = pg.delta_bic(s, [], None, periods)
    P = periods[np.argmax(d)]
    thr = ef.efap(d)['threshold']
    check('Gaia-4 b (S/N 38), jitter refitted at every period: peak within 3% of 571.3 d and above the plain eFAP '
          'threshold, baseline fitted', abs(P/571.3 - 1) < 0.03 and d.max() > thr,
          f'peak {P:.1f} d, dBIC {d.max():.0f}, threshold {thr:.0f}')
    dfix = pg.delta_bic(s, [], np.zeros(1), periods)
    thr_fix = ef.efap(dfix)['threshold']
    check('same data with the jitter held fixed: the peak is NOT above the eFAP threshold (why the jitter is refitted)',
          dfix.max() < thr_fix, f'peak {dfix.max():.0f}, threshold {thr_fix:.0f}')


def test_refit_known():
    """refit_known=True (RVSearch convention: known companions' P, e, M0
    re-optimised at every trial period)."""
    g = data.load_gaia_dr4(GAIA4)
    s = model.Stack([g])
    fill(s, [(300.0, 0.3, 1.0)], [[0.25, -0.15, 0.2, 0.3]], 5)
    periods = pg.period_grid(s.span, 20.0, 2*s.span)[::60]
    a = pg.delta_bic(s, [], None, periods)
    b = pg.delta_bic(s, [], None, periods, refit_known=True)
    check('refit_known with no known companions is identical to the default', np.allclose(a, b))
    known = [(300.0, 0.3, 1.0)]
    fixed = pg.delta_bic(s, known, None, periods)
    free = pg.delta_bic(s, known, None, periods, refit_known=True)
    check('refit_known with one known companion: finite, and within 15 of the fixed-known values',
          np.all(np.isfinite(free)) and np.max(np.abs(free - fixed)) < 15,
          f'max |free - fixed| = {np.max(np.abs(free - fixed)):.1f}')
    try:
        pg.delta_bic(s, known, np.zeros(1), periods, refit_known=True)
        ok = False
    except ValueError:
        ok = True
    check('refit_known with a fixed jitter is refused', ok)


def test_simthreshold():
    """Simulated per-step threshold: maxima finite, noise-only maxima in the
    white-noise range, a clear signal far above."""
    import simthreshold as st
    g = data.load_gaia_dr4(GAIA4)
    s = fill(model.Stack([g]), [], [], 3)
    r = st.simulated_threshold(s, [], M=16, seed=2)
    check('simulated threshold (noise only, M = 16): maxima finite, median in the white-noise range (-16 to -3), '
          'FAP 1e-3 threshold above FAP 1e-2', np.all(np.isfinite(r['maxima'])) and -16 < np.median(r['maxima']) < -3
          and r['thr'][1e-3] > r['thr'][1e-2], f"median {np.median(r['maxima']):.1f}, thresholds {r['thr']}")
    periods = pg.period_grid(s.span, 20.0, 2*s.span)
    sig = fill(model.Stack([g]), [(300.0, 0.0, 0.0)], [[0.25, -0.15, 0.2, 0.3]], 3)
    peak = pg.delta_bic(sig, [], None, periods).max()
    r2 = st.simulated_threshold(sig, [], M=16, seed=2)
    check('a clear companion (peak %.0f) is far above its own simulated threshold, which stays at the noise level' % peak,
          peak > r2['thr'][1e-3] + 50 and r2['thr'][1e-3] < 40, f"threshold {r2['thr'][1e-3]:.1f}")


def test_search():
    """The iterative search with the simulated threshold: one injected companion
    found, no spurious one."""
    import search
    g = data.load_gaia_dr4(GAIA4)
    s = fill(model.Stack([g]), [(300.0, 0.0, 0.0)], [[0.25, -0.15, 0.2, 0.3]], 3)
    r = search.run_search(s, M=8, max_companions=2, truth=[300.0], keep=0, verbose=False)
    sm = search.summarize(r)
    check('search loop: injected 300 d companion found and matched, no spurious companion, stops after one',
          sm['n_found'] == 1 and sm['n_matched'] == 1 and sm['n_spurious'] == 0, str(sm))


def test_injection_map_pieces():
    """Injection-map scaffold: parameter draws, Kepler/mass inversion round
    trip, kernel completeness (no search is run)."""
    import injection_map as im
    try:
        im.selftest()
        ok, msg = True, ''
    except AssertionError as err:
        ok, msg = False, str(err)[:200]
    check('injection_map pieces: draws in range, mass inversion round trip, amplitudes -> a*, kernel completeness recovers a known curve', ok, msg)


def test_refit_circular():
    """fit.refit_circular: on noise-free data a 3600 d circular companion is
    recovered from a start 26% away."""
    import fit
    g = data.load_gaia_dr4(GAIA4)
    s = model.Stack([g])
    s.y = s.companion_columns(3600.0, 0.0, 0.0)@np.asarray([0.2, 0.1, -0.15, 0.05])
    P, _ = fit.refit_circular(s, [], 2680.0, 0.5)
    check('circular refit recovers a noise-free 3600 d period from a start at 2680 d (window factor 1.5)', abs(P/3600.0 - 1) < 1e-3, f'{P:.1f} d')


def test_calibration(n=100):
    """Noise-only periodograms on the real Gaia-4 scan geometry: how often does
    the eFAP threshold (nominal FAP 1e-3 per periodogram) fire?"""
    g = data.load_gaia_dr4(GAIA4)
    hits, t0 = 0, time.time()
    for k in range(n):
        s = model.Stack([g])
        fill(s, [], [], 1000 + k)
        periods, d, r = detect(s, None, 20.0, 2*s.span)
        hits += d.max() > r['threshold']
    check('noise only (Gaia-4 geometry): false alarms at most 2 in 100 at nominal FAP 1e-3', hits <= 2,
          f'{hits} of {n}; {time.time() - t0:.0f} s')


if __name__ == '__main__':
    for fn in (test_orbit, test_injection, test_data, test_model, test_efap, test_periodogram,
               test_gaia4b_injection, test_refit_known, test_simthreshold, test_search, test_injection_map_pieces, test_refit_circular):
        fn()
    if '--quick' not in sys.argv:
        test_calibration()
    print(f'\n{len(failures)} failed' if failures else '\nall passed')
    sys.exit(1 if failures else 0)
