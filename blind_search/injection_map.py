"""Injection/recovery map on the Gaia-4 epochs and errors (planned 2026-11-12;
run on 2026-11-13).

Each injection is a different "system" around Gaia-4 (host mass and parallax
fixed at 0.644 M_sun and 13.58 mas) with its own noise draw: data = one injected
planet + white Gaussian noise at the formal errors (no reference motion, no
jitter, no real measurements). Planet: mass log-uniform in 0.5-20 M_J; physical
semi-major axis a log-uniform in 0.5-5 AU (the period follows from Kepler's
third law); CIRCULAR ORBITS ONLY for now (ECC_MODE = 'circular'; the e ~
Beta(0.867, 3.03) draw of Kipping 2013, capped at 0.95, is implemented: set
ECC_MODE = 'beta' to restore it); omega, Omega, M0 uniform; cos i uniform
(isotropic). The period grid runs from 20 d to the longest period in the range
(a = 5 AU, 0.5 M_J: about 5090 d), so every injected planet lies inside the
searched periods. No MCMC or Octofitter fits: the search statistic is the
Delta-BIC periodogram (a linear maximum-likelihood fit with a circular template
and free jitter at each trial period), and a detected peak is then optimised as
usual by a CIRCULAR refit of its period (fit.refit_circular: scan plus polish
over a factor 1.5 either way, wide enough to reach the neighbouring grid
periods; TO REVISIT: generalise to eccentric orbits). Per injection:
  1. detection: the highest peak must exceed the simulated threshold of THAT
     dataset (noise-sampled datasets from the baseline-only model fitted to it,
     white noise at sqrt(err^2 + s^2), M noise-sampled datasets, Gumbel fit, FAP
     1e-3), as in the iterative search at step 1;
  2. period: the refitted period (the peak period if not detected); mass:
     derived from the amplitudes (A, B, F, G) fitted at that period
     (Thiele-Innes -> astrometric semi-major axis), Kepler's third law and the
     fixed host mass and parallax;
  3. recovered = detected and period AND mass both within 25% of the injected
     values.
Injections run one per process (the noise-sampled datasets inside an injection
run serially); results are appended to a JSON-lines file, so a run can be
stopped and resumed. The map is the kernel-smoothed recovered fraction in (log
a, log m) with all injections shown.

Usage:
  python -I blind_search/injection_map.py --test                       checks of
      the pieces (draws, mass inversion, map machinery); no search is run
  python -I blind_search/injection_map.py run [--n 2000] [--M 50] [--method
      sim|efap] [--out PATH]    run (or resume) the injections
      (--method efap: same injections and noise, RVSearch's empirical FAP
      threshold; no noise-sampled datasets, M unused)
  python -I blind_search/injection_map.py plot [--out PATH]            draw the
      map from the results file
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np  # noqa: E402

import data  # noqa: E402
import fit  # noqa: E402
import injection  # noqa: E402
import model  # noqa: E402
import orbit  # noqa: E402
import periodogram as pg  # noqa: E402
import simthreshold as st  # noqa: E402

GAIA4_CSV = os.path.join(ROOT, 'data_conversion/GaiaDR4_data_Gaia4.csv')
MS, PLX = injection.GAIA4B['mstar'], injection.GAIA4B['plx']
MSUN_PER_MJ = injection.GM_JUP/injection.GM_SUN
M_RANGE, A_RANGE = (0.5, 20.0), (0.5, 5.0)          # M_J, AU
E_BETA = (0.867, 3.03)
E_MAX = 0.95
ECC_MODE = 'circular'        # TO REVISIT: 'beta' = e ~ Beta(0.867, 3.03) (Kipping 2013)
TOL = 0.25
REFIT_WINDOW = 0.5          # the refit may move the period by a factor 1.5 either way: adjacent grid periods differ by at most 1.45 (3510 -> 5087 d)
FAP = 1e-3
DEFAULT_OUT = os.path.join(HERE, 'results/injection_map.jsonl')
P_MAX = 365.25*np.sqrt(A_RANGE[1]**3/(MS + M_RANGE[0]*MSUN_PER_MJ))     # longest period in the injection range (d); the periodogram reaches it


def draw_injection(i, seed):
    """Injection number i: parameters from its own random stream (seed, i), so
    the set does not depend on how many are drawn."""
    rng = np.random.default_rng([seed, i])
    m = float(np.exp(rng.uniform(np.log(M_RANGE[0]), np.log(M_RANGE[1]))))
    a = float(np.exp(rng.uniform(np.log(A_RANGE[0]), np.log(A_RANGE[1]))))
    e = 0.0 if ECC_MODE == 'circular' else float(min(rng.beta(*E_BETA), E_MAX))
    om, Om, M0 = (float(x) for x in rng.uniform(0, 2*np.pi, 3))
    inc = float(np.arccos(rng.uniform(-1, 1)))
    P = float(365.25*np.sqrt(a**3/(MS + m*MSUN_PER_MJ)))
    return dict(i=i, mass_mj=m, a_au=a, P=P, e=e, omega=om, Omega=Om, M0=M0, inc=inc)


def draw_injections(n, seed=2026):
    return [draw_injection(i, seed) for i in range(n)]


def astar_from_mass(mass_mj, P):
    """Astrometric semi-major axis (mas) of the star's orbit for a planet of
    `mass_mj` and period P (d), as in injection.physical_companion."""
    m = mass_mj*MSUN_PER_MJ
    a_au = ((MS + m)*(P/365.25)**2)**(1/3)
    return a_au*PLX*m/(MS + m)


def mass_from_astar(astar_mas, P):
    """Planet mass (M_J) from the astrometric semi-major axis (mas) and period
    (d): the inverse of `astar_from_mass` (monotonic in the mass)."""
    from scipy.optimize import brentq
    lo, hi = 1e-3, 1e4
    f = lambda mj: astar_from_mass(mj, P) - astar_mas
    if f(lo) > 0:
        return lo
    if f(hi) < 0:
        return hi
    return float(brentq(f, lo, hi, xtol=1e-12, rtol=1e-12))


_G = {}


def _init():
    if _G:
        return
    g = data.load_gaia_dr4(GAIA4_CSV)
    s0 = model.Stack([g])
    _G.update(g=g, per=pg.period_grid(s0.span, 20.0, P_MAX), span=s0.span)


def run_one(args):
    """One injection end to end (detection, period, mass). Runs in a worker
    process."""
    inj, M, seed = args[:3]
    method = args[3] if len(args) > 3 else 'sim'          # 'sim' = simulated threshold, 'efap' = RVSearch's empirical FAP threshold
    _init()
    t0 = time.time()
    s = model.Stack([_G['g']])
    rng = np.random.default_rng(seed)
    c = injection.physical_companion(inj['P'], inj['e'], inj['omega'], inj['inc'], inj['Omega'], inj['mass_mj'], MS, PLX,
                                     s.tref - inj['M0']*inj['P']/(2*np.pi), s.tref)
    c['planet'] = (inj['P'], inj['e'], inj['M0'])
    signal = s.companion_columns(*c['planet'])@np.asarray(c['astro'])
    z = np.zeros(1)
    t = model.Stack([_G['g']])
    t.y = signal
    snr = float(np.sqrt(max(t.fit([], z)[2] - t.fit([c['planet']], z)[2], 0.0)))      # noise-free signal, formal errors, baseline refitted
    s.y = signal + rng.normal(0.0, s.err)
    d = pg.delta_bic(s, [], None, _G['per'])
    k = int(np.argmax(d))
    Pk, peak = float(_G['per'][k]), float(d[k])
    if method == 'efap':
        import efap
        ef = efap.efap(d, fap=FAP)                         # threshold from the periodogram itself; no noise-sampled datasets
        thr, r = float(ef['threshold']), dict(jitter2=[np.nan])
    else:
        r = st.simulated_threshold(s, [], M=M, seed=seed + 1, executor='serial', periods=_G['per'], boot=2)
        thr = float(r['thr'][FAP])
    detected = peak > thr
    P_peak = Pk
    if detected:                                           # optimise the detected signal as usual (circular refit of the period)
        Pk = fit.refit_circular(s, [], Pk, REFIT_WINDOW)[0]
    X = np.hstack([s.baseline, s.companion_columns(Pk, 0.0, 0.0)])
    amps = s.fit_free_jitter(X)[0][s.k_base:s.k_base + 4]
    m_rec = mass_from_astar(orbit.elements_from_ti(*amps)[0], Pk)
    p_ok = abs(Pk/inj['P'] - 1) < TOL
    m_ok = abs(m_rec/inj['mass_mj'] - 1) < TOL
    return dict(inj, M=M, method=method, snr=snr, astar_uas=float(c['a_star_mas']*1e3), peak=peak, P_peak=P_peak, P_rec=Pk, m_rec=float(m_rec), thr=thr,
                jitter=float(np.sqrt(r['jitter2'][0])), detected=bool(detected), p_ok=bool(p_ok), m_ok=bool(m_ok),
                recovered=bool(detected and p_ok and m_ok), time=time.time() - t0)


def run(n, M, out, seed=2026, method='sim'):
    """Run (or resume) n injections; results are appended to `out` (JSON lines).
    `method` 'efap' uses the same injections and noise (same seed) with
    RVSearch's empirical FAP threshold instead of the simulated one (M is then
    unused)."""
    from concurrent.futures import ProcessPoolExecutor
    injs = draw_injections(n, seed)
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)['i'] for l in open(out)}
    todo = [(inj, M, seed + 1000 + inj['i'], method) for inj in injs if inj['i'] not in done]
    print(f'{len(done)} injections done, {len(todo)} to run ({method}' + (f', M = {M} noise-sampled datasets each)' if method == 'sim' else ')'), flush=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    t0, nd = time.time(), 0
    with ProcessPoolExecutor(max(1, (os.cpu_count() or 2) - 1)) as ex, open(out, 'a') as f:
        for r in ex.map(run_one, todo, chunksize=1):
            f.write(json.dumps(r) + '\n')
            f.flush()
            nd += 1
            if nd % 20 == 0:
                print(f'  {nd}/{len(todo)} done, {time.time() - t0:.0f} s elapsed', flush=True)
    print(f'done in {time.time() - t0:.0f} s', flush=True)


def kernel_completeness(la, lm, rec, ga, gm, bw=0.12):
    """Gaussian-kernel-smoothed recovered fraction at grid points (ga, gm)
    (log10 a, log10 m); bandwidth bw in dex. Returns (fraction, effective number
    of injections)."""
    w = np.exp(-0.5*(((ga[..., None] - la)/bw)**2 + ((gm[..., None] - lm)/bw)**2))
    sw = w.sum(-1)
    return (w*rec).sum(-1)/np.maximum(sw, 1e-12), sw


def plot_map(rows, out, bw=0.12):
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    la = np.log10([r['a_au'] for r in rows])
    lm = np.log10([r['mass_mj'] for r in rows])
    rec = np.array([r['recovered'] for r in rows], float)
    ga, gm = np.meshgrid(np.linspace(*np.log10(A_RANGE), 90), np.linspace(*np.log10(M_RANGE), 90))
    frac, neff = kernel_completeness(la, lm, rec, ga, gm, bw)
    fig, axes = plt.subplots(1, 2, figsize=(22, 9.5))
    span = _G.get('span', 1805.0)
    a_base = (MS*(span/365.25)**2)**(1/3)
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    plain = FuncFormatter(lambda v, _: f'{v:g}')
    thr_name = 'eFAP threshold' if rows[0].get('method') == 'efap' else 'simulated threshold'
    for ax in axes:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(*A_RANGE)
        ax.set_ylim(*M_RANGE)
        ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 3, 5]))          # plain numbers instead of 6x10^-1, 2x10^0, ...
        ax.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20]))
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(plain)
            axis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.axvline(a_base, color='0.5', ls=':', lw=2, label=f'P = baseline ({span:.0f} d): a = {a_base:.1f} AU')
        ax.set_xlabel('planet semi-major axis a (AU)')
        ax.set_ylabel('planet mass (M$_J$)')
    ax = axes[0]
    ok = rec > 0
    ax.plot(10**la[~ok], 10**lm[~ok], 'o', color='#d55e00', ms=4, alpha=0.6, label=f'not recovered ({int((~ok).sum())})')
    ax.plot(10**la[ok], 10**lm[ok], 'o', color='#0072b2', ms=4, alpha=0.6, label=f'recovered ({int(ok.sum())})')
    fm = np.where(neff > 3, frac, np.nan)                       # no contours where fewer than 3 effective injections
    cs = ax.contour(10**ga, 10**gm, fm, levels=[0.1, 0.5, 0.9], colors='k', linewidths=2.5)
    if cs.allsegs and any(len(sg) for sg in cs.allsegs):
        ax.clabel(cs, fmt='%.1f', fontsize=15)
    ax.set_title(f'{len(rows)} injections, each with its own noise and its own {thr_name}; contours: smoothed recovered fraction', fontsize=14)
    ax.legend(loc='lower left', fontsize=15, frameon=True, framealpha=0.7, facecolor='white', edgecolor='0.3')
    ax = axes[1]
    im = ax.pcolormesh(10**ga, 10**gm, np.where(neff > 3, frac, np.nan), vmin=0, vmax=1, cmap='viridis', shading='auto')
    fig.colorbar(im, ax=ax, label='recovered fraction (kernel %.2f dex)' % bw)
    ax.contour(10**ga, 10**gm, fm, levels=[0.1, 0.5, 0.9], colors='w', linewidths=2)
    ax.set_title('smoothed completeness (recovered = detected AND period and mass within 25%); blank where fewer than 3 effective injections', fontsize=13)
    ax.legend(loc='lower left', fontsize=15, frameon=True, framealpha=0.7, facecolor='white', edgecolor='0.3')
    fig.suptitle('Injection/recovery map: Gaia-4 epochs and errors, host 0.644 M$_\\odot$, 13.58 mas; Delta-BIC periodogram, no full fits', fontsize=18, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=110)
    return out


def selftest():
    """Checks of the pieces; no search is run on injected data."""
    rng = np.random.default_rng(0)
    injs = draw_injections(4000, 1)
    assert draw_injections(10, 1) == injs[:10], 'draws must not depend on how many are drawn'
    assert abs(max(x['P'] for x in injs) - P_MAX) < 0.02*P_MAX and max(x['P'] for x in injs) <= P_MAX*1.0001
    m = np.array([x['mass_mj'] for x in injs])
    a = np.array([x['a_au'] for x in injs])
    e = np.array([x['e'] for x in injs])
    ci = np.cos([x['inc'] for x in injs])
    assert m.min() >= M_RANGE[0] and m.max() <= M_RANGE[1] and a.min() >= A_RANGE[0] and a.max() <= A_RANGE[1]
    assert abs(np.mean(np.log(m)) - np.log(M_RANGE[0]*M_RANGE[1])/2) < 0.05 and abs(np.mean(np.log(a)) - np.log(A_RANGE[0]*A_RANGE[1])/2) < 0.05
    assert (np.all(e == 0.0) if ECC_MODE == 'circular' else (abs(np.mean(e) - E_BETA[0]/sum(E_BETA)) < 0.02 and e.max() <= E_MAX)) and abs(np.mean(ci)) < 0.05
    # Kepler's third law and the mass inversion round trip
    for x in injs[:50]:
        P_k = 365.25*np.sqrt(x['a_au']**3/(MS + x['mass_mj']*MSUN_PER_MJ))
        assert abs(P_k/x['P'] - 1) < 1e-12
        assert abs(mass_from_astar(astar_from_mass(x['mass_mj'], x['P']), x['P'])/x['mass_mj'] - 1) < 1e-6
        c = injection.physical_companion(x['P'], x['e'], x['omega'], x['inc'], x['Omega'], x['mass_mj'], MS, PLX, 2457000.0 - 2400000.5, 57400.0)
        assert abs(c['a_star_mas']/astar_from_mass(x['mass_mj'], x['P']) - 1) < 1e-9
        assert abs(orbit.elements_from_ti(*c['astro'])[0]/c['a_star_mas'] - 1) < 1e-6        # amplitudes -> a* recovers the injected a*
    # map machinery on outcomes with a known smooth completeness
    la, lm = np.log10(a), np.log10(m)
    truth = lambda lm_: 1/(1 + np.exp(-(lm_ - 0.5)/0.1))
    rec = (rng.uniform(size=len(m)) < truth(lm)).astype(float)
    ga, gm = np.meshgrid(np.array([0.2, 0.4, 0.6]), np.array([0.2, 0.5, 0.8]))
    frac, neff = kernel_completeness(la, lm, rec, ga, gm, 0.12)
    assert np.all(np.abs(frac - truth(gm)) < 0.12) and np.all(neff > 20), (frac, truth(gm), neff)
    print('injection_map selftest passed: draws, Kepler/mass inversion round trip, map machinery (no search run)')


if __name__ == '__main__':
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else DEFAULT_OUT
    if '--test' in sys.argv:
        selftest()
    elif 'run' in sys.argv:
        n = int(sys.argv[sys.argv.index('--n') + 1]) if '--n' in sys.argv else 2000
        M = int(sys.argv[sys.argv.index('--M') + 1]) if '--M' in sys.argv else 50
        method = sys.argv[sys.argv.index('--method') + 1] if '--method' in sys.argv else 'sim'
        run(n, M, out, method=method)
    elif 'plot' in sys.argv:
        _init()
        rows = [json.loads(l) for l in open(out)]
        stem = os.path.splitext(os.path.basename(out))[0]
        print('wrote', plot_map(rows, os.path.join(HERE, 'figures_search/%s.png' % stem)))
    else:
        print(__doc__)
