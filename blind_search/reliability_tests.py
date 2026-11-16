"""Reliability tests of the blind search: does it rule out noise and find
signals?

All tests use the Gaia-4 epochs, scan angles and formal errors (93 transits) and
the simulated detection threshold of `simthreshold.py` (false-alarm probability
FAP). Each test prints PASS or FAIL against a stated criterion; results go to
blind_search/results/reliability_tests.json.

  T1  Rule out noise (calibration). N "model+noise" datasets: the zero-planet
      model (no signal) plus white noise at the formal errors, so any peak is
      noise. Each goes through the step-1 decision: its highest peak against its
      own simulated threshold, built from M "noise-sampled" datasets (sampled
      from the model fitted to that model+noise dataset, so any signal in them
      arises from noise only). If the threshold is calibrated, the probability
      that noise alone exceeds it equals the FAP. Two checks: the numbers of
      false alarms at FAP 0.1 and 0.01 must be consistent with the nominal rates
      (binomial 95% interval; with N = 45 the FAP 0.01 check is weak), and the
      model+noise dataset's p-value (the Gumbel fit to its own noise-sampled
      datasets evaluated at its highest peak) must be uniform on [0, 1]
      (Kolmogorov-Smirnov p > 0.01), which uses every dataset rather than only
      the rare exceedances. FAP 1e-3 cannot be checked with this few datasets.
  T2  Find signals. (a) Detection and recovery fractions versus SNR from the
      stored injection map (results/injection_map_M50.jsonl): about 0 below SNR
      4 and 1 above SNR 12 are expected. (b) Regression: the first
      `N_REGRESSION` injections of that map are rerun and must reproduce the
      stored peaks, thresholds and outcomes exactly (guards against code changes
      altering behaviour).
  T3  Iterative search on synthetic data (`search.run_search`, M noise-sampled
      datasets), several noise draws each: (a) no companion: none found; (b)
      Gaia-4 b alone: found, then the search stops; (c) Gaia-4 b plus a 200 d
      companion: both found, then it stops; (d) a faint companion (SNR 3): not
      detected and nothing spurious.
  T4  The real Gaia-4 pre-release data: Gaia-4 b found at 571.3 d (within 3%),
      then nothing more.

Usage: python -I blind_search/reliability_tests.py [--test] (--test: tiny sizes,
writes to /tmp, checks that everything runs.)
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, 'tests'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np  # noqa: E402

import data  # noqa: E402
import injection  # noqa: E402
import injection_map as im  # noqa: E402
import model  # noqa: E402
import periodogram as pg  # noqa: E402
import search  # noqa: E402
import simthreshold as st  # noqa: E402

N_NULL, M_NULL = 45, 50
N_REGRESSION = 5
N_DRAWS, M_SEARCH = 3, 50
FAPS = (0.1, 1e-2)
SNR_FAINT = 3.0
MAP_FILE = os.path.join(HERE, 'results/injection_map_M50.jsonl')
results = []


def check(name, ok, detail=''):
    """Print and record one PASS/FAIL line."""
    results.append(dict(name=name, passed=bool(ok), detail=detail))
    print(('PASS  ' if ok else 'FAIL  ') + name + (f'   [{detail}]' if detail else ''), flush=True)


def binomial_interval(k, n, level=0.95):
    """Clopper-Pearson interval for k successes in n trials."""
    from scipy.stats import beta
    a = (1 - level)/2
    lo = 0.0 if k == 0 else beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)


# ---- T1: noise only --------------------------------------------------------

def null_one(args):
    """One model+noise dataset (no signal): its highest peak, p-value and
    thresholds, built from M noise-sampled datasets."""
    i, M = args
    im._init()
    s = model.Stack([im._G['g']])
    rng = np.random.default_rng([4242, i])
    s.y = rng.normal(0.0, s.err)
    d = pg.delta_bic(s, [], None, im._G['per'])
    r = st.simulated_threshold(s, [], M=M, seed=900000 + i, executor='serial',
                               periods=im._G['per'], boot=2, faps=FAPS)
    peak = float(d.max())
    p = float(-np.expm1(-np.exp(-(peak - r['loc'])/r['scale'])))
    return dict(i=i, peak=peak, p=p, thr={str(f): r['thr'][f] for f in FAPS})


def test_null(n, M, ex):
    """T1: calibration of the threshold on noise-only datasets."""
    from scipy.stats import kstest
    t0 = time.time()
    rows = list(ex.map(null_one, [(i, M) for i in range(n)], chunksize=1))
    out = {}
    for f in FAPS:
        k = sum(r['peak'] > r['thr'][str(f)] for r in rows)
        lo, hi = binomial_interval(k, n)
        out[str(f)] = dict(false_alarms=k, n=n, rate=k/n, ci95=(lo, hi))
        check(f'T1 noise only: false-alarm rate consistent with FAP {f:g}',
              lo <= f <= hi,
              f'{k} of {n} = {k/n:.3f} (expected {f*n:.1f}), 95% interval '
              f'[{lo:.3f}, {hi:.3f}]')
    ks = kstest([r['p'] for r in rows], 'uniform')
    out['ks'] = dict(statistic=float(ks.statistic), pvalue=float(ks.pvalue))
    check('T1 noise only: the p-values of the highest peaks are uniform',
          ks.pvalue > 0.01,
          f'KS statistic {ks.statistic:.3f}, p = {ks.pvalue:.2f}, '
          f'{time.time() - t0:.0f} s')
    return rows, out


# ---- T2: detection versus SNR, regression ----------------------------------

def test_detection_curve():
    """T2a: detection and recovery fractions in SNR bins from the stored map."""
    rows = [json.loads(line) for line in open(MAP_FILE)]
    snr = np.array([r['snr'] for r in rows])
    det = np.array([r['detected'] for r in rows])
    rec = np.array([r['recovered'] for r in rows])
    edges = [0, 2, 4, 6, 8, 12, 20, 40, np.inf]
    table = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (snr >= lo) & (snr < hi)
        table.append(dict(snr_lo=lo, snr_hi=hi, n=int(s.sum()),
                          detected=float(det[s].mean()),
                          recovered=float(rec[s].mean())))
        print(f'   SNR {lo:g}-{hi:g}: n = {s.sum()}, detected '
              f'{det[s].mean():.3f}, recovered {rec[s].mean():.3f}')
    low, high = snr < 4, snr >= 12
    check('T2a signals: at most 2% of injections with SNR < 4 are detected',
          det[low].mean() <= 0.02,
          f'{det[low].sum()} of {low.sum()}')
    check('T2a signals: every injection with SNR >= 12 is detected',
          det[high].all(), f'{det[high].sum()} of {high.sum()}')
    check('T2a signals: at least 80% of injections with SNR >= 12 are '
          'recovered (period and mass within 25%)',
          rec[high].mean() >= 0.80, f'{rec[high].mean():.3f}')
    return table, rows


def test_regression(stored, n, ex):
    """T2b: rerun the first n injections; they must match the stored rows."""
    injs = im.draw_injections(n)
    args = [(inj, stored[0]['M'], 2026 + 1000 + inj['i'], 'sim')
            for inj in injs]
    new = list(ex.map(im.run_one, args, chunksize=1))
    old = {r['i']: r for r in stored}
    bad = [r['i'] for r in new
           if not (np.isclose(r['peak'], old[r['i']]['peak'], rtol=1e-9)
                   and np.isclose(r['thr'], old[r['i']]['thr'], rtol=1e-9)
                   and r['detected'] == old[r['i']]['detected']
                   and r['recovered'] == old[r['i']]['recovered']
                   and np.isclose(r['P_rec'], old[r['i']]['P_rec'], rtol=1e-6))]
    check(f'T2b regression: the first {n} injections reproduce the stored '
          'map results (peak, threshold, detection, recovery)', not bad,
          f'mismatches: {bad}' if bad else 'identical')


# ---- T3: iterative search on synthetic data ---------------------------------

def second_companion(tref):
    """An illustrative second companion (not a real planet): 8 M_J around the
    Gaia-4 host at 200 d, e = 0.1, i = 70 deg, omega (star) = 1.0 rad, Omega =
    0.7 rad, mean anomaly 1.0 rad at the reference time `tref`."""
    g = injection.GAIA4B
    P, e, M0 = 200.0, 0.1, 1.0
    c = injection.physical_companion(P, e, 1.0, np.radians(70.0), 0.7, 8.0,
                                     g['mstar'], g['plx'],
                                     tref - M0*P/(2*np.pi), tref)
    return c


def snr_of(stack, comp):
    """SNR of one companion: square root of the chi^2 gain of its noise-free
    signal at the formal errors, baseline refitted."""
    z = np.zeros(len(stack.channels))
    s = model.Stack(stack.channels)
    s.y = s.companion_columns(*comp['planet'])@np.asarray(comp['astro'])
    gain = s.fit([], z)[2] - s.fit([comp['planet']], z)[2]
    return float(np.sqrt(max(gain, 0.0)))


def scenario_stack(kind, seed):
    """Synthetic Gaia-4 data: reference motion + white noise + companions."""
    import test_blind_search as T
    g = data.load_gaia_dr4(T.GAIA4)
    s0 = model.Stack([g])
    base = np.linalg.lstsq(s0.baseline/g.err[:, None], g.y/g.err,
                           rcond=None)[0]
    b = injection.gaia4b(s0.tref)
    c2 = second_companion(s0.tref)
    s = model.Stack([g])
    if kind == 'none':
        T.fill(s, [], [], seed, base=base)
        return s, []
    if kind == 'b':
        T.fill(s, [b['planet']], [b['astro']], seed, base=base)
        return s, [571.3]
    if kind == 'b+200d':
        T.fill(s, [b['planet'], c2['planet']], [b['astro'], c2['astro']],
               seed, base=base)
        return s, [571.3, 200.0]
    if kind == 'faint':
        f = SNR_FAINT/snr_of(s, c2)
        c2 = dict(c2, astro=np.asarray(c2['astro'])*f)
        T.fill(s, [c2['planet']], [c2['astro']], seed, base=base)
        return s, [200.0]
    raise ValueError(kind)


def test_scenarios(n_draws, M, ex):
    """T3: the iterative search on four synthetic cases."""
    expect = {'none': 0, 'b': 1, 'b+200d': 2, 'faint': 0}
    text = {'none': 'no companion: nothing found',
            'b': 'Gaia-4 b alone: found, then stops',
            'b+200d': 'Gaia-4 b + 200 d companion: both found, then stops',
            'faint': f'faint companion (SNR {SNR_FAINT:g}): not detected'}
    out = {}
    for kind, n_exp in expect.items():
        sums = []
        for k in range(n_draws):
            s, truth = scenario_stack(kind, 5000 + 31*k)
            r = search.run_search(s, ex, M=M, truth=truth, keep=0,
                                  seed=6000 + 31*k, verbose=False)
            sums.append(search.summarize(r))
        found = [x['n_found'] for x in sums]
        ok = all(x['n_found'] == n_exp and x['n_spurious'] == 0
                 and x['n_matched'] == n_exp for x in sums)
        out[kind] = sums
        check(f'T3 {text[kind]} (all {n_draws} noise draws)', ok,
              f'companions found per draw: {found}; spurious: '
              f'{[x["n_spurious"] for x in sums]}')
    return out


# ---- T4: real data -----------------------------------------------------------

def test_real(M, ex):
    """T4: the real Gaia-4 data."""
    import test_blind_search as T
    g = data.load_gaia_dr4(T.GAIA4)
    r = search.run_search(model.Stack([g]), ex, M=M, truth=[571.3], keep=0,
                          seed=77, verbose=False)
    sm = search.summarize(r)
    check('T4 real Gaia-4 data: Gaia-4 b found at 571.3 d (within 3%), '
          'nothing more',
          sm['n_found'] == 1 and sm['n_matched'] == 1,
          f'{sm}; peaks {[round(x["peak"], 1) for x in r["steps"]]} vs '
          f'thresholds {[round(x["thr"], 1) for x in r["steps"]]}')
    return dict(summary=sm, steps=[dict(peak=x['peak'], P=x['P'], thr=x['thr'])
                                   for x in r['steps']])


def main():
    from concurrent.futures import ProcessPoolExecutor
    test = '--test' in sys.argv
    n_null, M_null, n_reg, n_draws, M_s = (
        (6, 10, 1, 1, 10) if test else
        (N_NULL, M_NULL, N_REGRESSION, N_DRAWS, M_SEARCH))
    out_json = ('/tmp/reliability_tests.json' if test else
                os.path.join(HERE, 'results/reliability_tests.json'))
    t0 = time.time()
    with ProcessPoolExecutor(max(1, (os.cpu_count() or 2) - 1)) as ex:
        null_rows, null_out = test_null(n_null, M_null, ex)
        if os.path.exists(MAP_FILE):
            table, rows = test_detection_curve()
            test_regression(rows, n_reg, ex)
        else:  # results/ is not in git: make the map with injection_map.py
            table = None
            print(f'T2 skipped: no stored map at {MAP_FILE}')
        scen = test_scenarios(n_draws, M_s, ex)
        real = test_real(M_s, ex)
    json.dump(dict(results=results, null=null_out, null_datasets=null_rows,
                   detection_table=table,
                   scenarios=scen, real=real, seconds=time.time() - t0,
                   settings=dict(n_null=n_null, M_null=M_null,
                                 n_regression=n_reg, n_draws=n_draws,
                                 M_search=M_s)),
              open(out_json, 'w'), indent=1, default=float)
    n_fail = sum(not r['passed'] for r in results)
    print(f'\n{len(results) - n_fail} passed, {n_fail} failed '
          f'({time.time() - t0:.0f} s); wrote {out_json}')
    if not test:
        sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
