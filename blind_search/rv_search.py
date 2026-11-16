"""Iterative blind search of radial velocities with the simulated threshold.

The same loop as `search.run_search`, with two additions for RVs: the model may
carry a linear or quadratic trend, and companions are refit with free
eccentricity (e up to `fit.E_MAX`), although the periodogram template is
circular.

At each step the current model has the known companions and a trend of order
`tau` (0, 1 or 2). Two kinds of candidate move compete:
  * add a companion at the highest peak of the Delta-BIC periodogram;
  * raise the trend order to tau + 1 (its Delta-BIC is 2 ln(L ratio) - ln N,
    the known companions held fixed, jitter free; not available at tau = 2).
The strength of the best move is S = max(highest peak, trend gain). The
threshold for S comes from M noise-sampled datasets drawn from the current
model (companions, trend, offsets and fitted jitter of each instrument, plus
white noise): each is refit like the real data, S is computed for it, and a
Gumbel distribution fitted to the M values gives the (1 - FAP) quantile. If S
exceeds the threshold the stronger move is made (a planet is added and all
companions refit together, or the trend order is raised and the companions
re-polished) and the loop repeats; otherwise the search stops. The trend order
is never lowered.

No RV noise model: the noise-sampled datasets carry white noise only.

Test mode:  python -I blind_search/rv_search.py --test   (a short synthetic
check). The HD 191939 run is validation/rv_hd191939.py.
"""

import copy
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
import fit  # noqa: E402
import model  # noqa: E402
import periodogram as pg  # noqa: E402
import simthreshold as st  # noqa: E402

WINDOW = 0.1          # refits may move a period by +-10% of its start


def trend_gain(stack, stack_hi, known):
    """Delta-BIC of raising the trend order by one: `stack_hi` is the same data
    with one more trend term; the known companions' (P, e, M0) stay fixed and
    every linear parameter and the jitter are refit in both models."""
    l0 = stack.fit_free_jitter(stack.design(known))[2]
    sh = copy.copy(stack_hi)
    sh.y = stack.y
    l1 = sh.fit_free_jitter(sh.design(known))[2]
    return float(2*(l1 - l0) - np.log(stack.n))


def one_sample(args):
    """One noise-sampled dataset: refit the known companions, then return its
    highest periodogram peak, its trend gain and (if `keep_full`) the whole
    periodogram."""
    stack, stack_hi, known, mean, s2, seed, periods, keep_full = args
    rng = np.random.default_rng(seed)
    s = copy.copy(stack)
    s.y = mean + rng.normal(0.0, np.sqrt(stack.err**2 + s2[stack.chan]))
    kn = list(known)
    if kn:
        kn, _ = fit._polish(s, [], kn, WINDOW)
    d = pg.delta_bic(s, kn, None, periods)
    gain = trend_gain(s, stack_hi, kn) if stack_hi is not None else -np.inf
    return float(d.max()), float(gain), (d if keep_full else None)


def threshold(stack, stack_hi, known, periods, M, seed, executor,
              faps=(1e-2, 1e-3), keep=0):
    """Simulated threshold for S = max(highest peak, trend gain). Returns a dict
    with the keys of `simthreshold.simulated_threshold` (maxima, loc, scale,
    thr, thr_sd, jitter2, examples) plus the peaks and gains separately."""
    mean, s2 = st.null_model(stack, known)
    args = [(stack, stack_hi, known, mean, s2, 1000003*seed + i, periods,
             i < keep) for i in range(M)]
    res = list(executor.map(one_sample, args, chunksize=2))
    peaks = np.array([r[0] for r in res])
    gains = np.array([r[1] for r in res])
    maxima = np.maximum(peaks, gains)
    loc, scale = st.gumbel_fit(maxima)
    rng = np.random.default_rng(seed + 1)
    boot = [st.gumbel_fit(rng.choice(maxima, len(maxima))) for _ in range(100)]
    return dict(maxima=maxima, peaks=peaks, gains=gains, loc=loc, scale=scale,
                jitter2=s2, periods=periods,
                examples=[r[2] for r in res if r[2] is not None],
                thr={f: float(st.gumbel_quantile(loc, scale, f))
                     for f in faps},
                thr_sd={f: float(np.std([st.gumbel_quantile(a, b, f)
                                         for a, b in boot])) for f in faps})


def amplitudes(stack, known):
    """RV semi-amplitude K (m/s) of each companion in the joint linear fit."""
    coef = stack.fit_free_jitter(stack.design(known))[0]
    k0, kc = stack.k_base, stack.k_comp_rv
    return [float(np.hypot(*coef[k0 + kc*i:k0 + kc*(i + 1)]))
            for i in range(len(known))]


def run_rv_search(channels, executor, M=50, fap=1e-3, max_steps=12, pmin=2.0,
                  max_trend=2, seed=0, label='', out_dir=None, truth=(),
                  verbose=True):
    """Run the search. Returns dict(known, trend, steps, amplitudes, periods).
    `steps` has one dict per step with the periodogram peak and period, the
    trend gain, the threshold, the move made ('planet', 'trend' or None), and
    the model before and after. With `out_dir` a step figure is written for
    each step (`truth` periods are drawn on it)."""
    import plots
    tau, known, steps = 0, [], []
    stack = model.Stack(channels, rv_trend=tau)
    pmax = 2*stack.span
    per = pg.period_grid(stack.span, pmin, pmax)
    for k in range(max_steps):
        t0 = time.time()
        stack_hi = (model.Stack(channels, rv_trend=tau + 1)
                    if tau < max_trend else None)
        d = pg.delta_bic(stack, known, None, per)
        j = int(np.argmax(d))
        peak, P = float(d[j]), float(per[j])
        gain = trend_gain(stack, stack_hi, known) if stack_hi else -np.inf
        sim = threshold(stack, stack_hi, known, per, M, seed + 17*k, executor,
                        faps=(1e-2, fap), keep=20 if out_dir else 0)
        thr = sim['thr'][fap]
        S = max(peak, gain)
        move = None if S <= thr else ('planet' if peak >= gain else 'trend')
        rec = dict(step=k + 1, trend_before=tau, known_before=list(known),
                   peak=peak, P_peak=P, trend_gain=gain, thr=thr,
                   thr_sd=sim['thr_sd'][fap], S=S, move=move,
                   jitter=[float(np.sqrt(x)) for x in sim['jitter2']])
        models = [list(known)]
        if move == 'planet':
            known = fit.add_and_refit(stack, known, P)[0]
            models.append(list(known))
        elif move == 'trend':
            tau += 1
            stack = model.Stack(channels, rv_trend=tau)
            if known:
                known, _ = fit._polish(stack, [], known, WINDOW)
        rec.update(trend_after=tau, known_after=list(known),
                   time=time.time() - t0)
        steps.append(rec)
        if verbose:
            what = {'planet': f'ADD PLANET (refit period(s) '
                    f'{[round(p[0], 2) for p in known]})',
                    'trend': f'RAISE TREND to order {tau}',
                    None: 'stop'}[move]
            print(f"  step {k + 1}: peak {peak:.1f} at {P:.2f} d, trend gain "
                  f"{gain:.1f}, threshold {thr:.1f} -> {what}  "
                  f"[{rec['time']:.0f} s]", flush=True)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            extra = (f'; raising the trend order would gain {gain:.1f}'
                     if np.isfinite(gain) else '')
            plots.plot_step_sim(
                stack if move != 'trend' else model.Stack(channels, tau - 1),
                per, d, sim, rec['known_before'], pmin, pmax,
                f'{out_dir}/{label}_step{k + 1}.png',
                f'{label}: step {k + 1}, {len(rec["known_before"])} '
                f'companion(s), trend order {rec["trend_before"]}{extra}',
                truth=truth, truth_label='known planet period',
                residual_models=models if move != 'trend' else
                [rec['known_before']], fap=fap, stat=S)
        if move is None:
            break
    return dict(known=known, trend=tau, steps=steps,
                amplitudes=amplitudes(stack, known), periods=per)


def selftest():
    """Synthetic RV system: two eccentric-free planets plus a linear trend, and
    a trend-only system. Checks that the loop finds them and stops."""
    from concurrent.futures import ProcessPoolExecutor
    import test_blind_search as T
    with ProcessPoolExecutor(3) as ex:
        s = model.Stack(T.rv_channels(5, n=160, span=3000.0), rv_trend=1)
        base = np.zeros(s.k_base)
        base[-1] = 12.0                                   # trend term, m/s
        T.fill(s, [(40.0, 0.0, 0.0), (210.0, 0.0, 0.0)],
               [[9.0, 4.0], [6.0, -5.0]], 3, base=base)
        ch = [data.Channel('rv', c.name, c.t, s.y[s.chan == i], c.err)
              for i, c in enumerate(s.channels)]
        r = run_rv_search(ch, ex, M=10, pmin=5.0, max_steps=6, seed=1)
        Ps = sorted(p[0] for p in r['known'])
        ok = (len(Ps) == 2 and abs(Ps[0]/40 - 1) < 0.03
              and abs(Ps[1]/210 - 1) < 0.05 and r['trend'] >= 1)
        assert ok, (Ps, r['trend'])
        print('rv_search selftest passed: planets', [round(p, 1) for p in Ps],
              'trend order', r['trend'])


if __name__ == '__main__':
    if '--test' in sys.argv:
        selftest()
    else:
        print(__doc__)
