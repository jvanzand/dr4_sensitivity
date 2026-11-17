"""The iterative blind search with the simulated per-step threshold
(QUESTIONS.md, 2026-11-11/12).

Loop (k = 1, 2, ...): compute the Delta-BIC periodogram of the current (k -
1)-companion model (`periodogram.delta_bic`, jitter free); build the detection
threshold by simulation (`simthreshold.simulated_threshold`: M noise-sampled
datasets from the fitted (k - 1)-companion model, each refit and searched, the
threshold from the Gumbel fit to their highest peaks at the chosen false-alarm
probability); if the highest peak of the real data exceeds the threshold, add a
companion at that period and refit all companions (`fit.add_and_refit`) and
continue; otherwise stop. `truth` periods (injection tests) are only used to
label each detection as matching an injected period (within `tol`) or spurious.

Trends (`max_trend` > 0; trend.py, RVSearch's logic): before step 1 the trend
order is raised while each term is warranted at `fap` (RV slope then
curvature; astrometric acceleration then jerk). At each step the periodogram
lets the new companion replace trend terms; a detected companion is refit at
every trend order up to the current one and the lowest-BIC combination kept
(`add_with_trend_choice`), and the trend is then re-checked (a term no longer
warranted is dropped). The simulated threshold uses the
same statistic.

Test mode:  python -I blind_search/search.py --test
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit  # noqa: E402
import periodogram as pg  # noqa: E402
import simthreshold as st  # noqa: E402
import trend as tr  # noqa: E402


def run_search(stack, executor=None, M=50, fap=1e-3, max_companions=4,
               seed=0, keep=20, noise='fitted', truth=(), label='',
               out_dir=None, pmin=20.0, tol=0.03,
               truth_label='injected / published period', verbose=True,
               method='sim', save_dir=None, e_max=None, max_trend=0):
    """Run the search on `stack`. Returns dict(known, steps, periods, trend,
    trend_log, stack): `steps`
    is one dict per search step with the peak, its period, the threshold,
    detected, matched (the truth period it matches, or None), spurious, the
    number of noise-sampled datasets with a highest peak at or above the real
    one, the fitted jitter, the time, and the arrays (dbic, maxima) needed to
    redraw that step. With `out_dir` (and the simulated threshold) a summary
    figure, `<out_dir>/<label>_step<k>.png`, is written per step. `method` is
    the detection rule: 'sim' = simulated threshold with M noise-sampled
    datasets, 'efap' = RVSearch's empirical FAP threshold (no figure). A
    detection matches a truth period p if |P/p - 1| < tol. With `save_dir`
    every step's arrays are saved (`save_step`), so `plots.redraw_step` can
    redraw its figure without searching again. `e_max` sets the eccentricity
    cap of every refit (`fit.set_e_max`; default 0.7, `fit.E_MAX`).
    `max_trend` = highest trend order allowed (0 = no trend, the default of
    the injection map and reliability tests; 2 = up to curvature / jerk); the
    returned `stack` is at the final trend order (its coefficients follow
    `stack.baseline_names`), `trend_log` lists every trend test."""
    import plots
    if e_max is not None:
        fit.set_e_max(e_max)
    pmax = 2*stack.span
    per = pg.period_grid(stack.span, pmin, pmax)
    truth_plot = list(truth)
    known, steps, unmatched = [], [], list(truth)
    trend_log = []
    if max_trend:
        stack, trend_log = tr.raise_trend(stack.with_trend(0), [], fap, max_trend)
        if verbose:
            for x in trend_log:
                print(f"  trend test: {x['action']}: 2 dlnL = {x['gain_2dlnL']:.1f} ({x['dof']} dof, "
                      f"FAP {x['fap']:.1e}) -> {'adopted' if x['adopted'] else 'not warranted'}", flush=True)
    for k in range(max_companions + 1):
        t0 = time.time()
        d, orders = pg.delta_bic(stack, known, None, per, return_order=True)
        sim = None
        if method == 'sim':
            sim = st.simulated_threshold(
                stack, known, M=M, seed=seed + 17*k, executor=executor,
                periods=per, keep=keep if out_dir else 0, noise=noise)
            thr = sim['thr'][fap]
        elif method == 'efap':
            import efap
            thr = float(efap.efap(d, fap=fap)['threshold'])
        else:
            raise ValueError(f"method must be 'sim' or 'efap', not {method!r}")
        j = int(np.argmax(d))
        peak, P = float(d[j]), float(per[j])
        detected = peak > thr
        new_known, new_stack, moves, P_match = None, stack, [], P
        if detected:
            if stack.trend:
                new_stack, new_known, moves = add_with_trend_choice(stack, known, per, fap)
                P_match = new_known[-1][0]            # the new companion's refit period
            else:
                new_known = fit.add_and_refit(stack, known, P)[0]
            trend_log += moves
        match = None
        if detected:
            near = [p for p in unmatched if abs(P_match/p - 1) < tol]
            if near:
                match = min(near, key=lambda p: abs(1/P_match - 1/p))
                unmatched.remove(match)
        rec = dict(step=k + 1, method=method, peak=peak, P=P, P_fit=new_known[-1][0] if detected else None,
                   thr=thr, trend_before=stack.trend, trend_after=new_stack.trend, trend_moves=moves,
                   thr_sd=sim['thr_sd'][fap] if sim else None,
                   detected=detected, matched=match,
                   spurious=bool(detected and truth and match is None),
                   n_fake_ge=int(np.sum(sim['maxima'] >= peak)) if sim else None,
                   jitter=float(np.sqrt(sim['jitter2'][0])) if sim else None,
                   time=time.time() - t0, dbic=d,
                   maxima=sim['maxima'] if sim else None,
                   known_before=list(known), known_after=new_known)
        steps.append(rec)
        if save_dir:
            save_step(save_dir, rec, per, sim, dict(pmin=pmin, pmax=pmax, fap=fap, label=label,
                                                    truth=list(truth), truth_label=truth_label))
        if verbose:
            note = ((f" (refit P = {P_match:.1f} d)" if detected and P_match != P else '')
                    + (f" (matches injected {match:g} d)" if match else
                       ' (SPURIOUS: no injected companion)' if rec['spurious'] else ''))
            print(f"  step {k + 1}: highest peak {peak:.1f} at {P:.1f} d vs "
                  f"threshold {thr:.1f} -> "
                  f"{'DETECTED' if detected else 'not detected'}{note}"
                  + (f"; trend order {stack.trend} -> {new_stack.trend}" if new_stack.trend != stack.trend else '')
                  + f"  [{rec['time']:.0f} s]", flush=True)
        if out_dir and sim is not None:
            os.makedirs(out_dir, exist_ok=True)
            models = [known] + ([new_known] if detected else [])
            title = step_title(label, k + 1, len(known), stack.trend,
                               'rv' if stack.has_rv and not stack.has_astro else 'astro')
            plots.plot_step_sim(stack, per, d, sim, known, pmin, pmax,
                                f'{out_dir}/{label}_step{k + 1}.png', title,
                                truth=truth_plot, truth_label=truth_label,
                                residual_models=models)
        if not detected:
            break
        known, stack = new_known, new_stack
    return dict(known=known, steps=steps, periods=per, unmatched=unmatched, trend=stack.trend,
                trend_log=trend_log, stack=stack)


def add_with_trend_choice(stack, known, periods, fap):
    """Add the detected companion when the model has a trend (RVSearch's
    logic, decided at the eccentric fit): for each trend order tau' = 0 ..
    stack.trend, the companion is added at the highest peak of the periodogram
    at that order and refit (`fit.add_and_refit`); the combination with the
    lowest BIC is kept, so the companion may replace trend terms. (The circular
    periodogram alone favours keeping the trend when an eccentric orbit of
    period ~ the time span is present, since the trend terms help the circular
    template; e.g. pi Men b in DR4.) The trend is then re-checked
    (`trend.settle_trend`). Returns (stack, companions, trend moves)."""
    best = None
    for tau in range(stack.trend, -1, -1):
        s = stack.with_trend(tau)
        d = pg.delta_bic(s, known, None, periods, drop_trend=False)
        pl, lnL = fit.add_and_refit(s, known, float(periods[int(np.argmax(d))]))
        b = -2*lnL + s.n_params(len(pl))*np.log(s.n)
        if best is None or b < best[0]:
            best = (b, s, pl)
    _, s, pl = best
    moves = ([dict(action=f'companion replaces trend terms ({stack.trend} -> {s.trend})')]
             if s.trend < stack.trend else [])
    if s.trend:
        s, pl, log = tr.settle_trend(s, pl, fap)
        moves += log
    return s, pl, moves


TREND_NAMES = {'astro': ('no trend', 'acceleration', 'acceleration + jerk'),
               'rv': ('no trend', 'linear trend', 'quadratic trend')}


def step_title(label, step, n_known, trend=0, kind='astro'):
    t = f', {TREND_NAMES[kind][trend]}' if trend else ''
    return f'{label}: step {step}, {n_known} companion(s){t} in the model'


def save_step(save_dir, rec, periods, sim, meta):
    """One search step's products in `<save_dir>/step<k>.npz`: the period grid
    and real Delta-BIC periodogram; the step's decision (peak, period,
    threshold, detected, matched); the companions (P, e, M0) before and after;
    and, for the simulated threshold, the noise-sampled highest peaks
    (`maxima`), the example periodograms, the Gumbel fit, the jitter and the
    thresholds at every FAP. `meta` (pmin, pmax, fap, label, truth periods and
    label) goes in as JSON text."""
    import json
    os.makedirs(save_dir, exist_ok=True)
    arr = dict(periods=periods, dbic=rec['dbic'],
               known_before=np.array(rec['known_before'], float).reshape(-1, 3),
               known_after=np.array(rec['known_after'] or [], float).reshape(-1, 3),
               decision=json.dumps({k: rec.get(k) for k in ('step', 'method', 'peak', 'P', 'thr', 'thr_sd',
                                                            'detected', 'matched', 'spurious', 'n_fake_ge',
                                                            'jitter', 'time', 'trend_before', 'trend_after',
                                                            'trend_moves')}, default=float),
               meta=json.dumps(meta, default=float))
    if sim is not None:
        faps = sorted(sim['thr'])
        arr.update(maxima=sim['maxima'], loc=sim['loc'], scale=sim['scale'],
                   jitter2=np.atleast_1d(sim['jitter2']),
                   examples=np.array(sim['examples']).reshape(-1, len(periods)),
                   faps=np.array(faps), thr=np.array([sim['thr'][f] for f in faps]),
                   thr_sd=np.array([sim['thr_sd'][f] for f in faps]))
    path = os.path.join(save_dir, f'step{rec["step"]}.npz')
    np.savez_compressed(path, **arr)
    return path


def load_step(path):
    """A saved step (`save_step`): (decision dict, meta dict, arrays dict,
    sim dict or None) with the companions as lists of (P, e, M0)."""
    import json
    z = np.load(path)
    a = {k: z[k] for k in z.files}
    dec, meta = json.loads(str(a.pop('decision'))), json.loads(str(a.pop('meta')))
    a['known_before'] = [tuple(r) for r in a['known_before']]
    a['known_after'] = [tuple(r) for r in a['known_after']] or None
    sim = None
    if 'maxima' in a:
        faps = [float(f) for f in a['faps']]
        sim = dict(maxima=a['maxima'], loc=float(a['loc']), scale=float(a['scale']),
                   jitter2=a['jitter2'], examples=list(a['examples']), periods=a['periods'],
                   thr=dict(zip(faps, map(float, a['thr']))),
                   thr_sd=dict(zip(faps, map(float, a['thr_sd']))))
    return dec, meta, a, sim


def summarize(res):
    """One-line summary: companions found, of which matching injected, spurious,
    missed."""
    st_ = res['steps']
    found = [s for s in st_ if s['detected']]
    return dict(n_found=len(found), n_matched=sum(s['matched'] is not None for s in found),
                n_spurious=sum(s['spurious'] for s in found), missed=list(res['unmatched']))


if __name__ == '__main__':
    if '--test' in sys.argv:
        import data
        import model
        sys.path.insert(0, os.path.join(data.ROOT, 'tests'))
        import test_blind_search as T
        g = data.load_gaia_dr4(T.GAIA4)
        s = model.Stack([g])
        T.fill(s, [(300.0, 0.0, 0.0)], [[0.25, -0.15, 0.2, 0.3]], 3)
        r = run_search(s, M=10, max_companions=2, truth=[300.0], keep=0)
        sm = summarize(r)
        assert sm['n_found'] >= 1 and sm['n_matched'] == 1 and sm['n_spurious'] == 0, sm
        print('search.py test passed:', sm)
