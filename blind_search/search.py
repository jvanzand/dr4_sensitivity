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


def run_search(stack, executor=None, M=50, fap=1e-3, max_companions=4,
               seed=0, keep=20, noise='fitted', truth=(), label='',
               out_dir=None, pmin=20.0, tol=0.03,
               truth_label='injected / published period', verbose=True,
               method='sim', save_dir=None, e_max=None):
    """Run the search on `stack`. Returns dict(known, steps, periods): `steps`
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
    cap of every refit (`fit.set_e_max`; default 0.7, `fit.E_MAX`)."""
    import plots
    if e_max is not None:
        fit.set_e_max(e_max)
    pmax = 2*stack.span
    per = pg.period_grid(stack.span, pmin, pmax)
    truth_plot = list(truth)
    known, steps, unmatched = [], [], list(truth)
    for k in range(max_companions + 1):
        t0 = time.time()
        d = pg.delta_bic(stack, known, None, per)
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
        match = None
        if detected:
            near = [p for p in unmatched if abs(P/p - 1) < tol]
            if near:
                match = min(near, key=lambda p: abs(1/P - 1/p))
                unmatched.remove(match)
        new_known = fit.add_and_refit(stack, known, P)[0] if detected else None
        rec = dict(step=k + 1, method=method, peak=peak, P=P, thr=thr,
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
            note = (f" (matches injected {match:g} d)" if match else
                    (' (SPURIOUS: no injected companion)' if rec['spurious']
                     else ''))
            print(f"  step {k + 1}: highest peak {peak:.1f} at {P:.1f} d vs "
                  f"threshold {thr:.1f} -> "
                  f"{'DETECTED' if detected else 'not detected'}{note}"
                  f"  [{rec['time']:.0f} s]", flush=True)
        if out_dir and sim is not None:
            os.makedirs(out_dir, exist_ok=True)
            models = [known] + ([new_known] if detected else [])
            title = step_title(label, k + 1, len(known))
            plots.plot_step_sim(stack, per, d, sim, known, pmin, pmax,
                                f'{out_dir}/{label}_step{k + 1}.png', title,
                                truth=truth_plot, truth_label=truth_label,
                                residual_models=models)
        if not detected:
            break
        known = new_known
    return dict(known=known, steps=steps, periods=per, unmatched=unmatched)


def step_title(label, step, n_known):
    return f'{label}: step {step}, {n_known} companion(s) in the model'


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
               decision=json.dumps({k: rec[k] for k in ('step', 'method', 'peak', 'P', 'thr', 'thr_sd',
                                                        'detected', 'matched', 'spurious', 'n_fake_ge',
                                                        'jitter', 'time')}, default=float),
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
