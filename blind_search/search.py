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
               method='sim'):
    """Run the search on `stack`. Returns dict(known, steps, periods): `steps`
    is one dict per search step with the peak, its period, the threshold,
    detected, matched (the truth period it matches, or None), spurious, the
    number of noise-sampled datasets with a highest peak at or above the real
    one, the fitted jitter, the time, and the arrays (dbic, maxima) needed to
    redraw that step. With `out_dir` (and the simulated threshold) a summary
    figure, `<out_dir>/<label>_step<k>.png`, is written per step. `method` is
    the detection rule: 'sim' = simulated threshold with M noise-sampled
    datasets, 'efap' = RVSearch's empirical FAP threshold (no figure). A
    detection matches a truth period p if |P/p - 1| < tol."""
    import plots
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
            title = (f'{label}: step {k + 1}, {len(known)} companion(s) in '
                     'the model')
            plots.plot_step_sim(stack, per, d, sim, known, pmin, pmax,
                                f'{out_dir}/{label}_step{k + 1}.png', title,
                                truth=truth_plot, truth_label=truth_label,
                                residual_models=models)
        if not detected:
            break
        known = new_known
    return dict(known=known, steps=steps, periods=per, unmatched=unmatched)


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
