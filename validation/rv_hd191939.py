"""Blind RV search of HD 191939 (RVs only; QUESTIONS.md 2026-11-13): runs
`rv_search.run_rv_search` on the project RV table and compares the detections
with the published planets (census table), matching periods within 25%.

Outputs (not in git): validation/results/rv_hd191939.json, step figures and a
summary figure in validation/results/figures/.

Usage:  python -I validation/rv_hd191939.py [--M 50]
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'blind_search'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import data  # noqa: E402
import figtools  # noqa: E402
from rv_search import run_rv_search  # noqa: E402

FIG_DIR = os.path.join(HERE, 'results', 'figures')


def match_to_truth(known, amps, truth, tol=0.25):
    """Pair detections with known planets, nearest period first, within a
    fractional period difference `tol`. truth: list of dict(name, P, K).
    Returns (rows per known planet, unmatched detections)."""
    left = list(range(len(known)))
    rows = []
    for t in sorted(truth, key=lambda x: x['P']):
        cand = [i for i in left if abs(known[i][0]/t['P'] - 1) < tol]
        i = min(cand, key=lambda i: abs(known[i][0]/t['P'] - 1)) if cand else None
        if i is not None:
            left.remove(i)
        rows.append(dict(name=t['name'], P=t['P'], K=t['K'],
                         found=i is not None,
                         P_found=known[i][0] if i is not None else None,
                         e_found=known[i][1] if i is not None else None,
                         K_found=amps[i] if i is not None else None))
    return rows, [dict(P=known[i][0], e=known[i][1], K=amps[i]) for i in left]


def hd191939_truth():
    """The published planets of HD 191939 from the census table."""
    import csv
    path = os.path.join(ROOT, 'sample/census_table.csv')
    out = []
    for r in csv.DictReader(open(path)):
        if r['hostname'] == 'HD 191939':
            out.append(dict(name=r['pl_name'].split()[-1],
                            P=float(r['period_d']), K=float(r['K_ms'])))
    return sorted(out, key=lambda x: x['P'])


def summary_figure(res, truth, rows, extra, path, title):
    """Known planets (vertical lines, with K) against the detections."""
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(18, 8))
    for i, t in enumerate(truth):            # short ticks at the bottom: the legend stays free above
        ax.plot([t['P'], t['P']], [0, 0.1], color='C2', lw=4, transform=ax.get_xaxis_transform(),
                label='published planet (period)' if i == 0 else None)
        ax.text(t['P'], 1.02, f"{t['name']}\nK = {t['K']:g}", ha='center',
                va='bottom', fontsize=14, color='C2',
                transform=ax.get_xaxis_transform())
    for r in rows:
        if r['found']:
            ax.plot(r['P_found'], r['K_found'], 'o', color='C0', ms=13,
                    label='detected, matches a published planet' if r is
                    next(x for x in rows if x['found']) else None)
    for i, x in enumerate(extra):
        ax.plot(x['P'], x['K'], 's', color='C3', ms=13,
                label='detected, no published planet within 10%' if i == 0
                else None)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('period (d)')
    ax.set_ylabel('fitted semi-amplitude K (m/s)')
    ax.set_title(title, fontsize=15, pad=60)
    fig.tight_layout()
    figtools.legend_clear(ax, fontsize=14, frameon=True, framealpha=0.7)
    figtools.check_legends(fig)
    fig.savefig(path, dpi=100)
    plt.close(fig)
    return path


def main():
    from concurrent.futures import ProcessPoolExecutor
    M = int(sys.argv[sys.argv.index('--M') + 1]) if '--M' in sys.argv else 50
    truth = hd191939_truth()
    channels = data.load_rv_project(os.path.join(
        ROOT, 'rv_data/systems/2248126315275354496_HD191939.csv'))
    print('HD 191939:', [(c.name, len(c.t)) for c in channels], flush=True)
    fig_dir = FIG_DIR
    t0 = time.time()
    with ProcessPoolExecutor(max(1, (os.cpu_count() or 2) - 1)) as ex:
        res = run_rv_search(channels, ex, M=M, label='hd191939', out_dir=fig_dir,
                            truth=[t['P'] for t in truth], seed=191939)
    rows, extra = match_to_truth(res['known'], res['amplitudes'], truth)
    print(f'\nfinal trend order {res["trend"]}; {len(res["known"])} planets')
    for r in rows:
        print(f"  {r['name']}: published P = {r['P']:g} d, K = {r['K']:g} m/s -> "
              + (f"FOUND P = {r['P_found']:.2f} d, e = {r['e_found']:.2f}, "
                 f"K = {r['K_found']:.2f} m/s" if r['found'] else 'not found'))
    for x in extra:
        print(f"  extra detection: P = {x['P']:.2f} d, e = {x['e']:.2f}, "
              f"K = {x['K']:.2f} m/s")
    out = dict(M=M, seconds=time.time() - t0, trend=res['trend'], rows=rows,
               extra=extra, steps=[{k: v for k, v in s.items()}
                                   for s in res['steps']])
    json.dump(out, open(os.path.join(HERE, 'results', 'rv_hd191939.json'), 'w'),
              indent=1, default=float)
    title = (f'HD 191939 (RVs only, M = {M} noise-sampled datasets, FAP 0.001): '
             f'{sum(r["found"] for r in rows)} of {len(rows)} published '
             f'planets found, {len(extra)} extra')
    print('wrote', summary_figure(res, truth, rows, extra,
                                  os.path.join(fig_dir, 'hd191939_summary.png'),
                                  title))
    print(f'done in {time.time() - t0:.0f} s')



if __name__ == '__main__':
    main()
