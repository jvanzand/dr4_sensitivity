"""Benchmark runner for the blind astrometric search.

The benchmark set (validation/benchmarks.toml) is a list of known RV systems,
each with one or more frozen synthetic DR4 datasets (validation/benchmarks/,
made by synthetic_dr4.py). Every version of the search runs on the same data,
so results can be compared over time.

Modes:
  generate   make the listed datasets that do not exist yet (never overwrites);
             `--system NAME --seed K` makes one (then add K to the toml)
  run        search every dataset (or `--system NAME`) with the toml settings;
             saves every step's products, step figures and a summary.json in
             validation/results/<NAME>_seed<K>/ (not in git)
  redraw     redraw every figure from the saved products (no search): step
             figures and the synthetic-data process figure
  compare    table of the current results against the stored reference
             (validation/benchmarks/reference.json, in git);
             `--update-reference` replaces the reference with the current
             results

Matching: an injected companion is recovered if a detection's period is within
`match_tol` of its period (nearest first, one detection per companion); a
detection matching no injected companion is spurious.

Usage:  python -I validation/benchmark.py run [--system HD39091]
Test mode:  python -I validation/benchmark.py --test   (the toml parses, every
            listed dataset exists, a saved step reloads; no search)
"""

import datetime
import json
import os
import subprocess
import sys
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _d in ('blind_search', 'injection'):
    sys.path.insert(0, os.path.join(ROOT, _d))
sys.path.insert(0, HERE)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import data  # noqa: E402
import model  # noqa: E402
import orbit  # noqa: E402

CONFIG = os.path.join(HERE, 'benchmarks.toml')
BENCH = os.path.join(HERE, 'benchmarks')
RESULTS = os.path.join(HERE, 'results')
REFERENCE = os.path.join(BENCH, 'reference.json')


def config():
    return tomllib.load(open(CONFIG, 'rb'))


def datasets(cfg, system=None):
    """(name, label, seed, stem) of every listed dataset (or one system's)."""
    out = []
    for s in cfg['system']:
        if system and s['name'] != system:
            continue
        out += [(s['name'], s['label'], k, f"{s['name']}_seed{k}") for k in s['seeds']]
    return out


def load(stem):
    truth = json.load(open(os.path.join(BENCH, f'{stem}_truth.json')))
    stack = model.Stack([data.load_gaia_dr4(os.path.join(BENCH, f'{stem}.csv'))])
    return stack, truth


def code_version():
    try:
        r = subprocess.run(['git', '-C', ROOT, 'describe', '--always', '--dirty'],
                           capture_output=True, text=True)
        return r.stdout.strip() or 'unknown'
    except OSError:
        return 'unknown'


def found_orbits(stack, known):
    """Each detected companion's (P, e, M0) and its orbit from the fitted
    Thiele-Innes constants: star's semi-major axis (mas), i, omega, Omega."""
    if not known:
        return []
    coef = stack.fit_free_jitter(stack.design(known))[0]
    out = []
    for i, (P, e, M0) in enumerate(known):
        k = stack.k_base + stack.k_comp_astro*i
        A, B, F, G = coef[k:k + 4]
        a, inc, w, Om = orbit.elements_from_ti(A, B, F, G)
        out.append(dict(P_d=float(P), e=float(e), M0_deg=float(np.degrees(M0) % 360), a_star_mas=float(a),
                        inc_deg=float(np.degrees(inc)), omega_deg=float(np.degrees(w) % 360),
                        Omega_deg=float(np.degrees(Om) % 360)))
    return out


def match(found, truth, tol):
    """Pair injected companions with detections (nearest period first, within
    tol). Returns (rows per injected companion, unmatched detections)."""
    left = list(range(len(found)))
    rows = []
    for c in sorted(truth['companions'], key=lambda c: -c['derived']['a_star_mas']):
        P = c['rv_based']['P_d']
        cand = [i for i in left if abs(found[i]['P_d']/P - 1) < tol]
        i = min(cand, key=lambda i: abs(found[i]['P_d']/P - 1)) if cand else None
        if i is not None:
            left.remove(i)
        rows.append(dict(name=c['name'], P_true=P, e_true=c['rv_based']['e'],
                         a_star_true_mas=c['derived']['a_star_mas'], inc_true_deg=c['random']['inc_deg'],
                         found=found[i] if i is not None else None))
    return rows, [found[i] for i in left]


def run(system=None):
    import search
    cfg = config()
    s = cfg['search']
    for name, label, seed, stem in datasets(cfg, system):
        stack, truth = load(stem)
        out = os.path.join(RESULTS, stem)
        periods = [c['rv_based']['P_d'] for c in truth['companions']]
        print(f'{stem} ({label}): {stack.n} transits over {stack.span:.0f} d; injected periods '
              f'{[round(p, 2) for p in periods]}', flush=True)
        res = search.run_search(stack, M=s['M'], fap=s['fap'], pmin=s['pmin_d'],
                                max_companions=s['max_companions'], truth=periods, tol=s['match_tol'],
                                label=stem, out_dir=out, save_dir=out, truth_label='injected period')
        found = found_orbits(stack, res['known'])
        rows, spurious = match(found, truth, s['match_tol'])
        steps = [{k: v for k, v in st.items() if k not in ('dbic', 'maxima', 'known_before', 'known_after')}
                 for st in res['steps']]
        summary = dict(dataset=stem, system=name, label=label, seed=seed, settings=s,
                       code_version=code_version(), date=datetime.date.today().isoformat(),
                       n_transits=int(stack.n), steps=steps, found=found, companions=rows,
                       spurious=spurious)
        json.dump(summary, open(os.path.join(out, 'summary.json'), 'w'), indent=1, default=float)
        for r in rows:
            f = r['found']
            print(f"  {r['name']}: P = {r['P_true']:.2f} d, a_star = {r['a_star_true_mas']:.3g} mas -> "
                  + (f"FOUND P = {f['P_d']:.2f} d, e = {f['e']:.3f}, a_star = {f['a_star_mas']:.3g} mas"
                     if f else 'not found'))
        for f in spurious:
            print(f"  SPURIOUS: P = {f['P_d']:.2f} d, a_star = {f['a_star_mas']:.3g} mas")


def redraw(system=None):
    import plots
    import synthetic_dr4
    cfg = config()
    for name, label, seed, stem in datasets(cfg, system):
        stack, _ = load(stem)
        out = os.path.join(RESULTS, stem)
        k = 1
        while os.path.exists(os.path.join(out, f'step{k}.npz')):
            print('wrote', plots.redraw_step(stack, os.path.join(out, f'step{k}.npz'),
                                             os.path.join(out, f'{stem}_step{k}.png')))
            k += 1
        synthetic_dr4.redraw(name, seed)


def generate(system=None, seed=None):
    import synthetic_dr4
    cfg = config()
    todo = ([(system, None, seed, f'{system}_seed{seed}')] if system and seed is not None
            else datasets(cfg, system))
    for name, _, k, stem in todo:
        if os.path.exists(os.path.join(BENCH, f'{stem}_truth.json')):
            print(f'{stem}: exists (frozen), skipped')
            continue
        synthetic_dr4.make_system(name, k)
        if seed is not None and k not in next(s['seeds'] for s in cfg['system'] if s['name'] == name):
            print(f'add seed {k} to {name} in {CONFIG}')


def _summaries(cfg):
    out = {}
    for _, _, _, stem in datasets(cfg):
        p = os.path.join(RESULTS, stem, 'summary.json')
        if os.path.exists(p):
            out[stem] = json.load(open(p))
    return out


def compare(update=False):
    cfg = config()
    cur = _summaries(cfg)
    if update:
        json.dump(cur, open(REFERENCE, 'w'), indent=1, default=float)
        print(f'reference updated: {sorted(cur)} -> {REFERENCE}')
        return
    ref = json.load(open(REFERENCE)) if os.path.exists(REFERENCE) else {}

    def cell(f):
        return 'missed' if f is None else f"P {f['P_d']:.1f}, a {f['a_star_mas']:.3g}, e {f['e']:.2f}"
    for stem in sorted(set(cur) | set(ref)):
        c, r = cur.get(stem), ref.get(stem)
        head = (f"{stem}: reference {r['code_version']} ({r['date']})" if r else f'{stem}: no reference')
        head += f" | current {c['code_version']} ({c['date']})" if c else ' | not run'
        print(head)
        rows_c = {x['name']: x for x in c['companions']} if c else {}
        rows_r = {x['name']: x for x in r['companions']} if r else {}
        for nm in list(dict.fromkeys(list(rows_r) + list(rows_c))):
            x = rows_c.get(nm) or rows_r.get(nm)
            fr, fc = (rows_r.get(nm) or {}).get('found'), (rows_c.get(nm) or {}).get('found')
            changed = (fr is None) != (fc is None) or (fr and fc and abs(fc['P_d']/fr['P_d'] - 1) > 0.01)
            print(f"  {nm:12s} truth P {x['P_true']:.1f}, a {x['a_star_true_mas']:.3g} | "
                  f"ref {cell(fr) if r else '-':32s} | now {cell(fc) if c else '-':32s}"
                  + ('  CHANGED' if changed and r and c else ''))
        ns = (len(r['spurious']) if r else None, len(c['spurious']) if c else None)
        print(f'  spurious detections: reference {ns[0]}, current {ns[1]}'
              + ('  CHANGED' if None not in ns and ns[0] != ns[1] else ''))


def selftest():
    import search
    cfg = config()
    for name, _, k, stem in datasets(cfg):
        for suf in ('.csv', '_noplanet.csv', '_truth.json'):
            assert os.path.exists(os.path.join(BENCH, stem + suf)), stem + suf
    stack, truth = load(datasets(cfg)[0][3])
    assert stack.n > 0 and truth['companions']
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        rec = dict(step=1, method='sim', peak=1.0, P=100.0, thr=0.5, thr_sd=0.1, detected=True, matched=None,
                   spurious=False, n_fake_ge=0, jitter=0.0, time=0.0, dbic=np.zeros(5),
                   known_before=[], known_after=[(100.0, 0.1, 0.2)])
        sim = dict(maxima=np.arange(3.0), loc=0.0, scale=1.0, jitter2=np.zeros(1), examples=[np.ones(5)],
                   thr={1e-3: 2.0}, thr_sd={1e-3: 0.1})
        p = search.save_step(d, rec, np.arange(5.0), sim, dict(pmin=1, pmax=5, fap=1e-3, label='t', truth=[],
                                                                truth_label='x'))
        dec, meta, a, s = search.load_step(p)
        assert dec['detected'] and a['known_after'] == [(100.0, 0.1, 0.2)] and s['thr'][1e-3] == 2.0
    print(f'benchmark.py test passed: {len(datasets(cfg))} dataset(s) listed and present')


if __name__ == '__main__':
    args = sys.argv[1:]
    opt = lambda k: args[args.index(k) + 1] if k in args else None   # noqa: E731
    if '--test' in args:
        selftest()
    elif 'run' in args:
        run(opt('--system'))
    elif 'redraw' in args:
        redraw(opt('--system'))
    elif 'generate' in args:
        generate(opt('--system'), int(opt('--seed')) if opt('--seed') else None)
    elif 'compare' in args:
        compare('--update-reference' in args)
    else:
        print(__doc__)
