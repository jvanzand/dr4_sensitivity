"""Benchmark runner for the blind astrometric search.

The benchmark set (validation/benchmarks.toml) is a list of known RV systems,
each with one or more frozen synthetic DR4 datasets (validation/benchmarks/,
made by synthetic_dr4.py). Every version of the search runs on the same data,
so results can be compared over time.

Modes:
  generate   make the listed datasets that do not exist yet (never overwrites);
             `--system NAME --seed K` makes one (then add K to the toml)
  run        search every dataset (or `--system NAME`) with the toml settings;
             saves every step's products, step figures, a summary.json and the
             known-vs-recovered summary figure in validation/results/<NAME>_seed<K>/
             (not in git)
  redraw     redraw every figure from the saved products (no search): step
             figures, the known-vs-recovered summary figure and the
             synthetic-data process figure
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
    tol). The trend companion (role 'trend', which stands in for the RV
    trend) is matched too, but it is not a target: finding it is neither a
    recovery nor a false alarm. Returns (rows per injected companion,
    unmatched detections)."""
    left = list(range(len(found)))
    rows = []
    for c in sorted(truth['companions'], key=lambda c: -c['derived']['a_star_mas']):
        P = c['rv_based']['P_d']
        cand = [i for i in left if abs(found[i]['P_d']/P - 1) < tol]
        i = min(cand, key=lambda i: abs(found[i]['P_d']/P - 1)) if cand else None
        if i is not None:
            left.remove(i)
        rows.append(dict(name=c['name'], role=c.get('role', 'target'), P_true=P, e_true=c['rv_based']['e'],
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
        if os.path.isdir(out):                 # a previous run's steps must not outlive it (redraw reads them all)
            import glob
            for f in glob.glob(os.path.join(out, 'step*.npz')) + glob.glob(os.path.join(out, f'{stem}_step*.png')):
                os.remove(f)
        periods = [c['rv_based']['P_d'] for c in truth['companions']]
        print(f'{stem} ({label}): {stack.n} transits over {stack.span:.0f} d; injected periods '
              f'{[round(p, 2) for p in periods]}', flush=True)
        res = search.run_search(stack, M=s['M'], fap=s['fap'], pmin=s['pmin_d'],
                                max_companions=s['max_companions'], truth=periods, tol=s['match_tol'],
                                label=stem, out_dir=out, save_dir=out, truth_label='injected period',
                                e_max=s.get('e_max'), max_trend=s.get('max_trend', 0))
        found = found_orbits(res['stack'], res['known'])
        rows, spurious = match(found, truth, s['match_tol'])
        steps = [{k: v for k, v in st.items() if k not in ('dbic', 'maxima', 'known_before', 'known_after')}
                 for st in res['steps']]
        summary = dict(dataset=stem, system=name, label=label, seed=seed, settings=s,
                       code_version=code_version(), date=datetime.date.today().isoformat(),
                       n_transits=int(stack.n), steps=steps, found=found, companions=rows,
                       spurious=spurious, trend=res['trend'], trend_log=res['trend_log'])
        json.dump(summary, open(os.path.join(out, 'summary.json'), 'w'), indent=1, default=float)
        for r in rows:
            f = r['found']
            print(f"  {r['name']}: P = {r['P_true']:.2f} d, a_star = {r['a_star_true_mas']:.3g} mas -> "
                  + (f"FOUND P = {f['P_d']:.2f} d, e = {f['e']:.3f}, a_star = {f['a_star_mas']:.3g} mas"
                     if f else 'not found'))
        for f in spurious:
            print(f"  SPURIOUS: P = {f['P_d']:.2f} d, a_star = {f['a_star_mas']:.3g} mas")
        print('wrote', summary_figure(stem))


def _angle_diff(a, b):
    return (a - b + 180) % 360 - 180


def summary_figure(stem):
    """Known (injected, from the RV reference fit) versus blindly recovered
    companions of one benchmark dataset, from its saved summary.json (no
    search): (left) star's astrometric semi-major axis versus period, with
    the searched period range and the median per-transit error; (right) the
    star's sky orbit, injected versus recovered; (bottom) a table of the orbit
    elements. Written to validation/results/<stem>/<stem>_summary.png."""
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    import figtools
    import synthetic_dr4 as sd
    s = json.load(open(os.path.join(RESULTS, stem, 'summary.json')))
    stack, truth = load(stem)
    cfg = s['settings']
    pmin, pmax = cfg['pmin_d'], 2*stack.span
    sig = float(np.median(stack.err))
    comps = {c['name']: c for c in truth['companions']}
    rows = s['companions']
    targets = [r for r in rows if r.get('role', 'target') == 'target']
    n_found = sum(r['found'] is not None for r in targets)

    fig = plt.figure(figsize=(22, 17))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 0.42], hspace=0.3, wspace=0.25)

    # left: a_star versus period
    ax = fig.add_subplot(gs[0, 0])
    ax.axvspan(1e-3, pmin, color='0.85', zorder=0, label=f'not searched (P < {pmin:g} d)')
    ax.axvspan(pmax, 1e7, color='0.85', zorder=0)
    ax.axvline(stack.span, color='0.4', ls=':', lw=2, label=f'DR4 time span ({stack.span:.0f} d)')
    ax.axhline(sig, color='C2', ls='--', lw=2, label=f'median per-transit error ({sig*1e3:.0f} $\\mu$as)')
    for k, r in enumerate(rows):
        f = r['found']
        trend = r.get('role') == 'trend'
        if trend and r['P_true'] > pmax and not f:     # beyond the searched periods: listed, not drawn
            ax.plot([], [], 's', ms=14, mfc='none', mec='k', mew=2.5,
                    label=(f"trend companion (not a target): P = {r['P_true']/365.25:.0f} yr,\n"
                           f"beyond the searched periods (not drawn)"))
            continue
        ax.plot(r['P_true'], r['a_star_true_mas'], 's' if trend else 'o', ms=20, mfc='none', mec='k', mew=2.5,
                label=('trend companion (reproduces the RV trend; not a target)' if trend else
                       'known companion (injected; from the RV fit)' if r is targets[0] else None), zorder=4)
        ax.annotate(f"{sd._pretty(r['name'])}\n{('found' if trend else 'recovered') if f else ('not found' if trend else 'not recovered')}",
                    (r['P_true'], r['a_star_true_mas']), xytext=(16, 8), textcoords='offset points',
                    fontsize=17, color='C0' if f else 'C3', va='bottom')
        if f:
            ax.annotate('', (f['P_d'], f['a_star_mas']), (r['P_true'], r['a_star_true_mas']),
                        arrowprops=dict(arrowstyle='-', color='C0', lw=1.5))
    fx = [r['found'] for r in rows if r['found']]
    if fx:
        ax.plot([f['P_d'] for f in fx], [f['a_star_mas'] for f in fx], '*', color='C0', ms=22, zorder=5,
                label='blindly recovered from DR4 astrometry alone')
    if s['spurious']:
        ax.plot([f['P_d'] for f in s['spurious']], [f['a_star_mas'] for f in s['spurious']], 'X', color='C3',
                ms=18, zorder=5, label='spurious detection (no known companion)')
    drawn = [r for r in rows if not (r.get('role') == 'trend' and r['P_true'] > pmax and not r['found'])]
    allP = [r['P_true'] for r in drawn] + [f['P_d'] for f in fx + s['spurious']] + [pmin, pmax]
    alla = [r['a_star_true_mas'] for r in drawn] + [f['a_star_mas'] for f in fx + s['spurious']] + [sig]
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(min(allP)/2.5, max(allP)*2.5)
    ax.set_ylim(min(alla)/10, max(alla)*10)
    ax.set_xlabel('period (d)')
    ax.set_ylabel("star's astrometric semi-major axis (mas)")
    ax.set_title(f"Known vs recovered: {n_found} of {len(targets)} known companions recovered, "
                 f"{len(s['spurious'])} spurious", fontsize=20)
    legs = [(ax, dict(fontsize=14, markerscale=0.6))]

    # right: sky orbits of matched companions
    ax2 = fig.add_subplot(gs[0, 1])
    for k, r in enumerate(rows):
        f = r['found']
        if not f:
            continue
        c = comps[r['name']]
        tt = np.linspace(0, c['rv_based']['P_d'], 3000) + stack.t.min()
        x, y = sd.companion_sky(tt, c)
        ax2.plot(x, y, '--', color=f'C{k}', lw=2.5, label=f"{sd._pretty(r['name'])}: injected orbit")
        A, B, F, G = orbit.ti_from_elements(f['a_star_mas'], np.radians(f['inc_deg']),
                                            np.radians(f['omega_deg']), np.radians(f['Omega_deg']))
        tf = np.linspace(0, f['P_d'], 3000) + stack.t.min()
        X, Y = orbit.xy(orbit.mean_anomaly(tf, f['P_d'], np.radians(f['M0_deg']), stack.tref), f['e'])
        ax2.plot(B*X + G*Y, A*X + F*Y, '-', color=f'C{k}', lw=1.5, label=f"{sd._pretty(r['name'])}: recovered orbit")
    ax2.plot(0, 0, '+', color='k', ms=16, mew=2, label='barycentre')
    ax2.set_aspect('equal', 'datalim')
    ax2.set_xlabel(r'$\Delta\alpha*$ (mas)')
    ax2.set_ylabel(r'$\Delta\delta$ (mas)')
    ax2.set_title("The star's orbit on the sky", fontsize=20)
    legs.append((ax2, dict(fontsize=13)))

    # bottom: element table
    ax3 = fig.add_subplot(gs[1, :])
    ax3.axis('off')
    head = ['companion', 'P (d)', 'e', "star's a (mas)", 'i (deg)', r'$\omega_\star$ (deg)', r'$\Omega$ (deg)']
    cells, colors = [], []
    for r in rows:
        c, f = comps[r['name']], r['found']
        tr = [c['rv_based']['P_d'], c['rv_based']['e'], c['derived']['a_star_mas'], c['random']['inc_deg'],
              c['rv_based']['omega_star_deg'] % 360, c['random']['Omega_deg'] % 360]
        trend = r.get('role') == 'trend'
        cells.append(['trend companion' if trend else sd._pretty(r['name']) + ' (injected)', f'{tr[0]:.2f}', f'{tr[1]:.3f}', f'{tr[2]:.3g}',
                      f'{tr[3]:.1f}', f'{tr[4]:.1f}', f'{tr[5]:.1f}'])
        colors.append(['w']*7)
        if f:
            w, Om = f['omega_deg'], f['Omega_deg']
            if abs(_angle_diff(Om, tr[5])) > 90:          # the (omega, Omega) + 180 deg twin: same sky orbit
                w, Om = (w + 180) % 360, (Om + 180) % 360
            cells.append(['   recovered', f"{f['P_d']:.2f} ({f['P_d']/tr[0] - 1:+.1%})", f"{f['e']:.3f}",
                          f"{f['a_star_mas']:.3g} ({f['a_star_mas']/tr[2] - 1:+.1%})", f"{f['inc_deg']:.1f}",
                          f'{w:.1f} ({_angle_diff(w, tr[4]):+.1f})', f'{Om:.1f} ({_angle_diff(Om, tr[5]):+.1f})'])
            colors.append(['#dbe9f6']*7)
        elif trend:
            cells.append(['   not found (not a target)'] + ['-']*6)
            colors.append(['0.93']*7)
        else:
            cells.append(['   not recovered'] + ['-']*6)
            colors.append(['#f8dcdc']*7)
    for f in s['spurious']:
        cells.append(['spurious', f"{f['P_d']:.2f}", f"{f['e']:.3f}", f"{f['a_star_mas']:.3g}", f"{f['inc_deg']:.1f}",
                      f"{f['omega_deg']:.1f}", f"{f['Omega_deg']:.1f}"])
        colors.append(['#f8dcdc']*7)
    tab = ax3.table(cellText=cells, colLabels=head, cellColours=colors, loc='upper center', cellLoc='center')
    tab.auto_set_font_size(False)
    tab.set_fontsize(15)
    tab.scale(1, 2.2)
    ax3.set_title('Orbit elements: injected vs recovered (differences in brackets; recovered angles shown on the '
                  r'($\omega$, $\Omega$) or ($\omega$+180, $\Omega$+180) branch nearer the truth: same sky orbit)',
                  fontsize=15)

    for a_, kw in legs:
        figtools.legend_clear(a_, **kw)
    figtools.check_legends(fig)
    fig.suptitle(f"{s['label']} ({stem}): blind search of synthetic DR4 astrometry, {s['n_transits']} transits; "
                 f"M = {cfg['M']}, FAP {cfg['fap']:g}, match within {cfg['match_tol']:.0%} in period; "
                 f"final model: {('no trend', 'acceleration', 'acceleration + jerk')[s.get('trend', 0)]}; "
                 f"code {s['code_version']}", fontsize=19, y=0.95)
    out = os.path.join(RESULTS, stem, f'{stem}_summary.png')
    fig.savefig(out, dpi=80, bbox_inches='tight')
    plt.close(fig)
    return out


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
        if os.path.exists(os.path.join(out, 'summary.json')):
            print('wrote', summary_figure(stem))
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
