"""Study sample, stages 1 and 2 (PLAN.md, Step 1b), built from census_table.csv.

Stage 1, the snapshot catalog: every census system with at least one precisely
characterized companion (`precise` in census_table.csv: orbit and mass at >= 3 sigma,
where a missing error passes), with all of its census companions. S/N plays no part.

Stage 2, the fit set: the stage-1 systems in which any companion has forecast S/N >=
SNR_MIN. The stage-1 membership never changes with SNR_MIN or the resonance policy.

Resonance (2026-11-06): resonant companions are poorly described by independent Keplerians,
and adopting literature orbits for them would defeat the purpose of refitting. The check
here is a crude period-ratio proxy: two census companions of a system whose period ratio
lies within RES_TOL (fractional) of a low-order ratio (RES_RATIOS). It over-flags (chance
near-commensurabilities) and misses companions outside the census cuts, so it is a
second-tier check, to be consulted after a system gives cause for concern (e.g., a
Keplerian fit that will not converge). The default policy is therefore `flag`: nothing is
dropped, and every flagged companion is marked (`resonance_member`) and counted.
Policy `exclude` applies the intended rule: resonant companions stay in the model (so the
non-resonant ones fit well) but are not scored as recoveries (`scored` = False), and a
system whose companions are all resonant is not fit (`fit_set` = False). Every output
table carries the flags, and `stage1_resonance.png` shows the flagged companions.

Systems are keyed by Gaia DR3 source_id (host name when the host has none).

Usage (from the project root, after census.py):
    python -I sample/catalog.py                       # SNR_MIN = 0.1, resonance flagged only
    python -I sample/catalog.py --snr-min 0.5 --resonance-policy exclude
    python -I sample/catalog.py --selftest            # checks the resonance rule
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS_FILE = os.path.join(HERE, 'census_table.csv')
CATALOG_FILE = os.path.join(HERE, 'stage1_catalog.csv')
SYSTEMS_FILE = os.path.join(HERE, 'stage1_systems.csv')
RESONANCE_FILE = os.path.join(HERE, 'resonance_candidates.csv')
MISSING_FILE = os.path.join(HERE, 'missing_errors.csv')
PLOT_FILE = os.path.join(HERE, 'stage1_resonance.png')

SNR_MIN_DEFAULT = 0.1

## Low-order period ratios P_outer/P_inner: first-order (2:1, 3:2, 4:3, 5:4), second-order
## (5:3, 5:2), and 3:1
RES_RATIOS = {'2:1': 2.0, '3:2': 1.5, '4:3': 4/3, '5:4': 1.25, '5:3': 5/3, '5:2': 2.5,
              '3:1': 3.0}
RES_TOL = 0.01


def resonant_pairs(periods, tol=RES_TOL):
    """Pairs (i, j, label, offset) of companions whose period ratio is within a fractional
    tol of a ratio in RES_RATIOS; i is the inner companion."""
    out = []
    for i, j in itertools.combinations(range(len(periods)), 2):
        inner, outer = sorted((i, j), key=lambda k: periods[k])
        ratio = periods[outer]/periods[inner]
        for label, value in RES_RATIOS.items():
            if abs(ratio/value - 1) <= tol:
                out.append((inner, outer, label, ratio/value - 1))
    return out


def selftest():
    p = [100.0, 200.5, 331.0, 1000.0]
    pairs = resonant_pairs(p, 0.01)
    labels = {(i, j, lab) for i, j, lab, _ in pairs}
    assert (0, 1, '2:1') in labels, labels          # 2.005 within 1% of 2:1
    assert (1, 2, '3:2') not in labels              # 1.651: not within 1% of 3:2 or 5:3
    assert (0, 3, '3:1') not in labels              # ratio 10
    assert resonant_pairs([100.0, 200.0, 300.0], 0.01) == [
        (0, 1, '2:1', 0.0), (0, 2, '3:1', 0.0), (1, 2, '3:2', 0.0)]
    assert resonant_pairs([100.0, 250.0], 0.001)[0][2] == '5:2'
    assert resonant_pairs([100.0, 237.0], 0.01) == []   # 2.37, between 2:1 and 5:2
    print('selftest passed')




def build(snr_min=SNR_MIN_DEFAULT, policy='flag', res_tol=RES_TOL):
    c = pd.read_csv(CENSUS_FILE, dtype={'source_id': str}, low_memory=False)
    c['system_id'] = c.source_id.fillna('host:' + c.hostname)
    c['missing_errors'] = c.missing_errors.fillna('')     # blank cells read back as NaN

    ## Resonance candidates, by pair and by companion
    rows = []
    for sid, sys in c.groupby('system_id'):
        if len(sys) < 2:
            continue
        sys = sys.reset_index(drop=True)
        for i, j, label, off in resonant_pairs(sys.period_d.to_numpy(), res_tol):
            rows.append({'system_id': sid, 'hostname': sys.hostname[i], 'ratio': label,
                         'offset_frac': off, 'inner': sys.pl_name[i],
                         'inner_period_d': sys.period_d[i], 'outer': sys.pl_name[j],
                         'outer_period_d': sys.period_d[j],
                         'inner_mass_mj': sys.mass_mj[i], 'outer_mass_mj': sys.mass_mj[j]})
    res = pd.DataFrame(rows)
    members = set(res.inner) | set(res.outer) if len(res) else set()
    c['resonance_member'] = c.pl_name.isin(members)

    ## Stage 1
    n_precise = c.groupby('system_id').precise.sum()
    stage1_ids = n_precise.index[n_precise > 0]
    cat = c[c.system_id.isin(stage1_ids)].copy()

    ## One row per system
    g = cat.groupby('system_id')
    systems = pd.DataFrame({
        'hostname': g.hostname.first(),
        'n_companions': g.size(),
        'n_precise': g.precise.sum(),
        'n_precision_unverified': g.precision_unverified.sum(),
        'n_missing_errors': g.missing_errors.apply(lambda s: (s != '').sum()),
        'any_msini_err_from_K': g.msini_err_from_K.any(),
        'any_standin_st_mass_err': g.st_mass_err_source.apply(
            lambda s: s.str.startswith('stand-in').any()),
        'max_snr': g.snr.max(),
        'n_recoverable': g.recoverable.sum(),
        'n_resonant_members': g.resonance_member.sum(),
        'G': g.G.first(), 'dist_pc': g.dist_pc.first(),
        'gaia_match': g.gaia_match.first(),
    })
    systems['resonance_candidate'] = systems.n_resonant_members > 0
    systems['all_resonant'] = systems.n_resonant_members == systems.n_companions
    systems['snr_ok'] = systems.max_snr >= snr_min
    systems['fit_set'] = systems.snr_ok & ~((policy == 'exclude') & systems.all_resonant)
    systems = systems.reset_index()          # system_id becomes a column
    cat = cat.merge(systems[['system_id', 'resonance_candidate', 'all_resonant', 'fit_set']],
                    on='system_id', how='left')
    ## Every companion stays in the model; under `exclude`, resonant ones are not scored
    cat['in_model'] = True
    cat['scored'] = ~((policy == 'exclude') & cat.resonance_member)
    return c, cat, systems, res


def plot(cat, snr_min, policy):
    import matplotlib
    matplotlib.rc_file(os.path.join(os.path.dirname(HERE), 'matplotlibrc'))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 8))
    snr = cat.snr.clip(lower=1e-7)
    groups = [
        ('companion in a system below the S/N cut,\nnot flagged', ~cat.fit_set & ~cat.resonance_member,
         dict(c='0.75', marker='o', s=14)),
        ('companion in a system below the S/N cut,\nresonance-flagged', ~cat.fit_set & cat.resonance_member,
         dict(c='none', edgecolors='tab:red', marker='o', s=60)),
        ('companion in a fit-set system,\nnot flagged', cat.fit_set & ~cat.resonance_member,
         dict(c='tab:blue', marker='o', s=20)),
        ('companion in a fit-set system,\nresonance-flagged\n'
         f'(ratio within {RES_TOL:.0%} of a low-order one)', cat.fit_set & cat.resonance_member,
         dict(c='tab:red', marker='D', s=45)),
    ]
    for label, m, kw in groups:
        n_sys = cat.system_id[m].nunique()        # a system can appear in more than one group
        ax.scatter(cat.period_d[m], snr[m], **kw,
                   label=f'{label}\n({int(m.sum())} companions in {n_sys} systems)')
    ax.axhline(snr_min, c='k', ls='--', lw=1, label=f'SNR_MIN = {snr_min}')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('companion period [d]')
    ax.set_ylabel('forecast astrometric S/N')
    ax.set_title(f'Stage-1 catalog by system (resonance policy: {policy}):\n'
                 'resonance-flagged companions are marked, not removed')
    ## Invisible entry explaining why the system counts do not add up to the stage totals
    ax.plot([], [], ' ', label='A system with companions in two groups\nis counted in both.')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14), fontsize=12, frameon=False, ncol=2)
    fig.savefig(PLOT_FILE, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--snr-min', type=float, default=SNR_MIN_DEFAULT,
                        help='stage-2 forecast S/N threshold (default 0.1)')
    parser.add_argument('--resonance-policy', choices=['flag', 'exclude'], default='flag',
                        help="'flag' marks resonance candidates only (default); 'exclude' "
                             'leaves resonant companions unscored and skips all-resonant systems')
    parser.add_argument('--res-tol', type=float, default=RES_TOL,
                        help='fractional tolerance on low-order period ratios')
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        selftest()
        raise SystemExit

    c, cat, systems, res = build(args.snr_min, args.resonance_policy, args.res_tol)
    cat.sort_values(['system_id', 'period_d']).to_csv(CATALOG_FILE, index=False)
    systems.to_csv(SYSTEMS_FILE, index=False)
    res.to_csv(RESONANCE_FILE, index=False)
    ## Companions whose errors no route could supply (concerning: neither the quantity's
    ## error nor a related one's was published)
    c[c.missing_errors != ''][['pl_name', 'hostname', 'catalog', 'fit_type', 'mass_type',
                               'missing_errors', 'period_d', 'period_err_source',
                               'mass_mj', 'mass_err_source', 'K_ms', 'sma_au',
                               'sma_err_source']].to_csv(MISSING_FILE, index=False)
    plot(cat, args.snr_min, args.resonance_policy)

    fit = cat[cat.fit_set]
    mult = systems.n_companions > 1
    print(f'Census: {len(c)} companions, {c.system_id.nunique()} systems')
    print(f'  precisely characterized companions: {int(c.precise.sum())} '
          f'({int(c.precision_unverified.sum())} only because an error is missing)')
    print(f'  companions with an error that no route could supply: '
          f'{int((c.missing_errors != "").sum())} (see {os.path.relpath(MISSING_FILE)})')
    print(f'Stage 1 (>= 1 precise companion): {len(systems)} systems, {len(cat)} companions')
    print(f'  systems with a derived M sin i error from K: '
          f'{int(systems.any_msini_err_from_K.sum())}; with a stand-in stellar-mass error: '
          f'{int(systems.any_standin_st_mass_err.sum())}; qualifying on a missing error only: '
          f'{int(((systems.n_precise - systems.n_precision_unverified) == 0).sum())}')
    print(f'Resonance proxy (pairs within {args.res_tol:.0%} of 2:1, 3:2, 4:3, 5:4, 5:3, 5:2, 3:1):')
    print(f'  stage 1: {int(cat.resonance_member.sum())} flagged companions in '
          f'{int(systems.resonance_candidate.sum())} systems; {int(systems.all_resonant.sum())} '
          f'systems have only flagged companions (these would not be fit under "exclude")')
    print(f'  fit set: {int(fit.resonance_member.sum())} flagged companions in '
          f'{int(systems[systems.fit_set].resonance_candidate.sum())} systems; under "exclude" '
          f'{int(systems[systems.fit_set].all_resonant.sum())} fit-set systems would be skipped '
          f'and {int(fit.resonance_member.sum())} companions would stay in the model unscored')
    print(f'Stage 2 (SNR_MIN = {args.snr_min}, policy {args.resonance_policy}): '
          f'{int(systems.fit_set.sum())} systems, {len(fit)} companions, '
          f'{int(fit.scored.sum())} scored')
