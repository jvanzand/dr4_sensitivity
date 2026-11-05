"""How many census companions each study supplies, and how many are recoverable.

Each companion is credited to the study whose fit the census adopts: the reference of the
chosen NEA solution, or the literature source. A study that supersedes the discovery paper
therefore gets the credit. One row per study, with two bars: all census companions whose
fit comes from that study, and the plausibly recoverable subset.

- NEA studies are shown if they supply at least MIN_COMPANIONS companions or any
  recoverable one; the rest are summed in one "other NEA studies" row, so every
  recoverable companion appears in a named row.
- Literature sources are all shown, grouped by tier (the precedence tiers of
  lit_search/README.md).
- Bars longer than X_CAP are cut at X_CAP (hatched end) with the full count printed.
- Each bar is annotated with up to two example hosts, drawn at random (fixed seed) from
  that bar, as space allows.

The rows run down three columns: the NEA studies fill the first two, the literature
the third.

Usage (from the project root, after census.py):
    python -I -W ignore sample/lit_sources_plot.py
"""

import html
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census import NEA_PS_FILE, P_MAX_YR, SNR_THRESHOLD

matplotlib.rc_file(os.path.join(HERE, '..', 'matplotlibrc'))

## Source color families, as in provenance_plot.py: light for all, dark for recoverable
COLORS = {'NEA': ('#8db8ee', '#133f7a'), 'lit': ('#f5a37c', '#a8430f')}
SEED = 1  # for the example hosts, so rebuilds give the same figure
MIN_COMPANIONS = 15
X_CAP = 70

NAMES = {
    'CLS I (Rosenthal+2021)': 'CLS I (Rosenthal+2021)',
    'Kiefer+2019': 'Kiefer+2019 (SOPHIE)',
    'sahlmann_2011': 'Sahlmann+2011 (CORALIE)',
    'unger_2023': 'Unger+2023 (CORALIE, HARPS)',
    'barbato_2023': 'Barbato+2023 (CORALIE XIX)',
    'barbato_2026': 'Barbato+2026 (GAPS LXXV)',
    'wilson_2016': 'Wilson+2016 (SOPHIE)',
    'diaz_2012': 'Díaz+2012 (SOPHIE)',
    'dalal_2021': 'Dalal+2021 (SOPHIE)',
    'nidever_2002': 'Nidever+2002 (Keck, Lick)',
    'parc_2026': 'Parc+2026 (CORALIE)',
    'frensch_2023': 'Frensch+2023 (HARPS)',
    'jenkins_2015': 'Jenkins+2015 (AAPS)',
    'okayama_eapsnet': 'Omiya+2009 (Okayama)',
    'ptps': 'PTPS / TAPAS',
    'teng_2023': 'Teng+2023 (Okayama)',
    'wittenmyer_2009': 'Wittenmyer+2009 (HET, McDonald)',
    'xiao_2023': 'Xiao+2023',
    'ren_2026': 'Ren+2026',
    'cui_2026': 'Cui+2026',
    'troup_2016': 'Troup+2016 (APOGEE)',
    'grieves_2017': 'Grieves+2017 (MARVELS)',
    'marvels_bd': 'MARVELS discovery papers',
    'exoplanet_eu': 'exoplanet.eu (OEC mirror)',
    'stevenson_2023': 'Stevenson+2023',
    'ma_ge_2014': 'Ma & Ge 2014',
    'kiefer_2021': 'Kiefer+2021',
    'sb9': 'SB9',
    'latham_2002': 'Latham+2002',
    'sb9_tic': 'SB9 (TIC primary masses)',
}

## (tier in the census, heading), in precedence order
TIERS = [
    ('survey paper', 'Survey papers (full RV orbits)'),
    ('RV+astrometry paper', 'RV+astrometry papers'),
    ('low-precision survey', 'Low-precision surveys (50-200 m/s)'),
    ('compilation', 'Compilations'),
    ('mass-function M sin i', 'SB1 catalogs (M sin i from mass function)'),
]

census = pd.read_csv(os.path.join(HERE, 'census_table.csv'), dtype={'source_id': str})
from_nea = census.catalog == 'NEA'
census['family'] = np.where(from_nea, 'NEA', 'lit')
census['study'] = census.pl_refname
## Literature tiers, as set in literature.py (CLS I and Kiefer+2019 are survey papers)
census['tier'] = census.catalog_status.where(
    ~census.catalog.isin(['CLS I (Rosenthal+2021)', 'Kiefer+2019']), 'survey paper'
).where(~from_nea, 'NEA')

## NEA display names from the archive's reference links ('Stassun et al. 2017')
refs = pd.read_csv(NEA_PS_FILE).pl_refname.str.extract(r'refstr=(\S+).*?>(.*?)</a>')
refs = refs.dropna().drop_duplicates(0).set_index(0)[1]
nea_names = refs.map(lambda t: html.unescape(t).strip().replace(' et al. ', '+'))
names = {**NAMES, **nea_names.to_dict()}

lit = census[~from_nea]
unnamed = set(census.study) - set(names)
assert not unnamed, f'studies without a display name: {unnamed}'
untiered = set(lit.tier) - {t for t, _ in TIERS}
assert not untiered, f'unexpected tiers: {untiered}'

counts = census.groupby(['family', 'tier', 'study']).agg(
    n=('hostname', 'size'), rec=('recoverable', 'sum')).reset_index()
counts['rec'] = counts.rec.astype(int)

## Up to two example hosts per bar; the 'all' examples avoid the recoverable ones when
## the study has enough other hosts
rng = np.random.RandomState(SEED)
examples = {}
for study, g in census.groupby('study'):
    rec_hosts = sorted(set(g.hostname[g.recoverable]))
    shown_rec = list(rng.choice(rec_hosts, min(2, len(rec_hosts)), replace=False))
    others = sorted(set(g.hostname) - set(shown_rec))
    pool = others if len(others) >= 2 else sorted(set(g.hostname))
    shown_all = list(rng.choice(pool, min(2, len(pool)), replace=False))
    examples[study] = (shown_all, shown_rec)

## Rows: ('heading', text) or ('row', family, label, n, rec, examples)
nea = counts[counts.family == 'NEA']
shown = (nea.n >= MIN_COMPANIONS) | (nea.rec > 0)
rest = nea[~shown]
nea_rows = [('row', 'NEA', names[r.study], r.n, r.rec, examples[r.study])
            for r in nea[shown].sort_values(['n', 'rec'], ascending=False).itertuples()]
nea_rows.append(('row', 'NEA', f'Other NEA studies ({len(rest)})', rest.n.sum(),
                 rest.rec.sum(), ([], [])))
nea_heading = (f'NASA Exoplanet Archive: {nea.n.sum()} companions, {nea.rec.sum()} '
               f'recoverable, from {len(nea)} studies')
half = (len(nea_rows) + 1)//2
columns = [[('heading', nea_heading)] + nea_rows[:half],
           [('heading', 'NASA Exoplanet Archive (continued)')] + nea_rows[half:], []]
for tier, heading in TIERS:
    t = counts[(counts.family == 'lit') & (counts.tier == tier)]
    columns[2].append(('heading', f'{heading}: {t.n.sum()} companions, '
                                  f'{t.rec.sum()} recoverable'))
    columns[2] += [('row', 'lit', names[r.study], r.n, r.rec, examples[r.study])
                   for r in t.sort_values(['n', 'rec'], ascending=False).itertuples()]

## Row positions: headings get a little extra space above them
GAP = 0.4
positions = []
for col in columns:
    y, pos = 0, []
    for i, item in enumerate(col):
        if item[0] == 'heading' and i > 0:
            y += GAP
        pos.append(y)
        y += 1
    positions.append(pos)
y_max = max(p[-1] for p in positions)

FONTSIZE = 22
EX_FONTSIZE = FONTSIZE - 5
h = 0.36
x_max = 1.45*X_CAP
fig, axes = plt.subplots(1, 3, figsize=(54, 34), gridspec_kw=dict(wspace=0.62))
fig.subplots_adjust(top=0.92)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()


def annotate(ax, n, y, hosts):
    """Count at the bar's end, then '(host, host)' after it, with as many hosts as fit
    before the axes' right edge (none if even one doesn't fit)."""
    count = ax.text(min(n, X_CAP) + 0.8, y, str(n), va='center', ha='left',
                    fontsize=FONTSIZE - 3)
    if not hosts:
        return
    ax_right = ax.get_window_extent(renderer).x1
    start = count.get_window_extent(renderer).x1 + 10
    for k in range(len(hosts), 0, -1):
        text = ax.text(0, y, '(' + ', '.join(hosts[:k]) + ')', fontsize=EX_FONTSIZE,
                       va='center', ha='left', color='0.25')
        if start + text.get_window_extent(renderer).width <= ax_right - 10:
            text.set_x(ax.transData.inverted().transform((start, 0))[0])
            return
        text.remove()


def bar(ax, y, n, color):
    """A bar of length n, cut at X_CAP with a hatched end when longer."""
    ax.barh(y, min(n, X_CAP), h, color=color, edgecolor='white', linewidth=1.5, zorder=3)
    if n > X_CAP:
        ax.barh(y, 6, h, left=X_CAP - 6, color=color, hatch='//', edgecolor='white',
                linewidth=0, zorder=4)


for ax, col, pos in zip(axes, columns, positions):
    ax.set_xlim(0, x_max)
    ax.set_ylim(y_max + 0.7, -0.8)
    rows = [(y, item) for y, item in zip(pos, col) if item[0] == 'row']
    for y, (_, family, label, n, rec, (ex_all, ex_rec)) in rows:
        c_all, c_rec = COLORS[family]
        bar(ax, y - h/2, n, c_all)
        bar(ax, y + h/2, rec, c_rec)
        annotate(ax, n, y - h/2, ex_all)
        annotate(ax, rec, y + h/2, ex_rec)
    for y, item in zip(pos, col):
        if item[0] == 'heading':
            ## Headings start at the left edge of the row labels
            ax.annotate(item[1], (0, y + 0.1), xycoords=('axes fraction', 'data'),
                        xytext=(-0.62, 0), textcoords=('axes fraction', 'offset points'),
                        va='center', ha='left', fontsize=FONTSIZE + 1, fontweight='bold')
            if y > 0:
                ax.axhline(y - 0.7, xmin=-0.62, color='0.6', lw=1.5, zorder=1,
                           clip_on=False)
    ax.set_yticks([y for y, _ in rows])
    ax.set_yticklabels([item[2] for _, item in rows], fontsize=FONTSIZE)
    ## Headings run across the label column into the plot, so only the x axis is drawn
    for side in ['left', 'right', 'top']:
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='y', which='both', left=False, right=False)
    ax.tick_params(axis='x', which='both', top=False)
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xticks(np.arange(0, X_CAP + 1, 10))
    ax.tick_params(axis='x', labelsize=FONTSIZE)
    ax.set_xlabel('Number of census companions', fontsize=FONTSIZE + 2)

fig.legend(handles=[
    matplotlib.patches.Patch(color=COLORS['NEA'][0], label='NEA study: all companions'),
    matplotlib.patches.Patch(color=COLORS['NEA'][1], label='NEA study: plausibly recoverable'),
    matplotlib.patches.Patch(color=COLORS['lit'][0], label='Literature source: all companions'),
    matplotlib.patches.Patch(color=COLORS['lit'][1],
                             label='Literature source: plausibly recoverable'),
    matplotlib.patches.Patch(facecolor='0.6', hatch='//', edgecolor='white',
                             label=f'Bar cut at {X_CAP} (full count printed)'),
], loc='upper center', bbox_to_anchor=(0.5, 0.955), ncol=5, fontsize=FONTSIZE + 2,
   frameon=False, handlelength=1.5, columnspacing=2.5)
fig.suptitle(
    f'Studies supplying the adopted fit for each census companion ({len(census)} '
    f'companions, {int(census.recoverable.sum())} plausibly recoverable: S/N $\\geq$ '
    f'{SNR_THRESHOLD} and $P \\leq$ {P_MAX_YR} yr)\n'
    f'Example hosts in parentheses; NEA studies with < {MIN_COMPANIONS} companions and '
    'none recoverable are summed in "Other NEA studies"',
    y=0.99, fontsize=32)
fig.savefig(os.path.join(HERE, 'lit_sources.png'), dpi=50, bbox_inches='tight')

print(counts.sort_values(['family', 'n'], ascending=[True, False]).to_string(index=False))
