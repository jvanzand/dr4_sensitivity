"""Where each host and companion parameter in the census comes from.

One stacked bar per parameter; each segment is a provenance, labeled directly with its
count. Host parameters count hosts; companion parameters count companions. Colors mark
the source in every bar (NEA blue, Gaia DR3 green, SIMBAD/2MASS purple, literature
orange, missing gray), and segments stack in that order from the bottom.

Usage (from the project root, after census.py):
    python -I -W ignore sample/provenance_plot.py
"""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
matplotlib.rc_file(os.path.join(HERE, '..', 'matplotlibrc'))

## One color family per source, dark to light (dataviz reference hues); the family order
## is also the stacking order
FAMILIES = {
    'NEA': ('NASA Exoplanet Archive', ['#133f7a', '#1f66c2', '#4f8fe0', '#8db8ee',
                                       '#c4dbf7']),
    'Gaia': ('Gaia DR3', ['#0f7a54', '#1baf7a', '#7fd8b6']),
    'SIMBAD': ('SIMBAD / 2MASS', ['#6250d6']),
    'lit': ('Literature (papers, surveys, catalogs)',
            ['#8f370b', '#d4561f', '#ef8f4f', '#eda100', '#f5cd73']),
    'none': ('Missing', ['#b3b3b3']),
}

census = pd.read_csv(os.path.join(HERE, 'census_table.csv'), dtype={'source_id': str})
hosts = census.drop_duplicates('hostname')
from_nea = census.catalog == 'NEA'

## Literature tiers, as set in literature.py (CLS I and Kiefer+2019 are survey papers)
tier = census.catalog_status.where(~from_nea)
tier = tier.where(~census.catalog.isin(['CLS I (Rosenthal+2021)', 'Kiefer+2019']),
                  'survey paper')
orbit_source = np.select(
    [from_nea & (census.fit_type == 'RV-only'),
     from_nea & (census.fit_type == 'RV+astrometry'),
     from_nea & (census.fit_type == 'transit+RV'),
     tier == 'survey paper', tier == 'RV+astrometry paper', tier == 'low-precision survey',
     tier == 'compilation', tier == 'mass-function M sin i'],
    ['NEA RV-only', 'NEA RV+astrometry', 'NEA transit + RV', 'Survey paper', 'RV+astrometry paper',
     'Low-precision survey (APOGEE, MARVELS)', 'Compilation',
     'SB1 catalogs (M sin i from mass function)'], 'other')

## Tp segments are split by source as well as by time system
tp_system = np.select(
    [census.tp_flag == 'no Tp',
     census.tp_flag == 'BJD_TDB as published',
     census.tp_flag.str.startswith('BJD with unstated'),
     census.tp_flag.str.startswith('converted from'),
     census.tp_flag.str.startswith('time system not stated')],
    ['no Tp', 'BJD_TDB as published', 'BJD, time scale unstated',
     'converted from JD, HJD, or BJD_UTC', 'time system unstated'], 'other')
tp_source = from_nea.map({True: 'NEA: ', False: 'Literature: '}) + tp_system
tp_source = tp_source.where(tp_system != 'no Tp', 'No Tp given')

gaia_match = hosts.gaia_match.replace({
    'archive ID': 'NEA Gaia DR3 ID', 'SIMBAD': 'SIMBAD / 2MASS crossmatch',
    'position': 'Gaia DR3 cone search', 'position (no DR3 parallax)': 'Gaia DR3 cone search',
    'not in DR3': 'Absent from DR3'}).where(~hosts.gaia_match.str.contains('suspect'),
                                            'NEA ID flagged suspect')

## The host mass called 'orbit source' comes from the row or paper giving the orbit
st_mass_source = hosts.st_mass_source.replace({
    'NEA composite': 'NEA composite table', 'Gaia FLAME': 'Gaia DR3 FLAME'})
st_mass_source = st_mass_source.where(
    st_mass_source != 'orbit source',
    np.where(hosts.catalog == 'NEA', 'NEA orbit solution', 'Literature orbit paper'))

TP_SYSTEMS = ['BJD_TDB as published', 'BJD, time scale unstated',
              'converted from JD, HJD, or BJD_UTC', 'time system unstated']

## (bar label, values, [(provenance, source family), ...] from bottom to top); a
## segment's shade is its rank within its family in that bar
bars = [
    ('Host:\nGaia DR3 match', gaia_match,
     [('NEA Gaia DR3 ID', 'NEA'), ('NEA ID flagged suspect', 'NEA'),
      ('Gaia DR3 cone search', 'Gaia'), ('SIMBAD / 2MASS crossmatch', 'SIMBAD'),
      ('Absent from DR3', 'none')]),
    ('Host:\nG magnitude', hosts.G_source.replace({
        'Gaia DR3': 'Gaia DR3 G', 'archive G': 'NEA G', 'archive V': 'NEA V (stand-in)'}),
     [('NEA G', 'NEA'), ('NEA V (stand-in)', 'NEA'), ('Gaia DR3 G', 'Gaia')]),
    ('Host:\ndistance', hosts.dist_source.replace({
        'archive parallax': 'NEA parallax', 'archive distance': 'NEA distance'}),
     [('NEA parallax', 'NEA'), ('NEA distance', 'NEA'), ('Gaia DR3 parallax', 'Gaia')]),
    ('Host:\nmass', st_mass_source,
     [('NEA orbit solution', 'NEA'), ('NEA composite table', 'NEA'),
      ('Gaia DR3 FLAME', 'Gaia'), ('Literature orbit paper', 'lit')]),
    ('Companion:\norbit and mass', pd.Series(orbit_source),
     [('NEA RV-only', 'NEA'), ('NEA RV+astrometry', 'NEA'), ('NEA transit + RV', 'NEA'),
      ('Survey paper', 'lit'), ('RV+astrometry paper', 'lit'),
      ('Low-precision survey (APOGEE, MARVELS)', 'lit'), ('Compilation', 'lit'),
      ('SB1 catalogs (M sin i from mass function)', 'lit')]),
    ('Companion:\nwhich solution', pd.Series(np.where(
        from_nea, 'NEA: ' + census.solution_rule.astype(str), 'Literature (precedence by tier)')),
     [('NEA: newest paper for the system', 'NEA'), ('NEA: newest paper for the planet', 'NEA'),
      ('NEA: more precise P (>2x)', 'NEA'), ('Literature (precedence by tier)', 'lit')]),
    ('Companion:\nperiastron time', pd.Series(tp_source),
     [(f'NEA: {s}', 'NEA') for s in TP_SYSTEMS]
     + [(f'Literature: {s}', 'lit') for s in TP_SYSTEMS] + [('No Tp given', 'none')]),
]

fig, ax = plt.subplots(figsize=(46, 22))
width = 0.18
LABEL_DX = 0.09     # gap between a bar's right edge and its labels [data units]
LABEL_MARGIN = 0.03  # clearance kept before the next bar [data units]
FONTSIZE = 26
x_right = len(bars) - 1 + 0.9
ax.set_xlim(-0.35, x_right)
ax.set_ylim(0, 1.2*max(len(census), len(hosts)))
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
px_per_x = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]


def wrap_to_width(text, max_px):
    """Break text at spaces so its rendered width fits max_px."""
    probe = ax.text(0, 0, text, fontsize=FONTSIZE)
    fits = lambda t: (probe.set_text(t) or probe.get_window_extent(renderer).width <= max_px)
    try:
        if fits(text):
            return text
        words = text.split(' ')
        for n_chars in range(len(text), 0, -1):
            lines, line = [], ''
            for word in words:
                trial = f'{line} {word}'.strip()
                if line and len(trial) > n_chars:
                    lines.append(line)
                    line = word
                else:
                    line = trial
            lines.append(line)
            wrapped = '\n'.join(lines)
            if fits(wrapped):
                return wrapped
        return text
    finally:
        probe.remove()


for i, (label, values, order) in enumerate(bars):
    counts = values.value_counts()
    unexpected = set(counts.index) - {name for name, _ in order}
    assert not unexpected, f'{label}: unlabeled provenance {unexpected}'
    total = counts.sum()
    bottom = 0
    segments = []
    rank = {}
    for name, family in order:
        color = FAMILIES[family][1][rank.get(family, 0)]
        rank[family] = rank.get(family, 0) + 1
        n = int(counts.get(name, 0))
        if n == 0:
            continue
        ax.bar(i, n, width, bottom=bottom, color=color, edgecolor='white',
               linewidth=1.5, zorder=3)
        ## A no-break space keeps the count on the same line as the last word
        segments.append((bottom + n/2, f'{name}: {n}', color))
        bottom += n
    ## Room for labels: from just right of this bar to just left of the next one
    label_x = i + width/2 + LABEL_DX
    next_edge = (i + 1 - width/2 - LABEL_MARGIN) if i < len(bars) - 1 else x_right - 0.02
    max_px = (next_edge - label_x - 0.03)*px_per_x
    ## Labels in a column, spread evenly over the bar's height so small segments don't
    ## collide, each joined to its segment by a thin leader line
    label_y = (np.linspace(0.04, 0.96, len(segments))*total if len(segments) > 1
               else [segments[0][0]])
    for (seg_y, text, color), y in zip(segments, label_y):
        ax.plot([i + width/2, label_x - 0.045], [seg_y, y], color='0.5', lw=1, zorder=2)
        ax.plot(label_x - 0.035, y, 's', color=color, ms=18, zorder=4)
        ax.text(label_x, y, wrap_to_width(text, max_px), va='center', ha='left',
                fontsize=FONTSIZE, linespacing=1.05)
    unit = 'hosts' if label.startswith('Host') else 'companions'
    ax.text(i, total*1.015, f'{total} {unit}', ha='center', va='bottom', fontsize=FONTSIZE)

ax.legend(handles=[matplotlib.patches.Patch(color=shades[0], label=name)
                   for name, shades in FAMILIES.values()],
          title='Source (shades distinguish provenances from the same source)',
          loc='upper left', ncol=len(FAMILIES), fontsize=FONTSIZE, title_fontsize=FONTSIZE,
          frameon=False, handlelength=1.2, columnspacing=1.5)
ax.set_xticks(range(len(bars)))
ax.set_xticklabels([b[0] for b in bars], fontsize=32)
ax.tick_params(axis='x', which='both', top=False, bottom=False)
ax.tick_params(axis='y', labelsize=28)
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_ylabel('Number of hosts or companions', fontsize=32)
ax.set_title('Provenance of host and companion parameters in the census '
             f'({len(census)} companions around {len(hosts)} hosts)', pad=16, fontsize=36)
fig.savefig(os.path.join(HERE, 'provenance.png'), dpi=100, bbox_inches='tight')

for label, values, order in bars:
    print(label.replace('\n', ' '), values.value_counts().to_dict())
