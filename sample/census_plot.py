"""Plot the census against the expected Gaia DR4 sensitivity.

Figures:
    alpha_vs_period.png              minimum astrometric signature vs. period, with
                                     DR4 detection thresholds for several G
    snr_vs_period_by_distance.png    forecast S/N vs. period in three distance bins,
                                     colored by absolute G

Usage (from the project root, after census.py):
    python -I -W ignore sample/census_plot.py
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census import P_ACCURATE_YR, P_MAX_YR, SNR_THRESHOLD, sigma_fov_gaiamock

matplotlib.rc_file(os.path.join(HERE, '..', 'matplotlibrc'))

census = pd.read_csv(os.path.join(HERE, 'census_table.csv'), dtype={'source_id': str})
census['period_yr'] = census.period_d/365.25
in_gaia = census.G_source == 'Gaia DR3'
rec = census.recoverable

## Figure 1: alpha_min vs. period with thresholds
cmap = plt.get_cmap('viridis')
norm = matplotlib.colors.Normalize(vmin=4, vmax=16)

Handle = matplotlib.lines.Line2D

fig, ax = plt.subplots(figsize=(14, 11))
ax.axvspan(P_MAX_YR, 1e3, color='0.88', zorder=0)
ax.axvline(P_ACCURATE_YR, color='0.4', ls=':', lw=2, zorder=1)
ax.scatter(census.period_yr[~rec & in_gaia], census.alpha_min_uas[~rec & in_gaia],
           c=census.G[~rec & in_gaia], cmap=cmap, norm=norm, s=25, alpha=0.6,
           edgecolor='none', zorder=2)
sc = ax.scatter(census.period_yr[rec], census.alpha_min_uas[rec], c=census.G[rec],
                cmap=cmap, norm=norm, s=90, edgecolor='k', linewidth=1.2, zorder=4)
ax.scatter(census.period_yr[~in_gaia], census.alpha_min_uas[~in_gaia], marker='x',
           color='r', s=60, zorder=3)

## Legend 1: thresholds, which stop at the end of the DR4 time span
thresh_handles = []
for G in [10, 15, 5, 17]:
    thresh = SNR_THRESHOLD*sigma_fov_gaiamock(G)
    ax.plot([1e-3, P_MAX_YR], [thresh, thresh], color=cmap(norm(G)), lw=3, zorder=5)
    thresh_handles.append(Handle([], [], color=cmap(norm(G)), lw=3,
                                 label=fr'$G = {G}$: {thresh:.0f} $\mu$as'))
leg1 = ax.legend(handles=thresh_handles, loc='upper left', fontsize=16, frameon=False,
                 title=r'Detection threshold, $\alpha = 1.5\,\sigma_\mathrm{FOV}(G)$',
                 title_fontsize=16)
ax.add_artist(leg1)

## Legend 2: markers and period markings
marker_handles = [
    Handle([], [], marker='o', ls='', ms=10, mfc='0.6', mec='k', mew=1.2,
           label=f'Plausibly recoverable ({rec.sum()}): above its\n'
                 r'own threshold and $P \leq 5.5$ yr'),
    Handle([], [], marker='o', ls='', ms=6, mfc='0.6', mec='none', alpha=0.6,
           label=f'Below threshold or $P > 5.5$ yr ({(~rec).sum()})'),
    Handle([], [], marker='x', ls='', ms=6, mew=1.5, color='r',
           label=f'No trusted Gaia DR3 photometry ({(~in_gaia).sum()}):\n'
                 'G from NEA G or, if absent, NEA V'),
    matplotlib.patches.Patch(color='0.88', label='$P > 5.5$ yr: longer than the DR4\n'
                                                 'time span'),
    Handle([], [], color='0.4', ls=':', lw=2, label='4 yr: DR4 periods accurate to\n'
                                                     '20% (Lammers & Winn 2025)'),
]
leg2 = ax.legend(handles=marker_handles, loc='lower right', fontsize=15, frameon=True,
                 framealpha=1, facecolor='white', edgecolor='0.5')
leg2.set_zorder(10)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 20)
ax.set_ylim(1e-3, 1e5)
ax.set_xlabel('Orbital period [yr]')
ax.set_ylabel(r'$\alpha_\mathrm{min}$ [$\mu$as] (from $M\,\sin i$, else true mass)')
ax.set_title('Astrometric signals of RV companions (NEA + RV survey literature)\n'
             'vs. expected Gaia DR4 sensitivity',
             pad=12)
cb = fig.colorbar(sc, ax=ax)
cb.set_label(r'Host $G$ [mag]')
fig.savefig(os.path.join(HERE, 'alpha_vs_period.png'), dpi=150, bbox_inches='tight')

## Figure 2: S/N vs. period in distance bins, colored by absolute G
dist_bins = [(0, 30), (30, 100), (100, np.inf)]
norm_abs = matplotlib.colors.Normalize(vmin=-1, vmax=12)
fig, axes = plt.subplots(1, 3, figsize=(30, 10), sharey=True)
for ax, (d_lo, d_hi) in zip(axes, dist_bins):
    s = census[(census.dist_pc >= d_lo) & (census.dist_pc < d_hi)]
    box = matplotlib.patches.Rectangle((1e-3, SNR_THRESHOLD), P_MAX_YR - 1e-3, 1e4,
                                       color='C2', alpha=0.12, zorder=0)
    ax.add_patch(box)
    bright = s.G < 6
    sc = ax.scatter(s.period_yr[~bright], s.snr[~bright], c=s.abs_G[~bright],
                    cmap='plasma', norm=norm_abs, s=40, edgecolor='k', linewidth=0.3)
    ax.scatter(s.period_yr[bright], s.snr[bright], c=s.abs_G[bright], cmap='plasma',
               norm=norm_abs, s=110, marker='*', edgecolor='k', linewidth=0.5)
    ax.axhline(SNR_THRESHOLD, color='k', ls='--')
    ax.axvline(P_MAX_YR, color='k', ls='--')
    ax.axvline(P_ACCURATE_YR, color='0.4', ls=':', lw=2)
    title = f'{d_lo}-{d_hi} pc' if np.isfinite(d_hi) else f'> {d_lo} pc'
    ax.set_title(f'{title}  ({len(s)} companions, {int(s.recoverable.sum())} recoverable)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-3, 20)
    ax.set_ylim(1e-5, 1e2)
    ax.set_xlabel('Orbital period [yr]')
axes[0].set_ylabel(r'Forecast S/N per transit, $\alpha_\mathrm{min}\,/\,\sigma_\mathrm{FOV}(G)$')
leg = axes[2].legend(handles=[
    Handle([], [], marker='o', ls='', ms=8, mfc='0.6', mec='k', mew=0.5,
           label='Host $G \\geq 6$'),
    Handle([], [], marker='*', ls='', ms=15, mfc='0.6', mec='k', mew=0.5,
           label='Host $G < 6$ (higher noise)'),
    matplotlib.patches.Patch(color='C2', alpha=0.2,
                             label='Plausibly recoverable:\nS/N $\\geq$ 1.5 and $P \\leq$ 5.5 yr'),
    Handle([], [], color='k', ls='--', label='S/N = 1.5 and $P$ = 5.5 yr'),
    Handle([], [], color='0.4', ls=':', lw=2, label='$P$ = 4 yr (accurate periods)'),
], loc='lower right', fontsize=17, frameon=True, framealpha=1, facecolor='white',
   edgecolor='0.5')
leg.set_zorder(10)
fig.suptitle('Forecast Gaia DR4 signal-to-noise of RV companions (NEA + RV survey '
             'literature), by host distance',
             y=1.02, fontsize=30)
cb = fig.colorbar(sc, ax=axes, pad=0.01)
cb.set_label(r'Host absolute magnitude $M_G$')
fig.savefig(os.path.join(HERE, 'snr_vs_period_by_distance.png'), dpi=120,
            bbox_inches='tight')

## Numbers quoted in CENSUS.md
best_thresh = SNR_THRESHOLD*sigma_fov_gaiamock(10)
easy = (census.alpha_min_uas >= best_thresh) & (census.period_yr <= P_MAX_YR)
print(f'alpha_min above the best-case threshold ({best_thresh:.0f} uas) and P <= 5.5 yr: '
      f'{easy.sum()}; recoverable: {(easy & rec).sum()}')
missed = census[easy & ~rec]
## The gaiamock error curve is elevated for 6 < G < 8.4 (up to 0.23 mas per CCD at G = 8.2,
## vs. 0.13 mas at G = 8.4-13)
reason = np.select([~in_gaia[missed.index], missed.G < 6, missed.G < 8.4, missed.G > 13],
                   ['not in Gaia DR3', 'bright (G < 6)', 'G = 6-8.4 noise bump',
                    'faint (G > 13)'], 'other')
print('  not recoverable because:', pd.Series(reason).value_counts().to_dict())
print(missed.assign(reason=reason)[['pl_name', 'G', 'dist_pc', 'alpha_min_uas', 'snr',
                                    'reason']].to_string(index=False))
print('\nRecoverable:', rec.sum(), 'planets,', census[rec].hostname.nunique(), 'hosts;',
      'discovery', census[rec].discoverymethod.value_counts().to_dict())
print('  distance percentiles (10/50/90):',
      census[rec].dist_pc.quantile([0.1, 0.5, 0.9]).round(0).tolist())
print('  multi-planet hosts:', census[rec & (census.sy_pnum > 1)].hostname.nunique())
print('  S/N 0.5-1.5 with P <= 5.5 yr (threshold region):',
      ((census.snr >= 0.5) & (census.snr < 1.5) & (census.period_yr <= P_MAX_YR)).sum())
print('  S/N >= 1.5 with 5.5 yr < P <= 2 T_DR4:',
      ((census.snr >= 1.5) & (census.period_yr > P_MAX_YR)
       & (census.p_over_T <= 2)).sum())
