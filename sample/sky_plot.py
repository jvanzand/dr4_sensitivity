"""Sky positions of the census hosts, to check the census for gaps on the sky.

Figure:
    sky_positions.png    RA/Dec (Mollweide) of every census host, NEA vs. literature,
                         with hosts of plausibly recoverable companions outlined

Coordinates are Gaia DR3 (epoch 2016.0) where the host has a DR3 match, else the NEA's.

Usage (from the project root, after census.py):
    python -I -W ignore sample/sky_plot.py
"""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
matplotlib.rc_file(os.path.join(HERE, '..', 'matplotlibrc'))

## Source colors as in provenance_plot.py (NEA blue, literature orange)
COLOR_NEA, COLOR_LIT = '#1f66c2', '#d4561f'

census = pd.read_csv(os.path.join(HERE, 'census_table.csv'), dtype={'source_id': str})
gaia = pd.concat([pd.read_csv(os.path.join(HERE, 'raw', f), dtype={'source_id': str})
                  for f in ['gaia_dr3_hosts.csv', 'gaia_dr3_cone_no_id.csv']])
gaia = gaia.drop_duplicates('source_id')[['source_id', 'ra', 'dec']]
nea = pd.read_csv(os.path.join(HERE, 'raw', 'nea_pscomppars_rv.csv'))
nea = nea.drop_duplicates('hostname')[['hostname', 'ra', 'dec']]

hosts = (census.groupby('hostname', sort=False)
         .agg(catalog=('catalog', 'first'), source_id=('source_id', 'first'),
              recoverable=('recoverable', 'any'))
         .reset_index())
hosts = hosts.merge(gaia, on='source_id', how='left').merge(
    nea.rename(columns={'ra': 'nea_ra', 'dec': 'nea_dec'}), on='hostname', how='left')
from_gaia = hosts.ra.notna()
hosts['ra'] = hosts.ra.fillna(hosts.nea_ra)
hosts['dec'] = hosts.dec.fillna(hosts.nea_dec)
missing = hosts.ra.isna() | hosts.dec.isna()
assert not missing.any(), f'hosts without coordinates: {list(hosts.hostname[missing])}'


def to_xy(ra_deg, dec_deg):
    """Mollweide axes coordinates [rad], with RA increasing to the left and 0h at
    the center."""
    ra = np.remainder(np.asarray(ra_deg, float) + 180, 360) - 180
    return -np.radians(ra), np.radians(np.asarray(dec_deg, float))


fig = plt.figure(figsize=(20, 11.5))
ax = fig.add_subplot(111, projection='mollweide')
ax.grid(True, color='0.8', lw=1)

## Galactic plane and the ecliptic, each split where it wraps around RA = 12h
lines = [(SkyCoord(l=np.linspace(0, 360, 721)*u.deg, b=0*u.deg, frame='galactic').icrs,
          dict(color='0.35', ls='--', lw=2), 'Galactic plane'),
         (SkyCoord(lon=np.linspace(0, 360, 721)*u.deg, lat=0*u.deg,
                   frame='barycentrictrueecliptic').icrs,
          dict(color='0.35', ls=':', lw=2), 'Ecliptic')]
handles = []
for coord, style, label in lines:
    x, y = to_xy(coord.ra.deg, coord.dec.deg)
    breaks = np.where(np.abs(np.diff(x)) > np.pi)[0] + 1
    for xs, ys in zip(np.split(x, breaks), np.split(y, breaks)):
        ax.plot(xs, ys, zorder=1, **style)
    handles.append(matplotlib.lines.Line2D([], [], label=label, **style))

lit = hosts.catalog != 'NEA'
for sel, color, label in [(~lit, COLOR_NEA, 'NASA Exoplanet Archive'),
                          (lit, COLOR_LIT, 'Literature (papers, surveys, catalogs)')]:
    x, y = to_xy(hosts.ra[sel], hosts.dec[sel])
    ax.scatter(x, y, s=14, color=color, alpha=0.7, edgecolor='none', zorder=2)
    handles.insert(len(handles) - 2, matplotlib.lines.Line2D(
        [], [], marker='o', ls='', ms=7, color=color, label=f'{label} host ({sel.sum()})'))
rec = hosts.recoverable
x, y = to_xy(hosts.ra[rec], hosts.dec[rec])
ax.scatter(x, y, s=70, facecolor='none', edgecolor='k', linewidth=1.3, zorder=3)
handles.insert(2, matplotlib.lines.Line2D(
    [], [], marker='o', ls='', ms=10, mfc='none', mec='k', mew=1.3,
    label=f'Host of a plausibly recoverable companion ({rec.sum()})'))

## The Kepler field holds the densest clump of hosts
kx, ky = to_xy(290.7, 44.5)
ax.annotate('Kepler field', (kx, ky), xytext=(kx + 0.55, ky + 0.35), fontsize=18,
            arrowprops=dict(arrowstyle='->', color='k', lw=1.5), ha='center', zorder=6,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='0.5',
                      alpha=0.85))

## RA labels as text above the points (tick labels would sit beneath them)
ax.set_xticks(np.radians(np.arange(-150, 151, 30)))
ax.set_xticklabels([])
for ra_hr in range(0, 24, 2):
    if ra_hr == 12:
        continue  # at the edge
    x, _ = to_xy(15*ra_hr, 0)
    ax.text(x, np.radians(3), f'{ra_hr}h', ha='center', va='bottom', fontsize=16,
            zorder=5, bbox=dict(facecolor='white', edgecolor='none', pad=1, alpha=0.8))
ax.tick_params(axis='y', labelsize=16)
ax.set_xlabel('Right ascension (increasing to the left)', fontsize=20, labelpad=12)
ax.set_ylabel('Declination', fontsize=20)
ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2,
          fontsize=17, frameon=False)
ax.set_title(f'Sky positions of census hosts ({len(hosts)} hosts, {len(census)} '
             'companions), ICRS', pad=24, fontsize=24)
fig.savefig(os.path.join(HERE, 'sky_positions.png'), dpi=110, bbox_inches='tight')

## Host counts by declination and RA band, to spot gaps
print('Hosts per 15-deg declination band:')
print(pd.crosstab(pd.cut(hosts.dec, np.arange(-90, 91, 15)),
                  lit.map({False: 'NEA', True: 'literature'})).to_string())
print('\nHosts per 2h RA band:')
print(pd.crosstab(pd.cut(hosts.ra, np.arange(0, 361, 30)),
                  lit.map({False: 'NEA', True: 'literature'})).to_string())
print(f'\nCoordinates: Gaia DR3 for {from_gaia.sum()} hosts, NEA for {(~from_gaia).sum()}')
