"""Data channels for the blind search, and the choice of which ones to use.

A channel is one dataset: along-scan astrometry (Gaia DR4 transits; Hipparcos
IAD later) or the RVs of one instrument. Names: 'gaia_dr4', 'hip_iad',
'rv:<dataset>'. `select(channels, spec)` keeps a subset, e.g. 'gaia_dr4',
'gaia_dr4,rv', 'rv:hires_j,rv:harps', 'astro', 'all'.

Times are MJD (days) throughout; project RV tables carry BJD_TCB, converted on
load.

Test mode:  python -I blind_search/data.py --test
"""

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BJD_MJD = 2400000.5


@dataclass
class Channel:
    kind: str                       # 'astro' or 'rv'
    name: str
    t: np.ndarray                   # MJD
    y: np.ndarray                   # mas (astro, along-scan abscissa) or m/s (rv)
    err: np.ndarray                 # formal errors, same units
    sinpsi: np.ndarray = None       # astro: sin and cos of the scan position angle
    cospsi: np.ndarray = None
    fal: np.ndarray = None          # astro: along-scan parallax factor

    def __len__(self):
        return len(self.t)


def load_gaia_dr4(path, name='gaia_dr4', reduction='tutorial', seed=0, nboot=256):
    """Gaia DR4 epoch astrometry from the archive-format CSV (one row per CCD
    crossing), reduced to one row per transit. Rows without used_by_agis_al, or
    with outlier_flag != 0, are dropped.

    reduction='tutorial' (the PLAN.md default, decided 2026-11-06): position =
    median of the CCD centroids; error = standard deviation of `nboot` bootstrap
    medians, floored at median(ipd_error_al)/sqrt(n); a one-CCD transit takes
    its ipd_error_al; scan angle and parallax factor from the first CCD; epoch =
    mean CCD epoch (an assumption: the PLAN.md description does not say). `seed`
    makes the bootstrap reproducible (the tutorial seeds its RNG with the source
    id; numpy's generator differs, so errors agree only statistically).

    reduction='demo' (conservative sensitivity test), as in Judah's demo fits:
    after sorting the CCDs by formal error, take the median row (mean of the two
    middle rows for even n), column by column; its error is a single CCD's error
    (about 2.2x the tutorial's for Gaia-4)."""
    d = pd.read_csv(path, comment='#')
    d = d[d.used_by_agis_al.astype(bool) & (d.outlier_flag == 0)]
    rows = []
    rng = np.random.default_rng(seed)
    for _, g in d.groupby('transit_id'):
        n = len(g)
        if reduction == 'demo':
            g = g.sort_values('centroid_pos_error_al')
            lo, hi = (n + 1)//2 - 1, n//2
            r = g[['epoch', 'scan_pos_angle', 'parallax_factor_al', 'centroid_pos_al',
                   'centroid_pos_error_al']].iloc[[lo, hi]].mean()
        elif reduction == 'tutorial':
            c = g.centroid_pos_al.to_numpy()
            if n == 1:
                err = g.ipd_error_al.iloc[0]
            else:
                boots = np.median(rng.choice(c, (nboot, n)), axis=1)
                err = max(np.std(boots, ddof=1), np.median(g.ipd_error_al)/np.sqrt(n))
            first = g.iloc[0]
            r = pd.Series({'epoch': g.epoch.mean(), 'scan_pos_angle': first.scan_pos_angle,
                           'parallax_factor_al': first.parallax_factor_al,
                           'centroid_pos_al': np.median(c), 'centroid_pos_error_al': err})
        else:
            raise ValueError(reduction)
        rows.append(r)
    d = pd.DataFrame(rows).sort_values('epoch')
    psi = np.deg2rad(d.scan_pos_angle.to_numpy())
    return Channel('astro', name, d.epoch.to_numpy(), d.centroid_pos_al.to_numpy(),
                   d.centroid_pos_error_al.to_numpy(), np.sin(psi), np.cos(psi),
                   d.parallax_factor_al.to_numpy())


def _rv_channels(t, rv, err, label, prefix='rv:'):
    out = []
    for lab in pd.unique(label):
        m = label == lab
        o = np.argsort(t[m])
        out.append(Channel('rv', f'{prefix}{lab}', t[m][o], rv[m][o], err[m][o]))
    return out


def load_rv_project(path):
    """A project RV table (rv_data/systems/*.csv): non-junk rows, one channel
    per dataset."""
    t = pd.read_csv(path, low_memory=False)
    ok = ~t.isjunk.astype(bool) & t.rv_ms.notna() & t.rv_err_ms.notna()
    t = t[ok]
    return _rv_channels(t.bjd_tcb.to_numpy() - BJD_MJD, t.rv_ms.to_numpy(),
                        t.rv_err_ms.to_numpy(), t.dataset.to_numpy())


def load_rv_simple(path, epoch_col='epoch', rv_col='rv', err_col='σ_rv', tel_col='tel'):
    """A small RV table (data_conversion/RV_data_*.csv; epochs already MJD)."""
    t = pd.read_csv(path)
    lab = t[tel_col].to_numpy() if tel_col in t else np.array(['rv']*len(t))
    return _rv_channels(t[epoch_col].to_numpy(float), t[rv_col].to_numpy(float),
                        t[err_col].to_numpy(float), lab)


def select(channels, spec='all'):
    """Keep the channels named in the comma-separated `spec`: 'all', 'astro',
    'rv' (every RV channel), or exact names ('gaia_dr4', 'rv:harps')."""
    keys = [s.strip() for s in spec.split(',') if s.strip()]
    keep = [c for c in channels if 'all' in keys or c.kind in keys or c.name in keys]
    if not keep:
        raise ValueError(f'spec {spec!r} selects nothing from {[c.name for c in channels]}')
    return keep


def span(channels):
    t = np.concatenate([c.t for c in channels])
    return t.max() - t.min()


if __name__ == '__main__':
    import sys
    if '--test' in sys.argv:
        g = load_gaia_dr4(os.path.join(ROOT, 'data_conversion/GaiaDR4_data_Gaia4.csv'))
        assert len(g) == 93 and g.kind == 'astro', len(g)
        r = load_rv_simple(os.path.join(ROOT, 'data_conversion/RV_data_HD114762.csv'))
        assert len(select(r + [g], 'gaia_dr4')) == 1 and len(select(r + [g], 'astro')) == 1
        assert len(select(r + [g], 'all')) == len(r) + 1
        print(f'data.py test passed: Gaia-4 {len(g)} transits over {span([g])/365.25:.2f} yr')
