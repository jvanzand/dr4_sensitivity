"""Build troup_2016_orbits.csv from Troup et al. 2016, AJ 151, 85 (APOGEE-1 companions; arXiv:1601.00688).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_1601.00688.gz (arXiv source, unpacked into src_1601.00688/), which contains goldOrbit-dr12.fits:
        the 382 gold-sample orbits (same content as VizieR J/AJ/151/85 table4, whose table4.dat URL returned 404).
        ReadMe is the VizieR description of J/AJ/151/85 (documentation only, not parsed).
Keeps orbits with 5 <= M sin i <= 85 MJ, in FITS row order.
Conversions: MSINI (Msun) -> MJ with 1047.57 MJ/Msun; OMEGA radians -> degrees. TPERI is used as given (JD; the
time standard is not stated). SEMIAMP_ERR = 0 means "missing" (ReadMe: ?=0) but is printed as e_K=0 in the notes.
"""
import csv
import gzip
import os
import re
import shutil
import sys
import tarfile
import urllib.request

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'troup_2016_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_1601.00688.gz': 'https://arxiv.org/e-print/1601.00688',
           'ReadMe': 'https://cdsarc.cds.unistra.fr/ftp/J/AJ/151/85/ReadMe'}
ARCHIVES = {'src_1601.00688.gz': 'src_1601.00688'}

COLS = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp_bjd', 'msini_mj',
        'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes']


# ---------------------------------------------------------------- inputs
def get_inputs():
    """Check (or with --download fetch) every raw input, then unpack archives that are not unpacked yet."""
    download = '--download' in sys.argv[1:]
    for fn, url in SOURCES.items():
        path = os.path.join(HERE, fn)
        if os.path.exists(path):
            continue
        if not download:
            sys.exit('missing input %s (from %s); rerun with --download' % (fn, url))
        print('downloading', url)
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (build.py)'})
        with urllib.request.urlopen(req, timeout=300) as r, open(path + '.part', 'wb') as f:
            shutil.copyfileobj(r, f)
        os.replace(path + '.part', path)
    for fn, dest in ARCHIVES.items():
        unpack(os.path.join(HERE, fn), os.path.join(HERE, dest))


def unpack(archive, dest):
    """Unpack an arXiv e-print (tar, usually gzipped; or a single gzipped file) into dest, unless dest exists.
    Only regular files and directories that land inside dest are extracted; nothing is executed."""
    if os.path.isdir(dest):
        return
    os.makedirs(dest)
    root = os.path.realpath(dest)
    try:
        tf = tarfile.open(archive, 'r:*')
    except tarfile.ReadError:
        with gzip.open(archive) as g, open(os.path.join(dest, 'main.tex'), 'wb') as f:
            shutil.copyfileobj(g, f)
        return
    with tf:
        for m in tf.getmembers():
            target = os.path.realpath(os.path.join(dest, m.name))
            if (m.isfile() or m.isdir()) and target.startswith(root + os.sep):
                tf.extract(m, dest)


# ---------------------------------------------------------------- LaTeX cell helpers
NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def val(cell):
    """Central value of a LaTeX cell like '$1718\\pm11$' or '$2.6^{+0.1}_{-0.2}$', as the printed string."""
    s = cell.replace('$', '').replace('{', '').replace('}', '').replace('~', '')
    s = re.sub(r'\s+', '', s.replace('\\,', '')).replace('--', '')
    m = re.match(NUM, s)
    return m.group(0) if m else ''


def err(cell):
    """Symmetric error from a LaTeX cell: the +/- value, or the mean of |+| and |-| if asymmetric."""
    s = re.sub(r'\s+', '', cell.replace('$', '').replace('{', '').replace('}', ''))
    m = re.search(r'\\pm(' + NUM + ')', s)
    if m:
        return m.group(1)
    m = re.search(r'\^\+?(' + NUM + r')_-?(' + NUM + ')', s) or re.search(r'_-?(' + NUM + r')\^\+?(' + NUM + ')', s)
    if m:
        return '%.6g' % ((abs(float(m.group(1))) + abs(float(m.group(2)))) / 2)
    return ''


def cells(line):
    return [c.strip() for c in line.split('\\\\')[0].split('&')]


def add_offset(cell, off):
    """Reduced date + offset (e.g. BJD-2400000 -> BJD), keeping the number of printed decimals."""
    v = val(cell)
    nd = len(v.split('.')[1]) if '.' in v else 0
    return '%.*f' % (nd, float(v) + off)


def bd(name):
    """'BD+73 0275' -> SIMBAD-style 'BD+73 275'; other names unchanged."""
    m = re.match(r'BD\s*([+-])\s*(\d{2})\s*(\d+)$', name)
    return 'BD%s%s %d' % (m.group(1), m.group(2), int(m.group(3))) if m else name


def star_name(cell):
    """'HD\\,134113$^\\dagger$' -> 'HD 134113'; 'BD+73~0275B' -> 'BD+73 0275B'."""
    return cell.split('$')[0].replace('\\,', ' ').replace('~', ' ').strip()


def block(lines, start_pat, end_pat):
    """Lines after the first line matching start_pat, up to the next line matching end_pat."""
    i0 = next(i for i, l in enumerate(lines) if re.search(start_pat, l))
    out = []
    for l in lines[i0 + 1:]:
        if re.search(end_pat, l):
            break
        out.append(l)
    return out


# ---------------------------------------------------------------- per-row text
TABLE_REF = 'arXiv:1601.00688 goldOrbit-dr12.fits (= Table 4, VizieR J/AJ/151/85)'
NOTE = ('APOGEE-1 DR12 gold sample; class=%s; NRV=%d over %.0f d baseline; sigma_RV~%.0f m/s; e_K=%.3g m/s; e_e=%.3g; '
        'PUI=%.2f VUI=%.2f; slope=%.3f m/s/d. Tp is JD as given (time standard not stated). msini converted from Msun '
        'with 1047.57 MJ/Msun. Few-epoch, low-precision orbit: treat as candidate; period aliases possible.')
MJ_PER_MSUN = 1047.57   # Jupiter masses per solar mass, used to convert MSINI (Msun) to MJ
MSINI_RANGE = (5, 85)   # MJ


def main():
    get_inputs()
    with fits.open(os.path.join(HERE, 'src_1601.00688', 'goldOrbit-dr12.fits')) as hdul:
        d = hdul[1].data
        rows = []
        for r in d:
            msini = r['MSINI'] * MJ_PER_MSUN  # Msun -> MJ
            if not (MSINI_RANGE[0] <= msini <= MSINI_RANGE[1]):
                continue
            aid = r['APOGEE_ID'].strip()
            host = '2MASS J' + aid[2:]  # APOGEE_ID is '2M' + 2MASS designation
            rows.append(dict(
                host=host, companion='b', P_d='%.6g' % r['PERIOD'], e_P_d='%.3g' % r['PERIOD_ERR'],
                K_ms='%.1f' % r['SEMIAMP'], e='%.4f' % r['ECC'],
                omega_deg='%.2f' % np.degrees(r['OMEGA']),  # radians -> degrees
                Tp_bjd='%.5f' % r['TPERI'], msini_mj='%.2f' % msini, a_au='%.4f' % r['SEMIMAJ'],
                mstar_msun='%.3f' % r['MSTAR'], fit_type='RV-only', table_ref=TABLE_REF,
                notes=NOTE % (r['STARTYPE'].strip(), r['NVISITS'], r['BASELINE'], r['SIGMA_V'], r['SEMIAMP_ERR'],
                              r['ECC_ERR'], r['PUI'], r['VUI'], r['SLOPE'])))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
