"""Build marvels_bd_orbits.csv from the SDSS-III MARVELS brown-dwarf discovery papers.

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs (arXiv LaTeX sources, each unpacked into src_<id>/):
  1006.1627  Fleming+2010, ApJ 718, 1186   TYC 2949-557-1 b  2949-00557.tex: finalresults:orbit (eccentric case),
                                                             finalresults:mags (M_primary)
  1206.5514  Fleming+2012, AJ 144, 72      TYC 2930-872-1 b  tyc2930.tex: orbprops, stellarprops; M2,min from the
                                                             Sect. 4.2 equation (eqn:m2min)
  1211.6140  Ma+2013, AJ 145, 20           TYC 2087-255-1 b  MC12_arxiv.tex: tab:companionparams (MARVELS+TNG
                                                             column), tab:stellarparams
  1307.3308  Jiang+2013, AJ 146, 65        HIP 67526 b       ms.tex: tbl-2 (orbit), tbl-3 (first M_star row, Torres
                                                             relation), tbl-4 (sin i = 1 mass; same 65.0 MJ as the abstract)
  1304.2815  De Lee+2013, AJ 145, 155      GSC 03546-01452 b ms.tex: orbittab, hosttab
  1608.03597 Ma+2016, AJ 152, 112          HD 87646 b, c     tyc1415_201608_aj.tex: tbl:companionparams (no-trend
                                                             fit), tab:stellarparams
The first transcription typed these values by hand; they are now parsed from the LaTeX tables above.
Unit conversions: K in km/s -> m/s (Fleming+2010, Fleming+2012, Ma+2013); omega in radians -> degrees
(Fleming+2010, Ma+2013, De Lee+2013, Ma+2016); De Lee+2013 Tp from BJD_TDB-2450000 to full BJD_TDB.
Tp is left blank where the paper gives only a time of conjunction (Fleming+2010, Fleming+2012).
"""
import csv
import gzip
import math
import os
import re
import shutil
import sys
import tarfile
import urllib.request
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'marvels_bd_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {
    'src_1006.1627.gz': 'https://arxiv.org/e-print/1006.1627',
    'src_1206.5514.gz': 'https://arxiv.org/e-print/1206.5514',
    'src_1211.6140.gz': 'https://arxiv.org/e-print/1211.6140',
    'src_1307.3308.gz': 'https://arxiv.org/e-print/1307.3308',
    'src_1304.2815.gz': 'https://arxiv.org/e-print/1304.2815',
    'src_1608.03597.gz': 'https://arxiv.org/e-print/1608.03597',
}
ARCHIVES = {fn: fn[:-3] for fn in SOURCES}

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


# ---------------------------------------------------------------- source-specific helpers
def read(src, fn):
    with open(os.path.join(HERE, src, fn), encoding='utf-8') as f:
        return f.read().split('\n')


def table(lines, label, end=r'\\enddata|\\end\{tabular\}'):
    """{first cell: all cells} for the '&' rows of the table whose \\label is `label`."""
    out = {}
    for l in block(lines, r'\\label\{%s\}' % re.escape(label), end):
        if '&' in l and not l.strip().startswith('%'):
            c = cells(l)
            out.setdefault(c[0], c)
    return out


def get(tab, prefix):
    """Row of tab whose first cell starts with prefix."""
    return next(c for k, c in tab.items() if k.startswith(prefix))


def kms_to_ms(cell):
    """K in km/s -> m/s (x1000), keeping the printed precision (e.g. '6.109' -> '6109')."""
    v = val(cell)
    nd = len(v.split('.')[1]) if '.' in v else 0
    return '%.*f' % (max(nd - 3, 0), Decimal(v) * 1000)


def rad_to_deg(cell):
    """omega in radians -> degrees, rounded to 0.1 deg."""
    return '%.1f' % math.degrees(float(val(cell)))


def lo_hi_err(lo, hi):
    """Error from separate sigma(low)/sigma(high) columns: the printed value if equal, else their mean."""
    return lo if lo == hi else '%.6g' % ((float(lo) + float(hi)) / 2)


# ---------------------------------------------------------------- per-row text (hand-written)
# Numbers quoted in these notes come from the same papers (tables and text); they are not used as data.
REF = {
    'f10': 'arXiv:1006.1627 (Fleming+2010, ApJ 718, 1186) Table 4 eccentric case (Mstar Table 5)',
    'f12': 'arXiv:1206.5514 (Fleming+2012, AJ 144, 72) orbprops table (Mstar stellarprops table); msini from Sect. 4.2 text',
    'm13': 'arXiv:1211.6140 (Ma+2013, AJ 145, 20) companionparams table, MARVELS+TNG column (Mstar stellarparams table)',
    'j13': 'arXiv:1307.3308 (Jiang+2013, AJ 146, 65) Table 2 (Mstar Table 3); msini from abstract',
    'd13': 'arXiv:1304.2815 (De Lee+2013, AJ 145, 155) orbittab table (Mstar hosttab table)',
    'm16': 'arXiv:1608.03597 (Ma+2016, AJ 152, 112) companionparams table (no-trend fit; Mstar stellarparams table)',
}
NOTES = {
    'f10': 'MARVELS pilot (KeckET) + KPNO ET + HET/HRS RVs. Only Tc (time of conjunction) given: Tc = 2454868.9878 BJD_TDB; '
           'Tp left blank (e~0, omega 4.69+/-1.58 rad unconstrained). Circular fit essentially identical.',
    'f12': 'MARVELS-2b. MARVELS + TNG/SARG RVs. Only Tc given: Tc = HJD_UTC 2454842.2640; Tp blank. Long-term trend from '
           'tertiary stellar companion (P>2000 d). Paper argues true mass likely stellar (median M2 0.34-0.42 Msun if '
           'synchronized/coplanar).',
    'm13': 'MARVELS-4b. Tp is HJD_UTC. omega given in radians (4.086+/-0.041), converted. Active subgiant host.',
    'j13': 'MARVELS-5b. MARVELS + TNG/SARG RVs. Tp is HJD. msini = 65.0+/-2.9 MJ (sin i=1); with isotropic/ mass-ratio '
           'priors median mass 76-244 MJ, so may be stellar. Wilson+2016 compilation lists 62.6 MJ.',
    'd13': 'MARVELS-6b. MARVELS + TNG/SARG RVs; Tp in BJD_TDB; omega given in radians (1.998), converted. Residual linear '
           'slope 0.180 m/s/d.',
    'm16b': 'MARVELS-7b; massive planet. Two-Keplerian fit (b,c). HD 87646 is a close visual binary (~22 AU), jitter 240 m/s; '
            'Tp in JD_UTC; omega in radians (5.20+/-0.52), converted. Probably in NASA Exoplanet Archive.',
    'm16c': 'MARVELS-7c BD candidate. Host is in a close visual binary (~22 AU separation): Gaia astrometry likely '
            'complicated. Tp JD_UTC; omega in radians (1.95), converted.',
}


def row(host, comp, key, **kw):
    r = dict(host=host, companion=comp, P_d='', e_P_d='', K_ms='', e='', omega_deg='', Tp_bjd='', msini_mj='',
             a_au='', mstar_msun='', fit_type='RV-only', table_ref=REF[key[:3]], notes=NOTES[key])
    r.update(kw)
    return r


def main():
    get_inputs()
    rows = []

    # Fleming+2010: columns Parameter, value, error (eccentric case), value, error (circular case).
    L = read('src_1006.1627', '2949-00557.tex')
    t = table(L, 'finalresults:orbit')
    st = table(L, 'finalresults:mags')
    P = get(t, 'Period')
    rows.append(row('TYC 2949-557-1', 'b', 'f10', P_d=P[1], e_P_d=P[2],
                    K_ms=kms_to_ms(get(t, 'K (km')[1]),            # km/s -> m/s
                    e=get(t, '$e$')[1],
                    omega_deg=rad_to_deg(get(t, '$\\omega$ (rad)')[1]),  # rad -> deg
                    msini_mj=get(t, '$m_{\\rm{min}}')[1], mstar_msun=get(st, '$M_{primary}$')[1]))

    # Fleming+2012: columns Parameter, median, sigma(low), sigma(high). omega has no unit (degrees).
    L = read('src_1206.5514', 'tyc2930.tex')
    t = table(L, 'orbprops')
    st = table(L, 'stellarprops')
    P = get(t, '$P$')
    # M2,min appears only in the Sect. 4.2 equation "M_{2,\rm min}= 68.1 \pm 3.0 ~ \rm{M_{Jup}} = ..." (eqn:m2min).
    msini = re.search(r'M_\{2,\\rm min\}=\s*(' + NUM + ')', '\n'.join(L)).group(1)
    rows.append(row('TYC 2930-872-1', 'b', 'f12', P_d=P[1], e_P_d=lo_hi_err(P[2], P[3]),
                    K_ms=kms_to_ms(get(t, '$K$')[1]),              # km/s -> m/s
                    e=get(t, '$e$')[1], omega_deg=get(t, '$\\omega$')[1], msini_mj=msini,
                    mstar_msun=get(st, '$M_{*}$')[1]))

    # Ma+2013: columns Parameter, MARVELS+TNG (used), MARVELS. Tp is full HJD_UTC.
    L = read('src_1211.6140', 'MC12_arxiv.tex')
    t = table(L, 'tab:companionparams')
    st = table(L, 'tab:stellarparams')
    P = get(t, '$P$')[1]
    rows.append(row('TYC 2087-255-1', 'b', 'm13', P_d=val(P), e_P_d=err(P),
                    K_ms=kms_to_ms(get(t, '$K$')[1]),              # km/s -> m/s
                    e=val(get(t, '$e$')[1]),
                    omega_deg=rad_to_deg(get(t, '$w$')[1]),         # rad -> deg
                    Tp_bjd=val(get(t, '$T_{\\rm periastron}')[1]), msini_mj=val(get(t, 'Minimum Mass')[1]),
                    a_au=val(get(t, '$a$')[1]), mstar_msun=val(get(st, 'Mass')[1])))

    # Jiang+2013: columns Parameter, Units, Value. Tp is full HJD; omega in degrees.
    L = read('src_1307.3308', 'ms.tex')
    t = table(L, 'tbl-2')
    st = table(L, 'tbl-3')          # first M_star row = Torres et al. (2010) relation, 1.10 Msun (as in the abstract)
    mt = table(L, 'tbl-4')          # "None (sin i = 1)" prior row: 65.0 MJ
    P = get(t, '$P$')[2]
    rows.append(row('HIP 67526', 'b', 'j13', P_d=val(P), e_P_d=err(P), K_ms=val(get(t, '$K$')[2]),
                    e=val(get(t, '$e$')[2]), omega_deg=val(get(t, '$\\omega$')[2]),
                    Tp_bjd=val(get(t, '$T_0$')[2]), msini_mj=val(get(mt, 'None')[1]),
                    mstar_msun=val(get(st, '$M_{\\star}$')[1])))

    # De Lee+2013: columns Parameter, Units, Value. Tp in BJD_TDB-2450000; omega in radians.
    L = read('src_1304.2815', 'ms.tex')
    t = table(L, 'orbittab')
    st = table(L, 'hosttab')
    P = get(t, '$P$')[2]
    rows.append(row('GSC 03546-01452', 'b', 'd13', P_d=val(P), e_P_d=err(P), K_ms=val(get(t, '$K$')[2]),
                    e=val(get(t, '$e$')[2]),
                    omega_deg=rad_to_deg(get(t, '$\\omega$')[2]),   # rad -> deg
                    Tp_bjd=add_offset(get(t, '$T_P$')[2], 2450000),  # BJD_TDB - 2450000 -> BJD_TDB
                    msini_mj=val(get(t, '$M_b')[2]), mstar_msun=val(get(st, 'M$_{*}$')[1])))

    # Ma+2016: columns Parameter, HD 87646 b, HD 87646 c. Tp is full JD_UTC; omega in radians.
    L = read('src_1608.03597', 'tyc1415_201608_aj.tex')
    t = table(L, 'tbl:companionparams')
    st = table(L, 'tab:stellarparams')
    for k, comp in ((1, 'b'), (2, 'c')):
        P = get(t, '$P$')[k]
        rows.append(row('HD 87646', comp, 'm16' + comp, P_d=val(P), e_P_d=err(P), K_ms=val(get(t, '$K$')[k]),
                        e=val(get(t, '$e$')[k]),
                        omega_deg=rad_to_deg(get(t, '$\\omega$ (radians)')[k]),  # rad -> deg
                        Tp_bjd=val(get(t, '$T_{\\rm periastron}')[k]), msini_mj=val(get(t, 'Minimum Mass')[k]),
                        a_au=val(get(t, '$a$')[k]), mstar_msun=val(get(st, 'Mass')[1])))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
