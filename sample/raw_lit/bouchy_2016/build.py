"""Build bouchy_2016_orbits.csv from Bouchy et al. 2016, A&A 585, A46 (SOPHIE VIII; arXiv:1511.08397).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_1511.08397.gz (arXiv LaTeX source), unpacked into src_1511.08397/.
  - Orbits: table "fit" ("Table 4", transposed: one column per star) in 26347.tex. P is given in years.
  - Host masses: table "spectro" ("Table 2") in 26347.tex.
Only companions with a <= 7 AU are kept (HD 14348, HD 72946, HD 209262; all beyond 5 AU, kept as borderline).
"""
import csv
import gzip
import os
import re
import shutil
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'bouchy_2016_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_1511.08397.gz': 'https://arxiv.org/e-print/1511.08397'}
ARCHIVES = {'src_1511.08397.gz': 'src_1511.08397'}

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
TABLE_REF = 'arXiv:1511.08397 Table 4 (Mstar Table 2)'
NOTE = ('OUTSIDE a<=5 AU target range (a=%s AU, P~%.0f yr); kept as borderline. ELODIE+SOPHIE RVs, Yorbit+MCMC; '
        'P given in years in paper (converted with 365.25 d/yr). No Hipparcos detection; Hipparcos upper limit %s Msun.')
A_MAX = 7.0  # AU; HD 10844 (10.2 AU) and HD 18757 (22.2 AU) are dropped
DAYS_PER_YEAR = 365.25


def obj(cell):
    """'\\object{HD14348}' -> 'HD14348'."""
    return cell.replace('\\object{', '').replace('}', '').strip()


def main():
    get_inputs()
    with open(os.path.join(HERE, 'src_1511.08397', '26347.tex'), encoding='utf-8') as f:
        L = f.read().split('\n')

    # Table 4 is transposed: "<parameter> & <star1> & ... & <star5>".
    tab = {}
    for l in block(L, r'\\label\{fit\}', r'\\end\{tabular\}'):
        c = cells(l)
        if len(c) == 6:
            tab[c[0]] = c[1:]

    def row(prefix):
        return next(v for k, v in tab.items() if k.startswith(prefix))

    names = [obj(x) for x in tab['Parameters']]
    Py, K, ecc, om, Tp = row('P [years]'), row('K ['), tab['e'], row('$\\omega$'), row('$T_p$')
    a, msini, mup = row('$a$ [AU]'), row('Mc\\,\\sini'), row('Mc,upper')

    # Table 2 (stellar parameters): M_star [M_sun] is the 6th column.
    ms = {}
    for l in block(L, r'\\label\{spectro\}', r'\\end\{tabular\}'):
        c = cells(l)
        if len(c) == 8 and 'object' in c[0]:
            ms[obj(c[0])] = val(c[5])

    rows = []
    for i, n in enumerate(names):
        if float(val(a[i])) > A_MAX:
            continue
        # P and its (mean asymmetric) error are given in years: convert to days with 365.25 d/yr.
        P = '%.2f' % (float(val(Py[i])) * DAYS_PER_YEAR)
        eP = '%.2f' % (float(err(Py[i])) * DAYS_PER_YEAR)
        rows.append(dict(host='HD ' + n[2:], companion='b', P_d=P, e_P_d=eP, K_ms=val(K[i]), e=val(ecc[i]),
                         omega_deg=val(om[i]),
                         Tp_bjd=add_offset(Tp[i], 2450000),  # Tp is given as (BJD) - 2450000
                         msini_mj=val(msini[i]), a_au=val(a[i]), mstar_msun=ms[n], fit_type='RV-only',
                         table_ref=TABLE_REF, notes=NOTE % (val(a[i]), float(val(Py[i])), val(mup[i]))))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
