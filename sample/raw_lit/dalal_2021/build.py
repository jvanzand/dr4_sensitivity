"""Build dalal_2021_orbits.csv from Dalal et al. 2021, A&A 651, A11 (SOPHIE XVII; arXiv:2105.09741).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_2105.09741.gz (arXiv LaTeX source), unpacked into src_2105.09741/.
  - Orbits: table tab:BDorbit ("Table 2", transposed: one column per BD) in aa.tex.
  - Host masses: table tab:stellar ("Table A.2", M_Torres column) in aa.tex.
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
OUT = os.path.join(HERE, 'dalal_2021_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_2105.09741.gz': 'https://arxiv.org/e-print/2105.09741'}
ARCHIVES = {'src_2105.09741.gz': 'src_2105.09741'}

COLS = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp_bjd', 'msini_mj',
        'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes', 'e_K_ms', 'e_msini_mj', 'e_e', 'e_a_au']


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


def err_nz(cell):
    """err() but a published error of 0 becomes blank."""
    e = err(cell)
    return e if e and float(e) != 0 else ''


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


# ---------------------------------------------------------------- per-row text (hand-written)
NOTE = 'SOPHIE/SOPHIE+ RVs; Yorbit single-Keplerian fit; Tp is periastron (BJD-2450000 in paper, converted).'
NOTES_EXTRA = {  # from the paper's Sect. on HD 205521 (Hipparcos and Gaia DR1 analyses)
    'HD 205521': ' Hipparcos astrometric orbit (>3 sigma) gives Mc=0.48+/-0.04 Msun and Gaia DR1 excess-noise '
                 'analysis 0.20 Msun: likely stellar, not a BD.'}
NOTE_ERR = (' e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au: published 1-sigma errors from Table 2 (those on M sin i and '
            'a include the host-mass uncertainty, per the table note).')
TABLE_REF = 'arXiv:2105.09741 Table 2 (BD orbits; Mstar Table A.2)'


def main():
    get_inputs()
    with open(os.path.join(HERE, 'src_2105.09741', 'aa.tex'), encoding='utf-8') as f:
        L = f.read().split('\n')

    # Table 2 is transposed: each row is "<parameter> & <unit> & <BD1> & <BD2> & <BD3> & <BD4>".
    tab = {}
    for l in block(L, r'\\label\{tab:BDorbit\}', r'\\end\{tabular\}'):
        c = cells(l)
        if len(c) == 6:
            tab[c[0]] = c[2:]

    def row(pattern):
        """Values of the first table row whose parameter cell matches pattern."""
        return next(v for k, v in tab.items() if re.match(pattern, k))

    names = [star_name(n) for n in row(r'Name$')]
    P, K, ecc, om = row(r'P$'), row(r'K$'), row(r'\$e\$'), row(r'\$\\omega\$')
    Tp, msini, a = row(r'\$T_p'), row(r'\\msini'), row(r'\$a\^')

    # Table A.2 (stellar parameters): M_Torres [M_sun] is the 7th column.
    mstar = {}
    for l in block(L, r'\\label\{tab:stellar\}', r'\\end\{tabular\}'):
        c = cells(l)
        if len(c) >= 7 and '\\pm' in l:
            mstar[star_name(c[0])] = val(c[6])

    rows = []
    for k, n in enumerate(names):
        rows.append(dict(host=bd(n), companion='b', P_d=val(P[k]), e_P_d=err(P[k]), K_ms=val(K[k]), e=val(ecc[k]),
                         omega_deg=val(om[k]),
                         Tp_bjd=add_offset(Tp[k], 2450000),  # Tp is given as BJD-2450000
                         msini_mj=val(msini[k]), a_au=val(a[k]), mstar_msun=mstar.get(n, ''),
                         fit_type='RV-only', table_ref=TABLE_REF, notes=NOTE + NOTES_EXTRA.get(n, '') + NOTE_ERR,
                         e_K_ms=err_nz(K[k]), e_msini_mj=err_nz(msini[k]), e_e=err_nz(ecc[k]),
                         e_a_au=err_nz(a[k])))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
