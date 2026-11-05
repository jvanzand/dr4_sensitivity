"""Build grieves_2017_orbits.csv from Grieves et al. 2017, MNRAS 467, 4264 (SDSS-III MARVELS; arXiv:1702.01784).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_1702.01784.gz (arXiv LaTeX source), unpacked into src_1702.01784/.
  - Orbits: table tab:compar ("Table 2", two stacked tabulars: M sin i, P, e, K; then a, Tp, omega, gamma,
    dv/dt, fit rms, chi2) in grieves2017_arxiv.tex.
  - Host masses: table tab:starpar ("Table 1") in grieves2017_arxiv.tex.
MARVELS-14b (100.2 MJ) and MARVELS-15b (97.6 MJ) are excluded (M sin i > 85 MJ).
(arxiv_1702.01784.src and arxiv_src/ in this directory are an identical earlier download of the same e-print;
they are not used.)
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
OUT = os.path.join(HERE, 'grieves_2017_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_1702.01784.gz': 'https://arxiv.org/e-print/1702.01784'}
ARCHIVES = {'src_1702.01784.gz': 'src_1702.01784'}

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
TABLE_REF = 'arXiv:1702.01784 Table 2 (Mstar Table 1)'
NOTE = ('%s. MARVELS UF2D RVs (+ any follow-up as in paper), RVLIN fit, bootstrap errors. Tp JD-2450000 in paper '
        '(converted; JD not BJD). fit rms %s m/s; dv/dt %s m/s/d.')
NOTE_END = ' Candidate: paper notes some fits have weak phase coverage/aliases.'
MSINI_MAX = 85.0  # MJ


def tyc(name):
    """'TYC 0123-00549-1' -> 'TYC 123-549-1'; other names unchanged."""
    m = re.match(r'TYC\s*0*(\d+)-0*(\d+)-(\d+)', name.strip())
    return 'TYC %s-%s-%s' % m.groups() if m else name.strip()


def main():
    get_inputs()
    with open(os.path.join(HERE, 'src_1702.01784', 'grieves2017_arxiv.tex'), encoding='utf-8') as f:
        L = f.read().split('\n')

    # Table 1 (host stars): name (with optional dagger and alias in parentheses), Teff, [Fe/H], logg, M*, R*, V.
    ms = {}
    for l in block(L, r'\\label\{tab:starpar\}', r'\\end\{tabular\}'):
        c = cells(l)
        if len(c) == 7 and '\\pm' in l:
            ms[re.sub(r'\$.*?\$', '', c[0]).split('(')[0].strip()] = val(c[4])

    # Table 2, part A: star, designation, M sin i, P, e, K.  Part B: designation, a, Tp, omega, gamma, dv/dt, rms, chi2.
    A, B = {}, {}
    for l in block(L, r'\\label\{tab:compar\}', r'\\end \{table\*\}'):
        c = cells(l)
        if len(c) == 6 and 'MARVELS-' in c[1]:
            A[c[1]] = (c[0].split('(')[0].strip(), c[2:])
        if len(c) == 8 and c[0].startswith('MARVELS-'):
            B[c[0]] = c[1:]

    rows = []
    for d, (star, (msini, P, ecc, K)) in A.items():
        a, Tp, w, gam, dvdt, rms, chi = B[d]
        m = float(val(msini))
        if m > MSINI_MAX:
            continue
        note = NOTE % (d, val(rms), val(dvdt))
        if 13 < m < 14:
            note += ' msini %s+/-%s near deuterium limit.' % (val(msini), err(msini))
        rows.append(dict(host=tyc(star), companion='b', P_d=val(P), e_P_d=err(P), K_ms=val(K), e=val(ecc),
                         omega_deg=val(w),
                         Tp_bjd=add_offset(Tp, 2450000),  # Tp is given as JD-2450000
                         msini_mj=val(msini), a_au=val(a), mstar_msun=ms[star], fit_type='RV-only',
                         table_ref=TABLE_REF, notes=note + NOTE_END))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
