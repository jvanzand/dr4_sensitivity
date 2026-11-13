"""Build diaz_2012_orbits.csv from Diaz et al. 2012, A&A 538, A113 (SOPHIE; arXiv:1111.1168).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_1111.1168.gz (arXiv LaTeX source), unpacked into src_1111.1168/.
  - Orbits: TableParameters_v3.tex ("Table 5", transposed, two blocks of four stars).
  - Host masses: TableStellarParameters_v1.tex ("Table 3", M_s row).
HD 167215 (Mc sin i = 92 MJ) is excluded.
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
OUT = os.path.join(HERE, 'diaz_2012_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_1111.1168.gz': 'https://arxiv.org/e-print/1111.1168'}
ARCHIVES = {'src_1111.1168.gz': 'src_1111.1168'}

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
NOTE = 'SOPHIE RVs, single Keplerian; T0 = periastron (BJD-2454000 in paper, converted). '
NOTES_EXTRA = {  # from the paper's table footnotes and per-object discussion
    'HD22781': 'Likely in NASA Exoplanet Archive (msini 13.65). A second companion probably present (linear+quadratic drift fitted); parameters may change.',
    'HD156279': 'Massive planet (msini 9.71); likely in NASA Exoplanet Archive.',
    'HD14651': 'Hipparcos 3-sigma upper limit on Mc 0.28-0.40 Msun.',
    'HD30246': 'Poorly constrained orbit (K +1600/-560 m/s, e 0.84+/-0.08); msini +20.3/-8.2 MJ.',
    'HD92320': 'Hipparcos upper limit 0.26-0.34 Msun.',
    'HD137510': 'Fit combines SOPHIE with McDonald, TLS and Lick RVs (Endl+2004 discovery). Hipparcos constrains Mc to 0.019-0.057 Msun (99.7%); so a BD.',
    'HD16702': 'Hipparcos astrometric orbit: Mc = 0.35-0.45 Msun (99.7%): low-mass star, not BD.',
}
NOTE_ERR = (' e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au: published 1-sigma errors from Table 5; where the table gives '
            'asymmetric errors, the two sides are averaged.')
TABLE_REF = 'arXiv:1111.1168 Table 5 (Mstar Table 3)'
SKIP = {'HD167215'}  # Mc sin i = 92 MJ


def transposed(lines, keys):
    """Parse a transposed table made of blocks: a header row naming four stars, then '<parameter> & v1 & ... & v4'
    rows. Returns {star: {key: cell}} for parameter cells starting with one of keys, plus the star order."""
    out, order, names = {}, [], None
    for l in lines:
        if l.strip().startswith('%'):
            continue
        c = cells(l)
        if len(c) == 5 and c[0] in ('', 'Parameters') and c[1].startswith('HD'):
            names = [re.sub(r'\\tablefootmark\{[^}]*\}', '', x) for x in c[1:]]
            order += names
            continue
        for k in keys:
            if names and c[0].startswith(k) and len(c) >= 5:
                for n, x in zip(names, c[1:5]):
                    out.setdefault(n, {})[k] = x
    return out, order


def main():
    get_inputs()
    src = os.path.join(HERE, 'src_1111.1168')
    with open(os.path.join(src, 'TableParameters_v3.tex'), encoding='utf-8') as f:
        T = f.read().split('\n')
    with open(os.path.join(src, 'TableStellarParameters_v1.tex'), encoding='utf-8') as f:
        S = f.read().split('\n')
    kP, ke, kw, kK, kT0, ka, kM = '$P$', '$e$', '$\\omega$', '$K$', '$T_0$', '$a$', '$M_\\textrm{c} \\sin i'
    P, order = transposed(T, [kP, ke, kw, kK, kT0, ka, kM])
    M, _ = transposed(S, ['$M_\\textrm{s}$'])

    rows = []
    for n in order:
        if n in SKIP:
            continue
        p = P[n]
        rows.append(dict(host='HD ' + n[2:], companion='b', P_d=val(p[kP]), e_P_d=err(p[kP]), K_ms=val(p[kK]),
                         e=val(p[ke]), omega_deg=val(p[kw]),
                         Tp_bjd=add_offset(p[kT0], 2454000),  # T0 is given as BJD-2454000
                         msini_mj=val(p[kM]), a_au=val(p[ka]), mstar_msun=val(M[n]['$M_\\textrm{s}$']),
                         fit_type='RV-only', table_ref=TABLE_REF, notes=NOTE + NOTES_EXTRA.get(n, '') + NOTE_ERR,
                         e_K_ms=err_nz(p[kK]), e_msini_mj=err_nz(p[kM]), e_e=err_nz(p[ke]),
                         e_a_au=err_nz(p[ka])))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
