"""Build teng_2023_orbits.csv from Teng et al. 2023, PASJ 75, 1030 (Okayama Planet Search revisit; arXiv:2308.05343).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_2308.05343.gz (arXiv LaTeX source), unpacked into src_2308.05343/.
  - Orbits: table tab:orbpar_orbits ("Orbital parameters", RV-only fits) in Revisiting.tex.
  - Host masses: table tab:star_pars ("Stellar parameters"), column M_star.
Only the six companions with M sin i above ~8 MJ (the brown-dwarf candidates and the most massive planets) are kept
(KEEP below). omega is given in radians (converted to degrees) and Tp in JD-2450000 (converted to JD).
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

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'teng_2023_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_2308.05343.gz': 'https://arxiv.org/e-print/2308.05343'}
ARCHIVES = {'src_2308.05343.gz': 'src_2308.05343'}

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
    """err() but a published error of 0 (e.g. '-0.000, +0.000') becomes blank."""
    e = err(cell)
    return e if e and float(e) != 0 else ''


def cells(line):
    return [c.strip() for c in line.split('\\\\')[0].split('&')]


def add_offset(cell, off):
    """Reduced date + offset (e.g. JD-2450000 -> JD), keeping the number of printed decimals."""
    v = val(cell)
    nd = len(v.split('.')[1]) if '.' in v else 0
    return '%.*f' % (nd, float(v) + off)


def name(cell):
    """'$\\nu$ Oph b' -> 'nu Oph b'; '\\textbf{75 Cet c}' -> '75 Cet c'."""
    s = re.sub(r'\$\\(\w+)\$', r'\1', cell).replace('\\textbf{', '').replace('}', '')
    return re.sub(r'\s+', ' ', s).strip()


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
NOTE = ('Okayama Planet Search (HIDES) + literature RVs, RV-only refit. omega given in radians, converted; '
        'Tp JD-2450000 in paper (JD). Likely in NASA Exoplanet Archive if <30 MJ. ')
KEEP = {  # companion -> extra note (from the paper's discussion and cross-checks with other sources)
    '11 Com b': 'Liu+2008 discovery (Xinglong); also Lick BD list.',
    'nu Oph b': 'Also in Quirrenbach+2019 (lick_giants).',
    'nu Oph c': 'a=5.83 AU slightly >5 AU. Also in Quirrenbach+2019.',
    'HD 175679 b': 'Wang+2012 discovery. Joint RV+HGCA fit (table "Joint fitting results") gives true mass 61.8 '
                   '(+15/-9.4) MJ, i=39 or 152 deg.',
    'HD 180314 b': 'Sato+2010 discovery.',
    'HD 14067 b': 'msini 8.8 MJ (massive planet) but joint RV+HGCA fit gives true mass 15.3 (+6.2/-4.5) MJ; '
                  'a=5.11 AU borderline.',
}
NOTE_ERR = (' e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au: published 1-sigma errors from the orbit table; asymmetric '
            'errors averaged; blank where the table prints 0.')
TABLE_REF = 'arXiv:2308.05343 (Teng+2023, PASJ) tab:orbpar_orbits (Mstar tab:star_pars)'


def main():
    get_inputs()
    with open(os.path.join(HERE, 'src_2308.05343', 'Revisiting.tex'), encoding='utf-8') as f:
        L = f.read().split('\n')

    # Stellar parameters: name, parallax, V, B-V, SpT, Teff_sp, [Fe/H]_sp, logg_sp, Teff, [Fe/H], logg, L, M, R.
    mstar = {}
    for l in block(L, r'\\tbl\{Stellar parameters', r'^Source &'):
        c = cells(l)
        if len(c) >= 14 and '$' in c[1]:
            mstar[name(c[0])] = val(c[12])

    # Orbital parameters (the table body precedes its \label): name, P [d], K [m/s], e, omega [rad],
    # Tp [JD-2450000], Mp sin i [MJ], a [au], trends.
    orb = {}
    for l in block(L, r'\\tbl\{Orbital parameters', r'\\label\{tab:orbpar_orbits\}'):
        c = cells(l)
        if len(c) >= 8 and '_{' in l:
            orb[name(c[0])] = c[1:]

    rows = []
    for comp, extra in KEEP.items():
        P, K, e, w, Tp, m, a = orb[comp][:7]
        host, letter = comp.rsplit(' ', 1)
        rows.append(dict(host=host, companion=letter, P_d=val(P), e_P_d=err(P), K_ms=val(K), e=val(e),
                         omega_deg='%.2f' % math.degrees(float(val(w))),
                         Tp_bjd=add_offset(Tp, 2450000),  # JD (not BJD) despite the column name
                         msini_mj=val(m), a_au=val(a), mstar_msun=mstar[host], fit_type='RV-only',
                         table_ref=TABLE_REF, notes=NOTE + extra + NOTE_ERR,
                         e_K_ms=err_nz(K), e_msini_mj=err_nz(m), e_e=err_nz(e), e_a_au=err_nz(a)))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
