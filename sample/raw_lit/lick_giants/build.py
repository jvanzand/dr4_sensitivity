"""Build lick_giants_orbits.csv: brown-dwarf companions of Lick giants (Quirrenbach+2019; Mitchell+2013).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs (arXiv LaTeX sources, each unpacked into src_<id>/):
  1904.03557  Quirrenbach+2019, A&A 624, A18 (nu Oph b, c): hip88048_arXiv.tex, table:orb_par_stable (Keplerian
              block, the first of three fits), table:phys_param (HB-model mass), and the fit epoch from Sect. 3 text.
  1305.5107   Mitchell+2013, A&A 555, A87 (tau Gem b): planets.tex, tables "orbit" (Best Fit Orbital Parameters) and
              "param" (Stellar Parameters).
  1412.4634   Reffert+2015 (Lick giants VII, occurrence rates): kept as the Lick-sample reference; not parsed.
Derived values: nu Oph Tp from the mean anomaly M0 at the fit epoch, Tp = epoch - M0/360 * P (JD);
tau Gem Tp from JD-2450000 to full JD.
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
OUT = os.path.join(HERE, 'lick_giants_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {
    'src_1904.03557.gz': 'https://arxiv.org/e-print/1904.03557',
    'src_1305.5107.gz': 'https://arxiv.org/e-print/1305.5107',
    'src_1412.4634.gz': 'https://arxiv.org/e-print/1412.4634',
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


def table(lines, start, end):
    """{first cell: all cells} for the '&' rows between the line matching start and the next line matching end."""
    out = {}
    for l in block(lines, start, end):
        if '&' in l and not l.strip().startswith('%'):
            c = cells(l)
            out.setdefault(c[0], c)
    return out


def get(tab, prefix):
    return next(c for k, c in tab.items() if k.startswith(prefix))


# ---------------------------------------------------------------- per-row text (hand-written)
REF_Q19 = ('arXiv:1904.03557 (Quirrenbach+2019, A&A 624, A18) table:orb_par_stable, Keplerian block '
           '(Mstar HB model, table:phys_param)')
REF_M13 = 'arXiv:1305.5107 (Mitchell+2013, A&A) Table "Best Fit Orbital Parameters" (Mstar Table "Stellar Parameters")'
NOTES = {
    ('nu Oph', 'b'): 'Lick + OAO RVs, double-Keplerian. Paper gives mean anomaly M0=235.78 deg at first epoch JD 2451853.595; '
                     'Tp derived here as epoch - M0/360*P (JD). Likely in NASA Exoplanet Archive (<30 MJ). 6:1 MMR with c.',
    ('nu Oph', 'c'): 'a=5.93 AU slightly beyond 5 AU (P=8.7 yr). M0=222.75 deg at JD 2451853.595; Tp derived (JD). Likely in '
                     'NASA Exoplanet Archive (<30 MJ).',
    ('tau Gem', 'b'): 'Lick RVs; T0 = periastron time, JD-2450000 in paper (JD, not BJD); poorly constrained (+/-14.8 d) since '
                      'e~0.03. Quirrenbach+2019 lists msini 22.1 with Stock+2018 Mstar 2.47. Likely in NASA Exoplanet '
                      'Archive (<30 MJ).',
}


def main():
    get_inputs()
    rows = []

    # Quirrenbach+2019: Keplerian block of table:orb_par_stable (columns: parameter, nu Oph b, nu Oph c).
    L = read('src_1904.03557', 'hip88048_arXiv.tex')
    t = table(L, r'\\label\{table:orb_par_stable\}', r'N-body')
    ms = get(table(L, r'\\label\{table:phys_param\}', r'\\end\{tabular\}'), 'Mass')[2]  # HB model (>99% probability)
    # The elements (incl. M0) refer to "the first observational epoch, which in our case is always JD 2451853.595".
    epoch = float(re.search(r'first observational epoch, which in our case is always JD\\,(' + NUM + ')',
                            ' '.join(L)).group(1))
    for k, comp in ((1, 'b'), (2, 'c')):
        P, M0 = get(t, '$P$')[k], get(t, '$M_0$')[k]
        # Tp from the mean anomaly: the star passed periastron M0/360 of a period before the epoch.
        tp = '%.2f' % (epoch - float(val(M0)) / 360.0 * float(val(P)))
        rows.append(dict(host='nu Oph', companion=comp, P_d=val(P), e_P_d=err(P), K_ms=val(get(t, '$K$')[k]),
                         e=val(get(t, '$e$')[k]), omega_deg=val(get(t, '$\\omega$')[k]), Tp_bjd=tp,
                         msini_mj=val(get(t, '$m \\sin i$')[k]), a_au=val(get(t, '$a$')[k]), mstar_msun=val(ms),
                         fit_type='RV-only', table_ref=REF_Q19, notes=NOTES[('nu Oph', comp)]))

    # Mitchell+2013: table "orbit" (columns: parameter, tau Gem, 91 Aqr); T0 in JD-2450000.
    L = read('src_1305.5107', 'planets.tex')
    t = table(L, r'\\label\{orbit\}', r'\\end\{tabular\}')
    ms = get(table(L, r'\\label\{param\}', r'\\end\{tabular\}'), 'Mass')[1]
    P = get(t, '$P$')[1]
    rows.append(dict(host='tau Gem', companion='b', P_d=val(P), e_P_d=err(P), K_ms=val(get(t, '$K_1$')[1]),
                     e=val(get(t, '$e$')[1]), omega_deg=val(get(t, '$\\omega$')[1]),
                     Tp_bjd=add_offset(get(t, '$T_0$')[1], 2450000),  # JD-2450000 -> JD
                     msini_mj=val(get(t, '$m_{2}\\sin(i)$')[1]), a_au=val(get(t, '$a_2$')[1]), mstar_msun=val(ms),
                     fit_type='RV-only', table_ref=REF_M13, notes=NOTES[('tau Gem', 'b')]))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
