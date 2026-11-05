"""Build express_jones_orbits.csv from two EXPRESS papers (Jones et al.; FEROS + CHIRON, southern survey).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs (arXiv LaTeX sources, unpacked into src_<id>/):
  - src_1409.7429.gz: Jones et al. 2015, A&A 573, A3. HIP 97233 b from the "Orbital parameters" table
    (T_P in JD-2455000) and host mass from "Stellar parameters" in new_paper4_ref.tex. HIP 67851 b (1.4 MJ) is
    not kept.
  - src_1612.06252.gz: Jones et al. 2017, A&A 602, A58. HIP 67537 b from "Stellar properties and orbital
    elements" in paper7_v5.tex (T_P in JD-2450000).
T_P is converted to JD. The first transcription was typed by hand from these tables; this script parses them
and reproduces it.
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
OUT = os.path.join(HERE, 'express_jones_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_1409.7429.gz': 'https://arxiv.org/e-print/1409.7429',
           'src_1612.06252.gz': 'https://arxiv.org/e-print/1612.06252'}
ARCHIVES = {'src_1409.7429.gz': 'src_1409.7429', 'src_1612.06252.gz': 'src_1612.06252'}

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
    """Reduced date + offset (e.g. JD-2450000 -> JD), keeping the number of printed decimals."""
    v = val(cell)
    nd = len(v.split('.')[1]) if '.' in v else 0
    return '%.*f' % (nd, float(v) + off)


def block(lines, start_pat, end_pat):
    """Lines after the first line matching start_pat, up to the next line matching end_pat."""
    i0 = next(i for i, l in enumerate(lines) if re.search(start_pat, l))
    out = []
    for l in lines[i0 + 1:]:
        if re.search(end_pat, l):
            break
        out.append(l)
    return out


def param_rows(text):
    """Rows of a 'Parameter & value(s)' table as {label: [cells]}; rows may span lines. Comments, rules and
    vertical spacing are dropped."""
    text = '\n'.join(re.sub(r'(?<!\\)%.*', '', l) for l in text.split('\n'))
    out = {}
    for r in text.split('\\\\'):
        r = re.sub(r'\\(vspace|hspace)\{[^}]*\}|\\begin\{tabular\}\{[^}]*\}|\\renewcommand\{[^}]*\}\{[^}]*\}'
                   r'|\\hline|\\normalsize', '', r).strip()
        if '&' in r:
            c = [x.strip() for x in r.split('&')]
            out[c[0]] = c[1:]
    return out


def pick(rows, pat):
    """Cells of the first row whose label matches the regex pat."""
    return next(v for k, v in rows.items() if re.search(pat, k))


def table(text, start_pat, end_pat=r'\\end\{tabular\}'):
    """Text between the first match of start_pat and the next match of end_pat."""
    i0 = re.search(start_pat, text).end()
    return text[i0:i0 + re.search(end_pat, text[i0:]).start()]


def write(rows):
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


# ---------------------------------------------------------------- per-row text (hand-written)
REF15 = ('arXiv:1409.7429 (Jones+2015, A&A 573, A3) Table "Orbital parameters of HIP 67851 b and HIP 97233 b" '
         '(Mstar Table "Stellar parameters")')
REF17 = 'arXiv:1612.06252 (Jones+2017, A&A 602, A58) Table "Stellar properties and orbital elements"'
NOTE15 = ('EXPRESS (FEROS + CHIRON; southern survey). Tp JD-2455000 in paper (JD). '
          'Likely in NASA Exoplanet Archive.')
NOTE17 = ('EXPRESS (FEROS + CHIRON). Massive planet/BD boundary (11.1 +0.4/-1.1 MJ); a=4.91 AU near the 5 AU '
          'limit. Likely in NASA Exoplanet Archive.')


def orbit(host, rows, col, tp_off, mstar, ref, note):
    c = lambda pat: pick(rows, pat)[col]
    return dict(host=host, companion='b', P_d=val(c(r'^P \(days\)')), e_P_d=err(c(r'^P \(days\)')),
                K_ms=val(c(r'^K ')), e=val(c(r'^\$?e\$?$')), omega_deg=val(c(r'omega')),
                Tp_bjd=add_offset(c(r'^T\$_'), tp_off),  # JD despite the column name
                msini_mj=val(c(r'^M\\,sin\$i\$|^\\mplan')), a_au=val(c(r'^\$?a\$? \(AU\)')), mstar_msun=mstar,
                fit_type='RV-only', table_ref=ref, notes=note)


def main():
    get_inputs()
    rows = []
    with open(os.path.join(HERE, 'src_1409.7429', 'new_paper4_ref.tex'), encoding='utf-8') as f:
        T = f.read()
    orb = param_rows(table(T, r'\\label\{HIP67851_orb_par\}\}'))
    star = param_rows(table(T, r'\\label\{stellar_par\}\}'))
    col = [re.sub(r'\\,|\{\\it b\}|\s', '', c) for c in orb['']].index('HIP97233')
    rows.append(orbit('HIP 97233', orb, col, 2455000, val(pick(star, r'^M\$_\\star')[col]), REF15, NOTE15))

    with open(os.path.join(HERE, 'src_1612.06252', 'paper7_v5.tex'), encoding='utf-8') as f:
        T = f.read()
    tab = param_rows(table(T, r'\\label\{orb_par\}\}'))
    rows.append(orbit('HIP 67537', tab, 0, 2450000, val(pick(tab, r'^M\$_\\star')[0]), REF17, NOTE17))
    write(rows)


if __name__ == '__main__':
    main()
