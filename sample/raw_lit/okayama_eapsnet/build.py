"""Build okayama_eapsnet_orbits.csv from Omiya et al. 2009, PASJ 61, 825 (HD 119445 b; arXiv:0906.3762).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_0906.3762.gz (arXiv LaTeX source), unpacked into src_0906.3762/.
  - Orbit: table "Orbital parameters of HD 119445 b" in paper.tex (T is periastron time in JD).
  - Host mass: table "Stellar parameters of the host star HD 119445".
src_1101.4108.gz (Omiya et al. 2011, the program overview) is kept for reference only; it repeats the same orbit.
The first transcription was typed by hand from these tables; this script parses them and reproduces it.
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
OUT = os.path.join(HERE, 'okayama_eapsnet_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_0906.3762.gz': 'https://arxiv.org/e-print/0906.3762'}
ARCHIVES = {'src_0906.3762.gz': 'src_0906.3762'}

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


def cells(line):
    return [c.strip() for c in line.split('\\\\')[0].split('&')]


def err_nz(cell):
    """err() but a published error of 0 becomes blank."""
    e = err(cell)
    return e if e and float(e) != 0 else ''


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
TABLE_REF = ('arXiv:0906.3762 (Omiya+2009, PASJ 61, 825) Table "Orbital parameters of HD 119445 b" '
             '(Mstar Table "Stellar parameters")')
NOTE = ('BOAO + OAO (Korean-Japanese EAPS-Net) RVs, 36 obs. T is periastron time in JD (not BJD). '
        'Massive (3.9 Msun) G6III host at 289 pc. Not in Teng+2023 revisit list.')

NOTE_ERR = ' e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au: published 1-sigma errors from the orbit table.'


def main():
    get_inputs()
    with open(os.path.join(HERE, 'src_0906.3762', 'paper.tex'), encoding='latin-1') as f:
        T = f.read()
    orb = param_rows(table(T, r'\\caption\{Orbital parameters of HD 119445 b\}'))
    star = param_rows(table(T, r'\\caption\{Stellar parameters of the host star HD 119445\}'))
    P = pick(orb, r'^\$P\$')[0]
    write([dict(host='HD 119445', companion='b', P_d=val(P), e_P_d=err(P), K_ms=val(pick(orb, r'^\$K_\{1\}\$')[0]),
                e=val(pick(orb, r'^\$e\$')[0]), omega_deg=val(pick(orb, r'^\$\\omega\$')[0]),
                Tp_bjd=val(pick(orb, r'^\$T\$ \(JD\)')[0]),  # JD despite the column name
                msini_mj=val(pick(orb, r'^\$M_\{\\mathrm\{2\}\}')[0]), a_au=val(pick(orb, r'^\$a\$ \(AU\)')[0]),
                mstar_msun=val(pick(star, r'^\$M\$')[0]), fit_type='RV-only', table_ref=TABLE_REF, notes=NOTE + NOTE_ERR,
                e_K_ms=err_nz(pick(orb, r'^\$K_\{1\}\$')[0]), e_msini_mj=err_nz(pick(orb, r'^\$M_\{\\mathrm\{2\}\}')[0]),
                e_e=err_nz(pick(orb, r'^\$e\$')[0]), e_a_au=err_nz(pick(orb, r'^\$a\$ \(AU\)')[0]))])


if __name__ == '__main__':
    main()
