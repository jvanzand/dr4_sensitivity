"""Build ptps_orbits.csv from two Penn State-Torun Planet Search (PTPS/TAPAS) papers.

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs (arXiv LaTeX sources, unpacked into src_<id>/):
  - src_0906.1804.gz: Niedzielski et al. 2009, ApJ 707, 768 (BD+20 2457 b, c). Orbits in tab4.tex
    (T0 in MJD); the host mass is only a range (1.3-4.3 Msun), so mstar_msun is left empty.
  - src_2101.10410.gz: Niedzielski et al. 2021, A&A 648, A58 (TAPAS VII). Orbits in table_master-new.tex
    (T0 in MJD = JD - 2400000.5); host masses in the "Basic parameters" table of TAPAS-7-final.tex.
    HD 96992 b (1.14 MJ) is not kept.
T0 is converted from MJD to JD. The first transcription was typed by hand from these tables; this script parses
them and reproduces it.
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
OUT = os.path.join(HERE, 'ptps_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_0906.1804.gz': 'https://arxiv.org/e-print/0906.1804',
           'src_2101.10410.gz': 'https://arxiv.org/e-print/2101.10410'}
ARCHIVES = {'src_0906.1804.gz': 'src_0906.1804', 'src_2101.10410.gz': 'src_2101.10410'}

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


def tyc(name):
    """'TYC 0434-04538-1' -> 'TYC 434-4538-1'; other names unchanged."""
    m = re.match(r'TYC\s*0*(\d+)-0*(\d+)-(\d+)$', name)
    return 'TYC %s-%s-%s' % m.groups() if m else name


def mjd_to_jd(cell):
    return '%.2f' % (float(val(cell)) + 2400000.5)


# ---------------------------------------------------------------- per-row text (hand-written)
REF09 = 'arXiv:0906.1804 (Niedzielski+2009, ApJ 707, 768) Table 4 (stellar Table 1)'
REF21 = ('arXiv:2101.10410 (Niedzielski+2021, A&A 648, A58; TAPAS VII) Table TableKeplerian '
         '(Mstar: basic parameters table)')
T0NOTE = (' T0 given in MJD (converted to JD); the paper does not define T0 explicitly (PTPS convention is '
          'periastron passage): verify before use.')
NOTES = {
    'BD+20 2457 b': 'HET/HRS RVs, two-Keplerian fit. Host mass only given as range 1.3-4.3 Msun (metal-poor K2 II, '
                    '[Fe/H]=-1.0), so msini is model-dependent. Quirrenbach+2019 notes no stable configuration found '
                    '(Horner+2014, Trifonov+2014): reality doubtful. Likely in NASA Exoplanet Archive.',
    'BD+20 2457 c': 'Massive planet; see BD+20 2457 b notes (dynamical stability doubtful).',
    'HD 4760 b': 'HET/HRS + TNG/HARPS-N RVs (35 epochs); metal-poor giant ([Fe/H]=-0.91), jitter 64 m/s.',
    'BD+02 3313 b': 'HET/HRS + TNG/HARPS-N RVs (29 epochs). msini > 30 MJ so not in NASA Exoplanet Archive scope.',
    'TYC 434-4538-1 b': 'Massive planet (6.1 MJ); likely in NASA Exoplanet Archive.',
}
NOTE_ERR = (' e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au: published 1-sigma errors from the orbit table; asymmetric '
            'errors averaged; blank where none is printed (Niedzielski+2009 gives no error on M sin i or a).')
KEEP21 = ['HD 4760', 'BD+02 3313', 'TYC 434-4538-1']


def orbit(host, comp, col, rows, mstar, ref):
    c = lambda pat: pick(rows, pat)[col]
    return dict(host=host, companion=comp, P_d=val(c(r'^\$?P\$?[~ ]')), e_P_d=err(c(r'^\$?P\$?[~ ]')), K_ms=val(c(r'^\$?K')),
                e=val(c(r'^\$?e\$?$')), omega_deg=val(c(r'omega')), Tp_bjd=mjd_to_jd(c(r'T_0')),
                msini_mj=val(c(r'm_2\\sin')), a_au=val(c(r'^\$a\$')), mstar_msun=mstar, fit_type='RV-only',
                table_ref=ref, notes=NOTES['%s %s' % (host, comp)] + T0NOTE + NOTE_ERR,
                e_K_ms=err_nz(c(r'^\$?K')), e_msini_mj=err_nz(c(r'm_2\\sin')), e_e=err_nz(c(r'^\$?e\$?$')),
                e_a_au=err_nz(c(r'^\$a\$')))


def main():
    get_inputs()
    rows = []
    # Niedzielski+2009 Table 4: columns BD+20 2457 b, c.
    with open(os.path.join(HERE, 'src_0906.1804', 'tab4.tex'), encoding='utf-8') as f:
        t4 = param_rows(table(f.read(), r'\\begin\{tabular\}'))
    for col, comp in enumerate('bc'):
        rows.append(orbit('BD+20 2457', comp, col, t4, '', REF09))

    # TAPAS VII: column macros \starA, \starB, ... are defined in the main file.
    d = os.path.join(HERE, 'src_2101.10410')
    with open(os.path.join(d, 'TAPAS-7-final.tex'), encoding='utf-8') as f:
        main_tex = f.read()
    macro = {m.group(1): tyc(m.group(2).replace('~', ' ').strip())
             for m in re.finditer(r'\\def\\(star\w)\{(.*?)\\thinspace\}', main_tex)}
    with open(os.path.join(d, 'table_master-new.tex'), encoding='utf-8') as f:
        kep = param_rows(f.read())
    hosts = [macro[re.sub(r'\\scriptsize|\\|\s', '', c)] for c in pick(kep, r'^Parameter')]
    mstar = {tyc(k.strip()): val(v[4]) for k, v in
             param_rows(table(main_tex, r'\\caption\{Basic parameters of the program stars\. \}')).items()
             if k.strip() != 'Star'}
    for host in KEEP21:
        rows.append(orbit(host, 'b', hosts.index(host), kep, mstar[host], REF21))
    write(rows)


if __name__ == '__main__':
    main()
