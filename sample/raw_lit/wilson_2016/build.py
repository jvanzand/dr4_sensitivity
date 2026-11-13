"""Build wilson_2016_orbits.csv from Wilson et al. 2016, A&A 588, A144 (SOPHIE IX; arXiv:1602.02749).

Run:  python -I build.py [--download]   (from anywhere; paths are relative to this file)
Inputs: src_1602.02749.gz (arXiv LaTeX source), unpacked into src_1602.02749/.
  - Orbits: deluxetable tab:companion_properties ("Table 3") in rv.tex.
  - Host masses: table tab:stellar_properties ("Table 2") in rv.tex.
HD 50761 B (m2 sin i = 145 MJ) is excluded.

The last two table rows (HD 160508 b, BD+24 4697 b) were missed by the first transcription (its script read a
fixed line window, rv.tex lines 245-270, that stops just before them) and were appended by hand after an audit.
Error columns (appended; published 1-sigma from the same Table 3 cell as the value, asymmetric errors averaged,
a published 0 -> blank): e_K_ms, e_msini_mj (m2 sin i in M_Jup as printed), e_e, e_a_au.
They are parsed here like the other rows; only their table_ref and notes are hand-written (AUDIT_ROWS).
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
OUT = os.path.join(HERE, 'wilson_2016_orbits.csv')

# Raw inputs (local file name -> URL) and the directory each archive is unpacked into.
SOURCES = {'src_1602.02749.gz': 'https://arxiv.org/e-print/1602.02749'}
ARCHIVES = {'src_1602.02749.gz': 'src_1602.02749'}

COLS = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp_bjd', 'msini_mj',
        'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes',
        'e_K_ms', 'e_msini_mj', 'e_e', 'e_a_au']


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


def perr(cell):
    """err(), but a published error of 0 becomes blank."""
    e = err(cell)
    return e if e and float(e) > 0 else ''


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
NOTE = ('SOPHIE RVs; single-Keplerian fit; T0 is periastron time (BJD-2400000 in paper, converted). '
        'Separate Hipparcos IAD analysis does not alter RV parameters.')
NOTES_EXTRA = {  # from the paper's discussion of each object (rv.tex Sects. 5-6)
    'HD 6664': 'Wilson+2016 class: low-mass binary; Hipparcos astrometric orbit gives M2=0.35+/-0.03 Msun (stellar). msini 80+/-5.',
    'BD+73 0275': 'Wilson+2016 class: low-mass binary; Hipparcos astrometric orbit gives M2=0.20+/-0.03 Msun (stellar).',
    'HD 51813': 'Wilson+2016 class: low-mass binary; Hipparcos astrometric orbit gives M2=0.27+/-0.07 Msun (stellar).',
    'HD 56709': 'Wilson+2016 class: low-mass binary candidate; tentative Hipparcos orbit M2=0.54+/-0.12 Msun. msini 72+/-4.',
    'HD 134169': 'Wilson+2016 class: low-mass binary; msini 83+/-4 (just above 80 MJ); no Hipparcos detection (m2,max 1.03 Msun).',
    'HD 77065': 'Includes Latham+2002 RVs (as G9-42).',
    'HD 134113': 'Includes Latham+2002 RVs (as G66-65). Mstar from Sousa+2011.',
}
TABLE_REF = 'arXiv:1602.02749 Table 3 (Mstar Table 2)'
# Rows appended by hand after an audit (missed in the first transcription). Values are parsed from Table 3
# (rv.tex lines 271-272) like the rest; these table_ref/notes strings are kept exactly as hand-entered.
# (HD 160508's Mstar 1.25 is in both Table 2 and the appendix compilation Table A.1, hence its table_ref.)
AUDIT_ROWS = {
    'HD 160508': dict(table_ref='arXiv:1602.02749 Table 3 (Mstar Tables 2 and A.1)',
                      notes='Added after audit (missed in first transcription). SOPHIE RVs; T0 is periastron time '
                            '(BJD-2400000 in paper, converted). msini 48+/-3.'),
    'BD+24 4697': dict(table_ref='arXiv:1602.02749 Table 3 (Mstar Table 2)',
                       notes='Added after audit (missed in first transcription). SOPHIE RVs; T0 is periastron time '
                             '(BJD-2400000 in paper, converted). msini 53+/-3.'),
}
SKIP = {'HD 50761'}  # m2 sin i = 145 MJ: stellar


def main():
    get_inputs()
    with open(os.path.join(HERE, 'src_1602.02749', 'rv.tex'), encoding='utf-8') as f:
        L = f.read().split('\n')

    # Table 2 (stellar properties): mass in M_sun is the 10th column.
    mstar = {}
    for l in block(L, r'\\label\{tab:stellar_properties\}', r'\\end\{tabular\}'):
        c = cells(l)
        if len(c) >= 10 and '\\pm' in l:
            mstar[star_name(c[0])] = val(c[9])

    # Table 3 (companion properties). Columns: name, P [d], e, omega [deg], K [m/s], T0 [BJD-2400000],
    # gamma, gamma+, O-C, O-C+, a [AU], m2 sin i [MJup].
    rows = []
    for l in block(L, r'\\label\{tab:companion_properties\}', r'\\enddata'):
        if l.strip().startswith('%') or '\\pm' not in l:
            continue
        c = cells(l)
        name = star_name(c[0])
        host, comp = name[:-1].strip(), name[-1]
        if host in SKIP:
            continue
        if host in AUDIT_ROWS:
            ref, note = AUDIT_ROWS[host]['table_ref'], AUDIT_ROWS[host]['notes']
        else:
            ref, note = TABLE_REF, NOTE + (' ' + NOTES_EXTRA[host] if host in NOTES_EXTRA else '')
        rows.append(dict(host=bd(host), companion=comp, P_d=val(c[1]), e_P_d=err(c[1]), K_ms=val(c[4]), e=val(c[2]),
                         omega_deg=val(c[3]),
                         Tp_bjd=add_offset(c[5], 2400000),  # T0 is given as BJD-2400000
                         msini_mj=val(c[11]), a_au=val(c[10]), mstar_msun=mstar.get(host, ''),
                         fit_type='RV-only', table_ref=ref, notes=note,
                         e_K_ms=perr(c[4]), e_msini_mj=perr(c[11]), e_e=perr(c[2]), e_a_au=perr(c[10])))

    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
