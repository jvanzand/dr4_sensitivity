"""Build sahlmann_2011_orbits.csv from the arXiv LaTeX source of Sahlmann et al. 2011, A&A 525, A95
("Search for brown-dwarf companions of stars", arXiv:1009.5991).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_1009.5991.src  arXiv e-print (gzipped tar). The LaTeX file 15427.tex is read from the
                       extracted copy arxiv_src/15427.tex if present, else directly from the tarball.
Tables used (labels in 15427.tex):
  tab:KeplerOrbits  (Table 3: P, e, K, T0, omega)
  tab:KeplerParams  (Table 4: a_rel, M2 sin i)
  tab:stellarChar   (Table 2: M1)
  tab:mass2sigma    (Table 7: Hipparcos-astrometry true masses)
  tab:target2sigma  (Table 8: spectroscopic elements of literature targets)
Output: sahlmann_2011_orbits.csv
Error columns (appended; published 1-sigma, asymmetric errors averaged, a published 0 -> blank):
e_K_ms, e_e from Table 3 (Table 8 for the literature orbits); e_msini_mj, e_a_au from Table 4 (a_rel).
Table 8 gives no M2 sin i error, so e_msini_mj is blank for those rows.
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_1009.5991.src': 'https://arxiv.org/e-print/1009.5991'}
ARCHIVE = 'arxiv_1009.5991.src'
TEX = '15427.tex'
OUT = 'sahlmann_2011_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes",
        "e_K_ms", "e_msini_mj", "e_e", "e_a_au"]


# ---------------------------------------------------------------- input handling
def fetch(download):
    for fn, url in SOURCES.items():
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            continue
        if not download:
            sys.exit('missing input %s (from %s); rerun with --download' % (fn, url))
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as r:
            data = r.read()
        with open(p + '.part', 'wb') as f:
            f.write(data)
        os.replace(p + '.part', p)


def read_tex(member):
    """Text of one file of the arXiv source: extracted copy if present, else from the tarball."""
    p = os.path.join(HERE, 'arxiv_src', member)
    if os.path.exists(p):
        with open(p, 'rb') as f:
            raw = f.read()
    else:
        with tarfile.open(os.path.join(HERE, ARCHIVE)) as t:
            raw = t.extractfile(member).read()
    try:
        return raw.decode('utf-8')
    except UnicodeDecodeError:
        return raw.decode('latin-1')


# ---------------------------------------------------------------- LaTeX helpers
def clean(s):
    s = re.sub(r'\\tablefootmark\{[^}]*\}', '', s)
    s = re.sub(r'\\vspace\{[^}]*\}', '', s)
    s = re.sub(r'\\(?:begin|end)\{math\}', '', s)
    s = s.replace('$', '').replace('\\,', '').replace('~', ' ')
    s = s.replace('\\mathrm', '').replace('\\textbf', '')
    return s.strip()


NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def parse_val(s):
    """Return (value, symmetric_error) from a LaTeX cell; asymmetric errors are averaged."""
    s = clean(s)
    s = s.replace(' ', '').replace('{', '').replace('}', '')
    if s in ('', '---', '\\nodata', '\\cdots', '...', '-', '--'):
        return None, None
    m = re.match(r'^(%s)\^\+?(%s)_-?(%s)' % (NUM, NUM, NUM), s)
    if m:
        v, up, lo = m.groups()
        return float(v), (abs(float(up)) + abs(float(lo))) / 2
    m = re.match(r'^(%s)_-?(%s)\^\+?(%s)' % (NUM, NUM, NUM), s)
    if m:
        v, lo, up = m.groups()
        return float(v), (abs(float(up)) + abs(float(lo))) / 2
    m = re.match(r'^(%s)\\pm(%s)' % (NUM, NUM), s)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.match(r'^(%s)' % NUM, s)
    if m:
        return float(m.group(1)), None
    return None, None


def pos(x):
    """Published error -> keep only if > 0 (a published 0 is treated as no error)."""
    return x if x else None


def split_rows(body):
    rows = []
    for r in re.split(r'\\\\', body):
        r = '\n'.join(l for l in r.split('\n') if not l.strip().startswith('%'))
        r = re.sub(r'\\(hline|midrule|toprule|bottomrule)', '', r).strip()
        if r:
            rows.append([c.strip() for c in r.split('&')])
    return rows


def norm_name(n):
    n = clean(n).replace('\\object{', '').replace('{\\footnotesize HD}', 'HD').strip(' {}')
    n = re.sub(r'\s+', ' ', n)
    m = re.match(r'^(HD|HIP|HR|GJ|TYC)\s*(.+)$', n)
    if m:
        return '%s %s' % (m.group(1), m.group(2).strip())
    m = re.match(r'^(BD|CD)([+-])(\d{2})\s*0*(\d+)$', n)
    if m:
        return '%s%s%s %d' % (m.group(1), m.group(2), m.group(3), int(m.group(4)))
    return n


def fmt(x):
    if x is None:
        return ''
    if isinstance(x, float):
        return '%.10g' % x
    return str(x)


# ---------------------------------------------------------------- build
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--download', action='store_true', help='fetch missing inputs listed in SOURCES')
    args = ap.parse_args()
    fetch(args.download)
    src = read_tex(TEX)

    def table_body(label):
        """Rows of the tabular after the \\label, starting after the header's closing \\hline."""
        i = src.index('\\label{' + label + '}')
        j = src.index('\\hline', src.index('\\hline\\hline', i) + 12)
        k = src.index('\\end{tabular}', i)
        return src[j:k]

    orb = split_rows(table_body('tab:KeplerOrbits'))
    par = split_rows(table_body('tab:KeplerParams'))

    # Table 2: stellar mass M1 (column 8)
    mstar = {}
    for r in split_rows(table_body('tab:stellarChar')):
        if len(r) < 9:
            continue
        mstar[norm_name(r[0])] = parse_val(r[7])[0]

    # Table 7: Hipparcos true masses M2 (column 7) and orbit significance (column 12)
    m2 = {}
    for r in split_rows(table_body('tab:mass2sigma')):
        if len(r) < 12:
            continue
        fm = re.search(r'tablefootmark\{(\w)\}', r[0])
        v, e = parse_val(r[6])
        m2[(norm_name(clean(r[0])), fm.group(1) if fm else None)] = (v, e, clean(r[11]))

    def key(cell):
        fm = re.search(r'tablefootmark\{(\w)\}', cell)
        return norm_name(cell.split('\\tablefootmark')[0]), (fm.group(1) if fm else None)

    params = {}
    for r in par:
        if len(r) < 12:
            continue
        params[key(r[0])] = r

    rows = []
    for r in orb:
        if len(r) < 8:
            continue
        name, fm = key(r[0])
        # Table 3 and Table 4 use different footnote letters for the multiple solutions:
        # HD 43848 (Table 3 g/h -> Table 4 f/g); HD 4747 and HD 211847 (Table 3 d/e/f -> Table 4 b/c/d).
        pmap = {('HD 43848', 'g'): 'f', ('HD 43848', 'h'): 'g'}
        if name in ('HD 4747', 'HD 211847'):
            pmap = {(name, 'd'): 'b', (name, 'e'): 'c', (name, 'f'): 'd'}
        pr = params.get((name, pmap.get((name, fm), fm)))
        # Keep only the adopted solutions
        if name == 'HD 43848' and fm == 'h':
            continue  # formal high-e solution; mentioned in notes
        if name in ('HD 4747', 'HD 211847') and fm != 'e':
            continue  # keep the best-fit solution only (d/f are the P range limits)
        P, eP = parse_val(r[3])
        e, ee = parse_val(r[4])
        K, eK = parse_val(r[5])
        T0, _ = parse_val(r[6])
        om, _ = parse_val(r[7])
        a, ea = parse_val(pr[10]) if pr else (None, None)
        ms, ems = parse_val(pr[11]) if pr else (None, None)
        notes = []
        if name == 'HD 43848':
            # high-e values typed from Tables 3-4, rows HD 43848^h (Table 3) / HD 43848^g (Table 4)
            notes.append('adopted low-e solution; formal high-e alternative (e=0.885, K=1859 m/s, Msini=50.2 MJ) '
                         'also in Table 3/4')
        if name in ('HD 4747', 'HD 211847'):
            # P ranges typed from Table 3 (rows ^d and ^f give the lower/upper P limits)
            notes.append('period not covered; best-fit solution with P range given in paper '
                         '(HD 4747: 7900-29000 d; HD 211847: 3750-100000 d)')
        mk = m2.get((name, {'g': 'a', 'h': 'b'}.get(fm) if name == 'HD 43848' else None))
        if mk:
            notes.append('Hipparcos-astrometry true mass M2=%s +/- %.4g MJ (orbit significance %s%%, Sahlmann Table 7)'
                         % (mk[0], mk[1], mk[2]))
        if name == 'HIP 103019':
            notes.append('HARPS (not CORALIE) target')
        if name == 'HD 17289':
            notes.append('companion contaminates CCF (SB2-like RV bias); M2=0.52 Msun')
        if (a is not None and a > 5) or (ms is not None and ms > 80):
            notes.insert(0, 'OUT OF TARGET RANGE (a>5 AU)' if a and a > 5 else 'OUT OF TARGET RANGE (Msini>80 MJ)')
        rows.append(dict(host=name, companion=name + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                         # T0 in Table 3 is BJD - 2400000 -> full BJD
                         Tp_bjd=(T0 + 2400000.0) if T0 is not None else None,
                         e_K_ms=pos(eK), e_msini_mj=pos(ems), e_e=pos(ee), e_a_au=pos(ea),
                         msini_mj=ms, a_au=a, mstar_msun=mstar.get(name), fit_type='RV-only',
                         table_ref='arXiv:1009.5991 Tables 3-4 (stellar mass Table 2)',
                         notes='; '.join(notes)))

    # Table 8: literature targets with spectroscopic elements (reference numbers -> papers, Table 8 notes)
    refs = {'1': 'Minniti+2009', '2': 'Nidever+2002', '3': 'Perrier+2003', '4': 'Wright+2009', '5': 'Jenkins+2009'}
    for r in split_rows(table_body('tab:target2sigma')):
        if len(r) < 10:
            continue
        nm = clean(r[0])
        host = norm_name(re.sub(r'c$', '', nm))  # 'HD 168443c' -> host HD 168443, companion c
        P, eP = parse_val(r[4])
        e, ee = parse_val(r[5])
        K, eK = parse_val(r[6])
        T0, _ = parse_val(r[7])
        om, _ = parse_val(r[8])
        ms = parse_val(r[3])[0]
        notes = ['orbit compiled by Sahlmann+2011 from %s; not a new fit' % refs[clean(r[9])]]
        if ms < 5:
            notes.insert(0, 'BELOW 5 MJ (Msini=3.6 MJ; listed because Hipparcos mass is substellar)')
        mk = m2.get((host, None))
        if mk:
            notes.append('Hipparcos-astrometry M2=%s +/- %.4g MJ (significance %s%%)' % (mk[0], mk[1], mk[2]))
        rows.append(dict(host=host, companion=host + (' c' if nm.endswith('c') else ' b'), P_d=P, e_P_d=eP,
                         K_ms=K, e=e, omega_deg=om, e_K_ms=pos(eK), e_e=pos(ee),
                         Tp_bjd=T0 + 2400000.0,  # T0 in Table 8 is BJD - 2400000 -> full BJD
                         msini_mj=ms, mstar_msun=parse_val(r[2])[0],
                         fit_type='RV-only', table_ref='arXiv:1009.5991 Table 8 (literature elements)',
                         notes='; '.join(notes)))

    with open(os.path.join(HERE, OUT), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in COLS})
    print('wrote %s (%d rows)' % (OUT, len(rows)))


if __name__ == '__main__':
    main()
