"""Build rickman_2019_orbits.csv from the arXiv LaTeX source of Rickman et al. 2019, A&A 625, A71
("The CORALIE survey for southern extrasolar planets XVIII. Three new massive planets and two low
mass brown dwarfs at separation larger than 5 AU", arXiv:1904.01573).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_1904.01573.src  arXiv e-print (gzipped tar). main.tex is read from the extracted copy
                        arxiv_src/main.tex if present, else directly from the tarball.
Tables used (labels in main.tex; all are transposed: one column per star or companion):
  table:1, table:2  (stellar parameters: M*)
  table:3, table:4  (orbital solutions: P [yr], K, e, omega, Tp [BJD-2450000], M sin i, a)
Selection: companions with M sin i >= 5 MJ, ordered by decreasing M sin i (4 rows, all a > 5 AU).
Output: rickman_2019_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_1904.01573.src': 'https://arxiv.org/e-print/1904.01573'}
ARCHIVE = 'arxiv_1904.01573.src'
TEX = 'main.tex'
OUT = 'rickman_2019_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

DAYS_PER_YEAR = 365.25   # P is given in years in Tables 3-4
TP_OFFSET = 2450000.0    # Tables 3-4 footnote: "T_P is shown in BJD-2,450,000"
OUT_A = 'OUT OF TARGET RANGE (a>5 AU)'


# ---------------------------------------------------------------- input handling
def fetch(download):
    for fn, url in SOURCES.items():
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            continue
        if not os.path.isdir(os.path.dirname(p)):
            os.makedirs(os.path.dirname(p))
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


def write_csv(rows):
    with open(os.path.join(HERE, OUT), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in COLS})
    print('wrote %s (%d rows)' % (OUT, len(rows)))


# ---------------------------------------------------------------- build
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--download', action='store_true', help='fetch missing inputs listed in SOURCES')
    args = ap.parse_args()
    fetch(args.download)
    src = read_tex(TEX)

    def transposed(label):
        """(column names, {row label: cells}) for the tabular of the table holding \\label{label}."""
        i = src.index('\\label{' + label + '}')
        b = src[src.index('\\begin{tabular}', i):src.index('\\end{tabular}', i)]
        b = re.sub(r'(?<!\\)%.*', '', b)  # drop LaTeX comments
        out = {}
        for r in split_rows(b):
            out[clean(r[0].strip().split('\n')[-1])] = r
        names = [clean(c) for c in out['Parameters'][2:]]
        return names, out

    mstar = {}
    for lab in ('table:1', 'table:2'):
        names, t = transposed(lab)
        for k, n in enumerate(names):
            mstar[n] = parse_val(t['M_{*}'][2 + k])[0]

    def yr(p):
        return round(p * DAYS_PER_YEAR, 2)  # years -> days

    rows = []
    for lab, tref in (('table:3', 'arXiv:1904.01573 Table 3 (mstar Tables 1-2)'),
                      ('table:4', 'arXiv:1904.01573 Table 4 (mstar Table 1)')):
        names, t = transposed(lab)
        for k, n in enumerate(names):
            col = 2 + k
            m = re.match(r'^(.*?)([b-h])$', n)  # Table 4 names companions ('HD 98649b'); Table 3 names hosts
            host, letter = (m.group(1).strip(), m.group(2)) if m else (n, 'b')
            P, eP = parse_val(t['P'][col])
            K, _ = parse_val(t['K'][col])
            e, _ = parse_val(t['e'][col])
            om, _ = parse_val(t['\\omega'][col])
            Tp, _ = parse_val(t['T_p'][col])
            ms, _ = parse_val(t[[x for x in t if x.startswith('M.') and 'sin i' in x][0]][col])
            a, _ = parse_val(t['a'][col])
            if ms < 5:
                continue
            pnote = 'P given as %g yr' % P
            if host == 'HD 13724':
                # the original table wrote out the asymmetric errors for this one row only
                pnote = 'P given as %g (+%g/-%g) yr' % tuple(abs(float(x)) for x in re.findall(NUM, clean(t['P'][col])))
            notes = ([OUT_A] if a > 5 else []) + [pnote + ', converted with 365.25 d/yr']
            if host == 'HD 13724':
                notes[-1] += '; mode of posterior'
            rows.append(dict(host=host, companion=host + ' ' + letter, P_d=yr(P),
                             e_P_d=yr(eP),  # asymmetric errors averaged, then converted to days
                             K_ms=K, e=e, omega_deg=om, Tp_bjd=Tp + TP_OFFSET, msini_mj=ms, a_au=a,
                             mstar_msun=mstar[host], fit_type='RV-only', table_ref=tref, notes='; '.join(notes)))
    rows.sort(key=lambda r: -r['msini_mj'])
    write_csv(rows)


if __name__ == '__main__':
    main()
