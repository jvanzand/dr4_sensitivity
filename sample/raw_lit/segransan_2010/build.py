"""Build segransan_2010_orbits.csv from the arXiv LaTeX source of Segransan et al. 2010, A&A 511, A45
("The CORALIE survey for southern extrasolar planets XVI", arXiv:0908.1479).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_0908.1479.src  arXiv e-print (gzipped tar). coralieXVILongPeriodPlanets4AstroPh.tex is read from
                       the extracted copy in arxiv_src/ if present, else directly from the tarball.
Tables used (labels in the .tex; all three are transposed: one column per planet or star):
  table:stellarparameters                      (Table 1: M*)
  table:orbitalelementshd147018                (Table 2: HD 147018 b, c)
  table:orbitalelementshd171238andhd204313     (Table 3: HD 171238 b, HD 204313 b)
Selection: companions with m sin i >= 5 MJ (only HD 147018 c).
Output: segransan_2010_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_0908.1479.src': 'https://arxiv.org/e-print/0908.1479'}
ARCHIVE = 'arxiv_0908.1479.src'
TEX = 'coralieXVILongPeriodPlanets4AstroPh.tex'
OUT = 'segransan_2010_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

# Table 2 gives T0 = 55301 +/- 22 for HD 147018 c in "JD* = JD-2450000", which is inconsistent with
# that definition (HD 147018 b has T0 = 4459.49 in the same table; 55301 is probably 5301, i.e.
# JD 2455301). Tp is therefore left blank, and the published value is quoted in `notes`.
TP_BLANK = True


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
        """{row label: list of cells} for the tabular of the table holding \\label{label}."""
        i = src.index('\\label{' + label + '}')
        b = src[src.index('\\begin{tabular}', i):src.index('\\end{tabular}', i)]
        out = {}
        for r in split_rows(b):
            lab = re.sub(r'^\[[^\]]*\]', '', r[0].strip()).strip()  # drop '[1mm]' left from '\\[1mm]'
            out[clean(lab)] = r
        return out

    st = transposed('table:stellarparameters')
    # Table 1 columns: label, unit, HD 147018, HD 171238, HD 204313
    stars = ['HD 147018', 'HD 171238', 'HD 204313']
    mstar = {s: parse_val(st['M_{\\star}^{(1)}'][2 + k])[0] for k, s in enumerate(stars)}

    rows = []
    for label, comps in (('table:orbitalelementshd147018', ['HD 147018 b', 'HD 147018 c']),
                         ('table:orbitalelementshd171238andhd204313', ['HD 171238 b', 'HD 204313 b'])):
        t = transposed(label)
        for k, comp in enumerate(comps):
            col = 2 + k  # columns: label, unit, planet 1, planet 2
            host = comp[:-2]
            P, eP = parse_val(t['P'][col])
            K, _ = parse_val(t['K'][col])
            e, _ = parse_val(t['e'][col])
            om, _ = parse_val(t['\\omega'][col])
            T0, eT0 = parse_val(t['T_{0}'][col])
            ms, _ = parse_val(t['m_{p} \\sin{i}'][col])
            a, _ = parse_val(t['a'][col])
            if ms < 5:
                continue
            notes = ('Tp left blank: table gives T0=%g+/-%g in JD*=JD-2450000, inconsistent with that definition '
                     '(probably 2455301 i.e. a typo); likely already in NASA Exoplanet Archive. '
                     # descriptive note; the other three companions are parsed above and dropped (Msini < 5 MJ)
                     'Other companions in paper (HD 147018 b, HD 171238 b, HD 204313 b) have Msini 2.1-4.1 MJ'
                     % (T0, eT0))
            rows.append(dict(host=host, companion=comp, P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                             Tp_bjd=None if TP_BLANK else T0 + 2450000.0,  # JD* = JD - 2450000 (see TP_BLANK)
                             msini_mj=ms, a_au=a, mstar_msun=mstar[host], fit_type='RV-only',
                             table_ref='arXiv:0908.1479 Table 2 (mstar Table 1)', notes=notes))
    write_csv(rows)


if __name__ == '__main__':
    main()
