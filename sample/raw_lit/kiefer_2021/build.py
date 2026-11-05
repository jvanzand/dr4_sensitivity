"""Build kiefer_2021_orbits.csv from the arXiv LaTeX source of Kiefer et al. 2021, A&A 645, A7
("Determining the true mass of radial-velocity exoplanets with Gaia. 9 planet candidates in the
brown-dwarf/stellar regime and 27 confirmed planets", arXiv:2009.14164).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_2009.14164.src  arXiv e-print (gzipped tar). exoplanet_mass_gaia.tex is read from the extracted
                        copy arxiv_src/exoplanet_mass_gaia.tex if present, else directly from the tarball.
Tables used (labels in exoplanet_mass_gaia.tex):
  tab:list_planets         (Table 3: 29 detection-sample companions; RV elements compiled from exoplanets.org)
  tab:list_planets_nondet  (Table C.1 longtable: 227 non-detection-sample companions)
  tab:results_det          (Table 4: GASTON true masses for the detection sample, quoted in notes)
Selection: companions with m sin i >= 5 MJ (34 rows), detection sample first, in table order.
Output: kiefer_2021_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_2009.14164.src': 'https://arxiv.org/e-print/2009.14164'}
ARCHIVE = 'arxiv_2009.14164.src'
TEX = 'exoplanet_mass_gaia.tex'
OUT = 'kiefer_2021_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

TP_OFFSET = 2450000.0  # Tp column is labelled "-2,450,000" (JD scale as compiled from exoplanets.org)


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
def flags(ms, a=None, P=None):
    n = []
    if ms is not None and ms < 5:
        n.append('BELOW 5 MJ')
    if ms is not None and ms > 80:
        n.append('OUT OF TARGET RANGE (Msini>80 MJ)')
    if a is not None and a > 5:
        n.append('OUT OF TARGET RANGE (a>5 AU)')
    elif a is None and P is not None and P > 4000:
        n.append('LIKELY OUT OF TARGET RANGE (P>4000 d, a not given)')
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--download', action='store_true', help='fetch missing inputs listed in SOURCES')
    args = ap.parse_args()
    fetch(args.download)
    src = read_tex(TEX)

    rows = []
    # Columns: Name, P, m sin i, K, e, omega, Tp, M*, a_p, plx, a* sin i, Transit flag, Drift flag, ...
    for lab, tref, end in (('tab:list_planets}', 'Table 3 (detection sample)', '\\end{tabular}'),
                           ('tab:list_planets_nondet}', 'Table C.1 (non-detection sample)', '\\end{longtable}')):
        i = src.index('\\caption{\\label{' + lab)
        blk = src[i:src.index(end, i)]
        for r in split_rows(blk):
            if len(r) < 13 or '$P$' in r[1] or '(days)' in r[1] or 'caption' in r[0] or 'multicolumn' in r[0]:
                continue  # header / caption / continuation rows
            name = clean(r[0]).replace('\\endhead', '').replace('\\endfirsthead', '').strip()
            if not name:
                continue
            P, eP = parse_val(r[1])
            ms, _ = parse_val(r[2])
            K, _ = parse_val(r[3])
            e, _ = parse_val(r[4])
            om, _ = parse_val(r[5])
            Tp, _ = parse_val(r[6])
            Ms, _ = parse_val(r[7])
            a, _ = parse_val(r[8])
            if ms is None or ms < 5:
                continue
            m = re.match(r'^(.*?)\s+([b-h])$', name)
            host = norm_name(re.sub(r'\s+', ' ', m.group(1))) if m else name
            host = re.sub(r'^BD ([+-])(\d+) (\d+)$', r'BD\1\2 \3', host)  # 'BD +20 2457' -> 'BD+20 2457'
            # (the companion keeps the table's spelling, e.g. 'BD +20 2457 b')
            n = flags(ms, a, P) + ['RV elements compiled by Kiefer+2021 from exoplanets.org (not new fits); '
                                   'transit flag=' + clean(r[11])]
            rows.append(dict(host=host, companion=name, P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                             Tp_bjd=Tp + TP_OFFSET if Tp is not None else None,  # Tp - 2,450,000 -> full JD
                             msini_mj=ms, a_au=a, mstar_msun=Ms,
                             fit_type='RV-only', table_ref='arXiv:2009.14164 ' + tref, notes='; '.join(n)))

    # GASTON true masses (1-sigma column of M_c,true) for the detection sample (Table 4)
    res = src[src.index('\\caption{\\label{tab:results_det}'):]
    res = res[:res.index('\\end{tabular}')]
    gm = {}
    for r in split_rows(res):
        if len(r) >= 13 and re.search(r'\d', r[1] or ''):
            nm = clean(re.sub(r'\\tablefootmark\{[^}]*\}', '', r[0])).strip()
            gm[nm] = clean(r[9]).replace(' ', '')
    for row in rows:
        if row['companion'] in gm:
            row['notes'] += ('; GASTON (Gaia DR1 excess noise) mass %s MJ (1-sigma; Kiefer+2021 Table 4)'
                             % gm[row['companion']])
    write_csv(rows)


if __name__ == '__main__':
    main()
