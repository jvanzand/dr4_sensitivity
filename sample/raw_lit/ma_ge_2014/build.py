"""Build ma_ge_2014_orbits.csv from the arXiv LaTeX source of Ma & Ge 2014, MNRAS 439, 2781
("Statistical Properties of Brown Dwarf Companions: Implications for Different Formation Mechanisms",
arXiv:1303.6442; the arXiv version is the Feb 2013 preprint).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_1303.6442.src  arXiv e-print (gzipped tar). twopopulation_Feb18.tex is read from the extracted
                       copy arxiv_src/twopopulation_Feb18.tex if present, else directly from the tarball.
Table used: deluxetable tab:catalog (Table 1, \\startdata ... \\enddata) and its \\tablerefs
(reference numbers -> citation keys). Columns: Object, Mc, Mc sin i, Period, e, M*, Teff, log g,
[Fe/H], References. No K, omega or Tp. All rows are kept (compilation; a pointer list).
Uncertainties: the table prints +/- errors on Mc sin i, Period and e (asymmetric ones are averaged; a printed
0.000 or a missing error is blank). e_msini_mj is the error on the Mc sin i cell (the msini_mj value; for HD 30501 that
cell is the Hipparcos true mass, see below), e_e the error on e. No errors on K or a (not tabulated).
Output: ma_ge_2014_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_1303.6442.src': 'https://arxiv.org/e-print/1303.6442'}
ARCHIVE = 'arxiv_1303.6442.src'
TEX = 'twopopulation_Feb18.tex'
OUT = 'ma_ge_2014_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes", "e_K_ms", "e_msini_mj", "e_e", "e_a_au"]

# Hand-written host-name normalizations (paper spelling -> SIMBAD-style name)
HOST_FIX = {'11Com': '11 Com', 'BD+202457': 'BD+20 2457', 'TYC 3130-160': 'TYC 3130-160-1',
            'COROT-15': 'CoRoT-15', 'Corot-3': 'CoRoT-3', 'NGC 2423-3': 'NGC 2423 3',
            'NGC 4349-127': 'NGC 4349 127'}

# Known oddity (kept as published): HD 30501 has Mc sin i = 90 +/- 12 MJ in this table, but
# Sahlmann+2011 (ref. 1) give M2 sin i = 62.3 MJ; 90 +/- 12 is their Hipparcos true mass (89.6 +/- 12.4).


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

    blk = src[src.index('\\startdata') + 10:src.index('\\enddata')]
    refs = src[src.index('\\tablerefs{'):]
    refs = refs[:refs.index('\\end{deluxetable}')]
    refmap = dict(re.findall(r'(\d+)\s*\\citet\{(\w+)\}', refs))  # '1' -> 'sahlmann11a'

    rows = []
    for r in split_rows(blk):
        if len(r) < 10:
            continue
        nm = clean(r[0]).replace('$-$', '-').strip()
        Mc, _ = parse_val(r[1])
        ms, e_ms = parse_val(r[2])
        P, eP = parse_val(r[3])
        e, e_e = parse_val(r[4])
        Ms, _ = parse_val(r[5])
        rr = [x.strip() for x in clean(r[9]).split(',') if x.strip()]
        m = re.match(r'^(.*?)(\s?[bc])$', nm)  # 'KOI-423b' / 'HD 38529 c' -> host + letter
        comp = nm if m else nm + ' b'
        host = m.group(1).strip() if m else nm
        host = norm_name(HOST_FIX.get(host, host))
        if host == 'HD 168443':
            comp = 'HD 168443 c'
        elif comp == host + ' b' or comp == nm + ' b':
            comp = host + ' b'
        else:
            comp = host + ' ' + nm[-1]
        n = []
        if Mc is not None:
            n.append('true mass Mc=' + re.sub(r'\s+', '', clean(r[1])).replace('\\pm', '+/-') +
                     ' MJ listed by Ma & Ge (from astrometry, or RV+transit for transiting systems); orbit P/e are RV values')
        if 'HD 137510' in host:
            n.append('Mc given as range 20-60 MJ')  # Mc cell of the HD 137510 row reads '20.0-60.0'
        n += flags(ms if ms is not None else Mc, None, P)
        n.append('compilation; orig refs: ' + ', '.join(refmap.get(x, x) for x in rr))
        n.append('no K/omega/Tp in source; take orbit from original refs')
        if 'fixed' in r[4] or 'adopted' in r[4]:
            n.append('e fixed in source')
        rows.append(dict(host=host, companion=comp, P_d=P,
                         e_P_d=eP,  # asymmetric errors averaged
                         e=e, msini_mj=ms, e_msini_mj=e_ms if e_ms else None, e_e=e_e if e_e else None, mstar_msun=Ms, fit_type='RV-only',
                         table_ref='arXiv:1303.6442 Table 1 (compilation)', notes='; '.join(n)))
    write_csv(rows)


if __name__ == '__main__':
    main()
