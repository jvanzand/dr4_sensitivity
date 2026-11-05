"""Build parc_2026_orbits.csv from the arXiv LaTeX source of Parc et al. 2026, accepted in A&A
("Companions of Hot Jupiters: two decades of radial velocity monitoring with CORALIE", arXiv:2610.25152).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_2610.25152.src  arXiv e-print (gzipped tar). aa62811-26.tex is read from the extracted copy
                        arxiv_src/aa62811-26.tex if present, else directly from the tarball.
Tables used (labels in aa62811-26.tex; transposed, one column per system):
  tab:orbital_params, tab:orbital_params_cont  ("Planet c (outer companion)" block)
  tab:wasp34                                    (alternative WASP-34 models, quoted in notes)
Selection: outer companions with M sin i >= 5 MJ (7 rows).
Output: parc_2026_orbits.csv
"""
import argparse
import csv
import math
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_2610.25152.src': 'https://arxiv.org/e-print/2610.25152'}
ARCHIVE = 'arxiv_2610.25152.src'
TEX = 'aa62811-26.tex'
OUT = 'parc_2026_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

# The table gives omega in radians (converted to degrees below) and the mean anomaly M0 at an epoch
# the text does not state, so Tp is left blank (M0 is quoted in notes). M* is only in the online table.
OUT_A = 'OUT OF TARGET RANGE (a>5 AU)'
TABLE_REF = 'arXiv:2610.25152 appendix table tab:orbital_params'


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

    def tab_rows(label):
        i = src.index('\\label{' + label + '}')
        b = src[src.index('\\begin{tabular}', i):src.index('\\end{tabular}', i)]
        b = re.sub(r'(?<!\\)%.*', '', b)  # drop LaTeX comments
        return split_rows(b)

    rows = []
    for label in ('tab:orbital_params', 'tab:orbital_params_cont'):
        names, t, block = [], {}, None
        for r in tab_rows(label):
            lab = r[0].strip().split('\n')[-1].strip()
            if lab == 'Parameter':
                names = [clean(c).strip('{}') for c in r[1:]]  # '\textbf{WASP-8}' -> 'WASP-8'
            elif 'Planet c' in lab:
                block = 'c'
            elif 'Planet b' in lab:
                block = 'b'
            elif block == 'c':
                t[lab.split('[')[0].strip()] = r
        for k, host in enumerate(names):
            col = 1 + k
            P, eP = parse_val(t['$P_c$'][col])
            M0, _ = parse_val(t['$M_{0,c}$'][col])
            K, _ = parse_val(t['$K_c$'][col])
            e, _ = parse_val(t['$e_c$'][col])
            om, _ = parse_val(t['$\\omega_c$'][col])
            ms, _ = parse_val(t['$M_c\\sin i_c$'][col])
            a, _ = parse_val(t['$a_c$'][col])
            if ms < 5:
                continue
            n = ([OUT_A] if a > 5 else []) + [
                'outer companion to hot Jupiter %s b; CORALIE (kima); omega given as %s rad (converted); mean anomaly '
                'M0=%s rad at an epoch not stated in text, so Tp left blank; mstar only in online table'
                % (host, om, M0)]
            if host == 'WASP-34':
                # M sin i of companion c in the '1 Keplerian + trend' and '2 Keplerians' blocks of tab:wasp34
                alt = [parse_val(r[1])[0] for r in tab_rows('tab:wasp34') if r[0].strip().startswith('$M\\sin i$')]
                n.append('paper table tab:wasp34 gives alternative models (1 Kep + trend or 2 Kep) with Msini %s or '
                         '%s MJ; 1-Keplerian row listed here' % (alt[1], alt[2]))
            rows.append(dict(host=host, companion=host + ' c', P_d=P,
                             e_P_d=eP,  # asymmetric errors averaged
                             K_ms=K, e=e,
                             omega_deg=round(math.degrees(om), 2),  # rad -> deg
                             msini_mj=ms, a_au=a, fit_type='RV-only', table_ref=TABLE_REF, notes='; '.join(n)))
    write_csv(rows)


if __name__ == '__main__':
    main()
