"""Build jenkins_2009_orbits.csv from the arXiv LaTeX source of Jenkins et al. 2009, MNRAS 398, 911
("First results from the Calan-Hertfordshire Extrasolar Planet Search: exoplanets and the discovery of
an eccentric brown dwarf in the desert", arXiv:0905.2985).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_0905.2985.src  arXiv e-print (gzipped tar). Jenkins0509.tex is read from the extracted copy
                       arxiv_src/Jenkins0509.tex if present, else directly from the tarball.
Tables used (labels in Jenkins0509.tex; both transposed):
  tab:params      (Table 1: stellar parameters of HD 191760; M* from the row "M*/Msun - Here")
  tab:mag_params  (Table 5: orbital parameters of HD 191760 b, HD 48265 b, HD 143361 b, HD 154672 b)
Selection: companions with M sin i >= 5 MJ (HD 191760 b and HD 154672 b).
Output: jenkins_2009_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_0905.2985.src': 'https://arxiv.org/e-print/0905.2985'}
ARCHIVE = 'arxiv_0905.2985.src'
TEX = 'Jenkins0509.tex'
OUT = 'jenkins_2009_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

# Table 5 labels the T0 row "JD-2,450,000"; it is taken to be the time of periastron.
TP_OFFSET = 2450000.0

# Hand-written notes (descriptive, not values parsed from this paper). The HD 191760 Hipparcos mass
# is from Sahlmann+2011 (arXiv:1009.5991) Table 7: M2 = 185.3 (+94.1/-93.4) MJ, i.e. about 185 +/- 94.
NOTES = {
    'HD 191760': 'CHEPS/HARPS; T0 column labelled JD-2,450,000 (assumed time of periastron; paper does not say BJD); '
                 'Sahlmann+2011 Hipparcos M2=185 +/- 94 MJ (low significance)',
    'HD 154672': 'Magellan + HARPS update of Lopez-Morales+2008 planet; T0 in JD-2,450,000; mstar not given in this '
                 'paper; likely already in NASA Exoplanet Archive',
}


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
        out = {}
        for r in split_rows(b):
            out[r[0].strip().split('\n')[-1].strip()] = r
        hdr = [k for k in out if 'Parameter' in k][0]
        names = [clean(re.sub(r'\\multicolumn\{\d+\}\{\w+\}\{(.*)\}$', r'\1', c.strip())) for c in out[hdr][1:]]
        return names, out

    def row(t, prefix):
        return [v for k, v in t.items() if k.startswith(prefix)][0]

    snames, st = transposed('tab:params')
    mstar = {norm_name(snames[0]): parse_val(st['$M_{\\rm{\\star}}$/$M_{\\odot}$ - Here'][1])[0]}

    names, t = transposed('tab:mag_params')
    rows = []
    for k, n in enumerate(names):
        col = 1 + k
        host = norm_name(n[:-1])  # 'HD191760b' -> 'HD 191760'
        P, eP = parse_val(row(t, 'Orbital period')[col])
        K, _ = parse_val(row(t, 'Velocity amplitude')[col])
        e, _ = parse_val(row(t, 'Eccentricity')[col])
        om, _ = parse_val(row(t, '$\\omega$')[col])
        T0, _ = parse_val(row(t, '$T$$_{\\rm{0}}$')[col])
        ms, _ = parse_val(row(t, 'M~sin~$i$')[col])
        a, _ = parse_val(row(t, 'Semimajor axis')[col])
        if ms < 5:
            continue
        has_m = host in mstar
        rows.append(dict(host=host, companion=host + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                         Tp_bjd=T0 + TP_OFFSET,  # JD-2,450,000 -> full JD
                         msini_mj=ms, a_au=a, mstar_msun=mstar.get(host), fit_type='RV-only',
                         table_ref='arXiv:0905.2985 Table 5' + (' (mstar Table 1)' if has_m else ''),
                         notes=NOTES[host]))
    write_csv(rows)


if __name__ == '__main__':
    main()
