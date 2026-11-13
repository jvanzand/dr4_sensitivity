"""Build jenkins_2015_orbits.csv from the arXiv LaTeX source of Jenkins et al. 2015, MNRAS 453, 1439
("The Observed Distribution of Spectroscopic Binaries from the Anglo-Australian Planet Search",
arXiv:1507.04749).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_1507.04749.src  arXiv e-print (gzipped tar). Jenkins2015.tex is read from the extracted copy
                        arxiv_src/Jenkins2015.tex if present, else directly from the tarball.
Parsed from Jenkins2015.tex:
  Table 1 (label primaries): stellar mass
  Table 2(i) (label orbit) and Table 2(ii) (label orbit2): P, omega, e, Tp [JD-2450000], K, a
  Sect. 3 per-star paragraphs ("{\\it HD\\,39213}: ..."): M sin i in Msun ("minimum mass of X +/- Y Msun");
  the paper has no M sin i column.
Output: jenkins_2015_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_1507.04749.src': 'https://arxiv.org/e-print/1507.04749'}
ARCHIVE = 'arxiv_1507.04749.src'
TEX = 'Jenkins2015.tex'
OUT = 'jenkins_2015_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes", "e_K_ms", "e_msini_mj", "e_e", "e_a_au"]

MJ_PER_MSUN = 1047.57    # M sin i in the text is in Msun -> MJ
TP_OFFSET = 2450000.0    # Table 2 Tp column is labelled "JD-2450000"

# Hand-chosen selection (reproduces the existing table): the two binaries whose text M sin i is below
# 80 MJ.  KNOWN OMISSION: HD 42024 (Table 2(i): P = 76.26 d, a = 0.38 au; text: "minimum mass of
# 0.066 +/- 0.003 Msun ... another possible brown dwarf companion", i.e. 69 MJ) also qualifies but was
# left out of the original table, which assumed all other 14 binaries have M sin i >= 0.1 Msun.
# (HD 150248's text value "0.10 +/- 0.02 MJ" is evidently a typo for Msun.)
HOSTS = [
    # host, Table 2 part, table_ref, notes (hand-written; values quoted from the HD 39213 / HD 145825
    # paragraphs of Sect. 3 and Table 2)
    ('HD 39213', 'orbit', 'arXiv:1507.04749 Table 2(i) (Msini from star notes in text; mstar Table 1)',
     'AAT/UCLES (AAPS) + HARPS; Msini given only in text as 0.07+/-0.01 Msun, converted at 1047.57 MJ/Msun; '
     'authors note poor fit (chi2_nu=1.97, rms 17 m/s); Tp column labelled JD-2450000; could be BD only if edge-on'),
    ('HD 145825', 'orbit2', 'arXiv:1507.04749 Table 2(ii) (Msini and mstar from star notes in text)',
     'OUT OF TARGET RANGE (a>5 AU); minimum-period solution of a partially covered orbit; Msini 0.06+/-0.01 Msun '
     'in text, converted at 1047.57 MJ/Msun; Tp column labelled JD-2450000'),
]
NOTE_ERR = ('e_P_d, e_K_ms, e_e, e_a_au: published 1-sigma errors from Table 2; e_msini_mj: the published error on '
            'M sin i from the text (Msun) converted to MJ at 1047.57 MJ/Msun like msini_mj')


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

    def table_rows(label):
        """{star: cells} for the data rows (first cell starting with 'HD') of the tabular holding \\label{label}."""
        k = src.index('\\label{' + label + '}')
        i = src.rindex('\\begin{tabular}', 0, k)
        out = {}
        for r in split_rows(src[i:src.index('\\end{tabular}', k)]):
            cells = [c for c in r]
            # table 1 starts with an index column ('4&HD39213&...')
            if not clean(cells[0]).startswith('HD') and len(cells) > 1 and clean(cells[1]).startswith('HD'):
                cells = cells[1:]
            if clean(cells[0]).startswith('HD'):
                out[norm_name(cells[0])] = cells
        return out

    prim = table_rows('primaries')   # Star, B-V, V, Parallax, MV, SpT, [Fe/H] Casa, [Fe/H] Bond, logR'HK, Mass
    orb = {'orbit': table_rows('orbit'), 'orbit2': table_rows('orbit2')}  # Star, P, omega, e, Tp, K, a, chi2, rms

    rows = []
    for host, part, tref, notes in HOSTS:
        r = orb[part][host]
        P, eP = parse_val(r[1])
        om, _ = parse_val(r[2])
        e, ee = parse_val(r[3])
        Tp, _ = parse_val(r[4])
        K, eK = parse_val(r[5])
        a, ea = parse_val(r[6])
        mstar = parse_val(prim[host][9])[0]
        # M sin i from the star's paragraph in Sect. 3, e.g. "{\it HD\,39213}: ... minimum mass of 0.07\,$\pm$\,0.01\,M$_{\sun}$"
        para = src[src.index('{\\it HD\\,%s}:' % host[3:]):]
        para = para[:para.index('{\\it HD', 10)]
        m = re.search(r'minimum mass of\s*([\d.]+)\\,\$\\pm\$\\,([\d.]+)\\,M\$_\{\\sun\}\$', para)
        msini_msun = float(m.group(1))
        e_msini_msun = float(m.group(2))  # published error on M sin i, in Msun
        rows.append(dict(host=host, companion=host + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                         Tp_bjd=Tp + TP_OFFSET,  # JD-2450000 -> full JD
                         msini_mj=round(msini_msun * MJ_PER_MSUN, 1),  # Msun -> MJ
                         a_au=a, mstar_msun=mstar, fit_type='RV-only', table_ref=tref,
                         notes=notes + '; ' + NOTE_ERR,
                         e_K_ms=eK or None, e_e=ee or None, e_a_au=ea or None,
                         e_msini_mj=round(e_msini_msun * MJ_PER_MSUN, 1) or None))  # Msun -> MJ, as msini_mj
    write_csv(rows)


if __name__ == '__main__':
    main()
