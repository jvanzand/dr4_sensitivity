"""Build unger_2023_orbits.csv from the arXiv LaTeX source of Unger et al. 2023, A&A 680, A16
("Exploring the brown dwarf desert with precision radial velocities and Gaia DR3 astrometric orbits",
arXiv:2310.02758).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_2310.02758.src  arXiv e-print (gzipped tar). main_full.tex is read from the extracted copy
                        arxiv_src/main_full.tex if present, else directly from the tarball.
Tables used (labels in main_full.tex):
  tab:stellar_parameters  (Table 1: M*)
  tab:rvonlyMCMC          (Table 3: RV-only MCMC; -> fit_type RV-only rows)
  tab:jointMCMC           (Table 4: joint Gaia DR3 + RV MCMC; -> fit_type RV+astrometry rows)
Output: unger_2023_orbits.csv
Error columns (appended; published 1-sigma, asymmetric errors averaged): e_P_d, e_K_ms, e_msini_mj, e_e,
e_a_au, all from the same table row as the value (Table 3 for RV-only rows, Table 4 for joint rows;
in Table 4 the msini error is that of the RV-only Msini column and e_a_au is the error on a_r, the
relative semi-major axis, which is the a_au value stored).
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_2310.02758.src': 'https://arxiv.org/e-print/2310.02758'}
ARCHIVE = 'arxiv_2310.02758.src'
TEX = 'main_full.tex'
OUT = 'unger_2023_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes",
        "e_K_ms", "e_msini_mj", "e_e", "e_a_au"]
# Optional published 1-sigma errors (columns appended after the original ones); blank = not published.
# Asymmetric errors are averaged (parse_val); no error is estimated or propagated here.

# Tp in Tables 3 and 4 is relative to the Gaia DR3 epoch 2016.0 = BJD 2457389.0 (stated in the
# Table 4 note, main_full.tex line ~1826). Table 3 states no epoch; the same reference was assumed
# and checked (HD 3277 Tp agrees with Sahlmann+2011 to 0.001 in phase).
REF = 2457389.0


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

    def body(label, start_marker='\\midrule'):
        """Tabular body of the table holding \\label{label}, after the first \\midrule."""
        j = src.index('\\label{' + label + '}')
        t = src.rindex('\\begin{table}', 0, j)
        i = src.index('\\begin{tabular}', t)
        k = src.index('\\end{tabular}', i)
        b = src[i:k]
        return b[b.index(start_marker) + len(start_marker):]

    def uname(n):
        n = clean(n).strip()
        m = re.match(r'^(BD|CD)([+-])(\d{2})(\d+)$', n)  # 'BD-1763' -> 'BD-17 63'
        if m:
            return '%s%s%s %d' % (m.group(1), m.group(2), m.group(3), int(m.group(4)))
        return norm_name(n)

    # Table 1: stellar mass (column 7)
    mstar = {}
    for r in split_rows(body('tab:stellar_parameters')):
        if len(r) >= 13:
            mstar[uname(r[0])] = parse_val(r[6])[0]

    rows = []
    rv = {}
    # Table 3: name, K, P, a, Tp, e, omega, Msini
    for r in split_rows(body('tab:rvonlyMCMC')):
        if len(r) < 8:
            continue
        h = uname(r[0])
        K, eK = parse_val(r[1])
        P, eP = parse_val(r[2])
        a, ea = parse_val(r[3])
        Tp, _ = parse_val(r[4])
        e, ee = parse_val(r[5])
        om, _ = parse_val(r[6])
        ms, ems = parse_val(r[7])
        rv[h] = ms
        notes = []
        if ms is not None and ms < 5:
            notes.append('BELOW 5 MJ')
        if ms is not None and ms > 80:
            notes.append('OUT OF TARGET RANGE (Msini>80 MJ)')
        if a is not None and a > 5:
            notes.append('OUT OF TARGET RANGE (a>5 AU)')
        notes.append('Tp converted from paper value relative to BJD 2457389.0')
        if h + 'A' in mstar and h not in mstar:
            notes.append('listed as ' + h + 'A in Tables 1 and 4 (component A of a binary)')
        rows.append(dict(host=h, companion=h + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                         Tp_bjd=Tp + REF if Tp is not None else None,  # relative Tp -> full BJD
                         msini_mj=ms, a_au=a,
                         e_K_ms=eK, e_msini_mj=ems, e_e=ee, e_a_au=ea,
                         mstar_msun=mstar.get(h, mstar.get(h + 'A')), fit_type='RV-only',
                         table_ref='arXiv:2310.02758 Table 3 (RV-only MCMC); mstar Table 1',
                         notes='; '.join(notes)))

    # Table 4: name, plx, M1, Msini(RV-only), M2, K, a, P, Tp, e, omega, Omega, i
    section = 'robust'
    for r in split_rows(body('tab:jointMCMC')):
        if 'Challenging' in ' '.join(r):
            section = 'challenging'  # rows after the "Challenging cases" sub-header
            continue
        if len(r) < 13 or 'multicolumn' in r[0]:
            continue
        h = uname(r[0])
        M1 = parse_val(r[2])[0]
        ms, ems = parse_val(r[3])
        M2, eM2 = parse_val(r[4])
        K, eK = parse_val(r[5])
        a, ea = parse_val(r[6])
        P, eP = parse_val(r[7])
        Tp, _ = parse_val(r[8])
        e, ee = parse_val(r[9])
        om, _ = parse_val(r[10])
        Om, _ = parse_val(r[11])
        inc, _ = parse_val(r[12])
        host = h[:-1] if h.endswith('A') and h[:-1] in rv else h  # HD 68638A -> HD 68638 (Table 3 name)
        notes = ['joint Gaia DR3 + RV fit; true mass M2=%g +/- %.3g MJ; inc=%g deg; Omega=%g deg; '
                 'Msini column is from RV-only fit' % (M2, eM2, inc, Om),
                 'Tp relative to BJD 2457389.0 in paper; converted']
        if section == 'challenging':
            notes.insert(0, 'CHALLENGING CASE (authors caution: Gaia and RV solutions discrepant)')
        if M2 > 80:
            notes.insert(0, 'true mass %g MJ is stellar' % M2)
        rows.append(dict(host=host, companion=host + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                         Tp_bjd=Tp + REF if Tp is not None else None,  # relative Tp -> full BJD
                         msini_mj=ms, a_au=a, mstar_msun=M1,
                         e_K_ms=eK, e_msini_mj=ems, e_e=ee, e_a_au=ea,
                         fit_type='RV+astrometry', table_ref='arXiv:2310.02758 Table 4 (joint Gaia DR3+RV MCMC)',
                         notes='; '.join(notes)))
    write_csv(rows)


if __name__ == '__main__':
    main()
