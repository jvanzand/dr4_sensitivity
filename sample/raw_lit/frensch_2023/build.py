"""Build frensch_2023_orbits.csv from the arXiv LaTeX source of Frensch et al. 2023, A&A 675, A173
("The HARPS search for southern extra-solar planets. XLVII. Five Jupiter-mass planets in long-period
orbits, one highly irradiated Neptune, one brown dwarf, and five stellar binaries", arXiv:2306.04420).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  arxiv_2306.04420.src  arXiv e-print (gzipped tar). 46203.tex is read from the extracted copy
                        arxiv_src/46203.tex if present, else directly from the tarball.
Tables used (labels in 46203.tex; both transposed, one column per star or companion):
  table:stellar_parameters  (Table 2, two blocks: M*)
  table:brown_dwarfs        (Table 5: the four brown-dwarf candidates; all rows kept)
Output: frensch_2023_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = {'arxiv_2306.04420.src': 'https://arxiv.org/e-print/2306.04420'}
ARCHIVE = 'arxiv_2306.04420.src'
TEX = '46203.tex'
OUT = 'frensch_2023_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

# T_P is given in "rBJD", which the paper does not define. Values such as 55738 only make sense as
# BJD - 2400000 (the DACE convention), so that offset is assumed.
TP_OFFSET = 2400000.0
TPN = ('Tp given in rBJD (undefined in paper; values only consistent with BJD-2400000, the DACE convention, '
       'so converted with +2400000)')
OUT_A = 'OUT OF TARGET RANGE (a>5 AU)'
TABLE_REF = 'arXiv:2306.04420 Table 5 (label table:brown_dwarfs); mstar Table 2'


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
        return [r for r in split_rows(b)]

    def label_of(r):
        return r[0].strip().split('\n')[-1].strip()

    # Table 2: blocks start with a '\multicolumn{2}{l}{Parameter} & star names' row
    mstar = {}
    names = []
    for r in tab_rows('table:stellar_parameters'):
        if 'Parameter' in label_of(r):
            names = [clean(c) for c in r[1:]]
        elif label_of(r) == '$M_{\\star}$':
            for n, c in zip(names, r[2:]):
                mstar[n] = parse_val(c)[0]

    t = {}
    for r in tab_rows('table:brown_dwarfs'):
        lab = label_of(r)
        if 'Parameter' in lab:
            comps = [clean(c) for c in r[1:]]
        else:
            t[lab] = r

    def central(cell):
        """Central value as printed, e.g. '{18.77}_{-0.63}^{+0.66}' -> '18.77'."""
        return re.match(r'\{?([-\d.]+)', clean(cell)).group(1)

    rows = []
    first = True
    for k, comp in enumerate(comps):
        col = 2 + k
        m = re.match(r'^(HD)\s*(\d+)\s*([A-Za-z])$', comp)  # 'HD62364 b', 'HD56380B', 'HD33473C'
        base, letter = m.group(1) + m.group(2), m.group(3)
        skey = base if base in mstar else base + 'A'      # HD 33473 C orbits HD33473A in Table 2
        host = norm_name(skey)
        P, eP = parse_val(t['$P$'][col])
        Tp, _ = parse_val(t['$T_P$'][col])
        e, _ = parse_val(t['$e$'][col])
        om, _ = parse_val(t['$\\omega$'][col])
        K, _ = parse_val(t['$K$'][col])
        ms, _ = parse_val(t['$m_2\\sin i$'][col])
        a, _ = parse_val(t['$a_2$'][col])
        m2cell = clean(t['$m_2$$^{(a)}$'][col]).replace(' ', '')
        mm = re.match(r'^\{([\d.]+)\}_\{-([\d.]+)\}\^\{\+([\d.]+)\}$', m2cell)
        m2, m2lo, m2hi = mm.group(1), mm.group(2), mm.group(3)
        inc = central(t['$i$$^{(a)}$'][col])
        drift = clean(t['Drift'][col])
        notes = [OUT_A] if a > 5 else []
        notes.append('orvara %strue mass %s +%s/-%s MJ%s, i=%s deg'
                     % ('(RV+Hipparcos-Gaia PMa) ' if first else '', m2, m2hi, m2lo,
                        ' (stellar)' if float(m2) > 80 else '', inc))
        if drift not in ('-', ''):
            notes.append('fit includes drift %g m/s/yr' % parse_val(drift)[0])
        notes.append(TPN)
        first = False
        rows.append(dict(host=host, companion=norm_name(base) + ' ' + letter, P_d=P,
                         e_P_d=eP,  # asymmetric errors averaged
                         K_ms=K, e=e, omega_deg=om,
                         Tp_bjd=Tp + TP_OFFSET,  # rBJD -> BJD (assumed offset, see TP_OFFSET)
                         msini_mj=ms, a_au=a, mstar_msun=mstar[skey], fit_type='RV-only',
                         table_ref=TABLE_REF, notes='; '.join(notes)))
    write_csv(rows)


if __name__ == '__main__':
    main()
