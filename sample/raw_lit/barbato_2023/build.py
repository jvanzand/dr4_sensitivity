"""Build barbato_2023_orbits.csv from VizieR J/A+A/674/A114 (Barbato et al. 2023, A&A 674, A114,
"The CORALIE survey for southern extrasolar planets XIX. Brown dwarfs and stellar companions
unveiled by radial velocity and astrometry", arXiv:2303.16717).

Run:  python -I build.py [--download]     (from anywhere; paths resolve relative to this file)
Inputs (in this directory):
  vizier/ReadMe, vizier/table1.dat (stellar parameters), vizier/table2.dat (orbital solutions)
  The arXiv e-print arxiv_2303.16717.src (https://arxiv.org/e-print/2303.16717; extracted in
  arxiv_src/) is kept for reference only: its Table 2 is truncated, so nothing is parsed from it.
Selection: every table2 record with RV-only M sin i <= 80 MJ (29 of 221). The rows use the RV-only
columns; the RV+PMa (orvara) solution goes in `notes`.
Output: barbato_2023_orbits.csv
"""
import argparse
import csv
import os
import re
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
VZ = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/674/A114/'
SOURCES = {'vizier/ReadMe': VZ + 'ReadMe',
           'vizier/table1.dat': VZ + 'table1.dat',
           'vizier/table2.dat': VZ + 'table2.dat'}
ARCHIVE = 'arxiv_2303.16717.src'  # not parsed (see above)
OUT = 'barbato_2023_orbits.csv'

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp_bjd", "msini_mj",
        "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

# Barbato+2023 give the mean longitude lambda0 (l0rv) instead of Tp and do not state its reference
# epoch. BJD 2455500 (the DACE convention) was inferred and checked: Tp computed from it matches 16
# independent Tp values in Sahlmann+2011 and Unger+2023 to within 0.005 in phase for well-sampled
# orbits (see sample/lit_search/notes_southern.md, section 3).
EPOCH = 2455500.0


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
# Whitespace-split token index -> ReadMe label for table2.dat (50 tokens per line; '---' = null):
#  0 Star | 1 Prv 2 E_Prv 3 e_Prv [d] | 4 Krv [km/s] | 7 erv | 10 l0rv [deg] | 13 omrv [deg]
#  16 Msini [MJ] | 22 arv [au] | 25 tspan/P | 27 Mtrue 28 E_Mtrue 29 e_Mtrue [MJ] | 33 arva [au]
#  36 Prva [d] | 39 erva | 42 irva [deg] | 49 pub (previous publication bibcode)
def f(x):
    return None if x == '---' else float(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--download', action='store_true', help='fetch missing inputs listed in SOURCES')
    args = ap.parse_args()
    fetch(args.download)

    mstar = {}
    with open(os.path.join(HERE, 'vizier', 'table1.dat')) as fh:
        for line in fh:
            mstar[line[0:9].strip()] = float(line[51:55])  # bytes 52-55: Mstar [Msun]
    recs = []
    with open(os.path.join(HERE, 'vizier', 'table2.dat')) as fh:
        for line in fh:
            t = line.split()
            assert len(t) == 50, (len(t), t[0])
            recs.append(t)

    rows = []
    seen = {}
    for t in recs:
        star = t[0]
        P, Pu, Pl = f(t[1]), f(t[2]), f(t[3])
        K = f(t[4]) * 1000  # km/s -> m/s
        e = f(t[7])
        l0 = f(t[10])
        om = f(t[13])
        ms = f(t[16])
        a = f(t[22])
        tspan = f(t[25])
        if ms > 80:
            continue
        n = seen.get(star, 0)
        seen[star] = n + 1
        host = norm_name(star)
        # Tp from mean longitude: lambda0 = M0 + omega at EPOCH, and M0 = 360 (EPOCH - Tp) / P
        # => Tp = EPOCH + (omega - lambda0) / 360 * P
        tp = EPOCH + (om - l0) / 360.0 * P
        notes = []
        if a > 5:
            notes.append('OUT OF TARGET RANGE (a>5 AU)')
        if ms < 5:
            notes.append('BELOW 5 MJ (planet in hierarchical stellar system)')
        if sum(1 for r in recs if r[0] == star) > 1:
            notes.append('host has >1 companion in Table 2 (hierarchical system)')
        notes.append('Tp_bjd computed as 2455500 + (omega - lambda0)/360*P from lambda0=%g deg (epoch 2455500 not '
                     'stated in paper; inferred and verified against Sahlmann+2011/Unger+2023 Tp, see notes file)' % l0)
        notes.append('tspan/P=%g' % tspan)
        Mt = f(t[27])
        if Mt is not None:
            notes.append('RV+PMa (orvara) true mass M=%g +%g/-%g MJ, i=%g deg; RV+PMa P=%g d, e=%g, a=%g AU'
                         % (Mt, f(t[28]), f(t[29]), f(t[42]), f(t[36]), f(t[39]), f(t[33])))
        if t[49] != '---':
            notes.append('previous publication ' + t[49])
        rows.append(dict(host=host, companion=host + ' ' + 'bcd'[n], P_d=P,
                         e_P_d=(Pu + Pl) / 2,  # asymmetric P errors averaged
                         K_ms=round(K, 3), e=e, omega_deg=om, Tp_bjd=round(tp, 3), msini_mj=ms, a_au=a,
                         mstar_msun=mstar.get(star), fit_type='RV-only',
                         table_ref='VizieR J/A+A/674/A114 table2.dat (RV-only columns); mstar table1.dat',
                         notes='; '.join(notes)))
    write_csv(rows)


if __name__ == '__main__':
    main()
