import re, csv

COLS = ["host","companion","P_d","e_P_d","K_ms","e","omega_deg","Tp_bjd","msini_mj",
        "a_au","mstar_msun","fit_type","table_ref","notes"]

def clean(s):
    s = re.sub(r'\\tablefootmark\{[^}]*\}', '', s)
    s = re.sub(r'\\vspace\{[^}]*\}', '', s)
    s = re.sub(r'\\(?:begin|end)\{math\}', '', s)
    s = s.replace('$', '').replace('\\,', '').replace('~', ' ')
    s = s.replace('\\mathrm', '').replace('\\textbf', '')
    return s.strip()

NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'

def parse_val(s):
    """Return (value, symmetric_error) from a LaTeX cell; None if absent."""
    s = clean(s)
    s = s.replace(' ', '').replace('{', '').replace('}', '')
    if s in ('', '---', '\\nodata', '\\cdots', '...', '-', '--'):
        return None, None
    m = re.match(rf'^({NUM})\^\+?({NUM})_-?({NUM})', s) or None
    if m:
        v, up, lo = m.groups()
        return float(v), (abs(float(up)) + abs(float(lo))) / 2
    m = re.match(rf'^({NUM})_-?({NUM})\^\+?({NUM})', s)
    if m:
        v, lo, up = m.groups()
        return float(v), (abs(float(up)) + abs(float(lo))) / 2
    m = re.match(rf'^({NUM})\\pm({NUM})', s)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.match(rf'^({NUM})', s)
    if m:
        return float(m.group(1)), None
    return None, None

def split_rows(tabular_body):
    rows = []
    for r in re.split(r'\\\\', tabular_body):
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
        return f"{m.group(1)} {m.group(2).strip()}"
    m = re.match(r'^(BD|CD)([+-])(\d{2})\s*0*(\d+)$', n)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)} {int(m.group(4))}"
    return n

def fmt(x):
    if x is None:
        return ''
    if isinstance(x, float):
        return f"{x:.10g}"
    return str(x)

def write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k)) for k in COLS})
