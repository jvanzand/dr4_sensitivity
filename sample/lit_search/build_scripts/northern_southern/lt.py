"""Helpers to turn LaTeX table cells into numbers and write the standard orbit CSV."""
import csv
import re

COLS = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp_bjd',
        'msini_mj', 'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes']

NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def val(cell):
    """Central value of a LaTeX cell like '$1718\\pm11$' or '$2.6^{+0.1}_{-0.2}$'."""
    if cell is None:
        return ''
    s = cell.replace('$', '').replace('{', '').replace('}', '').replace('~', '')
    s = re.sub(r'\s+', '', s.replace('\\,', '')).replace('−', '-').replace('--', '')
    m = re.match(NUM, s)
    return m.group(0) if m else ''


def err(cell):
    """Symmetric error (mean of +/- if asymmetric) from a LaTeX cell, or ''."""
    if cell is None:
        return ''
    s = re.sub(r'\s+', '', cell.replace('$', '').replace('{', '').replace('}', ''))
    m = re.search(r'\\pm(' + NUM + ')', s)
    if m:
        return m.group(1)
    m = re.search(r'\^\+?(' + NUM + r')_-?(' + NUM + ')', s) or re.search(r'_-?(' + NUM + r')\^\+?(' + NUM + ')', s)
    if m:
        return '%.6g' % ((abs(float(m.group(1))) + abs(float(m.group(2)))) / 2)
    m = re.search(r'\+/-(' + NUM + ')', s)
    if m:
        return m.group(1)
    return ''


def cells(line):
    line = line.split('\\\\')[0]
    return [c.strip() for c in line.split('&')]


def write(path, rows):
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in COLS})
    print('wrote', len(rows), 'rows to', path)


def bd(name):
    """Normalize 'BD+730275' / 'BD+73~0275' / 'BD-004475' to SIMBAD-style 'BD+73 275'."""
    s = name.replace('~', ' ').replace('\\,', ' ').strip()
    m = re.match(r'(BD|TYC)?\s*([+-])\s*(\d{2})\s*(\d+)$', s.replace('BD', '').strip()) if s.startswith('BD') else None
    if not m:
        return s
    return 'BD%s%s %d' % (m.group(2), m.group(3), int(m.group(4)))


def addbjd(x, off=2400000):
    v = val(x) if not isinstance(x, (int, float)) else x
    if v == '':
        return ''
    nd = len(str(v).split('.')[1]) if '.' in str(v) else 0
    return '%.*f' % (nd, float(v) + off)


def tyc(name):
    """'TYC 0123-00549-1' -> 'TYC 123-549-1'."""
    m = re.match(r'TYC\s*0*(\d+)-0*(\d+)-(\d+)', name.strip())
    return 'TYC %s-%s-%s' % m.groups() if m else name.strip()
