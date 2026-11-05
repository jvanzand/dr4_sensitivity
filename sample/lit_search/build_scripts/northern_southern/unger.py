import sys, re
sys.path.insert(0, sys.argv[1])
from texhelp import *
src = open(sys.argv[2]).read()
out = sys.argv[3]
REF = 2457389.0

def body(label, start_marker='\\midrule'):
    j = src.index('\\label{' + label + '}')
    t = src.rindex('\\begin{table}', 0, j)
    i = src.index('\\begin{tabular}', t)
    k = src.index('\\end{tabular}', i)
    b = src[i:k]
    return b[b.index(start_marker) + len(start_marker):]

def uname(n):
    n = clean(n).strip()
    m = re.match(r'^(BD|CD)([+-])(\d{2})(\d+)$', n)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)} {int(m.group(4))}"
    return norm_name(n)

# stellar masses
mstar = {}
for r in split_rows(body('tab:stellar_parameters')):
    if len(r) >= 13:
        mstar[uname(r[0])] = parse_val(r[6])[0]

rows = []
rv = {}
for r in split_rows(body('tab:rvonlyMCMC')):
    if len(r) < 8: continue
    h = uname(r[0])
    K, _ = parse_val(r[1]); P, eP = parse_val(r[2]); a, _ = parse_val(r[3]); Tp, _ = parse_val(r[4])
    e, _ = parse_val(r[5]); om, _ = parse_val(r[6]); ms, _ = parse_val(r[7])
    rv[h] = ms
    notes = []
    if ms is not None and ms < 5: notes.append('BELOW 5 MJ')
    if ms is not None and ms > 80: notes.append('OUT OF TARGET RANGE (Msini>80 MJ)')
    if a is not None and a > 5: notes.append('OUT OF TARGET RANGE (a>5 AU)')
    notes.append('Tp converted from paper value relative to BJD 2457389.0')
    if h + 'A' in mstar and h not in mstar: notes.append('listed as ' + h + 'A in Tables 1 and 4 (component A of a binary)')
    rows.append(dict(host=h, companion=h + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                     Tp_bjd=Tp + REF if Tp is not None else None, msini_mj=ms, a_au=a,
                     mstar_msun=mstar.get(h, mstar.get(h + 'A')), fit_type='RV-only',
                     table_ref='arXiv:2310.02758 Table 3 (RV-only MCMC); mstar Table 1', notes='; '.join(notes)))

jb = body('tab:jointMCMC')
section = 'robust'
for r in split_rows(jb):
    if 'Challenging' in ' '.join(r): section = 'challenging'; continue
    if len(r) < 13 or 'multicolumn' in r[0]: continue
    h = uname(r[0])
    plx = parse_val(r[1])[0]; M1 = parse_val(r[2])[0]; ms = parse_val(r[3])[0]
    M2, eM2 = parse_val(r[4]); K, _ = parse_val(r[5]); a, _ = parse_val(r[6]); P, eP = parse_val(r[7])
    Tp, _ = parse_val(r[8]); e, _ = parse_val(r[9]); om, _ = parse_val(r[10]); Om, _ = parse_val(r[11]); inc, _ = parse_val(r[12])
    host = h[:-1] if h.endswith('A') and h[:-1] in rv else h
    notes = [f'joint Gaia DR3 + RV fit; true mass M2={M2:g} +/- {eM2:.3g} MJ; inc={inc:g} deg; Omega={Om:g} deg; Msini column is from RV-only fit',
             'Tp relative to BJD 2457389.0 in paper; converted']
    if section == 'challenging':
        notes.insert(0, 'CHALLENGING CASE (authors caution: Gaia and RV solutions discrepant)')
    if M2 > 80: notes.insert(0, f'true mass {M2:g} MJ is stellar')
    rows.append(dict(host=host, companion=host + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                     Tp_bjd=Tp + REF if Tp is not None else None, msini_mj=ms, a_au=a, mstar_msun=M1,
                     fit_type='RV+astrometry', table_ref='arXiv:2310.02758 Table 4 (joint Gaia DR3+RV MCMC)',
                     notes='; '.join(notes)))
write_csv(out, rows)
print(len(rows), sum(1 for r in rows if r['fit_type']=='RV-only'))
