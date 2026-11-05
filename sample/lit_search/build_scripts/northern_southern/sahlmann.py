import sys, re
sys.path.insert(0, sys.argv[1])
from texhelp import *
src = open(sys.argv[2]).read()
out = sys.argv[3]

def table_body(label):
    i = src.index('\\label{' + label.replace('label{','') if not label.startswith('label{') else '\\' + label)
    j = src.index('\\hline', src.index('\\hline\\hline', i) + 12)  # after header hline
    # header ends at the first \hline after the units row
    k = src.index('\\end{tabular}', i)
    return src[j:k]

orb = split_rows(table_body('tab:KeplerOrbits'))
par = split_rows(table_body('tab:KeplerParams'))
stel = split_rows(table_body('label{tab:stellarChar}'))
mstar = {}
for r in stel:
    if len(r) < 9: continue
    mstar[norm_name(r[0])] = parse_val(r[7])[0]
mstar['HD 4747'] = mstar.get('HD 4747')

# true masses from Hipparcos (Table mass2sigma)
m2 = {}
for r in split_rows(table_body('tab:mass2sigma')):
    if len(r) < 12: continue
    nm = clean(r[0]); fm = re.search(r'tablefootmark\{(\w)\}', r[0])
    v, e = parse_val(r[6])
    m2[(norm_name(nm), fm.group(1) if fm else None)] = (v, e, clean(r[11]))

def key(cell):
    fm = re.search(r'tablefootmark\{(\w)\}', cell)
    return norm_name(cell.split('\\tablefootmark')[0]), (fm.group(1) if fm else None)

params = {}
for r in par:
    if len(r) < 12: continue
    params[key(r[0])] = r
# Params table uses different footnote letters for HD4747/HD211847 (b/c/d) and HD43848 (f/g)
rows = []
for r in orb:
    if len(r) < 8: continue
    name, fm = key(r[0])
    # map orbit-table footnotes to params-table footnotes
    pmap = {('HD 43848','g'):'f', ('HD 43848','h'):'g'}
    if name in ('HD 4747','HD 211847'):
        pmap = {(name,'d'):'b', (name,'e'):'c', (name,'f'):'d'}
    pf = pmap.get((name, fm), fm)
    pr = params.get((norm_name(name), pf)) or params.get((name.replace(' ',''), pf))
    if pr is None:
        # HD4747 in params table is written 'HD4747'
        pr = params.get(('HD 4747', pf)) if name == 'HD 4747' else None
    # Keep only adopted solutions
    if name == 'HD 43848' and fm == 'h':
        continue   # formal high-e solution; mentioned in notes
    if name in ('HD 4747','HD 211847') and fm != 'e':
        continue   # keep best-fit only
    P, eP = parse_val(r[3]); e, _ = parse_val(r[4]); K, _ = parse_val(r[5])
    T0, _ = parse_val(r[6]); om, _ = parse_val(r[7])
    a = parse_val(pr[10])[0] if pr else None
    ms = parse_val(pr[11])[0] if pr else None
    notes = []
    if name == 'HD 43848':
        notes.append('adopted low-e solution; formal high-e alternative (e=0.885, K=1859 m/s, Msini=50.2 MJ) also in Table 3/4')
    if name in ('HD 4747','HD 211847'):
        notes.append('period not covered; best-fit solution with P range given in paper (HD 4747: 7900-29000 d; HD 211847: 3750-100000 d)')
    mk = m2.get((name, {'g':'a','h':'b'}.get(fm) if name=='HD 43848' else None))
    if mk:
        notes.append(f'Hipparcos-astrometry true mass M2={mk[0]} +/- {mk[1]:.4g} MJ (orbit significance {mk[2]}%, Sahlmann Table 7)')
    if name == 'HIP 103019':
        notes.append('HARPS (not CORALIE) target')
    if name == 'HD 17289':
        notes.append('companion contaminates CCF (SB2-like RV bias); M2=0.52 Msun')
    if (a is not None and a > 5) or (ms is not None and ms > 80):
        notes.insert(0, 'OUT OF TARGET RANGE (a>5 AU)' if a and a > 5 else 'OUT OF TARGET RANGE (Msini>80 MJ)')
    rows.append(dict(host=name, companion=name + ' b', P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                     Tp_bjd=(T0 + 2400000.0) if T0 is not None else None, msini_mj=ms, a_au=a,
                     mstar_msun=mstar.get(name), fit_type='RV-only',
                     table_ref='arXiv:1009.5991 Tables 3-4 (stellar mass Table 2)',
                     notes='; '.join(notes)))

# Literature targets with spectroscopic elements (Table 8 'tab:target2sigma')
refs = {'1':'Minniti+2009','2':'Nidever+2002','3':'Perrier+2003','4':'Wright+2009','5':'Jenkins+2009'}
for r in split_rows(table_body('tab:target2sigma')):
    if len(r) < 10: continue
    nm = clean(r[0]); comp = nm + ('' if nm.endswith('c') else ' b')
    host = norm_name(re.sub(r'c$', '', nm))
    P, eP = parse_val(r[4]); e, _ = parse_val(r[5]); K, _ = parse_val(r[6]); T0, _ = parse_val(r[7]); om, _ = parse_val(r[8])
    ms = parse_val(r[3])[0]
    notes = [f'orbit compiled by Sahlmann+2011 from {refs[clean(r[9])]}; not a new fit']
    if ms < 5: notes.insert(0, 'BELOW 5 MJ (Msini=3.6 MJ; listed because Hipparcos mass is substellar)')
    mk = m2.get((host, None))
    if mk:
        notes.append(f'Hipparcos-astrometry M2={mk[0]} +/- {mk[1]:.4g} MJ (significance {mk[2]}%)')
    rows.append(dict(host=host, companion=host + (' c' if nm.endswith('c') else ' b'), P_d=P, e_P_d=eP, K_ms=K, e=e,
                     omega_deg=om, Tp_bjd=T0 + 2400000.0, msini_mj=ms, mstar_msun=parse_val(r[2])[0],
                     fit_type='RV-only', table_ref='arXiv:1009.5991 Table 8 (literature elements)',
                     notes='; '.join(notes)))
write_csv(out, rows)
print(len(rows))
