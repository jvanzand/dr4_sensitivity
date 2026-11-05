import sys, re
sys.path.insert(0, sys.argv[1])
from texhelp import *
L = sys.argv[2]
def rd(p): return open(f'{L}/{p}').read()
def flags(ms, a=None, P=None, mstar=None):
    n = []
    if ms is not None and ms < 5: n.append('BELOW 5 MJ')
    if ms is not None and ms > 80: n.append('OUT OF TARGET RANGE (Msini>80 MJ)')
    if a is not None and a > 5: n.append('OUT OF TARGET RANGE (a>5 AU)')
    elif a is None and P is not None and P > 4000: n.append('LIKELY OUT OF TARGET RANGE (P>4000 d, a not given)')
    return n

# ---------------- Kiefer 2021 ----------------
src = rd('kiefer_2021/arxiv_src/exoplanet_mass_gaia.tex')
rows = []
for lab, tref in (('tab:list_planets}', 'Table 3 (detection sample)'), ('tab:list_planets_nondet}', 'Table C.1 (non-detection sample)')):
    i = src.index('\\caption{\\label{' + lab)
    j = src.index('\\end{', src.index('\\hline', src.index('-2,450,000', src.rindex('-2,450,000', i, src.index('\\end{', i + 10) if False else i + 6000))))
    blk = src[i:src.index('\\end{longtable}' if 'nondet' in lab else '\\end{tabular}', i)]
    for r in split_rows(blk):
        if len(r) < 13 or '$P$' in r[1] or '(days)' in r[1] or 'caption' in r[0] or 'multicolumn' in r[0]: continue
        name = clean(r[0]).replace('\\endhead','').replace('\\endfirsthead','').strip()
        if not name: continue
        P, eP = parse_val(r[1]); ms, _ = parse_val(r[2]); K, _ = parse_val(r[3]); e, _ = parse_val(r[4])
        om, _ = parse_val(r[5]); Tp, _ = parse_val(r[6]); Ms, _ = parse_val(r[7]); a, _ = parse_val(r[8])
        if ms is None or ms < 5: continue
        m = re.match(r'^(.*?)\s+([b-h])$', name)
        host = norm_name(re.sub(r'\s+', ' ', m.group(1))) if m else name
        host = re.sub(r'^BD ([+-])(\d+) (\d+)$', r'BD\1\2 \3', host)
        n = flags(ms, a, P) + ['RV elements compiled by Kiefer+2021 from exoplanets.org (not new fits); transit flag=' + clean(r[11])]
        rows.append(dict(host=host, companion=name, P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om,
                         Tp_bjd=Tp + 2450000 if Tp is not None else None, msini_mj=ms, a_au=a, mstar_msun=Ms,
                         fit_type='RV-only', table_ref='arXiv:2009.14164 ' + tref, notes='; '.join(n)))
# GASTON true-mass results for detection sample
res = src[src.index('\\caption{\\label{tab:results_det}'):]
res = res[:res.index('\\end{tabular}')]
gm = {}
for r in split_rows(res):
    if len(r) >= 13 and re.search(r'\d', r[1] or ''):
        nm = clean(re.sub(r'\\tablefootmark\{[^}]*\}', '', r[0])).strip()
        gm[nm] = clean(r[9]).replace(' ', '')
for row in rows:
    if row['companion'] in gm:
        row['notes'] += f"; GASTON (Gaia DR1 excess noise) mass {gm[row['companion']]} MJ (1-sigma; Kiefer+2021 Table 4)"
write_csv(f'{L}/kiefer_2021/kiefer_2021_orbits.csv', rows); print('kiefer', len(rows))

# ---------------- Ma & Ge 2014 ----------------
src = rd('ma_ge_2014/arxiv_src/twopopulation_Feb18.tex')
blk = src[src.index('\\startdata') + 10:src.index('\\enddata')]
refs = src[src.index('\\tablerefs{'):]
refs = refs[:refs.index('\\end{deluxetable}')]
refmap = dict(re.findall(r'(\d+)\s*\\citet\{(\w+)\}', refs))
rows = []
for r in split_rows(blk):
    if len(r) < 10: continue
    nm = clean(r[0]).replace('$-$', '-').strip()
    Mc, eMc = parse_val(r[1]); ms, _ = parse_val(r[2]); P, eP = parse_val(r[3]); e, _ = parse_val(r[4]); Ms, _ = parse_val(r[5])
    rr = [x.strip() for x in clean(r[9]).split(',') if x.strip()]
    host = nm
    for suf in ('b', 'c', 'Bb'):
        pass
    m = re.match(r'^(.*?)(\s?[bc])$', nm)
    comp = nm if m else nm + ' b'
    host = m.group(1).strip() if m else nm
    host = {'11Com': '11 Com', 'BD+202457': 'BD+20 2457', 'TYC 3130-160': 'TYC 3130-160-1', 'COROT-15': 'CoRoT-15',
            'Corot-3': 'CoRoT-3', 'NGC 2423-3': 'NGC 2423 3', 'NGC 4349-127': 'NGC 4349 127'}.get(host, host)
    host = norm_name(host)
    if host == 'HD 168443': comp = 'HD 168443 c'
    elif comp == host + ' b' or comp == nm + ' b':
        comp = host + ' b'
    else:
        comp = host + ' ' + nm[-1]
    n = []
    if Mc is not None: n.append('true mass Mc=' + re.sub(r'\s+', '', clean(r[1])).replace('\\pm', '+/-') + ' MJ listed by Ma & Ge (from astrometry, or RV+transit for transiting systems); orbit P/e are RV values')
    if 'HD 137510' in host: n.append('Mc given as range 20-60 MJ')
    n += flags(ms if ms is not None else Mc, None, P)
    n.append('compilation; orig refs: ' + ', '.join(refmap.get(x, x) for x in rr))
    n.append('no K/omega/Tp in source; take orbit from original refs')
    if 'fixed' in r[4] or 'adopted' in r[4]: n.append('e fixed in source')
    rows.append(dict(host=host, companion=comp, P_d=P, e_P_d=eP, e=e, msini_mj=ms, mstar_msun=Ms,
                     fit_type='RV-only',
                     table_ref='arXiv:1303.6442 Table 1 (compilation)', notes='; '.join(n)))
write_csv(f'{L}/ma_ge_2014/ma_ge_2014_orbits.csv', rows); print('ma_ge', len(rows))

# ---------------- Grieves 2017 ----------------
src = rd('grieves_2017/arxiv_src/grieves2017_arxiv.tex')
sp = src[src.index('\\label{tab:starpar}'):]; sp = sp[:sp.index('\\end{tabular}')]
mst = {}
for r in split_rows(sp[sp.index('GSC 2.3') + 10:]):
    if len(r) >= 7:
        mst[clean(r[0]).replace('^{\\dagger}','').split('(')[0].strip(' \\\n\t')] = parse_val(r[4])[0]
cp = src[src.index('\\label{tab:compar}'):]; cp = cp[:cp.index('\\end {table*}') if '\\end {table*}' in cp else cp.index('\\end{table*}')]
t1, t2 = cp.split('\\begin{tabular}')[1:3]
r1 = [r for r in split_rows(t1[t1.index('(m s$^{-1}$)') + 12:]) if len(r) >= 6 and 'MARVELS' in r[1]]
r2 = {clean(r[0]): r for r in split_rows(t2) if len(r) >= 8 and clean(r[0]).startswith('MARVELS-')}
rows = []
for r in r1:
    star = clean(r[0]).split('(')[0].strip(); des = clean(r[1])
    ms, _ = parse_val(r[2]); P, eP = parse_val(r[3]); e, _ = parse_val(r[4]); K, _ = parse_val(r[5])
    s = r2[des]; a, _ = parse_val(s[1]); Tp, _ = parse_val(s[2]); om, _ = parse_val(s[3]); dvdt = re.sub(r'\s+', '', clean(s[5])).replace('\\pm', '+/-')
    host = ('TYC ' + '-'.join(str(int(x)) for x in star[4:].split('-'))) if star.startswith('TYC') else star
    n = flags(ms, a) + [f'{des}; MARVELS (SDSS-III, northern) candidate', f'linear trend dv/dt={dvdt} m/s/d']
    if ms > 80: n.append('authors classify as stellar')
    n.append('paper name ' + star)
    rows.append(dict(host=host, companion=des, P_d=P, e_P_d=eP, K_ms=K, e=e, omega_deg=om, Tp_bjd=Tp + 2450000,
                     msini_mj=ms, a_au=a, mstar_msun=mst.get(star), fit_type='RV-only',
                     table_ref='arXiv:1702.01784 Tables 4-5', notes='; '.join(n)))
write_csv(f'{L}/grieves_2017/grieves_2017_orbits.csv', rows); print('grieves', len(rows))
