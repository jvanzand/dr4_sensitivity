"""Recover SB9 SB1 orbits that the first pass set aside for lack of a primary mass.

Run:  python -I build.py   (network needed only on first run; results cached here)
Inputs (read-only, from the first SB9 pass in ../sb9/; VizieR B/sb9 via TAPVizieR):
  ../sb9/sb9_excluded_no_primary_mass.csv  (242 SB1 orbits, P <= 4000 d, M2 sin i(M1=1) < 200 MJ)
  ../sb9/tap_orbits.csv, ../sb9/tap_main.csv, ../sb9/tap_alias.csv
Primary mass: Gaia DR3 FLAME (source_id from SB9 aliases, else SIMBAD), else TIC v8 (IV/39/tic82,
TIC id from SIMBAD, else TIC row whose GAIA (DR2) id matches).
Outputs: simbad_cache.json, gaia_dr3_query.csv, tic_v8_query.csv, tic_v8_by_gaiadr2_query.csv (cached),
  sb9_tic_all.csv (all 242 with whatever mass was found), sb9_tic_orbits.csv (M2 sin i <= 80 MJ, a<=5 AU or P<=4000 d)
"""
import csv, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mf_common as mf  # noqa: E402

os.chdir(HERE)
SB9 = os.path.join('..', 'sb9')


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


excl = list(csv.DictReader(open(os.path.join(SB9, 'sb9_excluded_no_primary_mass.csv'))))
orb = {(r['Seq'], r['o']): r for r in csv.DictReader(open(os.path.join(SB9, 'tap_orbits.csv')))}
main = {r['Seq']: r for r in csv.DictReader(open(os.path.join(SB9, 'tap_main.csv')))}
alias = {}
for a in csv.DictReader(open(os.path.join(SB9, 'tap_alias.csv'))):
    alias.setdefault(a['Seq'], []).append(' '.join(a['Name'].split()))


def norm(n):
    n = ' '.join(n.split())
    n = re.sub(r'^GAIADR([23]) ', r'Gaia DR\1 ', n)
    return re.sub(r'^(BD|CD|CPD)\s+([+-])\s*', r'\1\2', n)


# candidate SIMBAD names: SB9 main name, then every SB9 alias
cand = {x['sb9_seq']: [norm(main[x['sb9_seq']]['Name'])] + [norm(a) for a in alias.get(x['sb9_seq'], [])]
        for x in excl}
res = mf.resolve([c[0] for c in cand.values()], 'simbad_cache.json')
retry = [n for c in cand.values() if not res[c[0]].get('resolved') for n in c[1:]]
if retry:
    res = mf.resolve([c[0] for c in cand.values()] + retry, 'simbad_cache.json')
names = {}
for seq, c in cand.items():
    ok = [n for n in c if res.get(n, {}).get('resolved')]
    names[seq] = ok[0] if ok else c[0]

for x in excl:
    seq = x['sb9_seq']
    r = res.get(names[seq], {})
    g3 = [a.split()[1] for a in alias.get(seq, []) if a.startswith('GAIADR3 ')]
    g2 = [a.split()[1] for a in alias.get(seq, []) if a.startswith('GAIADR2 ')]
    x['gaia_dr3'] = g3[0] if g3 else r.get('gaia_dr3', '')
    x['gaia_dr2'] = g2[0] if g2 else r.get('gaia_dr2', '')
    x['tic'] = r.get('tic', '')
    x['simbad'] = r

gaia = mf.gaia_query([x['gaia_dr3'] for x in excl], 'gaia_dr3_query.csv')
tic = mf.tic_query([x['tic'] for x in excl], 'tic_v8_query.csv')

# TIC fallback by Gaia DR2 id for stars without a TIC alias in SIMBAD
need = [x['gaia_dr2'] for x in excl if not x['tic'] and x['gaia_dr2']]
if need and not os.path.exists('tic_v8_by_gaiadr2_query.csv'):
    q = ('SELECT TIC, Mass, s_Mass, logg, Teff, Rad, LClass, Dist, GAIA FROM "IV/39/tic82" WHERE GAIA IN (' +
         ','.join(sorted(set(need))) + ')')
    t = mf._tap('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync', q)
    open('tic_v8_by_gaiadr2_query.csv', 'w').write(t)
if os.path.exists('tic_v8_by_gaiadr2_query.csv'):
    byg2 = {r['GAIA']: r for r in csv.DictReader(open('tic_v8_by_gaiadr2_query.csv'))}
    for x in excl:
        if not x['tic'] and x['gaia_dr2'] in byg2:
            x['tic'] = byg2[x['gaia_dr2']]['TIC']
            tic[x['tic']] = byg2[x['gaia_dr2']]

cols = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp', 'tp_system', 'msini_mj',
        'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes']
allrows, keep = [], []
for x in excl:
    seq, o = x['sb9_seq'], x['orbit']
    ob, mr, r = orb[(seq, o)], main[seq], x['simbad']
    P, K, e = fnum(ob['Per']), fnum(ob['K1']), fnum(ob['e']) or 0.0
    fm = mf.fm_from_orbit(P, K, e)
    ids = {'gaia_dr3': x['gaia_dr3'], 'tic': x['tic']}
    m1, msrc, flags = mf.primary_mass(ids, gaia, tic)
    sp = mr['Sp1'].strip()
    if re.search(r'(III|II|Ib|Ia)', sp.replace('IV', '')):
        flags.append('SB9 SpT %s is giant/supergiant (TIC/FLAME dwarf-based mass may be unreliable)' % sp)
    if re.match(r'^(D[ABOQZC]|sd[BO]|sbB|WD)', sp) or '+' in sp or r.get('otype') in ('HS*', 'WD*', 'CV*', 'Psr', 'XB*', 'LXB', 'HXB'):
        flags.append('SB9 SpT %s / SIMBAD otype %s suggests compact, hot-subdwarf or composite host; check' % (sp, r.get('otype')))
    if not r.get('resolved'):
        flags.append('SIMBAD did not resolve %s' % names[seq])
    ms = a = None
    if m1:
        ms = mf.msini_from_fm(fm, m1)
        a = mf.a_au(P, m1 + ms)
    tp, tps = ob['T0'], 'JD'
    tnote = 'Tp = SB9 T0 (VizieR B/sb9 labels it periastron time, JD; original paper may use HJD/BJD)'
    if fnum(ob['f_T0']):
        tp, tps = '', ''
        tnote = 'SB9 T0 flagged as not a periastron epoch (f_T0=%s); Tp left blank' % ob['f_T0']
    elif e == 0:
        tnote += '; e=0 so T0 is a reference epoch, not a true periastron'
    host = r['main_id'] if r.get('main_id') else names[seq]
    notes = ['SB9 Seq %s orbit %s (SB9 name %s); SB1; f(m)=1.0361e-7(1-e^2)^1.5 K1^3 P = %.3g Msun; M2 sin i solved '
             'from f(m)=(M2 sin i)^3/(M1+M2 sin i)^2 (minimum-mass approximation); first-pass exclusion reason: %s'
             % (seq, o, names[seq], fm, x['reason_excluded'])]
    notes.append('M1 = ' + msrc if m1 else 'no Gaia DR3 FLAME or TIC v8 mass found')
    notes += flags
    notes.append('SB9 grade=%s; Sp1=%s; orbit ref %s; %s' % (ob['Grade'] or 'none', sp, ob['Ref'], tnote))
    for fl, lab in [('f_Per', 'P'), ('f_e', 'e'), ('f_omega', 'omega'), ('f_K1', 'K1'), ('u_K1', 'K1')]:
        if ob.get(fl, '').strip():
            notes.append('SB9 flag %s on %s' % (ob[fl].strip(), lab))
    notes.append('a from Kepler III with M1+M2 sin i; SIMBAD %s (otype %s)' % (r.get('main_id', '-'), r.get('otype', '-')))
    row = dict(host=host, companion='%s (SB9 %s.%s)' % (host, seq, o), P_d=P, e_P_d=ob['e_Per'],
               K_ms=round(K * 1000, 2), e=ob['e'], omega_deg=ob['omega'], Tp=tp, tp_system=tps,
               msini_mj=round(ms * mf.MJ_PER_MSUN, 2) if ms else '', a_au=round(a, 4) if a else '',
               mstar_msun=round(m1, 3) if m1 else '', fit_type='RV-only',
               table_ref='SB9 (Pourbaix+2004; VizieR B/sb9 as cached in ../sb9/tap_orbits.csv); orbit ref ' + ob['Ref'],
               notes='; '.join(notes))
    allrows.append(row)
    if ms and ms * mf.MJ_PER_MSUN <= 80 and (a <= 5 or P <= 4000):
        keep.append(row)

for fn, rows in [('sb9_tic_all.csv', allrows), ('sb9_tic_orbits.csv', keep)]:
    with open(fn, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
print('input', len(excl), 'resolved', sum(1 for x in excl if x['simbad'].get('resolved')),
      'gaia_dr3', sum(1 for x in excl if x['gaia_dr3']), 'tic', sum(1 for x in excl if x['tic']),
      'with M1', sum(1 for r in allrows if r['mstar_msun'] != ''), 'kept', len(keep),
      'BD 13-80', sum(1 for r in keep if 13 <= r['msini_mj'] <= 80))
