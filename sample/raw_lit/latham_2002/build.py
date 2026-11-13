"""Build latham_2002_orbits.csv from VizieR J/AJ/124/1144 (Latham et al. 2002, AJ 124, 1144).

Run from anywhere:  python -I build.py   (network needed only on first run; results cached here)
Inputs (downloaded from https://cdsarc.cds.unistra.fr/ftp/J/AJ/124/1144/):
  ReadMe, table1.dat (names in Paper XII), table4.dat (156 robust orbits), table5.dat (15 preliminary)
Outputs:
  simbad_cache.json, gaia_dr3_query.csv, tic_v8_query.csv  (cached network results)
  latham_2002_all_msini.csv  (every orbit with M1 and M2 sin i, for provenance)
  latham_2002_orbits.csv     (rows with M2 sin i <= 80 MJ and (a <= 5 AU or P <= 4000 d))
"""
import csv, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mf_common as mf  # noqa: E402

os.chdir(HERE)


def fnum(s):
    s = s.strip()
    return float(s) if s else None


# Paper XII names from table1 (Name -> PXIIname)
px = {}
for line in open('table1.dat'):
    px[line[0:8].strip()] = line[54:64].strip()

orbits = []
for fn, kind in [('table4.dat', 'robust'), ('table5.dat', 'preliminary')]:
    for line in open(fn):
        line = line.rstrip('\n').ljust(206)
        o = dict(name=line[0:8].strip(), kind=kind,
                 P=fnum(line[30:43]), eP=fnum(line[43:55]), gamma=fnum(line[56:64]),
                 K=fnum(line[72:78]), eK=fnum(line[79:85]), e=fnum(line[86:92]),
                 om=fnum(line[100:106]), t0=fnum(line[114:123]), et0=fnum(line[124:132]),
                 asini=fnum(line[133:141]), fm=fnum(line[152:163]), efm=fnum(line[164:172]),
                 nobs=line[173:176].strip(), MA=fnum(line[196:200]), MBmin=fnum(line[201:206]))
        o['ee'] = fnum(line[93:99])  # e_e, published uncertainty on e
        o['ra'] = line[9:19].strip()
        o['dec'] = line[20:29].strip()
        orbits.append(o)


def simbad_name(n):
    """Paper XII name -> SIMBAD-style identifier (Giclas 'G130-32' -> 'G 130-32')."""
    import re
    m = re.match(r'^G(\d+)-(\d+[A-Z]?)$', n)
    if m:
        return 'G %s-%s' % m.groups()
    return n


cands = {}
for o in orbits:
    p = px.get(o['name'], o['name'])
    o['pxname'] = p
    cands[o['name']] = [simbad_name(p)]
    if simbad_name(p)[-1] in 'AB' and simbad_name(p)[-2].isdigit():
        cands[o['name']].append(simbad_name(p)[:-1])  # fall back to name without component letter
allnames = [n for v in cands.values() for n in v]
res = mf.resolve(allnames, 'simbad_cache.json')
for o in orbits:
    o['res'] = {}
    for n in cands[o['name']]:
        if res.get(n, {}).get('resolved'):
            o['res'] = res[n]
            o['query_name'] = n
            break

gaia = mf.gaia_query([o['res'].get('gaia_dr3') for o in orbits], 'gaia_dr3_query.csv')
tic = mf.tic_query([o['res'].get('tic') for o in orbits], 'tic_v8_query.csv')

cols = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp', 'tp_system', 'msini_mj',
        'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes', 'e_K_ms', 'e_e']
allrows, keep = [], []
for o in orbits:
    r = o['res']
    m1, msrc, flags = mf.primary_mass(r, gaia, tic)
    if m1 is None and o['MA'] is not None:  # third choice: the paper's own Table 6 primary mass
        m1, msrc = o['MA'], 'Latham+2002 Table 6 M(A) (no Gaia DR3 FLAME or TIC v8 mass available)'
    host = r['main_id'] if r.get('main_id') else simbad_name(o['pxname'])
    ms = a = None
    if m1:
        ms = mf.msini_from_fm(o['fm'], m1)
        a = mf.a_au(o['P'], m1 + ms)
    tab = 'Table 4 (robust orbits)' if o['kind'] == 'robust' else 'Table 5 (preliminary orbits)'
    notes = ['Latham+2002 star %s (Paper XII name %s)' % (o['name'], o['pxname'])]
    if o['kind'] == 'preliminary':
        notes.append('PRELIMINARY orbit (Table 5)')
    notes.append('SB1; f(M)=%s+-%s Msun as published; M2sini solved from f(M)=(M2 sini)^3/(M1+M2 sini)^2 '
                 '(minimum-mass approximation)' % (o['fm'], o['efm']))
    if m1:
        notes.append('M1 = ' + msrc)
    else:
        notes.append('no FLAME or TIC v8 mass found')
    if o['MA'] is not None:
        notes.append('paper Table 6 (substellar suspected): M(A)=%s Msun, M(B)min=%s Msun' % (o['MA'], o['MBmin']))
    notes += flags
    notes.append('Tp = HJD of periastron (VizieR t0 + 2400000); K and gamma on CfA native system; '
                 'Nobs=%s; SIMBAD %s (otype %s)' % (o['nobs'], r.get('main_id', 'unresolved'), r.get('otype', '')))
    notes.append('e_K_ms, e_e: published uncertainties on K (km/s converted to m/s) and e (blank when 0); '
                 'e_msini_mj and e_a_au are blank because M2 sin i and a here are solved from f(M) and the '
                 'adopted M1, not published (paper gives only e_f(M), quoted above)')
    if not r:
        notes.append('SIMBAD/Sesame did not resolve %s' % cands[o['name']])
    row = dict(host=host, companion=host + ' B (Latham+2002 SB1)', P_d=o['P'], e_P_d=o['eP'],
               K_ms=round(o['K'] * 1000, 1), e=o['e'], omega_deg=o['om'],
               Tp=round(o['t0'] + 2400000, 3) if o['t0'] is not None else '',
               tp_system='HJD', msini_mj=round(ms * mf.MJ_PER_MSUN, 2) if ms else '',
               a_au=round(a, 4) if a else '', mstar_msun=round(m1, 3) if m1 else '', fit_type='RV-only',
               table_ref='Latham+2002 AJ 124, 1144; VizieR J/AJ/124/1144 ' + tab, notes='; '.join(notes),
               e_K_ms=round(o['eK'] * 1000, 1) if o['eK'] else None,  # km/s -> m/s
               e_e=o['ee'] if o['ee'] else None)
    allrows.append(row)
    if ms and ms * mf.MJ_PER_MSUN <= 80 and (a <= 5 or o['P'] <= 4000):
        keep.append(row)

for fn, rows in [('latham_2002_all_msini.csv', allrows), ('latham_2002_orbits.csv', keep)]:
    with open(fn, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
print('orbits', len(orbits), 'resolved', sum(1 for o in orbits if o['res']),
      'with Gaia DR3 id', sum(1 for o in orbits if o['res'].get('gaia_dr3')),
      'with M1', sum(1 for r in allrows if r['mstar_msun'] != ''), 'kept', len(keep),
      'BD 13-80', sum(1 for r in keep if 13 <= r['msini_mj'] <= 80))
for o in orbits:
    if not o['res']:
        print('UNRESOLVED', o['name'], o['pxname'])
