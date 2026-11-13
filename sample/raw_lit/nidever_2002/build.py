"""Build nidever_2002_orbits.csv from VizieR J/ApJS/141/503 table3 (Nidever et al. 2002, ApJS 141, 503).

Run:  python -I build.py   (network needed only on first run; results cached here)
Inputs: table3.dat (15 SB1 orbits, from https://cdsarc.cds.unistra.fr/ftp/J/ApJS/141/503/),
        arxiv_src/ms.tex (astro-ph/0112477 LaTeX; used to check units/epochs, not parsed).
M2 sin i, M1 and a_min are the paper's own (Table 3; M1 from Prieto & Lambert 1999 or B-V + Allen 2000).
For comparison only, M2 sin i is recomputed from the published f(M) with Gaia DR3 FLAME (else TIC v8)
mass and reported in notes.
Keeps rows with M2,min <= 80 MJ and (a <= 5 AU or P <= 4000 d); all 15 go to nidever_2002_all.csv.
"""
import csv, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mf_common as mf  # noqa: E402

os.chdir(HERE)
REFS = {'1': 'Duquennoy & Mayor 1991', '2': 'Martin+1998', '3': 'Blazit+1987', '4': 'Kamper 1987',
        '5': 'Wagman 1949', '6': 'Beavers & Salzer 1983', '7': 'Kamper & Lyons 1981',
        '8': 'Halbwachs+2000', '9': 'Mazeh+1996', '10': 'Oetiker+2001', '11': 'Mayor+1997',
        '12': 'McCarthy 1986', '13': 'Barbieri+1996', '14': 'Henry & McCarthy 1993',
        '15': 'Marcy & Moore 1989', '16': 'McCarthy & Henry 1987', '17': 'Lippincott & Borgman 1978'}


def fnum(s):
    s = s.strip()
    return float(s) if s else None


orbits = []
for line in open('table3.dat'):
    line = line.rstrip('\n').ljust(125)
    orbits.append(dict(name=line[0:9].strip(), P=fnum(line[10:19]), eP=fnum(line[20:29]),
                       K=fnum(line[30:35]), eK=fnum(line[36:40]), e=fnum(line[41:45]),
                       ee=line[46:51].strip(), efix=line[51].strip(), om=fnum(line[53:56]),
                       T0=fnum(line[62:73]), eT0=fnum(line[74:82]), M1=fnum(line[83:87]),
                       M2=fnum(line[88:93]), amin=fnum(line[94:98]), fm=fnum(line[99:107]),
                       refs=line[108:125].strip()))

res = mf.resolve([o['name'] for o in orbits], 'simbad_cache.json')
gaia = mf.gaia_query([res[o['name']].get('gaia_dr3') for o in orbits], 'gaia_dr3_query.csv')
tic = mf.tic_query([res[o['name']].get('tic') for o in orbits], 'tic_v8_query.csv')

cols = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp', 'tp_system', 'msini_mj',
        'a_au', 'mstar_msun', 'fit_type', 'table_ref', 'notes', 'e_K_ms', 'e_e']
allrows, keep = [], []
for o in orbits:
    r = res[o['name']]
    m1g, msrc, flags = mf.primary_mass(r, gaia, tic)
    notes = ['SB1 (orbit previously unpublished per the paper; Keck/Lick/CAT velocities); M2 sin i, M1 (Prieto & Lambert 1999 or '
             'B-V with Allen 2000) and a (a_min) as published in Table 3; f(M)=%s Msun' % o['fm']]
    if o['efix'] == 'f':
        notes.append('e fixed in fit')
    if m1g:
        ms = mf.msini_from_fm(o['fm'], m1g) * mf.MJ_PER_MSUN
        notes.append('for comparison, M2 sin i = %.1f MJ from f(M) with M1 = %s' % (ms, msrc))
    notes += flags
    if o['refs']:
        notes.append('earlier work on this system (paper refs): ' +
                     ', '.join(REFS[k] for k in o['refs'].split(',')))
    notes.append('T0 listed as T0-2450000 without a stated time system (observation epochs in the paper are '
                 'labelled JD; velocities are barycentric); K in km/s converted to m/s')
    notes.append('SIMBAD %s (otype %s)' % (r.get('main_id', 'unresolved'), r.get('otype', '')))
    notes.append('e_K_ms, e_e: rms uncertainties on K (km/s converted to m/s) and e as published in Table 3 '
                 '(e_e blank where e was fixed); the paper gives no uncertainty on M2 sin i or a_min, so '
                 'e_msini_mj and e_a_au are not provided')
    row = dict(host=o['name'], companion=o['name'] + ' B (Nidever+2002 SB1)', P_d=o['P'], e_P_d=o['eP'],
               K_ms=round(o['K'] * 1000, 1), e=o['e'], omega_deg=o['om'], Tp=o['T0'], tp_system='unknown',
               msini_mj=o['M2'], a_au=o['amin'], mstar_msun=o['M1'], fit_type='RV-only',
               table_ref='Nidever+2002 ApJS 141, 503, Table 3 (VizieR J/ApJS/141/503 table3; arXiv:astro-ph/0112477)',
               notes='; '.join(notes),
               e_K_ms=round(o['eK'] * 1000, 1) if o['eK'] else None,  # km/s -> m/s
               e_e=float(o['ee']) if o['ee'] else None)
    allrows.append(row)
    if o['M2'] <= 80 and (o['amin'] <= 5 or o['P'] <= 4000):
        keep.append(row)

for fn, rows in [('nidever_2002_all.csv', allrows), ('nidever_2002_orbits.csv', keep)]:
    with open(fn, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
print('orbits', len(orbits), 'kept', len(keep), 'BD 13-80', sum(1 for r in keep if 13 <= r['msini_mj'] <= 80))
