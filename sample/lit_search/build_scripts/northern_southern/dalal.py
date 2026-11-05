import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/dalal_2021/'
L = open(R + 'src_2105.09741/aa.tex').read().split('\n')
i0 = [i for i, l in enumerate(L) if 'label{tab:BDorbit}' in l][0]
tab = {}
for l in L[i0:i0 + 30]:
    c = cells(l)
    if len(c) == 6:
        tab[c[0].strip()] = c[2:]
names = tab['Name']
mstar = {}
j0 = [i for i, l in enumerate(L) if 'label{tab:stellar}' in l][0]
for l in L[j0:j0 + 40]:
    c = cells(l)
    if len(c) >= 7 and '\\pm' in l:
        mstar[c[0].replace('~', ' ').strip()] = val(c[6].split('$')[0])
def get(key, k):
    for kk in tab:
        if kk.startswith(key):
            return tab[kk][k]
extra = {'HD 205521': ' Hipparcos astrometric orbit (>3 sigma) gives Mc=0.48+/-0.04 Msun and Gaia DR1 excess-noise analysis 0.20 Msun: likely stellar, not a BD.'}
rows = []
for k, n in enumerate(names):
    host = bd(n.replace('~', ' ').strip())
    note = 'SOPHIE/SOPHIE+ RVs; Yorbit single-Keplerian fit; Tp is periastron (BJD-2450000 in paper, converted).' + extra.get(n.replace('~', ' ').strip(), '')
    rows.append(dict(host=host, companion='b', P_d=val(get('P', k)), e_P_d=err(get('P', k)), K_ms=val(get('K', k)),
        e=val(get('$e$', k)), omega_deg=val(get('$\\omega$', k)), Tp_bjd=addbjd(get('$T_p', k), 2450000),
        msini_mj=val(get('\\msini', k)), a_au=val(get('$a^', k)), mstar_msun=mstar.get(n.replace('~', ' ').strip(), ''),
        fit_type='RV-only', table_ref='arXiv:2105.09741 Table 2 (BD orbits; Mstar Table A.2)', notes=note))
write(R + 'dalal_2021_orbits.csv', rows)
