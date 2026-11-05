import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
import math
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/teng_2023/'
L = open(R + 'src_2308.05343/Revisiting.tex').read().split('\n')
star = {}
for l in L[195:240]:
    c = cells(l)
    if len(c) >= 14:
        star[c[0].replace('$\\nu$', 'nu').strip()] = val(c[12])
orb = {}
i0 = [i for i, l in enumerate(L) if 'label{tab:orbpar_orbits}' in l][0]
for l in L[i0 - 60:i0]:
    c = cells(l)
    if len(c) >= 8 and '_{' in l:
        orb[c[0].replace('$\\nu$', 'nu').replace('\\textbf{', '').replace('}', '').strip()] = c[1:]
keep = {'11 Com b': 'Liu+2008 discovery (Xinglong); also Lick BD list.',
        'nu Oph b': 'Also in Quirrenbach+2019 (lick_giants).', 'nu Oph c': 'a=5.83 AU slightly >5 AU. Also in Quirrenbach+2019.',
        'HD 175679 b': 'Wang+2012 discovery. Joint RV+HGCA fit (table "Joint fitting results") gives true mass 61.8 (+15/-9.4) MJ, i=39 or 152 deg.',
        'HD 180314 b': 'Sato+2010 discovery.',
        'HD 14067 b': 'msini 8.8 MJ (massive planet) but joint RV+HGCA fit gives true mass 15.3 (+6.2/-4.5) MJ; a=5.11 AU borderline.'}
rows = []
for name, note in keep.items():
    P, K, e, w, Tp, m, a = orb[name][:7]
    host, comp = name.rsplit(' ', 1)
    rows.append(dict(host=host, companion=comp, P_d=val(P), e_P_d=err(P), K_ms=val(K), e=val(e),
        omega_deg='%.2f' % math.degrees(float(val(w))), Tp_bjd=addbjd(Tp, 2450000), msini_mj=val(m), a_au=val(a),
        mstar_msun=star[host], fit_type='RV-only',
        table_ref='arXiv:2308.05343 (Teng+2023, PASJ) tab:orbpar_orbits (Mstar tab:star_pars)',
        notes='Okayama Planet Search (HIDES) + literature RVs, RV-only refit. omega given in radians, converted; Tp JD-2450000 in paper (JD). Likely in NASA Exoplanet Archive if <30 MJ. ' + note))
write(R + 'teng_2023_orbits.csv', rows)
