import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
import sys
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/wilson_2016/'
L = open(R + 'src_1602.02749/rv.tex').read().split('\n')
mstar = {}
for l in L[176:192]:
    c = cells(l)
    if len(c) >= 10:
        mstar[c[0].split('$')[0].replace('\\,', ' ').replace('~', ' ').strip()] = val(c[9])
rows = []
skip = {'HD 50761'}
notes_extra = {
 'HD 6664': 'Wilson+2016 class: low-mass binary; Hipparcos astrometric orbit gives M2=0.35+/-0.03 Msun (stellar). msini 80+/-5.',
 'BD+73 0275': 'Wilson+2016 class: low-mass binary; Hipparcos astrometric orbit gives M2=0.20+/-0.03 Msun (stellar).',
 'HD 51813': 'Wilson+2016 class: low-mass binary; Hipparcos astrometric orbit gives M2=0.27+/-0.07 Msun (stellar).',
 'HD 56709': 'Wilson+2016 class: low-mass binary candidate; tentative Hipparcos orbit M2=0.54+/-0.12 Msun. msini 72+/-4.',
 'HD 134169': 'Wilson+2016 class: low-mass binary; msini 83+/-4 (just above 80 MJ); no Hipparcos detection (m2,max 1.03 Msun).',
 'HD 77065': 'Includes Latham+2002 RVs (as G9-42).',
 'HD 134113': 'Includes Latham+2002 RVs (as G66-65). Mstar from Sousa+2011.',
}
for l in L[244:270]:
    if '\\pm' not in l or l.strip().startswith('%'):
        continue
    c = cells(l)
    name = c[0].split('$')[0].replace('\\,', ' ').replace('~', ' ').strip()
    host, comp = name[:-1].strip(), name[-1]
    if host in skip:
        continue
    note = 'SOPHIE RVs; single-Keplerian fit; T0 is periastron time (BJD-2400000 in paper, converted). Separate Hipparcos IAD analysis does not alter RV parameters.'
    if host in notes_extra:
        note += ' ' + notes_extra[host]
    rows.append(dict(host=bd(host), companion=comp, P_d=val(c[1]), e_P_d=err(c[1]), K_ms=val(c[4]), e=val(c[2]),
        omega_deg=val(c[3]), Tp_bjd=addbjd(c[5]), msini_mj=val(c[11]), a_au=val(c[10]),
        mstar_msun=mstar.get(host, ''), fit_type='RV-only', table_ref='arXiv:1602.02749 Table 3 (Mstar Table 2)', notes=note))
write(R + 'wilson_2016_orbits.csv', rows)
