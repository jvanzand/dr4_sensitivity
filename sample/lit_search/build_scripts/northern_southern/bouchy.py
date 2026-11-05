import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/bouchy_2016/'
L = open(R + 'src_1511.08397/26347.tex').read().split('\n')
i0 = [i for i, l in enumerate(L) if '\\label{fit}' in l][0]
tab = {}
for l in L[i0:i0 + 25]:
    c = cells(l.replace('[+3pt]', ''))
    if len(c) == 6:
        tab[c[0]] = c[1:]
names = [x.replace('\\object{', '').replace('}', '') for x in tab['Parameters']]
j0 = [i for i, l in enumerate(L) if '\\label{spectro}' in l][0]
ms = {}
for l in L[j0:j0 + 15]:
    c = cells(l)
    if len(c) == 8 and 'object' in c[0]:
        ms[c[0].replace('\\object{', '').replace('}', '')] = val(c[5])
def g(k, i):
    for kk in tab:
        if kk.startswith(k):
            return tab[kk][i]
rows = []
for i, n in enumerate(names):
    a = float(val(g('$a$ [AU]', i)))
    if a > 7:
        continue
    Py = g('P [years]', i)
    P = '%.2f' % (float(val(Py)) * 365.25)
    eP = '%.2f' % (float(err(Py)) * 365.25)
    rows.append(dict(host='HD ' + n[2:], companion='b', P_d=P, e_P_d=eP, K_ms=val(g('K [', i)), e=val(tab['e'][i]),
        omega_deg=val(g('$\\omega$', i)), Tp_bjd=addbjd(g('$T_p$', i), 2450000), msini_mj=val(g('Mc\\,\\sini', i)),
        a_au=val(g('$a$ [AU]', i)), mstar_msun=ms[n], fit_type='RV-only', table_ref='arXiv:1511.08397 Table 4 (Mstar Table 2)',
        notes='OUTSIDE a<=5 AU target range (a=%s AU, P~%.0f yr); kept as borderline. ELODIE+SOPHIE RVs, Yorbit+MCMC; P given in years in paper (converted with 365.25 d/yr). No Hipparcos detection; Hipparcos upper limit %s Msun.' % (val(g('$a$ [AU]', i)), float(val(Py)), val(g('Mc,upper', i)))))
write(R + 'bouchy_2016_orbits.csv', rows)
