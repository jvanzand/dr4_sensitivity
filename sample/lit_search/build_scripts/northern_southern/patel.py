import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/patel_2007/'
t = {cells(l)[0]: cells(l)[1] for l in open(R + 'src_0704.3418/tab2.tex') if '&' in l}
def g(k): return [v for kk, v in t.items() if kk.startswith(k)][0]
rows = [dict(host='HD 167665', companion='b', P_d=val(g('$P$')), e_P_d=err(g('$P$')), K_ms=val(g('$K$')), e=val(g('$e$')),
    omega_deg=val(g('$\\omega$')), Tp_bjd=addbjd(g('$T_p$'), 2440000), msini_mj=val(g('\\msini')), a_au=val(g('$a_{\\rm min}$')),
    mstar_msun='1.09', fit_type='RV-only', table_ref='arXiv:0704.3418 Table 2 (Mstar Table 1)',
    notes='Keck/HIRES RVs (23 obs). Tp is JD (not BJD) in paper: JD-2440000 converted. a slightly above 5 AU (P~12 yr): borderline for target range. Possibly also in CLS I (Keck star); check overlap.')]
write(R + 'patel_2007_orbits.csv', rows)
