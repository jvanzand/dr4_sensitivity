import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/okayama_eapsnet/'
rows = [dict(host='HD 119445', companion='b', P_d='410.2', e_P_d='0.6', K_ms='413.5', e='0.082', omega_deg='160.5',
    Tp_bjd='2452873.9', msini_mj='37.6', a_au='1.71', mstar_msun='3.9', fit_type='RV-only',
    table_ref='arXiv:0906.3762 (Omiya+2009, PASJ 61, 825) Table "Orbital parameters of HD 119445 b" (Mstar Table "Stellar parameters")',
    notes='BOAO + OAO (Korean-Japanese EAPS-Net) RVs, 36 obs. T is periastron time in JD (not BJD). Massive (3.9 Msun) G6III host at 289 pc. Not in Teng+2023 revisit list.')]
write(R + 'okayama_eapsnet_orbits.csv', rows)
