import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/lick_giants/'
ep = 2451853.595
tp = lambda M0, P: '%.2f' % (ep - M0 / 360.0 * P)
rows = [
 dict(host='nu Oph', companion='b', P_d='530.02', e_P_d='0.11', K_ms='288.27', e='0.124', omega_deg='9.88', Tp_bjd=tp(235.78, 530.02),
      msini_mj='22.204', a_au='1.789', mstar_msun='2.7', fit_type='RV-only',
      table_ref='arXiv:1904.03557 (Quirrenbach+2019, A&A 624, A18) table:orb_par_stable, Keplerian block (Mstar HB model, table:phys_param)',
      notes='Lick + OAO RVs, double-Keplerian. Paper gives mean anomaly M0=235.78 deg at first epoch JD 2451853.595; Tp derived here as epoch - M0/360*P (JD). Likely in NASA Exoplanet Archive (<30 MJ). 6:1 MMR with c.'),
 dict(host='nu Oph', companion='c', P_d='3183.02', e_P_d='5.89', K_ms='176.78', e='0.180', omega_deg='8.57', Tp_bjd=tp(222.75, 3183.02),
      msini_mj='24.674', a_au='5.929', mstar_msun='2.7', fit_type='RV-only',
      table_ref='arXiv:1904.03557 (Quirrenbach+2019, A&A 624, A18) table:orb_par_stable, Keplerian block (Mstar HB model, table:phys_param)',
      notes='a=5.93 AU slightly beyond 5 AU (P=8.7 yr). M0=222.75 deg at JD 2451853.595; Tp derived (JD). Likely in NASA Exoplanet Archive (<30 MJ).'),
 dict(host='tau Gem', companion='b', P_d='305.5', e_P_d='0.1', K_ms='350.2', e='0.031', omega_deg='137.4', Tp_bjd='2453270.7',
      msini_mj='20.6', a_au='1.17', mstar_msun='2.3', fit_type='RV-only',
      table_ref='arXiv:1305.5107 (Mitchell+2013, A&A) Table "Best Fit Orbital Parameters" (Mstar Table "Stellar Parameters")',
      notes='Lick RVs; T0 = periastron time, JD-2450000 in paper (JD, not BJD); poorly constrained (+/-14.8 d) since e~0.03. Quirrenbach+2019 lists msini 22.1 with Stock+2018 Mstar 2.47. Likely in NASA Exoplanet Archive (<30 MJ).'),
]
write(R + 'lick_giants_orbits.csv', rows)
