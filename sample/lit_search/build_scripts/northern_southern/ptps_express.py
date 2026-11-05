import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
B = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/'
mjd = lambda x: '%.2f' % (float(x) + 2400000.5)
t0note = ' T0 given in MJD (converted to JD); the paper does not define T0 explicitly (PTPS convention is periastron passage): verify before use.'
ptps = [
 dict(host='BD+20 2457', companion='b', P_d='379.63', e_P_d='2.01', K_ms='322.35', e='0.15', omega_deg='207.64', Tp_bjd=mjd(54677.03),
      msini_mj='21.42', a_au='1.45', mstar_msun='', fit_type='RV-only',
      table_ref='arXiv:0906.1804 (Niedzielski+2009, ApJ 707, 768) Table 4 (stellar Table 1)',
      notes='HET/HRS RVs, two-Keplerian fit. Host mass only given as range 1.3-4.3 Msun (metal-poor K2 II, [Fe/H]=-1.0), so msini is model-dependent. Quirrenbach+2019 notes no stable configuration found (Horner+2014, Trifonov+2014): reality doubtful. Likely in NASA Exoplanet Archive.' + t0note),
 dict(host='BD+20 2457', companion='c', P_d='621.99', e_P_d='10.20', K_ms='160.03', e='0.18', omega_deg='126.02', Tp_bjd=mjd(53866.95),
      msini_mj='12.47', a_au='2.01', mstar_msun='', fit_type='RV-only',
      table_ref='arXiv:0906.1804 (Niedzielski+2009, ApJ 707, 768) Table 4 (stellar Table 1)',
      notes='Massive planet; see BD+20 2457 b notes (dynamical stability doubtful).' + t0note),
 dict(host='HD 4760', companion='b', P_d='434', e_P_d='3', K_ms='370', e='0.23', omega_deg='265', Tp_bjd=mjd(53955),
      msini_mj='13.9', a_au='1.14', mstar_msun='1.05', fit_type='RV-only',
      table_ref='arXiv:2101.10410 (Niedzielski+2021, A&A 648, A58; TAPAS VII) Table TableKeplerian (Mstar: basic parameters table)',
      notes='HET/HRS + TNG/HARPS-N RVs (35 epochs); metal-poor giant ([Fe/H]=-0.91), jitter 64 m/s.' + t0note),
 dict(host='BD+02 3313', companion='b', P_d='1393', e_P_d='3', K_ms='690', e='0.47', omega_deg='351.3', Tp_bjd=mjd(54982),
      msini_mj='34.1', a_au='2.47', mstar_msun='1.03', fit_type='RV-only',
      table_ref='arXiv:2101.10410 (Niedzielski+2021, A&A 648, A58; TAPAS VII) Table TableKeplerian (Mstar: basic parameters table)',
      notes='HET/HRS + TNG/HARPS-N RVs (29 epochs). msini > 30 MJ so not in NASA Exoplanet Archive scope.' + t0note),
 dict(host='TYC 434-4538-1', companion='b', P_d='193.2', e_P_d='0.4', K_ms='209', e='0.08', omega_deg='196', Tp_bjd=mjd(54829),
      msini_mj='6.1', a_au='0.66', mstar_msun='1.04', fit_type='RV-only',
      table_ref='arXiv:2101.10410 (Niedzielski+2021, A&A 648, A58; TAPAS VII) Table TableKeplerian (Mstar: basic parameters table)',
      notes='Massive planet (6.1 MJ); likely in NASA Exoplanet Archive.' + t0note),
]
write(B + 'ptps/ptps_orbits.csv', ptps)
jd = lambda x, off: '%.1f' % (float(x) + off)
express = [
 dict(host='HIP 97233', companion='b', P_d='1058.8', e_P_d='6.7', K_ms='320.1', e='0.61', omega_deg='249.3', Tp_bjd=jd(856.3, 2455000),
      msini_mj='20.0', a_au='2.55', mstar_msun='1.93', fit_type='RV-only',
      table_ref='arXiv:1409.7429 (Jones+2015, A&A 573, A3) Table "Orbital parameters of HIP 67851 b and HIP 97233 b" (Mstar Table "Stellar parameters")',
      notes='EXPRESS (FEROS + CHIRON; southern survey). Tp JD-2455000 in paper (JD). Likely in NASA Exoplanet Archive.'),
 dict(host='HIP 67537', companion='b', P_d='2556.5', e_P_d='97.05', K_ms='112.7', e='0.59', omega_deg='119.6', Tp_bjd=jd(6290.6, 2450000),
      msini_mj='11.1', a_au='4.91', mstar_msun='2.41', fit_type='RV-only',
      table_ref='arXiv:1612.06252 (Jones+2017, A&A 602, A58) Table "Stellar properties and orbital elements"',
      notes='EXPRESS (FEROS + CHIRON). Massive planet/BD boundary (11.1 +0.4/-1.1 MJ); a=4.91 AU near the 5 AU limit. Likely in NASA Exoplanet Archive.'),
]
write(B + 'express_jones/express_jones_orbits.csv', express)
