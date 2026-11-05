import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
import math
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/marvels_bd/'
deg = lambda r: '%.1f' % math.degrees(float(r))
# Values transcribed from the LaTeX tables printed in the session (sources in src_*/).
rows = [
 dict(host='TYC 2949-557-1', companion='b', P_d='5.69459', e_P_d='0.00029', K_ms='6109', e='0.0017', omega_deg=deg('4.69'),
      Tp_bjd='', msini_mj='64.3', a_au='', mstar_msun='1.25', fit_type='RV-only',
      table_ref='arXiv:1006.1627 (Fleming+2010, ApJ 718, 1186) Table 4 eccentric case (Mstar Table 5)',
      notes='MARVELS pilot (KeckET) + KPNO ET + HET/HRS RVs. Only Tc (time of conjunction) given: Tc = 2454868.9878 BJD_TDB; Tp left blank (e~0, omega 4.69+/-1.58 rad unconstrained). Circular fit essentially identical.'),
 dict(host='TYC 2930-872-1', companion='b', P_d='2.430420', e_P_d='0.000006', K_ms='8723', e='0.0066', omega_deg='142',
      Tp_bjd='', msini_mj='68.1', a_au='', mstar_msun='1.21', fit_type='RV-only',
      table_ref='arXiv:1206.5514 (Fleming+2012, AJ 144, 72) orbprops table (Mstar stellarprops table); msini from Sect. 4.2 text',
      notes='MARVELS-2b. MARVELS + TNG/SARG RVs. Only Tc given: Tc = HJD_UTC 2454842.2640; Tp blank. Long-term trend from tertiary stellar companion (P>2000 d). Paper argues true mass likely stellar (median M2 0.34-0.42 Msun if synchronized/coplanar).'),
 dict(host='TYC 2087-255-1', companion='b', P_d='9.0090', e_P_d='0.0004', K_ms='3571', e='0.226', omega_deg=deg('4.086'),
      Tp_bjd='2455552.797', msini_mj='40.0', a_au='0.090', mstar_msun='1.16', fit_type='RV-only',
      table_ref='arXiv:1211.6140 (Ma+2013, AJ 145, 20) companionparams table, MARVELS+TNG column (Mstar stellarparams table)',
      notes='MARVELS-4b. Tp is HJD_UTC. omega given in radians (4.086+/-0.041), converted. Active subgiant host.'),
 dict(host='HIP 67526', companion='b', P_d='90.2695', e_P_d='0.01875', K_ms='2948.14', e='0.4375', omega_deg='-140.91',
      Tp_bjd='2455563.73', msini_mj='65.0', a_au='', mstar_msun='1.10', fit_type='RV-only',
      table_ref='arXiv:1307.3308 (Jiang+2013, AJ 146, 65) Table 2 (Mstar Table 3); msini from abstract',
      notes='MARVELS-5b. MARVELS + TNG/SARG RVs. Tp is HJD. msini = 65.0+/-2.9 MJ (sin i=1); with isotropic/ mass-ratio priors median mass 76-244 MJ, so may be stellar. Wilson+2016 compilation lists 62.6 MJ.'),
 dict(host='GSC 03546-01452', companion='b', P_d='47.8929', e_P_d='0.00625', K_ms='1644', e='0.1442', omega_deg=deg('1.998'),
      Tp_bjd='2455025.81', msini_mj='31.7', a_au='', mstar_msun='1.11', fit_type='RV-only',
      table_ref='arXiv:1304.2815 (De Lee+2013, AJ 145, 155) orbittab table (Mstar hosttab table)',
      notes='MARVELS-6b. MARVELS + TNG/SARG RVs; Tp in BJD_TDB; omega given in radians (1.998), converted. Residual linear slope 0.180 m/s/d.'),
 dict(host='HD 87646', companion='b', P_d='13.481', e_P_d='0.001', K_ms='956', e='0.05', omega_deg=deg('5.20'),
      Tp_bjd='2454088.3', msini_mj='12.4', a_au='0.117', mstar_msun='1.12', fit_type='RV-only',
      table_ref='arXiv:1608.03597 (Ma+2016, AJ 152, 112) companionparams table (no-trend fit; Mstar stellarparams table)',
      notes='MARVELS-7b; massive planet. Two-Keplerian fit (b,c). HD 87646 is a close visual binary (~22 AU), jitter 240 m/s; Tp in JD_UTC; omega in radians (5.20+/-0.52), converted. Probably in NASA Exoplanet Archive.'),
 dict(host='HD 87646', companion='c', P_d='674', e_P_d='4', K_ms='1370', e='0.50', omega_deg=deg('1.95'),
      Tp_bjd='2453707', msini_mj='57.0', a_au='1.58', mstar_msun='1.12', fit_type='RV-only',
      table_ref='arXiv:1608.03597 (Ma+2016, AJ 152, 112) companionparams table (no-trend fit; Mstar stellarparams table)',
      notes='MARVELS-7c BD candidate. Host is in a close visual binary (~22 AU separation): Gaia astrometry likely complicated. Tp JD_UTC; omega in radians (1.95), converted.'),
]
write(R + 'marvels_bd_orbits.csv', rows)
