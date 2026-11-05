import sys
sys.path.insert(0, sys.argv[1])
from texhelp import write_csv
L = sys.argv[2]
MJ_PER_MSUN = 1047.57
rows = [
 dict(host='HD 39213', companion='HD 39213 b', P_d=1309, e_P_d=159, K_ms=1265, e=0.2, omega_deg=265, Tp_bjd=2451636.7,
      msini_mj=round(0.07 * MJ_PER_MSUN, 1), a_au=2.3, mstar_msun=0.93, fit_type='RV-only',
      table_ref='arXiv:1507.04749 Table 2(i) (Msini from Sect. 3 text; mstar Table 1)',
      notes='AAT/UCLES (AAPS) + HARPS; Msini given only in text as 0.07+/-0.01 Msun, converted at 1047.57 MJ/Msun; authors note poor fit (chi2_nu=1.97, rms 17 m/s); Tp column labelled JD-2450000; could be BD only if edge-on'),
 dict(host='HD 145825', companion='HD 145825 b', P_d=6024, e_P_d=163, K_ms=669, e=0.33, omega_deg=149, Tp_bjd=2451095,
      msini_mj=round(0.06 * MJ_PER_MSUN, 1), a_au=6.54, mstar_msun=1.03, fit_type='RV-only',
      table_ref='arXiv:1507.04749 Table 2(ii) (Msini and mstar from Sect. 3 text)',
      notes='OUT OF TARGET RANGE (a>5 AU); minimum-period solution of a partially covered orbit; Msini 0.06+/-0.01 Msun in text, converted at 1047.57 MJ/Msun; Tp column labelled JD-2450000'),
]
write_csv(f'{L}/jenkins_2015/jenkins_2015_orbits.csv', rows)
