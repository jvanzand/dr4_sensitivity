import sys, math
sys.path.insert(0, sys.argv[1])
from texhelp import write_csv
L = sys.argv[2]
OUT_A = 'OUT OF TARGET RANGE (a>5 AU)'
# Jenkins 2009 (arXiv:0905.2985 Table 5; T0 in JD-2450000; HD191760 mass from Table 1 'Here')
write_csv(f'{L}/jenkins_2009/jenkins_2009_orbits.csv', [
 dict(host='HD 191760', companion='HD 191760 b', P_d=505.65, e_P_d=0.42, K_ms=1047.83, e=0.63, omega_deg=200.37,
      Tp_bjd=2454835.65, msini_mj=38.17, a_au=1.35, mstar_msun=1.28, fit_type='RV-only',
      table_ref='arXiv:0905.2985 Table 5 (mstar Table 1)',
      notes='CHEPS/HARPS; T0 column labelled JD-2,450,000 (assumed time of periastron; paper does not say BJD); Sahlmann+2011 Hipparcos M2=185 +/- 94 MJ (low significance)'),
 dict(host='HD 154672', companion='HD 154672 b', P_d=163.9, e_P_d=0.1, K_ms=226.52, e=0.61, omega_deg=266,
      Tp_bjd=2454520, msini_mj=5.02, a_au=0.60, fit_type='RV-only', table_ref='arXiv:0905.2985 Table 5',
      notes='Magellan + HARPS update of Lopez-Morales+2008 planet; T0 in JD-2,450,000; mstar not given in this paper; likely already in NASA Exoplanet Archive'),
])
# Segransan 2010 (arXiv:0908.1479 Table 2); only companion >=5 MJ
write_csv(f'{L}/segransan_2010/segransan_2010_orbits.csv', [
 dict(host='HD 147018', companion='HD 147018 c', P_d=1008, e_P_d=18, K_ms=141.2, e=0.133, omega_deg=-133.1,
      msini_mj=6.56, a_au=1.922, mstar_msun=0.927, fit_type='RV-only',
      table_ref='arXiv:0908.1479 Table 2 (mstar Table 1)',
      notes='Tp left blank: table gives T0=55301+/-22 in JD*=JD-2450000, inconsistent with that definition (probably 2455301 i.e. a typo); likely already in NASA Exoplanet Archive. Other companions in paper (HD 147018 b, HD 171238 b, HD 204313 b) have Msini 2.1-4.1 MJ'),
])
# Rickman 2019 (arXiv:1904.01573 Tables 3-4; P in years -> days x365.25; Tp in BJD-2450000)
def yr(p): return round(p * 365.25, 2)
rk = 'arXiv:1904.01573 Table 3 (mstar Tables 1-2)'
write_csv(f'{L}/rickman_2019/rickman_2019_orbits.csv', [
 dict(host='HD 13724', companion='HD 13724 b', P_d=yr(40.42), e_P_d=yr((13.42 + 4.38) / 2), K_ms=214.3, e=0.34, omega_deg=187.5,
      Tp_bjd=2456189.1, msini_mj=26.77, a_au=12.40, mstar_msun=1.14, fit_type='RV-only', table_ref=rk,
      notes=OUT_A + '; P given as 40.42 (+13.42/-4.38) yr, converted with 365.25 d/yr; mode of posterior'),
 dict(host='HD 92987', companion='HD 92987 b', P_d=yr(28.35), e_P_d=yr((1.51 + 0.74) / 2), K_ms=152.7, e=0.21, omega_deg=195.1,
      Tp_bjd=2457889, msini_mj=16.88, a_au=9.62, mstar_msun=1.08, fit_type='RV-only', table_ref=rk,
      notes=OUT_A + '; P given as 28.35 yr, converted with 365.25 d/yr'),
 dict(host='HD 181234', companion='HD 181234 b', P_d=yr(20.43), e_P_d=yr(0.215), K_ms=126.8, e=0.73, omega_deg=93.3,
      Tp_bjd=2457668.7, msini_mj=8.37, a_au=7.52, mstar_msun=1.01, fit_type='RV-only', table_ref=rk,
      notes=OUT_A + '; P given as 20.43 yr, converted with 365.25 d/yr'),
 dict(host='HD 98649', companion='HD 98649 b', P_d=yr(16.49), e_P_d=yr((1.13 + 0.70) / 2), K_ms=140.1, e=0.86, omega_deg=252.61,
      Tp_bjd=2455121.7, msini_mj=6.79, a_au=6.57, mstar_msun=1.03, fit_type='RV-only', table_ref='arXiv:1904.01573 Table 4 (mstar Table 1)',
      notes=OUT_A + '; P given as 16.49 yr, converted with 365.25 d/yr'),
])
# Parc+2026 (arXiv:2610.25152 Table C.1); omega in rad -> deg; M0 given at unstated epoch -> Tp blank
pc = []
for host, P, ePu, ePl, M0, K, e, om, ms, a in [
    ('WASP-8', 3054, 13, 5, 1.429, 115.8, 0.147, 2.18, 8.33, 4.161),
    ('WASP-34', 2033.4, 3.6, 6.5, 4.490, 86.3, 0.214, 5.717, 5.29, 3.158),
    ('WASP-53', 3083.7, 5.5, 6.5, 5.528, 477, 0.8365, 3.512, 16.67, 3.910),
    ('WASP-68', 6212, 337, 252, 6.096, 245, 0.804, 1.992, 15.35, 7.11),
    ('WASP-81', 1306, 11, 9, 0.269, 1190, 0.5649, 5.595, 55.6, 2.400),
    ('WASP-99', 5795, 21, 239, 1.84, 189.6, 0.202, 0.15, 21.4, 7.14),
    ('WASP-132', 1741, 17, 23, 5.93, 113.2, 0.171, 2.28, 5.63, 2.619)]:
    n = ([OUT_A] if a > 5 else []) + [f'outer companion to hot Jupiter {host} b; CORALIE (kima); omega given as {om} rad (converted); mean anomaly M0={M0} rad at an epoch not stated in text, so Tp left blank; mstar only in online table']
    if host == 'WASP-34': n.append('paper Table 1 also gives alternative models for WASP-34 c')
    pc.append(dict(host=host, companion=host + ' c', P_d=P, e_P_d=(ePu + ePl) / 2, K_ms=K, e=e,
                   omega_deg=round(math.degrees(om), 2), msini_mj=ms, a_au=a, fit_type='RV-only',
                   table_ref='arXiv:2610.25152 Table C.1', notes='; '.join(n)))
write_csv(f'{L}/parc_2026/parc_2026_orbits.csv', pc)
print('ok')
