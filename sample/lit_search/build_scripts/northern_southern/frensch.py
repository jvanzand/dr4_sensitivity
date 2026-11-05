import sys
sys.path.insert(0, sys.argv[1])
from texhelp import write_csv
L = sys.argv[2]
ref = 'arXiv:2306.04420 Table 6 (BD candidates; mstar Table 2)'
tpn = 'Tp given in rBJD (undefined in paper; values only consistent with BJD-2400000, the DACE convention, so converted with +2400000)'
OUT = 'OUT OF TARGET RANGE (a>5 AU)'
rows = [
 dict(host='HD 62364', companion='HD 62364 b', P_d=5138, e_P_d=23, K_ms=170, e=0.607, omega_deg=1.2, Tp_bjd=2457575,
      msini_mj=12.7, a_au=6.15, mstar_msun=1.157, notes=f'{OUT}; orvara (RV+Hipparcos-Gaia PMa) true mass 18.77 +0.66/-0.63 MJ, i=137.3 deg; {tpn}'),
 dict(host='HD 56380', companion='HD 56380 B', P_d=3254.0, e_P_d=0.5, K_ms=650.5, e=0.6507, omega_deg=23.83, Tp_bjd=2455738.1,
      msini_mj=33.2, a_au=4.17, mstar_msun=0.768, notes=f'orvara true mass 375.3 +8.6/-8.4 MJ (stellar), i=5.96 deg; {tpn}'),
 dict(host='HD 221638', companion='HD 221638 B', P_d=6064, e_P_d=14.5, K_ms=758, e=0.694, omega_deg=200, Tp_bjd=2454013,
      msini_mj=53, a_au=6.87, mstar_msun=1.123, notes=f'{OUT}; orvara true mass 110.0 +3.9/-3.7 MJ (stellar), i=30.68 deg; fit includes drift 4 m/s/yr; {tpn}'),
 dict(host='HD 33473A', companion='HD 33473 C', P_d=18236, e_P_d=905, K_ms=308, e=0.56, omega_deg=108, Tp_bjd=2454859,
      msini_mj=38.3, a_au=14.7, mstar_msun=1.248, notes=f'{OUT}; orvara true mass 271.0 +3.9/-3.8 MJ (stellar), i=170.79 deg; {tpn}'),
]
for r in rows:
    r.update(fit_type='RV-only', table_ref=ref)
write_csv(f'{L}/frensch_2023/frensch_2023_orbits.csv', rows)
