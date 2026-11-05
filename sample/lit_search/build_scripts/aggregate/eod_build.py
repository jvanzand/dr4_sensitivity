import sys, csv, collections
sys.path.insert(0, sys.argv[4])
from common import *
rows_in=list(csv.DictReader(open(sys.argv[1], newline='', encoding='latin-1')))
nk=nea_pl_keys(sys.argv[3])
out=[]; st=collections.Counter()
for r in rows_in:
    if r['PLANETDISCMETH']!='RV' and not (fnum(r['MSINI']) or 0)>=13: continue
    ms=fnum(r['MSINI']); P=fnum(r['PER'])
    if ms is None or P is None or ms>80 or P>4000: continue
    names=[r['NAME']]+([r['OTHERNAME']] if r['OTHERNAME'] else [])
    hostn=[r['STAR']]+[x for x in (r['HD'],r['HIPP'],r['SIMBADNAME'],r['GL']) if x]
    innea=bool(pl_keys(names, hostn) & nk)
    if not (ms>=13 or not innea): continue
    st[(r['PLANETDISCMETH'], ms>=13, innea)]+=1
    host=r['SIMBADNAME'] or r['STAR']
    out.append(dict(host=host, companion=r['NAME'], P_d=r['PER'], e_P_d=r['UPER'], K_ms=r['K'], e=r['ECC'], omega_deg=r['OM'],
        Tp_bjd=r['T0'] if r['PLANETDISCMETH']=='RV' else '', msini_mj=r['MSINI'], a_au=r['A'], mstar_msun=r['MSTAR'], fit_type='RV-only',
        table_ref='Exoplanet Orbit Database (exoplanets.org/csv-files/exoplanets.csv; Han+2014), orbit ref: %s'%r['ORBREF'],
        notes='; '.join(x for x in ['EOD discovery method=%s'%r['PLANETDISCMETH'], 'planet in NASA Exoplanet Archive (name match): %s'%('yes' if innea else 'no'),
              'EOD flag=%s'%r['EOD'], 'T0 as listed in EOD (epoch of periastron for RV orbits; time system as published)' if r['PLANETDISCMETH']=='RV' else 'T0 not periastron for transiting planets; omitted',
              'EOD snapshot frozen ~2018'] if x)))
write(sys.argv[2], out)
print(len(out)); [print(k,v) for k,v in st.items()]
for o in out: print(o['companion'], o['msini_mj'], o['P_d'], o['notes'][:80])
