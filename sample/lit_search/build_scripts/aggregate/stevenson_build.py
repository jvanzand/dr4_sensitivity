import sys, re
sys.path.insert(0, sys.argv[2])
from common import *
D=sys.argv[1]
rows=[]
src=list(csv.DictReader(open(D+'/BD_catalogue_20-11-23.csv',encoding='utf-8')))
for r in src:
    meth=r['method'].strip()
    if meth in ('Gaia binary_masses','imaging'): continue
    ms=fnum(r['msini']); m=fnum(r['mass']); P=fnum(r['period'])
    use=ms if ms is not None else m
    if use is None or use>80 or P is None or P>4000: continue
    host=clean_name(r['host_star']); host=re.sub(r'^XO-3b$','XO-3',host)
    rvonly = meth=='rv'
    msini_out = ms if ms is not None else (m if meth=='transit' else None)
    note=['method=%s'%meth, 'true mass %s (+%s/-%s) MJ'%(r['mass'],r['p_mass'],r['m_mass']) if m is not None else '',
          'msini_mj filled with transit true mass (sin i ~ 1)' if (ms is None and meth=='transit') else '',
          'no msini given; true mass only' if (ms is None and meth!='transit') else '',
          'source DOI: %s'%r['source'], 'catalog host name: %s'%r['host_star'],
          'no K, omega, Tp or a in this catalogue', 'mstar from catalogue']
    rows.append(dict(host=host, companion=host+' (BD)', P_d=r['period'], e_P_d=r['p_period'], e=r['ecc'], msini_mj=fmt(msini_out,5) if msini_out is not None else '',
        mstar_msun=r['Mstar'], fit_type='RV-only' if meth in ('rv','transit') else 'RV+astrometry',
        table_ref='Stevenson et al. 2023 (MNRAS, accepted 2023; arXiv:2310.02695) BD catalogue, github adam-stevenson/brown-dwarf-desert BD_catalogue_20-11-23.csv',
        notes='; '.join(x for x in note if x)))
write(D+'/stevenson_2023_orbits.csv', rows)
import collections
print(len(src), len(rows), collections.Counter(r['notes'].split(';')[0] for r in rows))
print(sorted({r['host'] for r in rows})[:400])
