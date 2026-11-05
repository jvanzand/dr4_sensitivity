import sys
sys.path.insert(0, sys.argv[3])
from common import *
D=sys.argv[1]
def sl(l,a,b): return l[a-1:b].strip()
refs={}
for l in open(D+'/refs.dat'):
    refs[sl(l,1,3)]=sl(l,5,23)+' '+sl(l,25,48)
rows=[]; n=0
for l in open(D+'/table3.dat'):
    n+=1
    comp=sl(l,1,13); host=sl(l,17,27); P=fnum(sl(l,67,78)); eP=max(fnum(sl(l,80,90)) or 0, fnum(sl(l,92,102)) or 0)
    inc=sl(l,104,110); e=sl(l,127,131); om=sl(l,147,153); K=sl(l,171,178); mc=fnum(sl(l,220,226)); a=fnum(sl(l,242,247)); Tp=sl(l,262,272)
    ref=sl(l,304,314); hip=sl(l,316,321); im=sl(l,15,15); binf=sl(l,29,29)
    if mc is None or mc>80: continue  # filter on true mass (or mc sin i when P<1000 d)
    if not ((P is not None and P<=4000) or (a is not None and a<=5)): continue
    try: Tpb='%.3f'%(float(Tp)+2400000)
    except: Tpb=''
    note=['Feng+2022 Table 3 (VizieR J/ApJS/262/21 table3.dat)',
          ('mass is TRUE mass mc from RV+Hipparcos-Gaia fit (i=%s deg)'%inc) if inc else 'P<1000 d: mass is mc sin i (no inclination fitted)',
          'Tp converted from BJD-2400000', 'HIP %s'%hip if hip else '',
          'directly imaged (flag d)' if im else '', 'wide-binary host (flag *)' if binf else '',
          'refs: '+'; '.join(refs.get(x.strip(),x) for x in ref.split(',') if x.strip()) if ref else '',
          'Use as system list; orbit is a joint RV+astrometry fit']
    host=re.sub(r'^(BD|CD|CPD)([+-]\d{2})\s*(\d+)$', r'\1\2 \3', host)
    import math
    if inc:
        msini=mc*math.sin(math.radians(float(inc))); note.insert(2,'msini_mj computed here as mc*sin(i) = %.3f*sin(%s deg); true mass mc=%s MJ'%(mc,inc,fmt(mc,5)))
    else: msini=mc
    rows.append(dict(host=host, companion=comp, P_d=fmt(P,10), e_P_d=fmt(eP,4), K_ms=K, e=e, omega_deg=om, Tp_bjd=Tpb,
        msini_mj=fmt(msini,5), a_au=fmt(a,5), mstar_msun='', fit_type='RV+astrometry', table_ref='Feng et al. 2022, ApJS 262, 21, Table 3',
        notes='; '.join(x for x in note if x)))
write(D+'/feng_2022_orbits.csv', rows)
print(n, len(rows), '>13', sum(float(r['msini_mj'])>=13 for r in rows), '>30', sum(float(r['msini_mj'])>30 for r in rows))
for r in rows:
    if float(r['msini_mj'])>=13: print(r['companion'], r['host'], r['msini_mj'], r['P_d'], r['a_au'])
