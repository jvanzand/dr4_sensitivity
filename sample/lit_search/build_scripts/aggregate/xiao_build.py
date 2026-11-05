import sys, re
sys.path.insert(0, sys.argv[2])
from common import *
D=sys.argv[1]
tex=open(D+'/src/RAA-2022-0418.tex').read()
i0=tex.index('label{Tab:orbparas}'); i1=tex.index('end{tabular}',tex.index('label{Tabcc}'))
def v(s):
    s=s.replace('$','').strip()
    m=re.match(r'\{?([-0-9.]+)\}?(?:_\{?([-0-9.]+)\}?\^\{?\+?([-0-9.]+)\}?)?',s)
    if not m or not re.match(r"^-?[0-9]*\.?[0-9]+$", m.group(1)): return None,None
    val=float(m.group(1)); err=None
    if m.group(2): err=max(abs(float(m.group(2))),abs(float(m.group(3))))
    return val,err
rows=[]; n=0
for l in tex[i0:i1].splitlines():
    if not re.match(r'^(HD|GJ|HIP|BD|HR|[a-zA-Z0-9])',l) or '&' not in l or l.startswith('Name') or l.startswith('\\'): continue
    c=[x.strip() for x in l.rstrip('\\').split('&')]
    if len(c)<12: continue
    n+=1
    name=c[0].replace('\\,',' ').strip()
    Mp,_=v(c[1]); a,_=v(c[4]); e,_=v(c[5]); Pyr,ePyr=v(c[6]); om,_=v(c[8]); Tp,_=v(c[10]); msini,_=v(c[11]); i1_,_=v(c[2])
    if msini is None or msini>80: continue
    P=Pyr*365.25 if Pyr else None
    if not ((P is not None and P<=4000) or (a is not None and a<=5)): continue
    host=re.sub(r'\s+[b-z]$','',name); host=re.sub(r'\s+(B|C)$','',host) if re.search(r'\s+[BC]$',host) else host
    rows.append(dict(host=host, companion=name, P_d=fmt(P,6), e_P_d=fmt(ePyr*365.25,3) if ePyr else '', e=fmt(e,4), omega_deg=fmt(om,5),
        Tp_bjd=fmt(Tp+2450000,9) if Tp is not None else '', msini_mj=fmt(msini,4), a_au=fmt(a,4), fit_type='RV+astrometry',
        table_ref='Xiao et al. 2023, RAA (accepted; arXiv:2303.12409), Table "Posteriors of RV Companions" (LaTeX source)',
        notes='orvara joint fit of literature RVs + HGCA (Hipparcos-Gaia EDR3); true mass Mp=%s MJ, i=%s or %s deg (prograde/retrograde ambiguity); P converted from yr; Tp = (Tp-2450000)+2450000, time system as in paper; omega in orvara convention; no K tabulated'%(fmt(Mp,4), fmt(i1_,4), fmt(180-i1_,4) if i1_ is not None else '')))
write(D+'/xiao_2023_orbits.csv', rows)
print(n, len(rows), '>=13', sum(float(r['msini_mj'])>=13 for r in rows), '>30', sum(float(r['msini_mj'])>30 for r in rows))
for r in rows:
    if float(r['msini_mj'])>=13: print(r['companion'], r['msini_mj'], r['P_d'], r['a_au'])
