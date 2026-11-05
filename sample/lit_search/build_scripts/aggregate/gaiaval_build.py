import sys, re
sys.path.insert(0, sys.argv[2])
from common import *
D=sys.argv[1]
def sl(l,a,b): return l[a-1:b].strip()
def nm(n):
    n=n.strip()
    n=re.sub(r'^Gl\s*',r'GJ ',n)
    n=re.sub(r'^(HD|HIP|GJ|HR)\s*(\d)',r'\1 \2',n)
    n=re.sub(r'^(BD|CD)\s*([+-]\d{2})\s*(\d+)$',lambda m:'%s%s %d'%(m.group(1),m.group(2),int(m.group(3))),n)
    return n
rows=[]
# Arenou 2023 Table (tab:known_EP_BD_comp) from arXiv:2206.05595 source
tex=open(D+'/arenou2023/src/aa43782_arenou_Gaia.tex').read()
blk=tex[tex.index('label{tab:known_EP_BD_comp}'):tex.index('\\tablebib',tex.index('label{tab:known_EP_BD_comp}'))]
bib=tex[tex.index('\\tablebib',tex.index('label{tab:known_EP_BD_comp}')):]; bib=bib[:bib.index('}\n')+1]
refmap=dict(re.findall(r'\((\d+)\)[\s~]*\\cite[pt]?\{([^}]+)\}',bib))
arenou_names=set()
for l in blk.splitlines():
    if l.strip().startswith('%') or '\\object' not in l: continue
    c=[x.strip() for x in l.split('&')]
    gid=c[0]; name=re.search(r'\\object\{([^}]+)\}',c[1]).group(1).replace('\\,',' ')
    def val(s):
        m=re.match(r'\$?\s*([0-9.]+)',s.replace('$',''))
        return m.group(1) if m else ''
    msini,mc,inc,Pl,Pg,el,eg=[val(x) for x in c[2:9]]
    refs=c[10].replace('\\\\','').strip()
    rl=[]
    for part in refs.split(','):
        part=part.strip()
        if '-' in part:
            a,b=part.split('-'); rl+= [refmap.get(str(i),str(i)) for i in range(int(a),int(b)+1)]
        elif part: rl.append(refmap.get(part,part))
    host=re.sub(r'\s+[b-z]$','',name)
    host=re.sub(r'^(BD)\s*([+-])(\d+)\s+(\d+)$',lambda m:'%s%s%02d %d'%(m.group(1),m.group(2),int(m.group(3)),int(m.group(4))),host)
    arenou_names.add(norm(host))
    P=fnum(Pl)
    rows.append(dict(host=host, companion=name, P_d=Pl, e=el, msini_mj=msini, fit_type='RV-only',
        table_ref='Gaia Collaboration, Arenou et al. 2023, A&A 674, A34, Table "known substellar companions" (arXiv:2206.05595 LaTeX source); literature values',
        notes='Gaia DR3 %s; literature RV refs (bibtex keys): %s; Gaia DR3 astrometric true mass %s MJ (i=%s deg), P_Gaia=%s d, e_Gaia=%s; P/e/Msini are the literature RV values quoted in the table; K/omega/Tp not tabulated'%(gid, ', '.join(rl), mc, inc, Pg, eg)))
# Holl et al. 2023 refsols (VizieR J/A+A/674/A10)
for l in open(D+'/holl2023/refsols.dat'):
    idx=sl(l,144,144); ref=sl(l,172,281); name=nm(sl(l,147,169))
    if idx=='2' or ref.startswith('Gaia DR3'): continue
    if norm(name) in arenou_names: 
        for r in rows:
            if norm(r['host'])==norm(name): r['notes']+='; also in Holl+2023 refsols (ref %s)'%ref
        continue
    P=fnum(sl(l,54,64))
    if P is not None and P>4000: continue
    rows.append(dict(host=name, companion=name+' (Gaia DR3 %s)'%sl(l,1,19), P_d=sl(l,54,64), e_P_d=sl(l,67,76), e=sl(l,104,111),
        fit_type='RV-only', table_ref='Holl et al. 2023, A&A 674, A10, refsols.dat (VizieR J/A+A/674/A10); reference ref: %s'%ref,
        notes='Gaia DR3 %s (%s); Holl pseudo-mass index %s (%s, assuming 1 Msun primary); reference P/e only (no K, omega, Tp, Msini); Gaia DR3 P=%s d, e=%s; %s'%(
            sl(l,1,19), sl(l,22,51), idx, '<20 MJ' if idx=='0' else '20-120 MJ', sl(l,79,89), sl(l,124,131), sl(l,284,403))))
write(D+'/gaia_dr3_val_orbits.csv', rows)
print(len(rows))
for r in rows: print(r['host'],'|',r['P_d'],r['e'],r.get('msini_mj',''),'|',r['table_ref'][:30], r['notes'][:90])
