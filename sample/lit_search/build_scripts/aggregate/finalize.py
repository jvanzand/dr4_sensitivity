import sys, csv, json, re, gzip, subprocess, os
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(__file__))
from common import *
R, NEA, CACHE, mode = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
FILES=['sb9/sb9_orbits.csv','exoplanet_eu/exoplanet_eu_orbits.csv','exoplanets_org/exoplanets_org_orbits.csv','feng_2022/feng_2022_orbits.csv',
       'gaia_dr3_val/gaia_dr3_val_orbits.csv','stevenson_2023/stevenson_2023_orbits.csv','xiao_2023/xiao_2023_orbits.csv']
# exoplanet.eu alternative names
alt={}
root=ET.parse(gzip.open(R+'/oec/systems-exoplaneteu.xml.gz')).getroot()
for s in root.iter('system'):
    sn=[n.text for n in s.findall('name') if n.text]
    for st in s.iter('star'):
        stn=[n.text for n in st.findall('name') if n.text]
        if stn: alt.setdefault(stn[0],[]).extend(stn+sn)
def cands(h):
    c=[h, clean_name(h)]
    for x in list(c):
        c.append(re.sub(r'\s*\(?AB\)?$','',x)); c.append(re.sub(r'\s*A$','',x)); c.append(re.sub(r'\.0\d$','',x))
        c.append(re.sub(r'^TYC (\d+)-0*(\d+)-(\d)$',lambda m:'TYC %d-%d-%s'%(int(m.group(1)),int(m.group(2)),m.group(3)),x))
    for a in alt.get(h,[]): c+= [a, clean_name(a), re.sub(r'\s*A$','',a)]
    out=[]
    for x in c:
        x=x.strip()
        if x and x not in out: out.append(x)
    return out
cache=json.load(open(CACHE)) if os.path.exists(CACHE) else {}
if len(sys.argv)>5: FILES=sys.argv[5].split(',')
tables={f:list(csv.DictReader(open(R+'/'+f))) for f in FILES}
if mode=='collect':
    need=set()
    for f,rows in tables.items():
        for r in rows:
            cs=cands(r['host'])
            if any(cache.get(c,{}).get('resolved') for c in cs): continue
            need.update(c for c in cs if c not in cache)
    open(os.path.dirname(CACHE)+'/need.txt','w').write('\n'.join(sorted(need)))
    print('need', len(need)); sys.exit()
# NEA index by HD/HIP/Gaia
nea_by={}
for r in csv.DictReader(open(NEA)):
    for k in ('hd_name','hip_name','gaia_dr3_id','hostname'):
        if r[k]: nea_by.setdefault(r[k].strip(),[]).append(r)
lit_by={}
for r in csv.DictReader(open(R+'/literature_companions.csv')):
    for c in cands(r['hostname']):
        info=cache.get(c,{})
        if info.get('resolved'):
            for a in info['aliases']+[info['main_id']]: lit_by.setdefault(a,[]).append(r)
            break
    if r.get('source_id'): lit_by.setdefault('Gaia DR3 '+r['source_id'],[]).append(r)
OLD=[r'planet in NASA Exoplanet Archive( \(name match\))?: (yes|no)', r'host in NASA Exoplanet Archive: (yes|no)']
summary={}
for f,rows in tables.items():
    out=[]; st={'rows':0,'unres':0,'nea_pl':0,'nea_host':0,'lit':0,'dropped_nea_lt13':0}
    for r in rows:
        cs=cands(r['host']); hit=None
        for c in cs:
            if cache.get(c,{}).get('resolved'): hit=c; break
        notes=r['notes']
        for p in OLD: notes=re.sub(r';?\s*'+p,'',notes)
        if hit is None:
            st['unres']+=1; notes+='; host NOT resolved by SIMBAD/Sesame (name as given by source)'
            ids=[]
        else:
            info=cache[hit]
            if hit!=r['host']: notes+='; host renamed from "%s" for SIMBAD'%r['host']; r['host']=hit
            ids=info['aliases']+[info['main_id']]
            notes+='; SIMBAD main_id=%s (otype %s)'%(info['main_id'], info.get('otype'))
        nea=[p for i in ids for p in nea_by.get(i,[])]
        P=fnum(r['P_d'])
        plmatch=[p for p in nea if P and fnum(p['pl_orbper']) and abs(fnum(p['pl_orbper'])-P)/fnum(p['pl_orbper'])<0.1]
        lit=[x for i in ids for x in lit_by.get(i,[])]
        litm=[x for x in lit if P and fnum(x['pl_orbper']) and abs(fnum(x['pl_orbper'])-P)/fnum(x['pl_orbper'])<0.1]
        notes+='; NASA Exoplanet Archive overlap (SIMBAD HD/HIP/Gaia DR3 ID match, P within 10%%): %s'%(
            'companion YES (%s, mass %s MJ)'%(plmatch[0]['pl_name'], plmatch[0]['pl_bmassj']) if plmatch else ('host yes, companion no' if nea else 'no'))
        if litm: notes+='; also in project literature_companions.csv (%s)'%litm[0]['catalog']
        if mode=='write' and f.split('/')[0] in ('exoplanet_eu','exoplanets_org') and plmatch and (fnum(r['msini_mj']) or 0)<13:
            st['dropped_nea_lt13']+=1; continue
        st['rows']+=1; st['nea_pl']+=bool(plmatch); st['nea_host']+=bool(nea); st['lit']+=bool(litm)
        r['notes']=notes.lstrip('; '); out.append(r)
    summary[f]=st
    if mode=='write': write(R+'/'+f, out)
for f,s in summary.items(): print(f, s)
