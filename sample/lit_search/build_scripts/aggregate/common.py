import csv, re
COLS=['host','companion','P_d','e_P_d','K_ms','e','omega_deg','Tp_bjd','msini_mj','a_au','mstar_msun','fit_type','table_ref','notes']
MJ_MSUN=1/1047.57  # IAU nominal M_J/M_sun ratio (1.8982e27 kg / 1.98847e30 kg ~ 1/1047.5)
def norm(n):
    return re.sub(r'[\s_\-]','',n or '').upper()
def nea_sets(path):
    hosts=set(); pls={}
    for r in csv.DictReader(open(path)):
        for k in ('hostname','hd_name','hip_name'):
            if r[k]: hosts.add(norm(r[k]))
        pls[norm(r['pl_name'])]=r
    return hosts, pls
def m2_from_fm(fm, m1):
    lo,hi=0.,200.
    for _ in range(200):
        mid=(lo+hi)/2
        if mid**3/(m1+mid)**2>fm: hi=mid
        else: lo=mid
    return lo
def a_au(P_d, mtot):
    return (mtot*(P_d/365.25)**2)**(1/3)
def fnum(x):
    try:
        v=float(x)
        return v if v==v else None
    except: return None
def fmt(x, nd=6):
    if x is None: return ''
    if isinstance(x,float): return ('%.*g'%(nd,x))
    return str(x)
def write(path, rows):
    with open(path,'w',newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=COLS); w.writeheader()
        for r in rows: w.writerow({k:r.get(k,'') for k in COLS})
def nea_pl_keys(path):
    """Set of normalized (host-alias + letter) keys for NEA planets."""
    keys=set()
    for r in csv.DictReader(open(path)):
        letter=r['pl_name'].strip().split()[-1] if r['pl_name'].strip() else ''
        keys.add(norm(r['pl_name']))
        for k in ('hostname','hd_name','hip_name'):
            if r[k]:
                keys.add(norm(re.sub(r'\s+A$','',r[k].strip())+letter))
    return keys
def pl_keys(pnames, hostnames):
    keys=set()
    for n in pnames:
        n=n.strip(); keys.add(norm(n)); keys.add(norm(re.sub(r'\s+A\s*([b-z])$',r' \1',n)))
        m=re.search(r'\s+A?([b-z])$',n)
        if m:
            for h in hostnames: keys.add(norm(re.sub(r'\s+A$','',h.strip())+m.group(1)))
    return keys
GREEK={'alpha':'alf','beta':'bet','gamma':'gam','delta':'del','epsilon':'eps','zeta':'zet','eta':'eta','theta':'tet','iota':'iot','kappa':'kap','lambda':'lam','mu':'mu.','nu':'nu.','xi':'ksi','omicron':'omi','pi':'pi.','rho':'rho','sigma':'sig','tau':'tau','upsilon':'ups','phi':'phi','chi':'chi','psi':'psi','omega':'ome'}
def clean_name(n):
    n=(n or '').replace('−','-').replace('â\u0088\u0092','-').strip()
    n=re.sub(r'\$\\?([a-zA-Z]+)\$', lambda m: GREEK.get(m.group(1).lower(), m.group(1)), n)
    n=re.sub(r'^(HD|HIP|GJ|HR|LHS|Gl|NGC|HAT-P-|WASP-)(?=\d)', lambda m: m.group(1) if m.group(1).endswith('-') else m.group(1)+' ', n)
    n=re.sub(r'^Gl ', 'GJ ', n)
    n=re.sub(r'^TYC(?=\d)', 'TYC ', n)
    n=re.sub(r'^(BD|CD|CPD)\s*([+-])\s*(\d{2})\s*0*(\d+)$', lambda m:'%s%s%s %s'%(m.group(1),m.group(2),m.group(3),m.group(4)), n)
    n=re.sub(r'^(HD|HIP|GJ) (\d+)([A-C])$', r'\1 \2 \3', n)
    n=re.sub(r'^(\d+)([A-Z][a-z]{2})$', r'\1 \2', n)
    n=re.sub(r'^NGC (\d+)-(\d+)$', r'NGC \1 \2', n)
    return n
