import sys, csv, re
sys.path.insert(0, sys.argv[3])
from common import *
D=sys.argv[1]; NEA=sys.argv[2]
main={r['Seq']:r for r in csv.DictReader(open(D+'/tap_main.csv'))}
orb=list(csv.DictReader(open(D+'/tap_orbits.csv')))
alias={}
for a in csv.DictReader(open(D+'/tap_alias.csv')): alias.setdefault(a['Seq'],[]).append(a['Name'])
flame={r['Source']:r for r in csv.DictReader(open(D+'/gaia_dr3_flame_for_candidates.csv'))}
gast={r['Source']:r for r in csv.DictReader(open(D+'/gaia_dr3_astrometry_for_candidates.csv'))}
import math
# Mamajek dwarf table: SpT -> Msun
sptm={}; absmag={}
for i,l in enumerate(open(D+'/EEM_dwarf_UBVIJHK_colors_Teff.txt')):
    if 23<=i<140 and l[:1] in 'OBAFGKMLTY' and not l.startswith('#'):
        t=l.split(); 
        if fnum(t[-2]) is not None: sptm[t[0]]=float(t[-2]); absmag[t[0]]=(fnum(t[13]), fnum(t[7]))
def spt_mass(sp):
    s=(sp or '').strip()
    m=re.match(r'^([OBAFGKM])\s*(\d(?:\.\d)?)?',s)
    if not m: return None,'no parsable SpT'
    rest=s[m.end():]
    cls='V'
    if re.search(r'(III|II|Ia|Ib|I(?![V]))',rest.replace('IV','')): cls='giant'
    elif re.search(r'IV',rest): cls='IV'
    elif not re.search(r'V',rest): cls='none'
    if cls=='giant': return None,'giant/supergiant SpT (%s); no dwarf-table mass'%s
    sub=m.group(2) or '5'
    key=m.group(1)+sub+'V'
    if key not in sptm:
        key=m.group(1)+str(int(float(sub)))+'V'
    if key not in sptm: return None,'SpT %s not in table'%s
    note='SpT %s -> %s in Mamajek dwarf table'%(s,key)
    if not m.group(2): note+=' (no subtype; assumed 5)'
    if cls=='IV': note+=' (class IV; dwarf mass used, likely underestimate)'
    if cls=='none': note+=' (no luminosity class; dwarf assumed)'
    return sptm[key],note+'|'+key
hosts,_=nea_sets(NEA)
rows=[]; skipped=[]
for o in orb:
    if fnum(o['K2']) is not None: continue
    P=fnum(o['Per']);K=fnum(o['K1']);e=fnum(o['e']) or 0.0
    if P is None or K is None or P>4000: continue
    fm=1.0361e-7*(1-e*e)**1.5*K**3*P
    if m2_from_fm(fm,1.0)/MJ_MSUN>=200: continue
    seq=o['Seq']; mrow=main[seq]; al=alias.get(seq,[])
    gids=[a.split()[1] for a in al if a.startswith('GAIADR3 ')]
    m1=None; src=''
    for g in gids:
        fr=flame.get(g)
        if fr and fnum(fr['Mass-Flame']):
            m1=fnum(fr['Mass-Flame']); src='M1 = Gaia DR3 FLAME mass (Gaia DR3 %s; assumes single star)'%g; break
    if m1 is None:
        lg=[fnum(flame[g]['logg']) for g in gids if g in flame and fnum(flame[g]['logg']) is not None]
        if lg and lg[0]<3.8:
            src='Gaia DR3 GSP-Phot logg=%.2f < 3.8 (evolved) and no FLAME mass'%lg[0]
        else:
            m1,src=spt_mass(mrow['Sp1'])
            if m1 is not None:
                src,key=src.split('|')
                MG,MV=absmag[key]; chk=None
                ga=[gast[g] for g in gids if g in gast and fnum(gast[g]['Plx']) and fnum(gast[g]['Plx'])>0]
                if ga:
                    plx=fnum(ga[0]['Plx']); G=fnum(ga[0]['Gmag'])
                    if MG is not None and G is not None: chk=('M_G',G+5*math.log10(plx/100.),MG)
                    elif MV is not None and mrow['n_mag1'].strip()=='V' and fnum(mrow['mag1']) is not None: chk=('M_V',fnum(mrow['mag1'])+5*math.log10(plx/100.),MV)
                if chk and chk[1] < chk[2]-1.5:
                    src='SpT %s but Gaia DR3 parallax gives %s=%.2f vs dwarf %.2f (overluminous; evolved or misclassified)'%(mrow['Sp1'].strip(),chk[0],chk[1],chk[2]); m1=None
                elif chk: src+='; dwarf status consistent with Gaia DR3 parallax (%s=%.2f vs table %.2f)'%chk
                elif 'no luminosity class' in src:
                    src='SpT %s has no luminosity class and no Gaia DR3 parallax check is possible; excluded'%mrow['Sp1'].strip(); m1=None
                else: src+='; dwarf status NOT verified (no Gaia DR3 parallax/abs-mag check)'
            if m1 is not None and lg: src+='; Gaia DR3 GSP-Phot logg=%.2f'%lg[0]
            elif m1 is not None: src+='; no Gaia DR3 logg to confirm dwarf status'
    if m1 is None:
        skipped.append((seq,o['o'],mrow['Name'].strip(),mrow['Sp1'].strip(),P,K,e,round(m2_from_fm(fm,1.0)/MJ_MSUN,1),o['Ref'].strip(),src)); continue
    ms=m2_from_fm(fm,m1)/MJ_MSUN
    if ms>100: continue
    host=mrow['Name'].strip()
    host=re.sub(r'^(BD|CD|CPD)\s+([+-])',r'\1\2',host)
    innea=any(norm(a) in hosts for a in al+[host])
    flags=[]
    if o['f_Per'].strip()=='a': flags.append('P fixed')
    if o['f_e'].strip(): flags.append('e flag %s'%o['f_e'].strip())
    if o['u_K1'].strip(): flags.append('K1 uncertain (:)')
    if o['f_K1'].strip(): flags.append('K1 flag %s'%o['f_K1'].strip())
    if fnum(o['f_T0']): flags.append('T0 is not periastron (flag %s)'%o['f_T0'])
    note=('SB9 Seq %s orbit %s; SB1. f(m)=1.0361e-7(1-e^2)^1.5 K1^3 P = %.3g Msun; M2sini solved from f(m)=(M2 sini)^3/(M1+M2 sini)^2 (sin^3 i -> M2 sin i approximation, minimum mass); %s; ROUGH. SB9 grade=%s; Sp1=%s; host in NASA Exoplanet Archive: %s; T0 is JD as listed in SB9 (time system as published, not necessarily BJD); a from Kepler III with M1+M2sini'
          %(seq,o['o'],fm,src,o['Grade'] or 'none',mrow['Sp1'].strip() or 'none','yes' if innea else 'no'))
    if flags: note+='; flags: '+', '.join(flags)
    rows.append(dict(host=host, companion='%s (SB9 %s.%s)'%(host,seq,o['o']), P_d=fmt(P,10), e_P_d=o['e_Per'],
        K_ms=fmt(K*1000,6), e=o['e'], omega_deg=o['omega'], Tp_bjd=o['T0'], msini_mj=fmt(ms,4),
        a_au=fmt(a_au(P,m1+ms*MJ_MSUN),4), mstar_msun=fmt(m1,4), fit_type='RV-only',
        table_ref='SB9 (Pourbaix+2004; VizieR B/sb9, version 2024-04-22); orbit ref %s'%o['Ref'].strip(), notes=note))
best={}
for r in rows:
    seq=r['companion'].split('SB9 ')[1].split('.')[0]
    yr=int(re.search(r'orbit ref (\d{4})',r['table_ref']).group(1))
    gr=float(re.search(r'SB9 grade=([\d.]+|none)',r['notes']).group(1).replace('none','-1'))
    best.setdefault(seq,[]).append(((yr,gr),r))
rows=[]
for seq,lst in best.items():
    lst.sort(key=lambda t:t[0], reverse=True)
    r=lst[0][1]
    if len(lst)>1: r['notes']+='; %d SB9 orbits for this system in range, kept most recent (others: %s)'%(len(lst), ', '.join(re.search(r'orbit ref (\S+)',x[1]['table_ref']).group(1)+' msini=%s'%x[1]['msini_mj'] for x in lst[1:]))
    rows.append(r)
rows.sort(key=lambda r: float(r['msini_mj']))
write(D+'/sb9_orbits.csv', rows)
with open(D+'/sb9_excluded_no_primary_mass.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['sb9_seq','orbit','name','Sp1','P_d','K1_kms','e','msini_mj_if_M1_1msun','orbit_ref','reason_excluded'])
    for r in skipped: w.writerow(r)
print('rows',len(rows),'<=80:',sum(float(r['msini_mj'])<=80 for r in rows),'13-80:',sum(13<=float(r['msini_mj'])<=80 for r in rows), 'unique hosts', len({r['host'] for r in rows}))
print('M1 from FLAME:',sum('FLAME' in r['notes'] for r in rows))
print('in NEA:',sum('Archive: yes' in r['notes'] for r in rows))
print('skipped (no M1):',len(skipped))
for s in skipped[:15]: print('  ',s)
for r in rows: print(r['msini_mj'], r['host'], r['P_d'], r['K_ms'], r['mstar_msun'], 'NEA' if 'Archive: yes' in r['notes'] else '', r['table_ref'][-19:])
