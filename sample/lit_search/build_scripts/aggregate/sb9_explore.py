import sys, csv, math, collections
D=sys.argv[1]
main={r['Seq']:r for r in csv.DictReader(open(D+'/tap_main.csv'))}
orb=list(csv.DictReader(open(D+'/tap_orbits.csv')))
def f(x):
    try: return float(x)
    except: return None
MJ=1/1047.57
def m2sini(fm, m1):
    # solve (m2)^3/(m1+m2)^2 = fm for m2 (=M2 sin i assuming sin^3 i ~ sin i scaling; minimum mass)
    lo,hi=0.0,100.0
    for _ in range(200):
        mid=(lo+hi)/2
        if mid**3/(m1+mid)**2 > fm: hi=mid
        else: lo=mid
    return lo
c=0; out=[]
for o in orb:
    if f(o['K2']) is not None: continue
    P=f(o['Per']); K=f(o['K1']); e=f(o['e']) or 0.0
    if P is None or K is None: continue
    fm=1.0361e-7*(1-e*e)**1.5*K**3*P
    m=main[o['Seq']]
    ms=m2sini(fm,1.0)/MJ
    if ms<200 and P<=4000: out.append((round(ms,1),o['Seq'],m['Name'],m['Sp1'],P,K,e,o['Grade'],o['Ref']))
out.sort()
print(len(out))
for r in out: print(r)
