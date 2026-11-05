import sys, csv
sys.path.insert(0, sys.argv[1])
from texhelp import *
vz = sys.argv[2]; out = sys.argv[3]; sahl = sys.argv[4]; unger = sys.argv[5]
EPOCH = 2455500.0
def f(x): return None if x == '---' else float(x)
mstar = {}
for line in open(vz + '/table1.dat'):
    mstar[line[0:9].strip()] = float(line[51:55])
recs = []
for line in open(vz + '/table2.dat'):
    t = line.split()
    assert len(t) == 50, (len(t), t[0])
    recs.append(t)
# --- verify reference epoch of lambda0 against independent Tp values ---
ref = {}
for p in (sahl, unger):
    for r in csv.DictReader(open(p)):
        if r['fit_type'] == 'RV-only' and r['Tp_bjd']:
            ref.setdefault(r['host'].replace(' ', ''), []).append((float(r['P_d']), float(r['Tp_bjd']), r['table_ref'][:15]))
print('lambda0 epoch check (phase offset of Tp vs reference, in units of P):')
for t in recs:
    if t[0] in ref:
        P, l0, om = float(t[1]), float(t[10]), float(t[13])
        tp = EPOCH + (om - l0) / 360.0 * P
        for (Pr, Tr, src) in ref[t[0]]:
            d = ((tp - Tr) / P + 0.5) % 1 - 0.5
            print(f'  {t[0]:10s} {src} dphase={d:+.4f}')
rows = []; seen = {}
for t in recs:
    star = t[0]
    P, Pu, Pl = f(t[1]), f(t[2]), f(t[3]); K = f(t[4]) * 1000; e = f(t[7]); l0 = f(t[10]); om = f(t[13])
    ms = f(t[16]); a = f(t[22]); tspan = f(t[25])
    if ms > 80: continue
    n = seen.get(star, 0); seen[star] = n + 1
    host = norm_name(star)
    tp = EPOCH + (om - l0) / 360.0 * P
    notes = []
    if a > 5: notes.append('OUT OF TARGET RANGE (a>5 AU)')
    if ms < 5: notes.append('BELOW 5 MJ (planet in hierarchical stellar system)')
    if sum(1 for r in recs if r[0] == star) > 1: notes.append('host has >1 companion in Table 2 (hierarchical system)')
    notes.append(f'Tp_bjd computed as 2455500 + (omega - lambda0)/360*P from lambda0={l0:g} deg (epoch 2455500 not stated in paper; inferred and verified against Sahlmann+2011/Unger+2023 Tp, see notes file)')
    notes.append(f'tspan/P={tspan:g}')
    Mt = f(t[27])
    if Mt is not None:
        notes.append(f'RV+PMa (orvara) true mass M={Mt:g} +{f(t[28]):g}/-{f(t[29]):g} MJ, i={f(t[42]):g} deg; RV+PMa P={f(t[36]):g} d, e={f(t[39]):g}, a={f(t[33]):g} AU')
    if t[49] != '---': notes.append('previous publication ' + t[49])
    rows.append(dict(host=host, companion=host + ' ' + 'bcd'[n], P_d=P, e_P_d=(Pu + Pl) / 2, K_ms=round(K, 3), e=e, omega_deg=om,
                     Tp_bjd=round(tp, 3), msini_mj=ms, a_au=a, mstar_msun=mstar.get(star), fit_type='RV-only',
                     table_ref='VizieR J/A+A/674/A114 table2.dat (RV-only columns); mstar table1.dat', notes='; '.join(notes)))
write_csv(out, rows)
print(len(recs), 'records;', len(rows), 'with Msini<=80 MJ;', sum(1 for r in rows if r['a_au'] <= 5), 'with a<=5 AU')
