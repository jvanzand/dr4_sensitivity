import sys, csv, re
sys.path.insert(0, sys.argv[1])
from texhelp import write_csv, norm_name
inp, out = sys.argv[2], sys.argv[3]
def fl(x):
    try: return float(x)
    except: return None
def hname(n):
    n = n.strip().replace('\u2212', ' ')
    fix = {'LHS6343': 'LHS 6343', 'XO-3b': 'XO-3', 'NGC2423 3': 'NGC 2423 3', 'NGC4349 127': 'NGC 4349 127'}
    if n in fix: return fix[n]
    m = re.match(r'^COROT-(\d+)$', n, re.I)
    if m: return 'CoRoT-' + m.group(1)
    m = re.match(r'^TYC\s*(\d+)-(\d+)-(\d+)$', n)
    if m: return 'TYC ' + '-'.join(str(int(g)) for g in m.groups())
    m = re.match(r'^(BD|CD)([+-])(\d{2})(\d+)$', n.replace(' ', ''))
    if m: return f"{m.group(1)}{m.group(2)}{m.group(3)} {int(m.group(4))}"
    m = re.match(r'^(\d+)\s*(Com|UMi|Ari)\b(.*)$', n)
    if m: return f"{m.group(1)} {m.group(2)}{m.group(3)}"
    m = re.match(r'^(HD|HIP|HR|GJ|TYC|TOI|KOI|WASP|CoRoT|KELT|NLTT|LHS)[\s-]?(.+)$', n)
    if m and m.group(1) in ('HD','HIP','HR','GJ','TYC'):
        return f"{m.group(1)} {m.group(2).strip()}"
    return n
rows = []; seen = {}
for r in csv.DictReader(open(inp)):
    meth = r['method']
    if meth in ('Gaia binary_masses', 'imaging'): continue
    P = fl(r['period']); ms = fl(r['msini']); M = fl(r['mass'])
    host = hname(r['host_star'])
    k = seen.get(host, 0); seen[host] = k + 1
    n = [f"method={meth}"]
    if M is not None: n.append(f"mass={M:g} (+{r['p_mass']}/-{r['m_mass']}) MJ")
    if ms is None: n.append('no Msini listed (transiting or mass-only entry)')
    if (ms or M or 0) > 80: n.append('OUT OF TARGET RANGE (>80 MJ)')
    if P is not None and P > 4000: n.append('LIKELY OUT OF TARGET RANGE (P>4000 d; a not listed)')
    if meth not in ('rv', 'transit'): n.append('P/e may be from the joint RV+astrometry fit in the source; check source')
    if host != r['host_star']: n.append('catalogue name ' + r['host_star'])
    n.append('source doi ' + r['source'])
    n.append('compilation: no K/omega/Tp; take full orbit from source paper')
    ep = None
    if fl(r['p_period']) is not None and fl(r['m_period']) is not None:
        ep = (fl(r['p_period']) + fl(r['m_period'])) / 2
    rows.append(dict(host=host, companion=host + ' ' + 'bcd'[k], P_d=P, e_P_d=ep, e=fl(r['ecc']), msini_mj=ms,
                     mstar_msun=fl(r['Mstar']), fit_type='RV-only' if meth in ('rv', 'transit') else 'RV+astrometry',
                     table_ref='Stevenson+2023 (arXiv:2310.02695) catalogue github BD_catalogue_20-11-23.csv',
                     notes='; '.join(n)))
write_csv(out, rows)
inr = [x for x in rows if (x['msini_mj'] or 0) <= 80 and (x['P_d'] or 0) <= 4000]
print(len(rows), 'RV/transit rows;', len(inr), 'with Msini<=80 (or none) and P<=4000 d')
