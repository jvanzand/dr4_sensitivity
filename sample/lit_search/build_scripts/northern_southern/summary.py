import csv, glob, os, sys, re, collections
L = sys.argv[1]
lit = {re.sub(r'\s', '', r['hostname']).upper(): r['catalog'] for r in csv.DictReader(open(f'{L}/literature_companions.csv'))}
allin = collections.defaultdict(set)
for f in [f'{L}/{d}/{d}_orbits.csv' for d in 'sahlmann_2011 unger_2023 barbato_2023 rickman_2019 segransan_2010 kiefer_2021 ma_ge_2014 grieves_2017 jenkins_2009 parc_2026 stevenson_2023 frensch_2023'.split()]:
    sid = os.path.basename(os.path.dirname(f))
    rows = list(csv.DictReader(open(f)))
    inr = []
    for r in rows:
        ms = float(r['msini_mj']) if r['msini_mj'] else None
        if 'OUT OF' in r['notes'] or 'BELOW 5' in r['notes']: continue
        if ms is not None and not (5 <= ms <= 80): continue
        if ms is None and 'true mass' not in r['notes'] and 'mass=' not in r['notes']: continue
        inr.append(r)
    rvonly = [r for r in inr if r['fit_type'] == 'RV-only']
    hosts = sorted({r['host'] for r in inr})
    ov = sorted(h for h in hosts if re.sub(r'\s', '', h).upper() in lit)
    for r in inr: allin[r['host']].add(sid)
    print(f"{sid:16s} rows={len(rows):3d} in-range rows={len(inr):3d} (RV-only {len(rvonly):3d}) unique hosts={len(hosts):3d}; overlap w/ existing literature_companions.csv: {ov}")
print('unique in-range hosts across all sources:', len(allin))
only = [h for h, s in allin.items() if s <= {'stevenson_2023', 'ma_ge_2014', 'kiefer_2021'}]
print('hosts whose only source is a compilation without full orbit (stevenson/ma_ge/kiefer):', len(only))
full = [h for h, s in allin.items() if s - {'stevenson_2023', 'ma_ge_2014'}]
print('hosts with full orbital elements (K, omega, Tp) from at least one source:', len(full))
