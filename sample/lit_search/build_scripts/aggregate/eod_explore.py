import sys, csv, collections
rows=list(csv.DictReader(open(sys.argv[1], newline='', encoding='latin-1')))
print(len(rows)); print([k for k in rows[0].keys()][:400])
def f(x):
    try: return float(x)
    except: return None
print(collections.Counter(r['PLANETDISCMETH'] for r in rows))
print('DATE range', min(r['DATE'] for r in rows if r['DATE']), max(r['DATE'] for r in rows if r['DATE']))
big=[r for r in rows if (f(r['MSINI']) or 0)>13]
print('msini>13', len(big), 'EOD', sum(1 for r in big if r['EOD']=='1'))
for r in sorted(big,key=lambda r:f(r['MSINI'])): print(r['NAME'], r['MSINI'], r['MASS'], r['PER'], r['A'], r['K'], r['PLANETDISCMETH'], r['EOD'], r['ORBREF'][:40], r['DATE'])
