# Re-implements dace_query parse_parameters (read, not executed from package) to flatten JSON to CSV
import sys, json, csv
from collections import defaultdict
d=json.load(open(sys.argv[1]))
data=defaultdict(list); meta={}
def occ(vals, occs):
    out=[]
    for v,o in zip(vals,occs): out.extend([v]*o)
    return out
for p in d['parameters']:
    name=p['variableName']; meta[name]=(p.get('label'),p.get('unit'),p.get('group'))
    vals=next((p[k] for k in ['doubleValues','floatValues','intValues','stringValues','boolValues'] if p.get(k) is not None), [])
    o=p.get('occurrences')
    data[name].extend(occ(vals,o) if o else vals)
    ev=p.get('minErrorValues')
    if ev is not None: data[name+'_err'].extend(occ(ev,o) if o else ev)
lens={k:len(v) for k,v in data.items()}
n=max(lens.values())
bad=[k for k,v in lens.items() if v!=n]
print('nrows',n,'ncols',len(data),'ragged',bad[:10], file=sys.stderr)
for k,(l,u,g) in meta.items(): print(k,'|',l,'|',u,'|',g, file=sys.stderr)
keys=[k for k in data if lens[k]==n]
w=csv.writer(open(sys.argv[2],'w',newline=''))
w.writerow(keys)
for i in range(n): w.writerow([data[k][i] for k in keys])
