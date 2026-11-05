# Resolve names via CDS Sesame (SIMBAD); cache results in JSON
import sys, json, os, re, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import xml.etree.ElementTree as ET
cache_path=sys.argv[1]; names=[l.strip() for l in open(sys.argv[2]) if l.strip()]
cache=json.load(open(cache_path)) if os.path.exists(cache_path) else {}
def q(n):
    url='https://cds.unistra.fr/cgi-bin/nph-sesame/-oxI/S?'+urllib.parse.quote(n)
    for _ in range(3):
        try:
            t=urllib.request.urlopen(url, timeout=30).read().decode('utf-8','replace')
            break
        except Exception as e:
            t=None
    if t is None: return n, {'error':'network'}
    try: r=ET.fromstring(t)
    except Exception: return n, {'error':'parse'}
    res=r.find('.//Resolver')
    if res is None or res.find('oname') is None: return n, {'resolved':False}
    al=[a.text for a in res.findall('alias') if a.text]
    return n, {'resolved':True, 'main_id':' '.join(res.findtext('oname').split()), 'otype':res.findtext('otype'),
               'aliases':[a for a in al if re.match(r'^(HD|HIP|Gaia DR3|TIC|GJ|2MASS|TYC|BD|CD|KOI|Kepler|TOI|WASP|HAT|NAME|\*|V\*)',a)]}
todo=[n for n in names if n not in cache or 'error' in cache[n]]
with ThreadPoolExecutor(8) as ex:
    for n,r in ex.map(q, todo): cache[n]=r
json.dump(cache, open(cache_path,'w'), indent=0)
print('resolved', sum(1 for n in names if cache[n].get('resolved')), 'of', len(names), 'errors', sum(1 for n in names if 'error' in cache[n]))
for n in names:
    if not cache[n].get('resolved'): print('UNRESOLVED:', n, cache[n])
