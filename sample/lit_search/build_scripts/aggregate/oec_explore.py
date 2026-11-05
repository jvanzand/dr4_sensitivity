import sys, gzip, xml.etree.ElementTree as ET, collections
fn = sys.argv[1]
root = ET.parse(gzip.open(fn)).getroot()
n=0; meth=collections.Counter(); upd=collections.Counter(); rows=[]
for sysel in root.iter('system'):
    for star in sysel.iter('star'):
        for pl in star.findall('planet'):
            n+=1
            m = pl.findtext('mass'); P = pl.findtext('period'); a=pl.findtext('semimajoraxis')
            dm = pl.findtext('discoverymethod'); meth[dm]+=1
            lu = pl.findtext('lastupdate') or ''; upd[lu[:2]]+=1
            try: m=float(m)
            except: m=None
            rows.append((pl.findtext('name'), star.findtext('name'), m, P, a, dm, lu, pl.findtext('list')))
print(n, meth.most_common(), sorted(upd.items())[-6:])
big=[r for r in rows if r[2] and r[2]>13]
print('mass>13:', len(big), 'RV:', sum(1 for r in big if r[5]=='RV'), '>30:', sum(1 for r in big if r[2]>30), '>30 RV:', sum(1 for r in big if r[2]>30 and r[5]=='RV'))
print('max mass', max(r[2] for r in rows if r[2]))
for r in sorted([r for r in big if r[2]>30 and r[5]=='RV'], key=lambda r:r[2])[:200]: print(r)
