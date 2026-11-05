import sys, gzip, re, collections
import xml.etree.ElementTree as ET
sys.path.insert(0, sys.argv[4])
from common import *
src, out, NEA = sys.argv[1], sys.argv[2], sys.argv[3]
hosts, pls = nea_sets(NEA); nk = nea_pl_keys(NEA)
root = ET.parse(gzip.open(src)).getroot()
rows=[]; stats=collections.Counter()
def err(el):
    if el is None: return ''
    a=fnum(el.get('errorminus')); b=fnum(el.get('errorplus'))
    if a is None or b is None: return ''
    return fmt(max(a,b),4)
for sysel in root.iter('system'):
    sysnames=[n.text for n in sysel.findall('name') if n.text]
    for star in sysel.iter('star'):
        snames=[n.text for n in star.findall('name') if n.text]
        for pl in star.findall('planet'):
            dm=pl.findtext('discoverymethod')
            m=fnum(pl.findtext('mass')); P=fnum(pl.findtext('period')); a=fnum(pl.findtext('semimajoraxis'))
            if m is None or m>80: continue
            if not ((P is not None and P<=4000) or (P is None and a is not None and a<=5)): continue
            pnames=[n.text for n in pl.findall('name') if n.text]
            in_nea_pl=bool(pl_keys(pnames, snames+sysnames) & nk)
            in_nea_host=any(norm(n) in hosts for n in snames+sysnames)
            if dm=='RV':
                keep = m>=13 or not in_nea_pl
            elif dm=='transit':
                keep = m>=13 and not in_nea_pl
            else: keep=False
            if not keep: continue
            stats[(dm, m>=13, m>30, in_nea_pl)]+=1
            inc=pl.findtext('inclination')
            mnote=('mass column is exoplanet.eu "mass" (true mass when exoplanet.eu lists one, else M sin i)'
                   + ('; inclination %s deg listed, so this may be a TRUE mass'%inc if inc else '; no inclination listed, likely M sin i'))
            notes=['exoplanet.eu discovery method=%s'%dm, mnote,
                   'planet in NASA Exoplanet Archive: %s'%('yes' if in_nea_pl else 'no'),
                   'host in NASA Exoplanet Archive: %s'%('yes' if in_nea_host else 'no'),
                   'exoplanet.eu last update 20%s (yy/mm/dd)'%(pl.findtext('lastupdate') or '?'),
                   'alt names: %s'%'; '.join(pnames[1:5]) if len(pnames)>1 else '',
                   'no K or Tp in this mirror (fetch from exoplanet.eu or the discovery paper)',
                   'fit_type not given by exoplanet.eu; RV-only assumed']
            rows.append(dict(host=snames[0] if snames else sysnames[0], companion=pnames[0], P_d=fmt(P,10), e_P_d=err(pl.find('period')),
                K_ms='', e=fmt(fnum(pl.findtext('eccentricity')),6), omega_deg=fmt(fnum(pl.findtext('periastron')),6), Tp_bjd='',
                msini_mj=fmt(m,6), a_au=fmt(a,6), mstar_msun=fmt(fnum(star.findtext('mass')),4), fit_type='RV-only',
                table_ref='exoplanet.eu via OpenExoplanetCatalogue/oec_gzip systems-exoplaneteu.xml.gz (nightly mirror, downloaded 2026-10-24)',
                notes='; '.join(n for n in notes if n)))
rows.sort(key=lambda r:-float(r['msini_mj']))
write(out, rows)
print(len(rows))
for k,v in sorted(stats.items(), key=str): print(k,v)
