"""Shared helpers for mass-function SB1 sources (copied into latham_2002, nidever_2002, sb9_tic).

- resolve(names, cache): SIMBAD (via CDS Sesame) main id, otype, Gaia DR3 / TIC / HIP / HD aliases.
- gaia_query(ids): Gaia DR3 FLAME mass, GSP-Phot logg, RUWE, NSS flag, parallax, G.
- tic_query(ids): TIC v8 (VizieR IV/39/tic82) mass, logg, lumclass.
- msini_from_fm(f, m1): solve f = m2^3/(m1+m2)^2 for m2 (Msun).
All network results are cached as CSV/JSON in the calling directory.
"""
import csv, io, json, os, re, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

MJ_PER_MSUN = 1047.5654  # IAU 2015 nominal (GM_sun / GM_J)


def _get(url, data=None, timeout=60, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data)
            return urllib.request.urlopen(req, timeout=timeout).read().decode('utf-8', 'replace')
        except Exception as e:  # noqa
            err = e
            time.sleep(2)
    raise err


def _sesame(n):
    url = 'https://cds.unistra.fr/cgi-bin/nph-sesame/-oxI/S?' + urllib.parse.quote(n)
    try:
        t = _get(url, timeout=30)
        r = ET.fromstring(t)
    except Exception as e:  # noqa
        return n, {'error': str(e)}
    res = r.find('.//Resolver')
    if res is None or res.find('oname') is None:
        return n, {'resolved': False}
    al = [' '.join(a.text.split()) for a in res.findall('alias') if a.text]
    out = {'resolved': True, 'main_id': ' '.join(res.findtext('oname').split()),
           'otype': res.findtext('otype'), 'sptype': res.findtext('spType'),
           'ra': res.findtext('jradeg'), 'dec': res.findtext('jdedeg')}
    for key, pat in [('gaia_dr3', r'^Gaia DR3 (\d+)$'), ('gaia_dr2', r'^Gaia DR2 (\d+)$'),
                     ('tic', r'^TIC (\d+)$'), ('hip', r'^HIP (\d+)$'), ('hd', r'^HD (\d+[A-Z]?)$')]:
        m = [re.match(pat, a).group(1) for a in al if re.match(pat, a)]
        out[key] = m[0] if m else ''
    return n, out


def resolve(names, cache_path):
    cache = json.load(open(cache_path)) if os.path.exists(cache_path) else {}
    todo = [n for n in set(names) if n not in cache or 'error' in cache[n]]
    with ThreadPoolExecutor(8) as ex:
        for n, r in ex.map(_sesame, todo):
            cache[n] = r
    with open(cache_path, 'w') as f:
        json.dump(cache, f, indent=0, sort_keys=True)
    return cache


def _tap(url, query):
    data = urllib.parse.urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
                                   'QUERY': query}).encode()
    return _get(url, data=data, timeout=300)


def gaia_query(ids, cache_csv):
    """Gaia DR3 FLAME + GSP-Phot + basic astrometry for a list of source_ids (cached)."""
    if not os.path.exists(cache_csv):
        ids = sorted({str(i) for i in ids if i})
        rows, header = [], None
        for k in range(0, len(ids), 300):
            q = ('SELECT g.source_id, g.parallax, g.phot_g_mean_mag, g.ruwe, g.non_single_star, '
                 'g.radial_velocity, a.mass_flame, a.mass_flame_lower, a.mass_flame_upper, '
                 'a.flags_flame, a.evolstage_flame, a.logg_gspphot, a.teff_gspphot, a.mh_gspphot '
                 'FROM gaiadr3.gaia_source AS g LEFT OUTER JOIN gaiadr3.astrophysical_parameters AS a '
                 'ON g.source_id = a.source_id WHERE g.source_id IN (' + ','.join(ids[k:k + 300]) + ')')
            t = _tap('https://gea.esac.esa.int/tap-server/tap/sync', q)
            r = list(csv.reader(io.StringIO(t)))
            header = r[0]
            rows += r[1:]
        with open(cache_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
    return {r['source_id']: r for r in csv.DictReader(open(cache_csv))}


def tic_query(ids, cache_csv):
    """TIC v8 (VizieR IV/39/tic82) masses for a list of TIC ids (cached)."""
    if not os.path.exists(cache_csv):
        ids = sorted({str(i) for i in ids if i})
        rows, header = [], None
        for k in range(0, len(ids), 300):
            q = ('SELECT TIC, Mass, s_Mass, logg, Teff, Rad, LClass, Dist, GAIA FROM "IV/39/tic82" '
                 'WHERE TIC IN (' + ','.join(ids[k:k + 300]) + ')')
            t = _tap('https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync', q)
            r = list(csv.reader(io.StringIO(t)))
            header = r[0]
            rows += r[1:]
        with open(cache_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
    return {r['TIC']: r for r in csv.DictReader(open(cache_csv))}


def fm_from_orbit(P_d, K_kms, e):
    """Spectroscopic mass function in Msun."""
    return 1.0361e-7 * (1 - e * e) ** 1.5 * K_kms ** 3 * P_d


def msini_from_fm(f, m1):
    """Solve f = m2^3/(m1+m2)^2 for m2 by bisection (Msun); minimum-mass (sin i -> 1 in total mass)."""
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid ** 3 / (m1 + mid) ** 2 < f:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def a_au(P_d, mtot):
    return (P_d / 365.25) ** (2 / 3) * mtot ** (1 / 3)


def primary_mass(res, gaia, tic):
    """Pick M1: Gaia DR3 FLAME, else TIC v8. Returns (mass, source_text, flags list)."""
    flags = []
    g = gaia.get(res.get('gaia_dr3', ''), {}) if res else {}
    t = tic.get(res.get('tic', ''), {}) if res else {}
    if g:
        try:
            if float(g.get('ruwe') or 0) > 1.4:
                flags.append('Gaia RUWE=%.2f (>1.4; binary motion may bias FLAME/parallax)' % float(g['ruwe']))
        except ValueError:
            pass
        if g.get('non_single_star') not in ('', '0', None):
            flags.append('Gaia DR3 non_single_star=%s' % g['non_single_star'])
        if g.get('logg_gspphot'):
            lg = float(g['logg_gspphot'])
            flags.append('GSP-Phot logg=%.2f' % lg)
            if lg < 3.8:
                flags.append('EVOLVED (logg<3.8)')
    if g.get('mass_flame'):
        m = float(g['mass_flame'])
        src = 'Gaia DR3 FLAME mass (Gaia DR3 %s; %s-%s Msun; flags_flame=%s, evolstage=%s)' % (
            res['gaia_dr3'], g['mass_flame_lower'], g['mass_flame_upper'], g['flags_flame'],
            g['evolstage_flame'])
        try:
            if int(float(g['evolstage_flame'])) > 360:
                flags.append('EVOLVED (FLAME evolstage>360, post-MS)')
        except ValueError:
            pass
        return m, src, flags
    if t.get('Mass'):
        m = float(t['Mass'])
        src = 'TIC v8 mass (TIC %s; s_Mass=%s; TIC logg=%s; LClass=%s)' % (
            res['tic'], t.get('s_Mass', ''), t.get('logg', ''), t.get('LClass', ''))
        if t.get('LClass', '').strip().upper() == 'GIANT':
            flags.append('TIC LClass=GIANT (TIC mass unreliable for giants)')
        return m, src, flags
    return None, '', flags
