"""Build exoplanet_eu_orbits.csv: exoplanet.eu companions, read from the Open Exoplanet Catalogue's nightly mirror.

Run:  python -I build.py [--download] [--refresh]   (from anywhere; paths are relative to this file)

Raw input (exoplanet.eu and its TAP service timed out from this machine in 2026-10, so the mirror is used; the
OEC oec_continuity job converts http://exoplanet.eu/catalog/csv/ to OEC XML every night):
  ../oec/systems-exoplaneteu.xml.gz
      https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems-exoplaneteu.xml.gz
      (downloaded 2026-10-24; a fresh download is a newer snapshot and will change the table)
Uncertainties: the mirror carries errorminus/errorplus on period, eccentricity, semi-major axis and mass (no K, so
no e_K_ms). e_P_d is the larger side (original rule); new columns e_msini_mj (on the same "mass" as msini_mj),
e_e and e_a_au are the mean of the two sides; blank when a side is missing or 0. These are exoplanet.eu's
compiled values, taken from the paper cited there.
The mirror keeps P, e, omega, a, i, mass, star mass, discovery method and alternate names. It drops K, Tp,
mass_detection_type and the separate M sin i column. Its "mass" is exoplanet.eu's true mass when one exists, else
M sin i; msini_mj is that "mass", and the notes say "inclination X deg listed" when an inclination is given
(sample/literature.py parses that string and the "otype <type>" in the SIMBAD note).

Selection: mass <= 80 MJ and (P <= 4000 d, or a <= 5 AU when P is missing). Keep RV discoveries with mass >= 13 MJ
or not name-matched to the NASA Exoplanet Archive (NEA), and transit discoveries with mass >= 13 MJ that are not
name-matched to the NEA. After the SIMBAD step, rows below 13 MJ whose companion the SIMBAD-ID match finds in the
NEA are dropped. Rows are sorted by decreasing mass.

Lookups, cached in this directory (used offline; redone over the network only with --refresh):
  simbad_cache.csv                  CDS Sesame (SIMBAD) answers for the host-name candidates (main_id, otype, ids)
  nea_pscomppars_cache.csv          NEA pscomppars rows reachable by planet-name or host-ID match
  literature_companions_cache.csv   CLS I + Kiefer+2019 rows of ../literature_companions.csv, with SIMBAD ids
The caches were seeded from the original pass's own query results of 2026-10-24 (its Sesame cache, its pscomppars
download, and the 239-row literature_companions.csv of that day), which reproduce the table exactly. A --refresh
queries today's SIMBAD and NEA and reads today's ../literature_companions.csv, so notes (and, through the NEA
name match, the row selection) can change; that is expected.
"""
import csv
import gzip
import io
import os
import re
import sys
import tarfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

# Raw inputs (path relative to this directory -> URL)
SOURCES = {'../oec/systems-exoplaneteu.xml.gz':
           'https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems-exoplaneteu.xml.gz'}
OUT = 'exoplanet_eu_orbits.csv'
TABLE_REF = ('exoplanet.eu via OpenExoplanetCatalogue/oec_gzip systems-exoplaneteu.xml.gz '
             '(nightly mirror, downloaded 2026-10-24)')

# ----------------------------------------------------------------------------------------------------------
# Shared helpers. The same block is copied into each aggregate-catalog build.py (exoplanet_eu, exoplanets_org,
# stevenson_2023, gaia_dr3_val, feng_2022, xiao_2023, sb9) so that every script is self-contained; if you fix
# something here, fix it in the others too.
# ----------------------------------------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))

COLS = ['host', 'companion', 'P_d', 'e_P_d', 'K_ms', 'e', 'omega_deg', 'Tp_bjd', 'msini_mj', 'a_au',
        'mstar_msun', 'fit_type', 'table_ref', 'notes', 'e_K_ms', 'e_msini_mj', 'e_e', 'e_a_au']
MJ_MSUN = 1 / 1047.57  # M_Jup / M_sun (IAU nominal masses: 1.8982e27 kg / 1.98847e30 kg)


def here(*parts):
    return os.path.join(HERE, *parts)


def fnum(x):
    """float(x), or None for blanks, non-numbers and NaN."""
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def fmt(x, nd=6):
    """Floats with nd significant digits (%g); None -> ''; anything else via str()."""
    if x is None:
        return ''
    if isinstance(x, float):
        return '%.*g' % (nd, x)
    return str(x)


def norm(n):
    """Name key for loose matching: drop whitespace, '_' and '-', upper case."""
    return re.sub(r'[\s_\-]', '', n or '').upper()


GREEK = {'alpha': 'alf', 'beta': 'bet', 'gamma': 'gam', 'delta': 'del', 'epsilon': 'eps', 'zeta': 'zet',
         'eta': 'eta', 'theta': 'tet', 'iota': 'iot', 'kappa': 'kap', 'lambda': 'lam', 'mu': 'mu.', 'nu': 'nu.',
         'xi': 'ksi', 'omicron': 'omi', 'pi': 'pi.', 'rho': 'rho', 'sigma': 'sig', 'tau': 'tau', 'upsilon': 'ups',
         'phi': 'phi', 'chi': 'chi', 'psi': 'psi', 'omega': 'ome'}


def clean_name(n):
    """Catalogue host name -> SIMBAD-friendly form (spacing, LaTeX Greek letters, BD/CD zones)."""
    n = (n or '').replace('−', '-').replace('â\u0088\u0092', '-').strip()
    n = re.sub(r'\$\\?([a-zA-Z]+)\$', lambda m: GREEK.get(m.group(1).lower(), m.group(1)), n)
    n = re.sub(r'^(HD|HIP|GJ|HR|LHS|Gl|NGC|HAT-P-|WASP-)(?=\d)',
               lambda m: m.group(1) if m.group(1).endswith('-') else m.group(1) + ' ', n)
    n = re.sub(r'^Gl ', 'GJ ', n)
    n = re.sub(r'^TYC(?=\d)', 'TYC ', n)
    n = re.sub(r'^(BD|CD|CPD)\s*([+-])\s*(\d{2})\s*0*(\d+)$',
               lambda m: '%s%s%s %s' % (m.group(1), m.group(2), m.group(3), m.group(4)), n)
    n = re.sub(r'^(HD|HIP|GJ) (\d+)([A-C])$', r'\1 \2 \3', n)
    n = re.sub(r'^(\d+)([A-Z][a-z]{2})$', r'\1 \2', n)
    n = re.sub(r'^NGC (\d+)-(\d+)$', r'NGC \1 \2', n)
    return n


def read_csv(path, encoding='utf-8'):
    with open(path, newline='', encoding=encoding) as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, cols):
    """csv default dialect (CRLF line ends, minimal quoting), as the original tables were written."""
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in cols})


def http_get(url, data=None, timeout=300):
    # identify as a script (cdsarc.cds.unistra.fr answers browser-like User-Agents with a bot-check page)
    ua = 'build.py (gaia_dr4_rv_sensitivity; python-urllib)'
    req = urllib.request.Request(url, data=data, headers={'User-Agent': ua})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(3)


def get_inputs(sources):
    """Check (or with --download fetch) every raw input.

    sources: {path relative to this directory: URL, or (URL of an arXiv source tarball, member to extract)}.
    Downloads are only read as data (one tar member is copied out by name); nothing is unpacked or executed."""
    for fn, src in sources.items():
        url, member = src if isinstance(src, tuple) else (src, None)
        path = here(fn)
        if os.path.exists(path):
            continue
        if not DOWNLOAD:
            sys.exit('missing input %s (from %s%s); rerun with --download' % (
                fn, url, ', member %s' % member if member else ''))
        print('downloading', url)
        data = http_get(url)
        if data[:200].lstrip().lower().startswith((b'<!doctype html', b'<html')):
            sys.exit('%s returned an HTML page instead of data (bot check or error page?); fetch %s by hand'
                     % (url, fn))
        if member:
            with tarfile.open(fileobj=io.BytesIO(data), mode='r:*') as tf:
                data = tf.extractfile(tf.getmember(member)).read()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path + '.part', 'wb') as fh:
            fh.write(data)
        os.replace(path + '.part', path)


def need_caches(names):
    """Offline runs read the cache files; exit with a clear message if one is missing."""
    for fn in names:
        if not os.path.exists(here(fn)):
            sys.exit('missing cache %s; rerun with --refresh (network) to rebuild it' % fn)


DOWNLOAD = '--download' in sys.argv[1:]
REFRESH = '--refresh' in sys.argv[1:]

# --- SIMBAD (CDS Sesame) -----------------------------------------------------------------------------------
# simbad_cache.csv: one row per name sent to Sesame: name, resolved (1/0), main_id, otype, aliases ('|'-joined;
# only HD/HIP/Gaia DR3/TIC/GJ/2MASS/TYC/BD/CD/KOI/Kepler/TOI/WASP/HAT/NAME/* identifiers). A name that is not in
# the cache counts as unresolved.
SESAME = 'https://cds.unistra.fr/cgi-bin/nph-sesame/-oxI/S?'
SIMBAD_CACHE = 'simbad_cache.csv'
SIMBAD_COLS = ['name', 'resolved', 'main_id', 'otype', 'aliases']


def sesame(name):
    try:
        root = ET.fromstring(http_get(SESAME + urllib.parse.quote(name), timeout=30).decode('utf-8', 'replace'))
    except Exception as e:  # noqa: BLE001
        return {'error': str(e)}
    res = root.find('.//Resolver')
    if res is None or res.find('oname') is None:
        return {'resolved': False}
    al = [a.text for a in res.findall('alias') if a.text]
    keep = r'^(HD|HIP|Gaia DR3|TIC|GJ|2MASS|TYC|BD|CD|KOI|Kepler|TOI|WASP|HAT|NAME|\*|V\*)'
    return {'resolved': True, 'main_id': ' '.join(res.findtext('oname').split()), 'otype': res.findtext('otype'),
            'aliases': [a for a in al if re.match(keep, a)]}


def candidates(h, alt):
    """Names to try in SIMBAD for host h, in order: as given, cleaned, without a trailing (AB)/A or '.0N',
    TYC zero padding removed, then alternate names alt[h] (exoplanet.eu only)."""
    c = [h, clean_name(h)]
    for x in list(c):
        c.append(re.sub(r'\s*\(?AB\)?$', '', x))
        c.append(re.sub(r'\s*A$', '', x))
        c.append(re.sub(r'\.0\d$', '', x))
        c.append(re.sub(r'^TYC (\d+)-0*(\d+)-(\d)$',
                        lambda m: 'TYC %d-%d-%s' % (int(m.group(1)), int(m.group(2)), m.group(3)), x))
    for a in alt.get(h, []):
        c += [a, clean_name(a), re.sub(r'\s*A$', '', a)]
    out = []
    for x in c:
        x = x.strip()
        if x and x not in out:
            out.append(x)
    return out


def first_hit(h, cache, alt):
    """First candidate name for h that SIMBAD resolves (None if none does)."""
    for c in candidates(h, alt):
        if cache.get(c, {}).get('resolved'):
            return c
    return None


def load_simbad():
    cache = {}
    for r in read_csv(here(SIMBAD_CACHE)):
        cache[r['name']] = ({'resolved': True, 'main_id': r['main_id'], 'otype': r['otype'],
                             'aliases': r['aliases'].split('|') if r['aliases'] else []}
                            if r['resolved'] == '1' else {'resolved': False})
    return cache


def refresh_simbad(hosts, alt):
    """Query Sesame for each host's candidate names in order, stopping at the first that resolves."""
    def one(h):
        out = []
        for c in candidates(h, alt):
            r = sesame(c)
            out.append((c, r))
            if r.get('resolved'):
                break
        return out
    cache = {}
    with ThreadPoolExecutor(8) as ex:
        for out in ex.map(one, sorted(set(hosts))):
            cache.update(out)
    bad = sorted(n for n, v in cache.items() if 'error' in v)
    if bad:
        sys.exit('Sesame failed for %d names (e.g. %s); caches not written, try again' % (len(bad), bad[:3]))
    rows = [{'name': n, 'resolved': int(bool(v['resolved'])), 'main_id': v.get('main_id', ''),
             'otype': v.get('otype') or '', 'aliases': '|'.join(v.get('aliases', []))}
            for n, v in sorted(cache.items())]
    write_csv(here(SIMBAD_CACHE), rows, SIMBAD_COLS)
    return load_simbad()


def simbad_ids(h, cache, alt):
    hit = first_hit(h, cache, alt)
    return cache[hit]['aliases'] + [cache[hit]['main_id']] if hit else []


# --- NASA Exoplanet Archive (pscomppars) -------------------------------------------------------------------
# nea_pscomppars_cache.csv: the pscomppars rows this source can match (by host identifier, and for
# exoplanet_eu/exoplanets_org also by planet name), in archive order.
NEA_CACHE = 'nea_pscomppars_cache.csv'
NEA_COLS = ['pl_name', 'hostname', 'hd_name', 'hip_name', 'gaia_dr3_id', 'pl_bmassj', 'pl_orbper']
NEA_URL = ('https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query=' +
           urllib.parse.quote('select %s from pscomppars' % ','.join(NEA_COLS)))


def fetch_nea():
    print('querying NASA Exoplanet Archive pscomppars')
    return [{k: r[k] for k in NEA_COLS} for r in csv.DictReader(io.StringIO(http_get(NEA_URL).decode('utf-8')))]


def nea_id_match(r, ids):
    return any(r[k] and r[k].strip() in ids for k in ('hd_name', 'hip_name', 'gaia_dr3_id', 'hostname'))


# --- project literature list -------------------------------------------------------------------------------
# literature_companions_cache.csv: the CLS I (Rosenthal+2021) and Kiefer+2019 rows of ../literature_companions.csv
# (all that file held when these tables were first made), with each host's SIMBAD identifiers ('|'-joined).
LIT_CACHE = 'literature_companions_cache.csv'
LIT_COLS = ['pl_name', 'hostname', 'catalog', 'pl_orbper', 'source_id', 'simbad_ids']
LIT_CATALOGS = ('CLS I (Rosenthal+2021)', 'Kiefer+2019')
LIT_SOURCE = here('..', 'literature_companions.csv')  # read only by --refresh


def lookups(rows, alt, nea_full=None, name_match=None):
    """SIMBAD, NEA and literature-list lookups for the finalize step: (simbad cache, NEA rows, literature rows).

    Offline (default) everything comes from the cache files in this directory. With --refresh the lookups are
    redone over the network and the cache files rewritten. nea_full: the full pscomppars table if the caller
    already fetched it; name_match(nea_row) -> True keeps NEA rows that the caller's name matching can reach."""
    if not REFRESH:
        need_caches([SIMBAD_CACHE, NEA_CACHE, LIT_CACHE])
        lit = read_csv(here(LIT_CACHE))
        for x in lit:
            x['simbad_ids'] = x['simbad_ids'].split('|') if x['simbad_ids'] else []
        return load_simbad(), read_csv(here(NEA_CACHE)), lit
    lit = [{k: r[k] for k in LIT_COLS[:-1]} for r in read_csv(LIT_SOURCE)
           if r['catalog'] in LIT_CATALOGS]
    cache = refresh_simbad([r['host'] for r in rows] + [x['hostname'] for x in lit], alt)
    for x in lit:
        x['simbad_ids'] = simbad_ids(x['hostname'], cache, alt)
    write_csv(here(LIT_CACHE), [dict(x, simbad_ids='|'.join(x['simbad_ids'])) for x in lit], LIT_COLS)
    nea = nea_full if nea_full is not None else fetch_nea()
    ids = set(i for r in rows for i in simbad_ids(r['host'], cache, alt))
    nea = [r for r in nea if nea_id_match(r, ids) or (name_match is not None and name_match(r))]
    write_csv(here(NEA_CACHE), nea, NEA_COLS)
    return cache, nea, lit


def finalize(rows, cache, nea, lit, alt, drop_nea_lt13=False):
    """SIMBAD host resolution plus NEA and literature-list overlap notes (the original pass's finalize.py).

    - host: replaced by the first candidate name SIMBAD resolves, if that differs from the source's name
      (notes: 'host renamed from "..." for SIMBAD'); unresolved hosts are flagged and kept.
    - notes += 'SIMBAD main_id=... (otype ...)'  (downstream code parses the otype from this text).
    - NEA overlap: NEA planets whose hd_name/hip_name/gaia_dr3_id/hostname equals one of the host's SIMBAD
      identifiers ('host yes'), and of those, one whose period is within 10% of P_d ('companion YES').
    - drop_nea_lt13 (exoplanet_eu, exoplanets_org): drop rows below 13 MJ whose companion is in the NEA."""
    nea_by = {}
    for p in nea:
        for k in ('hd_name', 'hip_name', 'gaia_dr3_id', 'hostname'):
            if p[k]:
                nea_by.setdefault(p[k].strip(), []).append(p)
    lit_by = {}
    for x in lit:
        for i in x['simbad_ids']:
            lit_by.setdefault(i, []).append(x)
        if x.get('source_id'):
            lit_by.setdefault('Gaia DR3 ' + x['source_id'], []).append(x)

    def within10(x, P):
        Pr = fnum(x['pl_orbper'])
        return bool(P and Pr and abs(Pr - P) / Pr < 0.1)

    out = []
    for r in rows:
        hit = first_hit(r['host'], cache, alt)
        notes = r['notes']
        if hit is None:
            notes += '; host NOT resolved by SIMBAD/Sesame (name as given by source)'
            ids = []
        else:
            info = cache[hit]
            if hit != r['host']:
                notes += '; host renamed from "%s" for SIMBAD' % r['host']
                r['host'] = hit
            ids = info['aliases'] + [info['main_id']]
            notes += '; SIMBAD main_id=%s (otype %s)' % (info['main_id'], info['otype'])
        P = fnum(r['P_d'])
        hostnea = [p for i in ids for p in nea_by.get(i, [])]
        plmatch = [p for p in hostnea if within10(p, P)]
        litm = [x for i in ids for x in lit_by.get(i, []) if within10(x, P)]
        notes += '; NASA Exoplanet Archive overlap (SIMBAD HD/HIP/Gaia DR3 ID match, P within 10%%): %s' % (
            'companion YES (%s, mass %s MJ)' % (plmatch[0]['pl_name'], plmatch[0]['pl_bmassj']) if plmatch
            else ('host yes, companion no' if hostnea else 'no'))
        if litm:
            notes += '; also in project literature_companions.csv (%s)' % litm[0]['catalog']
        if drop_nea_lt13 and plmatch and (fnum(r['msini_mj']) or 0) < 13:
            continue
        r['notes'] = notes.lstrip('; ')
        out.append(r)
    return out
# ------------------------------------------------------------------------------- end of shared helpers ----


# ---------------------------------------------------------------- NEA planet-name matching
def nea_row_keys(p):
    """Normalised keys for one NEA planet: its name, and each host alias (minus a trailing ' A') + its letter."""
    name = p['pl_name'].strip()
    letter = name.split()[-1] if name else ''
    keys = {norm(p['pl_name'])}
    for k in ('hostname', 'hd_name', 'hip_name'):
        if p[k]:
            keys.add(norm(re.sub(r'\s+A$', '', p[k].strip()) + letter))
    return keys


def pl_keys(pnames, hostnames):
    """The same kind of keys for one exoplanet.eu planet (its names, and each host name + its letter)."""
    keys = set()
    for n in pnames:
        n = n.strip()
        keys.add(norm(n))
        keys.add(norm(re.sub(r'\s+A\s*([b-z])$', r' \1', n)))
        m = re.search(r'\s+A?([b-z])$', n)
        if m:
            for h in hostnames:
                keys.add(norm(re.sub(r'\s+A$', '', h.strip()) + m.group(1)))
    return keys


# ---------------------------------------------------------------- parse the mirror
def alt_names(root):
    """exoplanet.eu alternate names: first star name -> all star and system names (extra SIMBAD candidates)."""
    alt = {}
    for s in root.iter('system'):
        sn = [n.text for n in s.findall('name') if n.text]
        for st in s.iter('star'):
            stn = [n.text for n in st.findall('name') if n.text]
            if stn:
                alt.setdefault(stn[0], []).extend(stn + sn)
    return alt


def planets(root):
    """(planet element, star element, star names, system names) for every planet passing the mass/period cuts."""
    out = []
    for sysel in root.iter('system'):
        sysnames = [n.text for n in sysel.findall('name') if n.text]
        for star in sysel.iter('star'):
            snames = [n.text for n in star.findall('name') if n.text]
            for pl in star.findall('planet'):
                # OEC units: mass MJ, period d, semimajoraxis AU, star mass Msun (no conversions needed)
                m, P = fnum(pl.findtext('mass')), fnum(pl.findtext('period'))
                a = fnum(pl.findtext('semimajoraxis'))
                if m is None or m > 80:
                    continue
                if not ((P is not None and P <= 4000) or (P is None and a is not None and a <= 5)):
                    continue
                out.append((pl, star, snames, sysnames))
    return out


def build_rows(cands, nea):
    nk = set()
    for p in nea:
        nk |= nea_row_keys(p)

    def err(el):  # symmetric error = larger of the two one-sided errors (the original e_P_d column)
        if el is None:
            return ''
        a, b = fnum(el.get('errorminus')), fnum(el.get('errorplus'))
        if a is None or b is None:
            return ''
        return fmt(max(a, b), 4)

    def err_avg(el):  # new e_* columns: mean of the two one-sided errors; blank if either is missing or 0
        if el is None:
            return ''
        a, b = fnum(el.get('errorminus')), fnum(el.get('errorplus'))
        if a is None or b is None or a <= 0 or b <= 0:
            return ''
        return fmt((a + b) / 2, 4)

    rows = []
    for pl, star, snames, sysnames in cands:
        dm = pl.findtext('discoverymethod')
        m, P = fnum(pl.findtext('mass')), fnum(pl.findtext('period'))
        a = fnum(pl.findtext('semimajoraxis'))
        pnames = [n.text for n in pl.findall('name') if n.text]
        in_nea_pl = bool(pl_keys(pnames, snames + sysnames) & nk)
        if dm == 'RV':
            keep = m >= 13 or not in_nea_pl
        elif dm == 'transit':
            keep = m >= 13 and not in_nea_pl
        else:
            keep = False
        if not keep:
            continue
        inc = pl.findtext('inclination')
        mnote = ('mass column is exoplanet.eu "mass" (true mass when exoplanet.eu lists one, else M sin i)'
                 + ('; inclination %s deg listed, so this may be a TRUE mass' % inc if inc
                    else '; no inclination listed, likely M sin i'))
        # (the original pass also wrote name-match NEA flags here, which its finalize step then deleted again)
        notes = ['exoplanet.eu discovery method=%s' % dm, mnote,
                 'exoplanet.eu last update 20%s (yy/mm/dd)' % (pl.findtext('lastupdate') or '?'),
                 'alt names: %s' % '; '.join(pnames[1:5]) if len(pnames) > 1 else '',
                 'no K or Tp in this mirror (fetch from exoplanet.eu or the discovery paper)',
                 'fit_type not given by exoplanet.eu; RV-only assumed']
        # e_P_d: original rule (larger side); where that is blank (a side missing or 0), the mean-of-sides rule
        rows.append(dict(host=snames[0] if snames else sysnames[0], companion=pnames[0], P_d=fmt(P, 10),
                         e_P_d=err(pl.find('period')) or err_avg(pl.find('period')),
                         e_msini_mj=err_avg(pl.find('mass')), e_e=err_avg(pl.find('eccentricity')),
                         e_a_au=err_avg(pl.find('semimajoraxis')), K_ms='', e=fmt(fnum(pl.findtext('eccentricity')), 6),
                         omega_deg=fmt(fnum(pl.findtext('periastron')), 6), Tp_bjd='', msini_mj=fmt(m, 6),
                         a_au=fmt(a, 6), mstar_msun=fmt(fnum(star.findtext('mass')), 4), fit_type='RV-only',
                         table_ref=TABLE_REF, notes='; '.join(n for n in notes if n)))
    rows.sort(key=lambda r: -float(r['msini_mj']))  # stable: ties keep catalogue order
    return rows


def main():
    get_inputs(SOURCES)
    with gzip.open(here('../oec/systems-exoplaneteu.xml.gz')) as fh:
        root = ET.parse(fh).getroot()
    alt = alt_names(root)
    cands = planets(root)
    if REFRESH:
        nea_full = fetch_nea()
    else:
        need_caches([NEA_CACHE])
        nea_full = read_csv(here(NEA_CACHE))
    rows = build_rows(cands, nea_full)
    allkeys = set()
    for pl, star, snames, sysnames in cands:
        allkeys |= pl_keys([n.text for n in pl.findall('name') if n.text], snames + sysnames)
    cache, nea, lit = lookups(rows, alt, nea_full=nea_full, name_match=lambda p: bool(nea_row_keys(p) & allkeys))
    rows = finalize(rows, cache, nea, lit, alt, drop_nea_lt13=True)
    write_csv(here(OUT), rows, COLS)
    print('wrote %d rows to %s (%d with mass >= 13 MJ; %d hosts not resolved by SIMBAD)' % (
        len(rows), OUT, sum(float(r['msini_mj']) >= 13 for r in rows),
        sum('host NOT resolved' in r['notes'] for r in rows)))


if __name__ == '__main__':
    main()
