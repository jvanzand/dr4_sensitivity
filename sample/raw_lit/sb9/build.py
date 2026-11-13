"""Build sb9_orbits.csv (and sb9_excluded_no_primary_mass.csv) from SB9, the 9th Catalogue of Spectroscopic Binary
Orbits (Pourbaix et al. 2004; VizieR B/sb9, version 2024-04-22): single-lined orbits with a substellar-to-low-mass
companion, with a rough M2 sin i from the mass function and an estimated primary mass.

Run:  python -I build.py [--download] [--refresh]   (from anywhere; paths are relative to this file)

Raw inputs (VizieR TAP queries 'SELECT * FROM "B/sb9/<table>"', CSV; TAP_URL below). Re-downloading them gave
byte-identical files on 2026-10-26. raw_lit/sb9_tic/build.py also reads these three files.
  tap_main.csv     B/sb9/main    (systems: position, magnitudes, spectral types, name)
  tap_orbits.csv   B/sb9/orbits  (orbits: P, T0, e, omega, K1, K2, grade, reference, flags)
  tap_alias.csv    B/sb9/alias   (other designations, including 'GAIADR3 <source_id>')
  (The FTP copies main.dat, orbits.dat, alias.dat, notes.txt and ReadMe from https://cdsarc.cds.unistra.fr/ftp/
  B/sb9/ are kept for reference but not read: orbits.dat is pipe-delimited, not in the ReadMe's fixed-width layout.)

Method:
  1. SB1 orbits only (no K2), P <= 4000 d.
  2. f(m) = 1.0361e-7 (1-e^2)^1.5 K1^3 P  [Msun; K1 km/s, P d]; solve (M2 sin i)^3/(M1 + M2 sin i)^2 = f(m).
  3. Prefilter: M2 sin i < 200 MJ for M1 = 1 Msun.
  4. M1: (a) Gaia DR3 FLAME mass (Gaia DR3 id from the SB9 aliases); else (b) if Gaia DR3 GSP-Phot logg < 3.8,
     excluded; else (c) Mamajek dwarf mass for Sp1 (no subtype -> 5; class IV -> dwarf mass), rejected if
     giant/supergiant, unparsable, or overluminous by > 1.5 mag in Gaia M_G (or M_V) vs the dwarf value, or if Sp1
     has no luminosity class and no parallax check is possible. Excluded orbits go to
     sb9_excluded_no_primary_mass.csv with the reason and M2 sin i for M1 = 1 Msun.
  5. Keep M2 sin i <= 100 MJ; a from Kepler III with M1 + M2 sin i. When a system has several orbits in range,
     keep the most recent reference (then the highest grade) and list the others in the notes.

Lookups, cached in this directory (used offline; redone over the network only with --refresh):
  gaia_dr3_flame_for_candidates.csv       Gaia DR3 astrophysical parameters (VizieR I/355/paramp: Teff, logg,
                                          FLAME mass and its bounds, flags, evolution stage) for the Gaia DR3 ids
                                          of the prefiltered systems (ids not in I/355/paramp are absent)
  gaia_dr3_astrometry_for_candidates.csv  Gaia DR3 parallax, G, RUWE (VizieR I/355/gaiadr3) for the same ids
  EEM_dwarf_UBVIJHK_colors_Teff.txt       Mamajek dwarf table (version 2022.04.16; a refresh downloads the current
                                          version, which may change some M1 values); read by column name
  simbad_cache.csv                        CDS Sesame (SIMBAD) answers for the host names (main_id, otype, ids)
  nea_pscomppars_cache.csv                NASA Exoplanet Archive pscomppars rows reachable by host-ID match
  literature_companions_cache.csv         CLS I + Kiefer+2019 rows of ../literature_companions.csv, with SIMBAD ids
The Gaia and Mamajek files are the original pass's own downloads (2026-10-24). The SIMBAD/NEA/literature caches
were seeded from that pass's own query results (its Sesame cache, its pscomppars download, and the 239-row
literature_companions.csv of that day), which reproduce the table exactly; a --refresh can change the notes.

Outputs: sb9_orbits.csv (14-column table; header keeps the older name Tp_bjd; T0 is JD as published, not
necessarily BJD), sb9_excluded_no_primary_mass.csv (read by ../sb9_tic/build.py, so run this script first).
"""
import csv
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
import math

TAP_URL = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'


def tap_query_url(adql):
    return TAP_URL + '?' + urllib.parse.urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
                                                   'QUERY': adql})


# Raw inputs (path relative to this directory -> URL)
SOURCES = {'tap_%s.csv' % t: tap_query_url('SELECT * FROM "B/sb9/%s"' % t) for t in ('main', 'orbits', 'alias')}
OUT = 'sb9_orbits.csv'
EXCL = 'sb9_excluded_no_primary_mass.csv'
EXCL_COLS = ['sb9_seq', 'orbit', 'name', 'Sp1', 'P_d', 'K1_kms', 'e', 'msini_mj_if_M1_1msun', 'orbit_ref',
             'reason_excluded']
FLAME_CACHE = 'gaia_dr3_flame_for_candidates.csv'
ASTRO_CACHE = 'gaia_dr3_astrometry_for_candidates.csv'
EEM_CACHE = 'EEM_dwarf_UBVIJHK_colors_Teff.txt'
EEM_URL = 'https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt'

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


# ---------------------------------------------------------------- mass function
def m2_from_fm(fm, m1):
    """Solve m2^3/(m1+m2)^2 = fm for m2 (Msun) by bisection on [0, 200]."""
    lo, hi = 0., 200.
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid ** 3 / (m1 + mid) ** 2 > fm:
            hi = mid
        else:
            lo = mid
    return lo


def a_au(P_d, mtot):
    return (mtot * (P_d / 365.25) ** 2) ** (1 / 3)


# ---------------------------------------------------------------- primary-mass lookups (cached)
def refresh_gaia(gids):
    """Query VizieR for the Gaia DR3 FLAME parameters and astrometry of the candidate ids; rewrite the caches."""
    queries = {FLAME_CACHE: 'SELECT Source, Teff, logg, "Mass-Flame", "b_Mass-Flame", "B_Mass-Flame", '
                            '"Flags-Flame", Evol FROM "I/355/paramp" WHERE Source IN (%s)',
               ASTRO_CACHE: 'SELECT Source, Plx, e_Plx, Gmag, RUWE FROM "I/355/gaiadr3" WHERE Source IN (%s)'}
    for fn, q in queries.items():
        print('querying VizieR for', fn)
        rows, cols = [], None
        for i in range(0, len(gids), 200):
            data = urllib.parse.urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
                                           'QUERY': q % ', '.join(gids[i:i + 200])}).encode()
            rd = csv.DictReader(io.StringIO(http_get(TAP_URL, data=data).decode('utf-8')))
            rows += list(rd)
            cols = rd.fieldnames
        write_csv(here(fn), sorted(rows, key=lambda r: r['Source']), cols)


def load_eem():
    """Mamajek dwarf table -> {SpT: Msun}, {SpT: (M_G, M_V)} for rows with a mass (columns found by header name,
    between the first two '#SpT' header lines)."""
    with open(here(EEM_CACHE), encoding='utf-8') as fh:
        lines = fh.read().split('\n')
    h = [i for i, l in enumerate(lines) if l.startswith('#SpT')]
    cols = lines[h[0]][1:].split()
    sptm, absmag = {}, {}
    for l in lines[h[0] + 1:h[1]]:
        t = l.split()
        if len(t) != len(cols):
            continue
        r = dict(zip(cols, t))
        if fnum(r['Msun']) is not None:
            sptm[r['SpT']] = float(r['Msun'])
            absmag[r['SpT']] = (fnum(r['M_G']), fnum(r['Mv']))
    return sptm, absmag


def spt_mass(sp, sptm):
    """Dwarf mass for an SB9 spectral type: (mass or None, note['|' + table key])."""
    s = (sp or '').strip()
    m = re.match(r'^([OBAFGKM])\s*(\d(?:\.\d)?)?', s)
    if not m:
        return None, 'no parsable SpT'
    rest = s[m.end():]
    cls = 'V'
    if re.search(r'(III|II|Ia|Ib|I(?![V]))', rest.replace('IV', '')):
        cls = 'giant'
    elif re.search(r'IV', rest):
        cls = 'IV'
    elif not re.search(r'V', rest):
        cls = 'none'
    if cls == 'giant':
        return None, 'giant/supergiant SpT (%s); no dwarf-table mass' % s
    sub = m.group(2) or '5'
    key = m.group(1) + sub + 'V'
    if key not in sptm:
        key = m.group(1) + str(int(float(sub))) + 'V'
    if key not in sptm:
        return None, 'SpT %s not in table' % s
    note = 'SpT %s -> %s in Mamajek dwarf table' % (s, key)
    if not m.group(2):
        note += ' (no subtype; assumed 5)'
    if cls == 'IV':
        note += ' (class IV; dwarf mass used, likely underestimate)'
    if cls == 'none':
        note += ' (no luminosity class; dwarf assumed)'
    return sptm[key], note + '|' + key


def primary_mass(mrow, gids, flame, gast, sptm, absmag):
    """(M1 in Msun or None, provenance note / exclusion reason)."""
    for g in gids:
        fr = flame.get(g)
        if fr and fnum(fr['Mass-Flame']):
            return fnum(fr['Mass-Flame']), 'M1 = Gaia DR3 FLAME mass (Gaia DR3 %s; assumes single star)' % g
    lg = [fnum(flame[g]['logg']) for g in gids if g in flame and fnum(flame[g]['logg']) is not None]
    if lg and lg[0] < 3.8:
        return None, 'Gaia DR3 GSP-Phot logg=%.2f < 3.8 (evolved) and no FLAME mass' % lg[0]
    m1, src = spt_mass(mrow['Sp1'], sptm)
    if m1 is None:
        return None, src
    src, key = src.split('|')
    MG, MV = absmag[key]
    chk = None
    ga = [gast[g] for g in gids if g in gast and fnum(gast[g]['Plx']) and fnum(gast[g]['Plx']) > 0]
    if ga:
        plx, G = fnum(ga[0]['Plx']), fnum(ga[0]['Gmag'])
        if MG is not None and G is not None:
            chk = ('M_G', G + 5 * math.log10(plx / 100.), MG)
        elif MV is not None and mrow['n_mag1'].strip() == 'V' and fnum(mrow['mag1']) is not None:
            chk = ('M_V', fnum(mrow['mag1']) + 5 * math.log10(plx / 100.), MV)
    if chk and chk[1] < chk[2] - 1.5:
        return None, ('SpT %s but Gaia DR3 parallax gives %s=%.2f vs dwarf %.2f (overluminous; evolved or '
                      'misclassified)' % (mrow['Sp1'].strip(), chk[0], chk[1], chk[2]))
    elif chk:
        src += '; dwarf status consistent with Gaia DR3 parallax (%s=%.2f vs table %.2f)' % chk
    elif 'no luminosity class' in src:
        return None, ('SpT %s has no luminosity class and no Gaia DR3 parallax check is possible; excluded'
                      % mrow['Sp1'].strip())
    else:
        src += '; dwarf status NOT verified (no Gaia DR3 parallax/abs-mag check)'
    src += '; Gaia DR3 GSP-Phot logg=%.2f' % lg[0] if lg else '; no Gaia DR3 logg to confirm dwarf status'
    return m1, src


# ---------------------------------------------------------------- build
def candidates_in(orb):
    """(orbit row, P, K1, e, f(m)) for SB1 orbits with P <= 4000 d and M2 sin i(M1 = 1 Msun) < 200 MJ."""
    out = []
    for o in orb:
        if fnum(o['K2']) is not None:
            continue
        P, K, e = fnum(o['Per']), fnum(o['K1']), fnum(o['e']) or 0.0
        if P is None or K is None or P > 4000:
            continue
        fm = 1.0361e-7 * (1 - e * e) ** 1.5 * K ** 3 * P
        if m2_from_fm(fm, 1.0) / MJ_MSUN >= 200:
            continue
        out.append((o, P, K, e, fm))
    return out


def build_rows(cands, main, alias, flame, gast, sptm, absmag):
    rows, skipped = [], []
    for o, P, K, e, fm in cands:
        seq = o['Seq']
        mrow, al = main[seq], alias.get(seq, [])
        gids = [a.split()[1] for a in al if a.startswith('GAIADR3 ')]
        m1, src = primary_mass(mrow, gids, flame, gast, sptm, absmag)
        if m1 is None:
            skipped.append(dict(zip(EXCL_COLS, (seq, o['o'], mrow['Name'].strip(), mrow['Sp1'].strip(), P, K, e,
                                                round(m2_from_fm(fm, 1.0) / MJ_MSUN, 1), o['Ref'].strip(), src))))
            continue
        ms = m2_from_fm(fm, m1) / MJ_MSUN
        if ms > 100:
            continue
        host = re.sub(r'^(BD|CD|CPD)\s+([+-])', r'\1\2', mrow['Name'].strip())
        flags = []
        if o['f_Per'].strip() == 'a':
            flags.append('P fixed')
        if o['f_e'].strip():
            flags.append('e flag %s' % o['f_e'].strip())
        if o['u_K1'].strip():
            flags.append('K1 uncertain (:)')
        if o['f_K1'].strip():
            flags.append('K1 flag %s' % o['f_K1'].strip())
        if fnum(o['f_T0']):
            flags.append('T0 is not periastron (flag %s)' % o['f_T0'])
        # (the original pass also wrote a name-match "host in NASA Exoplanet Archive" flag after Sp1, which its
        # finalize step then deleted again)
        note = ('SB9 Seq %s orbit %s; SB1. f(m)=1.0361e-7(1-e^2)^1.5 K1^3 P = %.3g Msun; M2sini solved from '
                'f(m)=(M2 sini)^3/(M1+M2 sini)^2 (sin^3 i -> M2 sin i approximation, minimum mass); %s; ROUGH. SB9 '
                'grade=%s; Sp1=%s; T0 is JD as listed in SB9 (time system as published, not necessarily BJD); a from '
                'Kepler III with M1+M2sini'
                % (seq, o['o'], fm, src, o['Grade'] or 'none', mrow['Sp1'].strip() or 'none'))
        if flags:
            note += '; flags: ' + ', '.join(flags)
        # SB9 publishes 1-sigma errors on P (e_Per, d), e (e_e) and K1 (e_K1, km/s -> m/s; K1 is used as K); 0 or blank
        # means none. M2 sin i and a are computed here, not published, so e_msini_mj and e_a_au stay blank.
        def sb9_err(col, scale=1):
            v = fnum(o[col])
            return fmt(v * scale, 6) if v else ''
        rows.append(dict(host=host, companion='%s (SB9 %s.%s)' % (host, seq, o['o']), P_d=fmt(P, 10),
                         e_P_d=o['e_Per'], K_ms=fmt(K * 1000, 6), e_K_ms=sb9_err('e_K1', 1000), e_e=sb9_err('e_e'),
                         e=o['e'], omega_deg=o['omega'], Tp_bjd=o['T0'],
                         msini_mj=fmt(ms, 4), a_au=fmt(a_au(P, m1 + ms * MJ_MSUN), 4), mstar_msun=fmt(m1, 4),
                         fit_type='RV-only',
                         table_ref='SB9 (Pourbaix+2004; VizieR B/sb9, version 2024-04-22); orbit ref %s'
                                   % o['Ref'].strip(),
                         notes=note))
    # one row per system: most recent reference year, then highest grade
    best = {}
    for r in rows:
        seq = r['companion'].split('SB9 ')[1].split('.')[0]
        yr = int(re.search(r'orbit ref (\d{4})', r['table_ref']).group(1))
        gr = float(re.search(r'SB9 grade=([\d.]+|none)', r['notes']).group(1).replace('none', '-1'))
        best.setdefault(seq, []).append(((yr, gr), r))
    rows = []
    for seq, lst in best.items():
        lst.sort(key=lambda t: t[0], reverse=True)
        r = lst[0][1]
        if len(lst) > 1:
            r['notes'] += '; %d SB9 orbits for this system in range, kept most recent (others: %s)' % (
                len(lst), ', '.join(re.search(r'orbit ref (\S+)', x[1]['table_ref']).group(1) + ' msini=%s'
                                    % x[1]['msini_mj'] for x in lst[1:]))
        rows.append(r)
    rows.sort(key=lambda r: float(r['msini_mj']))
    return rows, skipped


def main():
    get_inputs(SOURCES)
    main_ = {r['Seq']: r for r in read_csv(here('tap_main.csv'))}
    orb = read_csv(here('tap_orbits.csv'))
    alias = {}
    for a in read_csv(here('tap_alias.csv')):
        alias.setdefault(a['Seq'], []).append(a['Name'])
    cands = candidates_in(orb)
    if REFRESH:
        gids = sorted({a.split()[1] for o in cands for a in alias.get(o[0]['Seq'], []) if a.startswith('GAIADR3 ')})
        refresh_gaia(gids)
        print('downloading', EEM_URL)
        with open(here(EEM_CACHE), 'wb') as fh:
            fh.write(http_get(EEM_URL))
    need_caches([FLAME_CACHE, ASTRO_CACHE, EEM_CACHE])
    flame = {r['Source']: r for r in read_csv(here(FLAME_CACHE))}
    gast = {r['Source']: r for r in read_csv(here(ASTRO_CACHE))}
    sptm, absmag = load_eem()
    rows, skipped = build_rows(cands, main_, alias, flame, gast, sptm, absmag)
    alt = {}  # no alternate host names in this source
    cache, nea, lit = lookups(rows, alt)
    rows = finalize(rows, cache, nea, lit, alt)
    write_csv(here(OUT), rows, COLS)
    with open(here(EXCL), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(EXCL_COLS)
        for r in skipped:
            w.writerow([r[k] for k in EXCL_COLS])
    print('wrote %d rows to %s (%d with M sin i <= 80 MJ; M1 from FLAME for %d) and %d excluded orbits to %s'
          % (len(rows), OUT, sum(float(r['msini_mj']) <= 80 for r in rows), sum('FLAME' in r['notes'] for r in rows),
             len(skipped), EXCL))


if __name__ == '__main__':
    main()
