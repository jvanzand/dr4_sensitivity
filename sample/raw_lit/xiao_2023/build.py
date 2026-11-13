"""Build xiao_2023_orbits.csv from Xiao et al. 2023, RAA (accepted; arXiv:2303.12409): orvara fits of literature
RVs plus the Hipparcos-Gaia catalogue of accelerations (HGCA) for 115 RV companions.

Run:  python -I build.py [--download] [--refresh]   (from anywhere; paths are relative to this file)

Raw input:
  src/RAA-2022-0418.tex   arXiv LaTeX source (https://arxiv.org/e-print/2303.12409, member RAA-2022-0418.tex);
      table "Posteriors of RV Companions" (labels Tab:orbparas ... Tabcc, two parts). (src/*.bbl is not read.)
Uncertainties: every posterior cell prints +/- (16th/84th percentile) errors. e_P_d keeps the original rule (larger
side); new columns e_msini_mj (on the tabulated M sin i), e_e and e_a_au are the mean of the two sides; no K error.
Columns used: Mp [MJ], i [deg], a [AU], e, P [yr], omega [deg], Tp [JD-2450000], M sin i [MJ]; no K is tabulated.

Selection: M sin i <= 80 MJ and (P <= 4000 d or a <= 5 AU). P is converted from years (x 365.25), Tp by adding
2450000 (time system as in the paper). fit_type RV+astrometry; omega is in orvara's convention, and i has a
prograde/retrograde ambiguity (both values are in the notes).

Lookups, cached in this directory (used offline; redone over the network only with --refresh):
  simbad_cache.csv                  CDS Sesame (SIMBAD) answers for the host-name candidates (main_id, otype, ids)
  nea_pscomppars_cache.csv          NASA Exoplanet Archive pscomppars rows reachable by host-ID match
  literature_companions_cache.csv   CLS I + Kiefer+2019 rows of ../literature_companions.csv, with SIMBAD ids
The caches were seeded from the original pass's own query results of 2026-10-24 (its Sesame cache, its pscomppars
download, and the 239-row literature_companions.csv of that day), which reproduce the table exactly. A --refresh
queries today's SIMBAD and NEA and reads today's ../literature_companions.csv, so the notes can change.
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

# Raw inputs (path relative to this directory -> (arXiv source tarball URL, member))
SOURCES = {'src/RAA-2022-0418.tex': ('https://arxiv.org/e-print/2303.12409', 'RAA-2022-0418.tex')}
OUT = 'xiao_2023_orbits.csv'

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


def v(s):
    """(value, symmetric error) of a LaTeX cell like '$12.3_{-0.4}^{+0.5}$'; error = larger of |-| and |+|."""
    s = s.replace('$', '').strip()
    m = re.match(r'\{?([-0-9.]+)\}?(?:_\{?([-0-9.]+)\}?\^\{?\+?([-0-9.]+)\}?)?', s)
    if not m or not re.match(r"^-?[0-9]*\.?[0-9]+$", m.group(1)):
        return None, None
    val, err = float(m.group(1)), None
    if m.group(2):
        err = max(abs(float(m.group(2))), abs(float(m.group(3))))
    return val, err


def err_mean(s):
    """Mean of the two printed sides of a '{x}_{-lo}^{+hi}' cell; '' unless both are > 0."""
    m = re.search(r'_\{?-?([0-9.]+)\}?\^\{?\+?([0-9.]+)\}?', s)
    if not m or float(m.group(1)) <= 0 or float(m.group(2)) <= 0:
        return ''
    return fmt((float(m.group(1)) + float(m.group(2))) / 2, 4)


def build_rows():
    with open(here('src/RAA-2022-0418.tex'), encoding='utf-8') as fh:
        tex = fh.read()
    i0 = tex.index('label{Tab:orbparas}')
    i1 = tex.index('end{tabular}', tex.index('label{Tabcc}'))
    rows = []
    for line in tex[i0:i1].splitlines():
        if (not re.match(r'^(HD|GJ|HIP|BD|HR|[a-zA-Z0-9])', line) or '&' not in line or line.startswith('Name')
                or line.startswith('\\')):
            continue
        c = [x.strip() for x in line.rstrip('\\').split('&')]
        if len(c) < 12:
            continue
        # columns: name, Mp [MJ], i<90 [deg], i>90, a [AU], e, P [yr], Omega, omega [deg], a_rel [mas],
        # Tp-2450000 [d], M sin i [MJ]
        name = c[0].replace('\\,', ' ').strip()
        Mp, a, e = v(c[1])[0], v(c[4])[0], v(c[5])[0]
        Pyr, ePyr = v(c[6])
        om, Tp, msini, inc = v(c[8])[0], v(c[10])[0], v(c[11])[0], v(c[2])[0]
        if msini is None or msini > 80:
            continue
        P = Pyr * 365.25 if Pyr else None
        if not ((P is not None and P <= 4000) or (a is not None and a <= 5)):
            continue
        host = re.sub(r'\s+[b-z]$', '', name)
        host = re.sub(r'\s+(B|C)$', '', host) if re.search(r'\s+[BC]$', host) else host
        rows.append(dict(
            host=host, companion=name, P_d=fmt(P, 6), e_P_d=fmt(ePyr * 365.25, 3) if ePyr else '', e=fmt(e, 4),
            omega_deg=fmt(om, 5), Tp_bjd=fmt(Tp + 2450000, 9) if Tp is not None else '', msini_mj=fmt(msini, 4),
            a_au=fmt(a, 4), fit_type='RV+astrometry', e_msini_mj=err_mean(c[11]), e_e=err_mean(c[5]),
            e_a_au=err_mean(c[4]),
            table_ref='Xiao et al. 2023, RAA (accepted; arXiv:2303.12409), Table "Posteriors of RV Companions" '
                      '(LaTeX source)',
            notes='orvara joint fit of literature RVs + HGCA (Hipparcos-Gaia EDR3); true mass Mp=%s MJ, i=%s or %s '
                  'deg (prograde/retrograde ambiguity); P converted from yr; Tp = (Tp-2450000)+2450000, time system '
                  'as in paper; omega in orvara convention; no K tabulated'
                  % (fmt(Mp, 4), fmt(inc, 4), fmt(180 - inc, 4) if inc is not None else '')))
    return rows


def main():
    get_inputs(SOURCES)
    rows = build_rows()
    alt = {}  # no alternate host names in this source
    cache, nea, lit = lookups(rows, alt)
    rows = finalize(rows, cache, nea, lit, alt)
    write_csv(here(OUT), rows, COLS)
    print('wrote %d rows to %s' % (len(rows), OUT))


if __name__ == '__main__':
    main()
