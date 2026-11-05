"""Identity table for the systems whose RVs we collect.

Each system is keyed by its host's Gaia DR3 source_id from the census, which every RV
source can be tied to. SIMBAD supplies the position, proper motion, and aliases
(resolved from the Gaia DR3 ID, so the census name's spelling doesn't matter). The
aliases are what the RV queries (Jump, DACE) match on; the position is what each matched
record is checked against.

Outputs (in rv_data/):
    targets.csv            system_id, target, simbad_main_id, ra, dec [deg, ICRS, J2000],
                           pmra, pmdec [mas/yr], aliases (semicolon-separated)
    targets_values.sql     the VALUES list for rv_data/jump_queries.sql (queries A and B)
    simbad_cache.csv       the raw SIMBAD rows, reused offline

Usage (from the project root, in the gaia_dr4_sensitivity conda environment):
    python -I rv_data/build_targets.py                  # the pilot systems
    python -I rv_data/build_targets.py "HD 191939" ...  # census host names
    python -I rv_data/build_targets.py --cached         # from simbad_cache.csv
"""

import argparse
import io
import os

import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
CENSUS = os.path.join(HERE, '..', 'sample', 'census_table.csv')
CACHE = os.path.join(HERE, 'simbad_cache.csv')

PILOT = ['TOI-1694', 'HD 111232', 'HD 191939', 'TOI-1736', 'HD 43691']
SIMBAD_HOSTS = ['https://simbad.cds.unistra.fr', 'https://simbad.harvard.edu']

## Identifier catalogs worth matching on (SIMBAD lists many more)
ALIAS_PREFIXES = ('HD ', 'HIP ', 'HR ', 'GJ ', 'Gl ', 'BD', 'CD', 'CPD', 'TOI-', 'TIC ', 'KOI-',
                  'Kepler-', 'K2-', 'EPIC ', 'WASP-', 'HAT-P-', 'KELT-', 'XO-', 'TrES-',
                  'CoRoT-', 'NGTS-', 'Qatar-', 'LHS ', 'LP ', 'Ross ', 'Wolf ', 'G ', 'NAME ',
                  '2MASS J', 'Gaia DR3 ', 'TYC ')


def simbad_query(adql):
    """Run a SIMBAD TAP query, falling back to the Harvard mirror."""
    for host in SIMBAD_HOSTS:
        try:
            resp = requests.post(f'{host}/simbad/sim-tap/sync', timeout=120,
                                 data={'request': 'doQuery', 'lang': 'adql',
                                       'format': 'csv', 'query': adql})
            resp.raise_for_status()
            return pd.read_csv(io.StringIO(resp.text))
        except requests.RequestException as err:
            print(f'{host}: {err.__class__.__name__}; trying the next host')
    raise RuntimeError('SIMBAD did not answer')


def fetch(gaia_ids):
    """Every SIMBAD identifier, plus position and proper motion, for these Gaia DR3 IDs."""
    rows = []
    for i in range(0, len(gaia_ids), 200):
        ids = ', '.join(f"'Gaia DR3 {g}'" for g in gaia_ids[i:i+200])
        rows.append(simbad_query(
            'SELECT q.id AS gaia_name, b.main_id, b.ra, b.dec, b.pmra, b.pmdec, a.id AS alias '
            'FROM ident AS q JOIN basic AS b ON b.oid = q.oidref '
            'JOIN ident AS a ON a.oidref = q.oidref '
            f'WHERE q.id IN ({ids})'))
    return pd.concat(rows, ignore_index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('hosts', nargs='*', help='census host names (default: the pilot)')
    parser.add_argument('--cached', action='store_true', help='use simbad_cache.csv')
    args = parser.parse_args()
    hosts = args.hosts or PILOT

    census = pd.read_csv(CENSUS, dtype={'source_id': str}).drop_duplicates('hostname')
    sel = census.set_index('hostname').reindex(hosts)
    missing = sel.index[sel.source_id.isna()].tolist()
    if missing:
        raise SystemExit(f'not in the census or without a Gaia DR3 ID: {missing}')
    gaia_ids = sel.source_id.tolist()

    if args.cached:
        simbad = pd.read_csv(CACHE, dtype=str)
    else:
        simbad = fetch(gaia_ids).astype(str)
        simbad.to_csv(CACHE, index=False)
    simbad['source_id'] = simbad.gaia_name.str.replace('Gaia DR3 ', '', regex=False)
    simbad['alias'] = simbad.alias.str.replace(r'\s+', ' ', regex=True).str.strip()

    rows = []
    for host, source_id in zip(hosts, gaia_ids):
        s = simbad[simbad.source_id == source_id]
        if s.empty:
            raise SystemExit(f'{host}: SIMBAD has no object for Gaia DR3 {source_id}')
        aliases = sorted({a for a in s.alias if a.startswith(ALIAS_PREFIXES)} | {host})
        first = s.iloc[0]
        rows.append({'system_id': source_id, 'target': host, 'simbad_main_id': first.main_id,
                     'ra': float(first.ra), 'dec': float(first.dec),
                     'pmra': pd.to_numeric(first.pmra, errors='coerce'),
                     'pmdec': pd.to_numeric(first.pmdec, errors='coerce'),
                     'aliases': '; '.join(aliases)})
    targets = pd.DataFrame(rows)
    targets.to_csv(os.path.join(HERE, 'targets.csv'), index=False)

    values = ',\n'.join("  ('{}', '{}')".format(r.target.replace("'", "''"),
                                                r.aliases.replace("'", "''"))
                        for r in targets.itertuples())
    with open(os.path.join(HERE, 'targets_values.sql'), 'w') as fh:
        fh.write('-- Paste in place of the VALUES list in queries A and B of jump_queries.sql\n'
                 f'  VALUES\n{values}\n')
    print(targets[['system_id', 'target', 'simbad_main_id']].to_string(index=False))


if __name__ == '__main__':
    main()
