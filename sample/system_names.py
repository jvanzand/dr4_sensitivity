"""Display names for the census systems, for figures, tables, and result folders.

Systems are keyed by the host's Gaia DR3 source_id, which is exact but unreadable. For
display, each system gets a recognizable name with spaces and special characters removed
(HD 191939 -> HD191939, TOI-1694 -> TOI1694). The key stays the source_id everywhere it
matters (RV table names, the catalog, fit inputs).

Two ways to pick the name (`--basis`):
  hierarchy  the first name in HIERARCHY that SIMBAD lists for the star, else SIMBAD's main
             identifier (the default; HD, GJ, TOI, Kepler, K2, WASP, HAT-P, HIP cover 88%
             of the census systems, and the rest fall back to the main id)
  simbad     SIMBAD's preferred (main) identifier alone, minus catalog-type prefixes such
             as "* " and "V* ". Its drawback: TOI-1736 is "BD+68 202" and 51 Peg is "51Peg",
             so a TOI host may get an unfamiliar catalog name.

SIMBAD identifiers are fetched once for every host with a Gaia DR3 ID and cached in
sample/raw/simbad_aliases.csv; later runs are offline. Hosts without a Gaia ID use their
census host name.

Output: sample/system_names.csv (system_id, display_name, name_basis, simbad_main_id,
census_hostname), and a coverage table on stdout: how many systems have an identifier in
each catalog, which is what decides a sensible hierarchy.

Usage (project root, gaia_dr4_sensitivity environment):
    python -I sample/system_names.py                     # fetch if needed; hierarchy names
    python -I sample/system_names.py --basis simbad      # SIMBAD main id only
"""

import argparse
import importlib.util
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CENSUS_FILE = os.path.join(HERE, 'census_table.csv')
CACHE = os.path.join(HERE, 'raw', 'simbad_aliases.csv')
OUT = os.path.join(HERE, 'system_names.csv')

## Catalogs whose identifiers are counted, and (in this order) preferred under
## --basis hierarchy
CATALOGS = {
    'HD': r'^HD \d',
    'HIP': r'^HIP \d',
    'GJ': r'^(GJ|Gl|NN) \d',
    'HR': r'^HR \d',
    'TOI': r'^TOI-\d',
    'Kepler': r'^Kepler-\d',
    'KOI': r'^KOI-\d',
    'K2': r'^K2-\d',
    'WASP': r'^WASP-\d',
    'HAT-P': r'^HAT-P-\d',
    'TIC': r'^TIC \d',
    'TYC': r'^TYC ',
    '2MASS': r'^2MASS J',
    'Gaia DR3': r'^Gaia DR3 \d',
}
## Order decided 2026-11-06 from the coverage table (cumulative: HD 36%, +GJ 41%, +TOI 81%,
## +Kepler 84%, +K2 86%, +WASP 86%, +HIP 88%); TIC, 2MASS, and Gaia DR3 cover nearly every
## star but are not recognizable names
HIERARCHY = ['HD', 'GJ', 'TOI', 'Kepler', 'K2', 'WASP', 'HAT-P', 'HIP']

## SIMBAD main identifiers carry object-type prefixes ("* 51 Peg", "V* AB Dor",
## "NAME Barnard's star"); they aren't part of the name
PREFIX = re.compile(r'^(\*\*|\*|V\*|NAME|EM\*|PM\*|Cl\*|\*iot|LTT)\s+', re.I)


def clean(name):
    """Remove the type prefix, then every character that is not a letter or digit."""
    name = re.sub(r'\s+', ' ', str(name)).strip()
    name = PREFIX.sub('', name)
    return re.sub(r'[^A-Za-z0-9]', '', name)


def load_fetch():
    """build_targets.fetch (the SIMBAD identifier query), loaded by path."""
    spec = importlib.util.spec_from_file_location('build_targets',
                                                  os.path.join(ROOT, 'rv_data', 'build_targets.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.fetch


def load_aliases(source_ids):
    """SIMBAD main id and identifiers for each Gaia DR3 ID, from the cache, fetching the
    IDs it lacks."""
    cache = (pd.read_csv(CACHE, dtype=str) if os.path.exists(CACHE)
             else pd.DataFrame(columns=['gaia_name', 'main_id', 'alias']))
    have = set(cache.gaia_name.str.replace('Gaia DR3 ', '', regex=False))
    missing = sorted(set(source_ids) - have)
    if missing:
        print(f'Fetching SIMBAD identifiers for {len(missing)} hosts')
        new = load_fetch()(missing).astype(str)[['gaia_name', 'main_id', 'alias']]
        cache = pd.concat([cache, new], ignore_index=True)
        cache.to_csv(CACHE, index=False)
    cache['alias'] = cache.alias.str.replace(r'\s+', ' ', regex=True).str.strip()
    cache['source_id'] = cache.gaia_name.str.replace('Gaia DR3 ', '', regex=False)
    return cache


def build(basis='hierarchy'):
    c = pd.read_csv(CENSUS_FILE, dtype={'source_id': str}, low_memory=False)
    c['system_id'] = c.source_id.fillna('host:' + c.hostname)
    systems = c.drop_duplicates('system_id')[['system_id', 'source_id', 'hostname']]
    al = load_aliases(systems.source_id.dropna().tolist())
    by_star = {sid: g for sid, g in al.groupby('source_id')}
    rows, has = [], {k: 0 for k in CATALOGS}
    n_with_simbad = 0
    for _, r in systems.iterrows():
        g = by_star.get(r.source_id)
        if g is None:
            rows.append({'system_id': r.system_id, 'display_name': clean(r.hostname),
                         'name_basis': 'census host name', 'simbad_main_id': '',
                         'census_hostname': r.hostname})
            continue
        n_with_simbad += 1
        main = g.main_id.iloc[0]
        names = g.alias.tolist()
        for k, pat in CATALOGS.items():
            has[k] += any(re.match(pat, a) for a in names)
        chosen, how = clean(main), 'SIMBAD main id'
        if basis == 'hierarchy':
            for k in HIERARCHY:
                hit = [a for a in names if re.match(CATALOGS[k], a)]
                if hit:
                    chosen, how = clean(sorted(hit)[0]), f'{k} (hierarchy)'
                    break
        rows.append({'system_id': r.system_id, 'display_name': chosen, 'name_basis': how,
                     'simbad_main_id': main, 'census_hostname': r.hostname})
    out = pd.DataFrame(rows)
    ## Names must be unique: a duplicate gets its census host name appended
    dup = out.display_name.duplicated(keep=False)
    out.loc[dup, 'display_name'] = (out.display_name[dup] + '_'
                                    + out.census_hostname[dup].map(clean))
    return out, has, n_with_simbad, len(systems)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--basis', choices=['simbad', 'hierarchy'], default='hierarchy')
    a = ap.parse_args()
    out, has, n_simbad, n = build(a.basis)
    out.to_csv(OUT, index=False)
    print(f'{n} systems, {n_simbad} with SIMBAD identifiers; names written to {OUT}')
    print('Systems with an identifier in each catalog (SIMBAD):')
    for k, v in has.items():
        print(f'  {k:9s} {v:5d}  ({v/n_simbad:.0%})')
    print('Name basis:'); print(out.name_basis.value_counts().to_string())
    print(out.sample(12, random_state=0)[['system_id', 'display_name', 'simbad_main_id']]
          .to_string(index=False))
