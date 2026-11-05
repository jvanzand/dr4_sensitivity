"""Rebuild the census from the cached downloads, in one command.

Steps:
    1. sample/raw_lit/<source>/build.py for every literature source, in precedence order
       (some sources read others' outputs, e.g., sb9_tic reads sb9)
    2. literature.py: VizieR catalogs, merge rules, SIMBAD host matching
    3. census.py: NEA solutions, Gaia DR3 matching, forecasts, cuts
    4. census_plot.py, provenance_plot.py, and lit_sources_plot.py

By default every step reads the cached downloads and lookup caches, so the result depends
only on the files in sample/raw/ and sample/raw_lit/ (dated in sample/raw/SNAPSHOT.json).
The network is used only to look up hosts that the caches don't cover yet (a newly added
source). --download and --refresh are passed only to the source scripts that accept them;
the second-pass scripts (barbato_2026, cui_2026, latham_2002, lee_2023, nidever_2002,
ren_2026, sb9_tic, wittenmyer_2009) take no arguments and ignore both flags. They read raw
files already in their directories (their docstrings name the sources), and latham_2002,
nidever_2002, and sb9_tic query SIMBAD, Gaia, or the TIC only for names missing from their
cache files (delete a cache file to re-query it).

Usage (from the project root):
    python -I sample/build_all.py              # rebuild from cached downloads
    python -I sample/build_all.py --download   # also fetch missing literature raw files
    python -I sample/build_all.py --refresh    # re-query the NEA, Gaia, VizieR, and the source
                                               # scripts' lookup caches (SIMBAD etc.; new snapshot)
"""

import argparse
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from literature import RAW_DIR, SOURCE_TIERS


def run(script, *args, log=None):
    cmd = [sys.executable, '-I', '-W', 'ignore', os.path.join(HERE, script), *args]
    print('>', ' '.join(os.path.relpath(c, os.getcwd()) if os.path.isabs(c) else c
                         for c in cmd[4:]))
    out = subprocess.run(cmd, capture_output=True, text=True)
    if log:
        with open(os.path.join(HERE, log), 'w') as fh:
            fh.write(out.stdout)
    if out.returncode != 0:
        sys.exit(f'{script} failed:\n{out.stderr[-3000:]}')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true',
                        help='let source build scripts fetch missing raw files')
    parser.add_argument('--refresh', action='store_true',
                        help='re-query the NEA, Gaia, VizieR, and source lookup caches')
    args = parser.parse_args()

    ## 1. Literature sources, in precedence order, then any not yet assigned a tier
    ordered = [s for tier in SOURCE_TIERS.values() for s in tier]
    found = sorted(os.path.basename(os.path.dirname(p))
                   for p in glob.glob(os.path.join(RAW_DIR, '*', 'build.py')))
    missing = [s for s in ordered if s not in found]
    if missing:
        sys.exit(f'Sources without a build.py: {missing}')
    # Flags go only to scripts that accept them (the second-pass scripts take none; see above)
    for source in ordered + [s for s in found if s not in ordered]:
        script = os.path.join('raw_lit', source, 'build.py')
        with open(os.path.join(HERE, script)) as fh:
            code = fh.read()
        flags = [f for f, on in (('--download', args.download), ('--refresh', args.refresh))
                 if on and f"'{f}'" in code]
        run(script, *flags)

    ## 2-4. Merge, census, figures
    cached = [] if args.refresh else ['--cached']
    run('literature.py', *cached, log='literature_output.txt')
    run('census.py', *cached, log='census_output.txt')
    run('census_plot.py', log='census_plot_output.txt')
    run('provenance_plot.py', log='provenance_output.txt')
    run('lit_sources_plot.py', log='lit_sources_output.txt')
    run('sky_plot.py', log='sky_output.txt')
    print('Done: sample/census_table.csv and figures rebuilt.')
