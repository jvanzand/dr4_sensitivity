"""Published 1-sigma uncertainties of the Kiefer+2019 companions (A&A 631, A125; VizieR J/A+A/631/A125, Table B5).

Reads the cached tableb5.dat with its ReadMe (no network) and writes kiefer_2019_errors.csv with columns
    pl_name_key, e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au
where pl_name_key is the pl_name that literature.load_kiefer() builds ('HD 217850 [Kiefer+2019]').
The catalog gives "rms uncertainty" columns, copied as published: e_Per (d), e_K (m/s), e_M2sini (M_Jup),
e_e, e_a2 (AU). A blank (masked) or non-positive value stays blank.

Run from the project root:  python -I sample/raw_lit/kiefer_2019/errors.py [--test]
--test checks the pl_name_key list against literature.load_kiefer().
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from astropy.io import ascii

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'kiefer_2019_errors.csv')

# output column -> ReadMe label of the published rms uncertainty
OUT_COLS = {'e_P_d': 'e_Per', 'e_K_ms': 'e_K', 'e_msini_mj': 'e_M2sini', 'e_e': 'e_e', 'e_a_au': 'e_a2'}


def simbad_name(raw):
    """Copy of literature.simbad_name (kept in sync by --test): catalog identifier -> host name."""
    raw = raw.strip()
    m = re.fullmatch(r'(\d+)([a-zA-Z]?)', raw)
    if m:
        return f'HD {m.group(1)}{m.group(2).upper()}'
    m = re.fullmatch(r'(?i)(hd|hip|gl|gj)\s*(\d+)([a-zA-Z]?)', raw)
    if m:
        prefix = {'gl': 'GJ'}.get(m.group(1).lower(), m.group(1).upper())
        return f'{prefix} {m.group(2)}{m.group(3).upper()}'
    m = re.fullmatch(r'(?i)bd([+-]\d{2})(\d+)', raw)
    if m:
        return f'BD{m.group(1)} {int(m.group(2))}'
    return raw


def build():
    t = ascii.read(os.path.join(HERE, 'tableb5.dat'), readme=os.path.join(HERE, 'ReadMe'),
                   format='cds').to_pandas()
    out = pd.DataFrame({'pl_name_key': t.Name.map(simbad_name) + ' [Kiefer+2019]'})
    for col, label in OUT_COLS.items():
        v = pd.to_numeric(t[label], errors='coerce')
        out[col] = v.where(v > 0).values
    return out


def main():
    out = build()
    if '--test' in sys.argv[1:]:
        sys.path.insert(0, os.path.join(HERE, '..', '..'))
        import literature
        assert list(literature.load_kiefer().pl_name) == list(out.pl_name_key), 'keys differ from load_kiefer()'
        print('test OK: %d keys match load_kiefer()' % len(out))
        return
    out.to_csv(OUT, index=False, float_format='%.10g')
    print('wrote %s (%d rows)' % (OUT, len(out)))
    print(out.drop(columns='pl_name_key').notna().sum().to_string())


if __name__ == '__main__':
    main()
