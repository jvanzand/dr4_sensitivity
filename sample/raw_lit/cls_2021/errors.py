"""Published 1-sigma uncertainties of the CLS I companions (Rosenthal+2021, ApJS 255, 8; VizieR J/ApJS/255/8).

Reads the cached planets.dat (no network) and writes cls_2021_errors.csv with columns
    pl_name_key, e_P_d, e_K_ms, e_msini_mj, e_e, e_a_au
where pl_name_key is the pl_name that literature.load_cls() builds ('HD 1234 [CLS 1]'); the same
rows are kept (Status K, C, or S).

planets.dat gives, for every parameter X, the max-likelihood value X, the MCMC median XMed, and the
"lower" and "upper" uncertainties e_X and E_X. These are the 16th and 84th percentile BOUNDS of the
posterior, not offsets (checked: e_X <= XMed <= E_X in all 330 rows for Msini, a, Per, K, e). The
1-sigma error written here is half the 68% interval, (E_X - e_X) / 2, which equals the average of the
upper offset (E_X - XMed) and lower offset (XMed - e_X), as for asymmetric errors elsewhere in the
census. Nothing is estimated or propagated; a width <= 0 would be left blank. Columns:
    e_P_d      Per   (d)
    e_K_ms     K     (m/s)
    e_msini_mj Msini (M_Jup; the catalog's Msini is the max-likelihood value, the interval is the
                      posterior's, so for poorly constrained signals (mostly Status S) the error can
                      exceed the value)
    e_e        e
    e_a_au     a     (AU)

Run from the project root:  python -I sample/raw_lit/cls_2021/errors.py [--test]
--test checks the pl_name_key set against literature.load_cls() and the interval convention.
"""
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cls_2021_errors.csv')

# planets.dat ('|' separated) column index -> label: (lower bound, upper bound) pairs per parameter
COLS = {1: 'CPS', 2: 'Status', 3: 'Npl',
        4: 'Msini', 5: 'MsiniMed', 6: 'lo_Msini', 7: 'hi_Msini',
        8: 'a', 9: 'aMed', 10: 'lo_a', 11: 'hi_a',
        12: 'Per', 13: 'PerMed', 14: 'lo_Per', 15: 'hi_Per',
        24: 'K', 25: 'KMed', 26: 'lo_K', 27: 'hi_K',
        28: 'e', 29: 'eMed', 30: 'lo_e', 31: 'hi_e'}
OUT_COLS = {'e_P_d': 'Per', 'e_K_ms': 'K', 'e_msini_mj': 'Msini', 'e_e': 'e', 'e_a_au': 'a'}


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
    pl = pd.read_csv(os.path.join(HERE, 'planets.dat'), sep='|', header=None)
    pl = pl[list(COLS)].rename(columns=COLS)
    pl['CPS'] = pl.CPS.str.strip()
    pl['Status'] = pl.Status.str.strip()
    pl = pl[pl.Status.isin(['K', 'C', 'S'])]
    out = pd.DataFrame({'pl_name_key': pl.CPS.map(simbad_name) + ' [CLS ' + pl.Npl.astype(int).astype(str) + ']'})
    for col, x in OUT_COLS.items():
        width = (pl['hi_' + x] - pl['lo_' + x]) / 2
        out[col] = width.where(width > 0).values
    return pl, out


def main():
    pl, out = build()
    if '--test' in sys.argv[1:]:
        sys.path.insert(0, os.path.join(HERE, '..', '..'))
        import literature
        ref = literature.load_cls().pl_name
        assert list(ref) == list(out.pl_name_key), 'pl_name_key differs from literature.load_cls()'
        for x in OUT_COLS.values():
            assert ((pl['lo_' + x] <= pl[x + 'Med']) & (pl[x + 'Med'] <= pl['hi_' + x])).all(), x
        assert out.pl_name_key.is_unique
        print('test OK: %d keys match load_cls()' % len(out))
        return
    out.to_csv(OUT, index=False, float_format='%.10g')
    print('wrote %s (%d rows)' % (OUT, len(out)))
    print(out.drop(columns='pl_name_key').notna().sum().to_string())


if __name__ == '__main__':
    main()
