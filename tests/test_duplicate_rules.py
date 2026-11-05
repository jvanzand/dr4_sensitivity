"""Tests of the duplicate-RV rules (rv_data/DUPLICATE_RULES.md) on synthetic systems, one
case per rule.

Run from the project root, in the gaia_dr4_sensitivity conda environment:
    python -I tests/test_duplicate_rules.py
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rv_data'))
from merge_rvs import resolve_duplicates  # noqa: E402

DAY = 86400.0
rows = []


def add(dataset, instrument, source, detail, t_s, rv, err=2.0, file_id=None,
        time_quality='exposure', within_rank=0, system='sys1'):
    """One synthetic RV; t_s is seconds after JD 2455000."""
    rows.append({
        'system_id': system, 'target': system, 'rv_source': source,
        'rv_source_detail': detail, 'study': None, 'source_record_name': system,
        'obs_id': f'{dataset}#{len(rows)}', 'dataset': dataset, 'instrument': instrument,
        'bjd_tdb': 2455000 + t_s/DAY, 'bjd_tcb': 2455000 + t_s/DAY, 'time_flag': '',
        'rv_ms': rv, 'rv_err_ms': err, 'isjunk': False, 'quality_flag': '', 's_index': None,
        'program': 'JUMP_LEGACY_PUBLIC', 'access': '', 'exptime_s': 300.0,
        'within_rank': within_rank, 'file_id': file_id, 'time_quality': time_quality})


rng = np.random.default_rng(1)

## Case A (E2, published window, check, remainder): Jump HIRES on 10 nights; a published
## HIRES table copies 8 of them with times off by up to 8 min and a +25 m/s offset, plus
## one pair whose RV is far from the offset
nights = np.arange(10)*DAY
for t in nights:
    add('hires_j', 'hires_j', 'Jump', 'CPS pipeline', t, rng.normal(0, 30))
jump_a = [r for r in rows if r['dataset'] == 'hires_j']
for k, r in enumerate(jump_a[:8]):
    t = (r['bjd_tdb'] - 2455000)*DAY + rng.uniform(-480, 480)
    add('pub:hires', 'HIRES', 'DACE', 'publication', t,
        r['rv_ms'] + 25 + (80 if k == 7 else rng.normal(0, 2)),
        time_quality='published', within_rank=1)

## Case B (E1): DACE's HARPS reduction and a Jump import of the same raw frames, with
## time stamps 5 min apart (outside the exposure window)
for k in range(4):
    t = 20*DAY + k*DAY
    stamp = f'2009-06-{28 + k % 3:02d}T0{k}:00:00'
    add('dace:HARPS03', 'HARPS03', 'DACE', 'DRS-3.5', t, 10.0*k,
        file_id=f'HARPS.{stamp}.123')
    add('harps03-import', 'harps03', 'Jump', 'dace archive', t + 300, 10.0*k + 3,
        file_id=f'r.HARPS.{stamp.replace(":", "-")}.123_ccf_A.fits')

## Case C (E2, exposure window): the same exposure 90 s apart is a copy; 300 s apart isn't
add('dace:HARPN', 'HARPN', 'DACE', 'DRS-2.3', 40*DAY, 5.0)
add('harpn-import', 'harpn', 'Jump', 'tng archive', 40*DAY + 90, 6.0)
add('dace:HARPN', 'HARPN', 'DACE', 'DRS-2.3', 41*DAY, 7.0)
add('harpn-import', 'harpn', 'Jump', 'tng archive', 41*DAY + 300, 8.0)

## Case D: back-to-back exposures in one dataset are all kept
for k in range(3):
    add('apf', 'apf', 'Jump', 'CPS pipeline', 50*DAY + 30*k, 1.0*k)

## Case E: different instruments at the same time are never copies
add('hires_j', 'hires_j', 'Jump', 'CPS pipeline', 60*DAY, 0.0)
add('apf', 'apf', 'Jump', 'CPS pipeline', 60*DAY, 0.0)

## Case F: an exact repeat inside one dataset is dropped
add('dace:SOPHIE', 'SOPHIE', 'DACE', 'archive reduction: HR', 70*DAY, 3.0)
add('dace:SOPHIE', 'SOPHIE', 'DACE', 'archive reduction: HR', 70*DAY, 3.0)

rv = pd.DataFrame(rows)
out, drop, checks, review, n_repeat = resolve_duplicates(rv)
dropped = set(drop.obs_id)
kept = set(out.obs_id)


def ids(dataset):
    return [r['obs_id'] for r in rows if r['dataset'] == dataset]


failures = []


def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        failures.append(name)


## A
pub = ids('pub:hires')
check('A: 7 published copies dropped by time', sum(i in dropped for i in pub) == 7)
check('A: the off-offset pair is rejected and both RVs kept',
      pub[7] in kept and len(review) == 1 and jump_a[7]['obs_id'] in kept)
row = checks[checks.dataset == 'pub:hires'].iloc[0]
check('A: offset measured near +25 m/s', abs(row.offset_ms - 25) < 3)
check('A: Jump HIRES all kept', all(i in kept for i in ids('hires_j')))
## B
check('B: Jump imports of the same raw frames dropped by exposure ID',
      all(i in dropped for i in ids('harps03-import'))
      and set(drop[drop.dataset == 'harps03-import'].evidence) == {'E1 exposure ID'})
## C
c = ids('harpn-import')
check('C: 90 s apart is a copy; 300 s apart is not', c[0] in dropped and c[1] in kept)
## D
check('D: back-to-back exposures in one dataset kept',
      sum(i in kept for i in ids('apf')) == 4)
## E
check('E: different instruments at the same time both kept',
      ids('hires_j')[-1] in kept and ids('apf')[-1] in kept)
## F
check('F: exact repeat dropped', n_repeat == 1
      and (out.dataset == 'dace:SOPHIE').sum() == 1)

print(f'\n{len(failures)} failed' if failures else '\nall passed')
sys.exit(1 if failures else 0)
