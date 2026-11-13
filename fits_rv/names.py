"""Display names for systems (sample/system_names.py writes sample/system_names.csv).

Result folders and figure titles use a readable name (HD191939, TOI1694), while the
Gaia DR3 source_id stays the key inside every file. Falls back to the source_id when the
table is missing or lacks the system."""

import os

import pandas as pd

TABLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     'sample', 'system_names.csv')
_cache = {}


def _table():
    if 'df' not in _cache:
        _cache['df'] = (pd.read_csv(TABLE, dtype=str) if os.path.exists(TABLE)
                        else pd.DataFrame(columns=['system_id', 'display_name']))
    return _cache['df']


def display_name(system_id):
    t = _table()
    hit = t.display_name[t.system_id == str(system_id)]
    return hit.iloc[0] if len(hit) else str(system_id)


def system_id(name_or_id):
    """The source_id for a display name (or the argument itself if it is already an id)."""
    t = _table()
    hit = t.system_id[t.display_name == str(name_or_id)]
    return hit.iloc[0] if len(hit) else str(name_or_id)
