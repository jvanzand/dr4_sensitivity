"""Check that transcribed orbit values appear verbatim (within rounding) in each source's raw files."""
import glob, os, re, sys
import numpy as np, pandas as pd

RAW = sys.argv[1]
SKIP = {'cls_2021', 'kiefer_2019'}
num = re.compile(r'(?<![\w.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?')

def numbers(d, own):
    vals = []
    for f in glob.glob(os.path.join(d, '**', '*'), recursive=True):
        if not os.path.isfile(f) or os.path.abspath(f) == own: continue
        if not re.search(r'\.(tex|dat|txt|csv|md|tsv|bbl)$|ReadMe$', f): continue
        txt = open(f, errors='ignore').read()
        txt = re.sub(r'(\d)\\[,!;:]\s*(\d)', r'\1\2', txt)   # 3\,508 -> 3508
        txt = re.sub(r'(\d)\{\}(\d)', r'\1\2', txt)
        vals += [float(x) for x in num.findall(txt)]
    return np.array(sorted(set(vals)))

def found(v, arr, rel=2e-3, ab=0):
    if pd.isna(v) or len(arr) == 0: return np.nan
    return bool(np.any(np.abs(arr - v) <= max(rel*abs(v), ab)))

rows = []
for path in sorted(glob.glob(os.path.join(RAW, '*', '*_orbits.csv'))):
    d = os.path.dirname(path); sid = os.path.basename(d)
    if sid in SKIP or os.path.basename(path) != f'{sid}_orbits.csv': continue
    arr = numbers(d, os.path.abspath(path))
    t = pd.read_csv(path)
    for _, r in t.iterrows():
        rows.append(dict(src=sid, host=r.host, comp=r.companion,
                         P=found(r.P_d, arr), K=found(r.K_ms, arr), e=found(r.e, arr, ab=1.5e-3),
                         msini=found(r.msini_mj, arr), Pv=r.P_d, msv=r.msini_mj))
a = pd.DataFrame(rows)
summ = a.groupby('src').agg(n=('P', 'size'), P_ok=('P', lambda s: f'{int(s.sum())}/{s.notna().sum()}'),
                            K_ok=('K', lambda s: f'{int(s.sum())}/{s.notna().sum()}'),
                            e_ok=('e', lambda s: f'{int(s.sum())}/{s.notna().sum()}'),
                            msini_ok=('msini', lambda s: f'{int(s.sum())}/{s.notna().sum()}'))
print(summ.to_string())
a.to_csv(os.path.join(os.path.dirname(__file__), 'audit_rows.csv'), index=False)
