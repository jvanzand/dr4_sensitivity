"""Empirical false-alarm threshold on a Delta-BIC periodogram: a port of
rvsearch's eFAP.

Source: rvsearch/periodogram.py, `Periodogram.eFAP` (J. B. Ruffio's 2022-02-17
version, after Howard & Fulton 2016; Rosenthal+2021, Sec. V.1). The
periodogram's own values set the noise level: sort them, keep the 50th-95th
percentile (the 95th cut keeps the real peaks out), shift by the lowest kept
value (the median) to give x, histogram x in 10 bins, and fit a straight line to
log10(count) against bin centre. With count(x) = B exp(-A x), A = -slope ln10 is
the decay rate per unit Delta-BIC. The chance that one period exceeds x is the
integral of the normalized decay, exp(-A x); to keep the false-alarm probability
FAP over N periods (union bound) each period may exceed it with probability
FAP/N, so
    threshold = median + ln(N/FAP)/A,    fap_min = N exp(-A (max - median)).
Detection (max > threshold) is the same as fap_min < FAP. B is unused. The only
deviation from rvsearch: if the fit is unusable (fewer than 2 populated bins or
A <= 0) the threshold is infinite (no detection) and `ok` is False, instead of a
meaningless number.

Test mode:  python -I blind_search/efap.py --test
"""

import numpy as np


def efap(dbic, fap=1e-3, lo=0.5, hi=0.95, nbins=10):
    dbic = np.asarray(dbic, float)
    n = len(dbic)
    s = np.sort(dbic)
    crop = s[int(lo*n):int(hi*n)]
    med = crop[0]
    hist, edge = np.histogram(crop - med, bins=nbins)
    cent = (edge[1:] + edge[:-1])/2
    with np.errstate(divide='ignore'):
        loghist = np.log10(hist)
    good = np.isfinite(loghist)
    out = dict(median=med, hist=hist, centres=cent, edges=edge, good=good, crop=(crop[0], crop[-1]),
               n_periods=n, fap=fap, ok=False, A=np.nan, B=np.nan, slope=np.nan, intercept=np.nan,
               threshold=np.inf, fap_min=np.nan)
    if good.sum() < 2:
        return out
    a, b = np.polyfit(cent[good], loghist[good], 1)
    A = -a*np.log(10)
    out.update(slope=a, intercept=b, A=A, B=10**b)
    if A <= 0:
        return out
    out.update(ok=True, threshold=np.log(fap/n)/(-A) + med,
               fap_min=np.exp(-A*(s[-1] - med))*n)
    return out


if __name__ == '__main__':
    import sys
    if '--test' in sys.argv:
        x = np.random.default_rng(0).exponential(2.0, 5000)
        r = efap(x)
        assert r['ok'] and abs(r['A'] - 0.5) < 0.1, r['A']
        print(f"efap.py test passed: A = {r['A']:.3f} (true 0.5), threshold {r['threshold']:.1f}")
