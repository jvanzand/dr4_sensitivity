"""Orbit geometry for the blind search: Kepler's equation and the linear
columns.

For a fixed (P, e, M0) both data types are linear in the remaining orbit
constants:
  * along-scan astrometry: eta = (B X + G Y) sin(psi) + (A X + F Y) cos(psi),
    with X = cos E - e, Y = sqrt(1 - e^2) sin E and the Thiele-Innes constants
    A, B, F, G [mas] of the STAR's own orbit about the barycentre (A, F: north;
    B, G: east);
  * RVs: v = (K cos w)(cos nu + e) + (K sin w)(-sin nu).
So omega, Omega, i and the semi-major axis never enter the nonlinear search; the
nonlinear parameters per companion are only (ln P, e, M0). M0 is the mean
anomaly at the reference time `tref` (it stays defined as e -> 0, unlike the
time of periastron).

Test mode:  python -I blind_search/orbit.py --test   (the full tests are
            tests/test_blind_search.py)
"""

import os
import sys

import numpy as np

# The project's compiled solver (fits_rv/kepler.c, Halley's method); numpy fallback if no C compiler.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fits_rv'))
try:
    import kepler_c
    C_SOLVER = kepler_c.AVAILABLE
except ImportError:
    kepler_c, C_SOLVER = None, False

TWO_PI = 2*np.pi


def kepler_E(M, e):
    """Eccentric anomaly from the mean anomaly, 0 <= e < 1. Circular orbits need
    no solve; otherwise the C solver when available (it is used from
    `kepler_c`), else the numpy Newton iteration."""
    if e == 0:
        return np.mod(np.asarray(M, float), TWO_PI)
    if C_SOLVER:
        M = np.asarray(M, float)
        return kepler_c.kepler_E(M.ravel(), e).reshape(M.shape)
    return kepler_E_numpy(M, e)


def kepler_E_numpy(M, e, tol=1e-13, itmax=60):
    """Pure-numpy Newton iteration (reference and fallback)."""
    M = np.mod(np.asarray(M, float), TWO_PI)
    E = M + 0.85*e*np.sign(np.sin(M))              # Danby's starting value
    for _ in range(itmax):
        dE = (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def mean_anomaly(t, P, M0, tref):
    return M0 + TWO_PI*(np.asarray(t, float) - tref)/P


def xy(M, e):
    """X = cos E - e and Y = sqrt(1 - e^2) sin E (the orbit in the orbital
    plane, units of a)."""
    E = kepler_E(M, e)
    return np.cos(E) - e, np.sqrt(1 - e**2)*np.sin(E)


def true_anomaly_terms(M, e):
    """(cos nu, sin nu)."""
    E = kepler_E(M, e)
    d = 1 - e*np.cos(E)
    return (np.cos(E) - e)/d, np.sqrt(1 - e**2)*np.sin(E)/d


def astro_columns(M, e, sinpsi, cospsi):
    """(n, 4) design columns for the coefficients (A, B, F, G)."""
    X, Y = xy(M, e)
    return np.column_stack([X*cospsi, X*sinpsi, Y*cospsi, Y*sinpsi])


def rv_columns(M, e):
    """(n, 2) design columns for the coefficients (K cos w, K sin w)."""
    cn, sn = true_anomaly_terms(M, e)
    return np.column_stack([cn + e, -sn])


def ti_from_elements(a, inc, omega, Omega):
    """Thiele-Innes constants (A, B, F, G) from the star's elements (a in mas,
    angles in rad)."""
    co, so, cO, sO, ci = np.cos(omega), np.sin(omega), np.cos(Omega), np.sin(Omega), np.cos(inc)
    return (a*(co*cO - so*sO*ci), a*(co*sO + so*cO*ci),
            a*(-so*cO - co*sO*ci), a*(-so*sO + co*cO*ci))


def elements_from_ti(A, B, F, G):
    """Inverse of `ti_from_elements`: (a, inc, omega, Omega). The pair (omega,
    Omega) is defined only up to a simultaneous shift of both by pi, which
    astrometry alone cannot break (the RV sign does); degenerate at inc = 0 or
    pi."""
    u = (A*A + B*B + F*F + G*G)/2
    v = A*G - B*F
    a2 = u + np.sqrt(max(u*u - v*v, 0.0))
    inc = np.arccos(np.clip(v/a2, -1, 1))
    wpO = np.arctan2(B - F, A + G)                  # omega + Omega
    wmO = np.arctan2(-(B + F), A - G)               # omega - Omega
    return np.sqrt(a2), inc, 0.5*(wpO + wmO), 0.5*(wpO - wmO)


if __name__ == '__main__':
    import sys
    if '--test' in sys.argv:
        t = np.linspace(0, 100, 7)
        X, Y = xy(mean_anomaly(t, 30.0, 0.4, 50.0), 0.0)
        assert np.allclose(X**2 + Y**2, 1), 'circular orbit is the unit circle'
        print('orbit.py test passed (full tests: tests/test_blind_search.py)')
