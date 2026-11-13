"""Physical companions for injection tests: from (P, e, omega, i, Omega, masses,
parallax) to the linear amplitudes the search model uses (Thiele-Innes constants
of the star's orbit, RV K).

Conventions are those of orbit.py (the star's own orbit; omega is the star's,
i.e. the companion's omega + pi). Whether they agree with Octofitter's angle
conventions is checked by V1 (not yet run), so the sky ORIENTATION of an
injected orbit is unverified; its size, period and shape are not.

GAIA4B: Gaia-4 b as published (Stefansson+2025; values from
sample/stage1_catalog.csv): P, e, omega, K, true mass, host mass and Tp. The
parallax is the one used in Judah's demo fit (13.58 mas, 73.6 pc). The catalog
has no inclination or Omega: the inclination follows from the published K and
true mass (sin i = m sin i / m, taking i > 90 deg as the demo fit does), and
Omega is the median of Judah's demo joint RV + pre-release-DR4 fit (2.96 rad);
both are noted where used.

Test mode:  python -I blind_search/injection.py --test
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orbit

GM_SUN = 1.32712440018e20          # m^3 s^-2 (IAU nominal)
GM_JUP = 1.26686534e17
AU_M = 1.495978707e11
DAY_S = 86400.0
MJD_BJD = 2400000.5

GAIA4B = dict(name='Gaia-4 b', P=571.3, e=0.338, K=368.0, omega_star=np.deg2rad(180.3),
              mass_mj=11.8, mstar=0.644, plx=13.58, tp_mjd=2457386.1 - MJD_BJD,
              Omega=2.96)


def k_from_mass(P, e, inc, mass_mj, mstar):
    """RV semi-amplitude (m/s) of the star."""
    m = mass_mj*GM_JUP
    M = mstar*GM_SUN + m
    return (2*np.pi/(P*DAY_S))**(1/3)*m*np.sin(inc)/M**(2/3)/np.sqrt(1 - e**2)


def inclination_from_k(P, e, K, mass_mj, mstar):
    """Inclination (> 90 deg branch) that gives K for the true mass."""
    m = mass_mj*GM_JUP
    M = mstar*GM_SUN + m
    s = K*np.sqrt(1 - e**2)*M**(2/3)/((2*np.pi/(P*DAY_S))**(1/3)*m)
    if s > 1:
        raise ValueError('K too large for this mass: sin i > 1')
    return np.pi - np.arcsin(s)


def physical_companion(P, e, omega_star, inc, Omega, mass_mj, mstar, plx, tp_mjd, tref):
    """The companion's nonlinear parameters (P, e, M0 at `tref`, MJD) and its
    linear amplitudes: astro = (A, B, F, G) in mas, rv = (K cos w, K sin w) in
    m/s. Also a_star (mas) and K."""
    m_sun = mass_mj*GM_JUP/GM_SUN
    a_au = ((mstar + m_sun)*(P/365.25)**2)**(1/3)           # relative semi-major axis
    a_star_mas = a_au*plx*m_sun/(mstar + m_sun)
    K = k_from_mass(P, e, inc, mass_mj, mstar)
    M0 = np.mod(2*np.pi*(tref - tp_mjd)/P, 2*np.pi)
    return dict(planet=(P, e, M0), astro=np.array(orbit.ti_from_elements(a_star_mas, inc, omega_star, Omega)),
                rv=np.array([K*np.cos(omega_star), K*np.sin(omega_star)]), a_star_mas=a_star_mas,
                a_au=a_au, K=K, inc=inc, distance_pc=1000/plx)


def gaia4b(tref):
    p = GAIA4B
    inc = inclination_from_k(p['P'], p['e'], p['K'], p['mass_mj'], p['mstar'])
    c = physical_companion(p['P'], p['e'], p['omega_star'], inc, p['Omega'], p['mass_mj'],
                           p['mstar'], p['plx'], p['tp_mjd'], tref)
    c['label'] = (f"{p['name']}: P = {p['P']} d, e = {p['e']}, M = {p['mass_mj']} M_J, host {p['mstar']} M_sun, "
                  f"d = {c['distance_pc']:.1f} pc ({p['plx']} mas), a = {c['a_au']:.2f} AU, "
                  f"astrometric semi-major axis {c['a_star_mas']*1e3:.0f} uas, i = {np.degrees(inc):.1f} deg, "
                  f"K = {c['K']:.0f} m/s")
    return c


if __name__ == '__main__':
    if '--test' in sys.argv:
        c = gaia4b(57500.0)
        assert abs(c['K'] - 368.0) < 0.5, c['K']
        print('injection.py test passed:', c['label'])
