/* Keplerian radial-velocity curve for companion_search.py.
 *
 * rv(t) = K [cos(nu + w) + e cos w],  M = 2 pi (t - tref)/P + M0,  M = E - e sin E.
 *
 * For each time the eccentric anomaly E comes from Halley's method on Kepler's equation,
 * started from a series guess that is accurate enough (error < ~0.1 rad even at e = 0.99)
 * for the cubic convergence to reach 1e-14 in 2-4 steps. cos(nu) and sin(nu) follow from
 * E without an arctangent:
 *     cos nu = (cos E - e)/(1 - e cos E),  sin nu = sqrt(1 - e^2) sin E/(1 - e cos E),
 * so rv = K [cos nu cos w - sin nu sin w + e cos w]. One pass over the times, no
 * temporary arrays. Build: see kepler_c.py (cc -O3 -shared -fPIC). */
#include <math.h>

static double solve_kepler(double M, double e)
{
    /* M in [0, 2 pi). Fold to [0, pi] using the symmetry E(2 pi - M) = 2 pi - E(M). */
    int flip = 0;
    if (M > M_PI) { M = 2.0*M_PI - M; flip = 1; }
    double E;
    if (e < 0.8) {
        E = M + e*sin(M);
    } else {
        /* near-parabolic: start at pi, or at a cube-root guess when M is small */
        E = (M < 0.5) ? cbrt(6.0*M) : M_PI;
        if (E > M_PI) E = M_PI;
    }
    for (int k = 0; k < 8; k++) {
        double s = sin(E), c = cos(E);
        double f = E - e*s - M;
        double fp = 1.0 - e*c;
        double fpp = e*s;
        /* Halley step: dE = -f / (f' - f f''/(2 f')) */
        double dE = -f/(fp - 0.5*f*fpp/fp);
        E += dE;
        if (fabs(dE) < 1e-14) break;
    }
    return flip ? 2.0*M_PI - E : E;
}

void keplerian_rv(const double *t, int n, double tref, double P, double K, double e,
                  double w, double M0, double *out)
{
    const double twopi = 2.0*M_PI;
    const double cw = cos(w), sw = sin(w), sq = sqrt(1.0 - e*e);
    const double ecw = e*cw;
    const double inv = twopi/P;
    for (int i = 0; i < n; i++) {
        double M = fmod((t[i] - tref)*inv + M0, twopi);
        if (M < 0.0) M += twopi;
        double E = solve_kepler(M, e);
        double sE = sin(E), cE = cos(E);
        double den = 1.0 - e*cE;
        double cnu = (cE - e)/den;
        double snu = sq*sE/den;
        out[i] += K*(cnu*cw - snu*sw + ecw);   /* accumulates: callers sum planets */
    }
}

/* Eccentric anomaly alone, for tests against the pure-Python solver */
void kepler_E(const double *M, int n, double e, double *E)
{
    for (int i = 0; i < n; i++) {
        double m = fmod(M[i], 2.0*M_PI);
        if (m < 0.0) m += 2.0*M_PI;
        E[i] = solve_kepler(m, e);
    }
}
