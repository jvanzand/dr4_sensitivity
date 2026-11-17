"""The linear model behind the search: stacked channels, baseline columns,
companion columns.

All selected channels are stacked into one weighted least-squares problem.
  * Baseline columns are channel-specific: astrometry gets the reference-point
    motion (position and proper motion along sin psi and cos psi, about the
    channel's mean epoch) and parallax; each RV channel gets an offset; an
    optional RV trend (order 1, 2) is shared by all RVs; an optional
    astrometric trend (order 1 = acceleration, 2 = + its rate of change, the
    "jerk": Gaia's 7- and 9-parameter solutions) is added per astrometric
    channel.
  * A companion with nonlinear parameters (P, e, M0) adds 4 columns to the
    astrometric rows (A, B, F, G) and 2 to the RV rows (K cos w, K sin w).
    Channels of the same type share the companion's amplitudes. In reality one
    orbit fixes both sets; the search treats them as independent (conservative:
    it pays for the extra parameters in the BIC) and the final Octofitter fit
    imposes the physical link. Parallax is fitted per astrometric channel here
    and shared in the final fit.
  * Each channel has a jitter added in quadrature to its formal errors (`jit`,
    one per channel).

Test mode:  python -I blind_search/model.py --test
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # python -I drops the script dir
import orbit


class Stack:
    """All data channels stacked into one weighted least-squares problem: the
    object every search step works on.

    A *channel* is one dataset (a `data.Channel`): along-scan Gaia astrometry
    (kind 'astro', y in mas) or the radial velocities of one instrument (kind
    'rv', y in m/s). The rows of all channels are concatenated in the order
    given, so every per-point array below has length N = `n`. The model for each
    row is

        y_i = baseline_i + sum over companions of (columns @ amplitudes)_i
              + noise,   noise variance = err_i^2 + jit^2,

    where jit is the jitter of the row's channel (extra white noise added in
    quadrature, in the channel's units). Given the nonlinear parameters (P, e,
    M0) of every companion, all other parameters (baseline coefficients and
    companion amplitudes) enter linearly and are solved exactly by weighted
    least squares; `fit_free_jitter` also fits the jitter. A companion is
    described by a `planets` entry (P, e, M0): period P in days, eccentricity e,
    and M0 = mean anomaly (rad) at the reference time `tref`.

    Astrometric convention: a transit measures only the offset along its scan
    direction,

        y_i = dRA* sin(psi_i) + dDec cos(psi_i) + parallax x f_i + orbit,

    with psi the scan position angle (from north through east) and f the
    along-scan parallax factor. The orbit term of a companion is (B X + G Y)
    sin(psi) + (A X + F Y) cos(psi), with X = cos E - e, Y = sqrt(1 - e^2) sin E
    (E from Kepler's equation) and A, B, F, G the Thiele-Innes constants in mas.
    RV rows get K cos(w) (cos(nu) + e) - K sin(w) sin(nu).

    Attributes (all set in __init__)
      channels   the list of channels.
      n          number of rows N.
      t, y, err  epoch (MJD, days), measurement (mas or m/s) and formal error
                 (same units) of every row.
      chan       index of each row's channel in `channels` (integer array).
      is_astro, is_rv, has_astro, has_rv
                 boolean masks per row, and whether any such rows exist.
      sinpsi, cospsi, fal
                 sin and cos of the scan angle and the along-scan parallax
                 factor of astrometric rows (0 for RV rows).
      tref       mid-baseline epoch (MJD), the reference time of M0.
      span       t_max - t_min (days).
      rv_trend   order (0, 1 or 2) of a polynomial in (t - tref)/span shared by
                 all RV rows (0 if there are no RVs).
      astro_trend
                 order (0, 1 or 2) of the astrometric trend: 1 adds an
                 acceleration (columns sin psi x dt^2, cos psi x dt^2), 2 also
                 its rate of change (dt^3); per astrometric channel, dt as for
                 the proper motion (0 if there is no astrometry).
      trend      the trend order of the model, max(rv_trend, astro_trend);
                 `with_trend` returns a copy at another order.
      baseline   (N, k_base) matrix of the baseline columns; baseline_names
                 their names.
      k_base     number of baseline columns.
      k_comp_astro, k_comp_rv
                 columns one companion adds to astrometric rows (4: A, B, F, G)
                 and RV rows (2: K cos w, K sin w); 0 if absent.

    Baseline columns. Per astrometric channel: position offsets 'ra0' (column
    sin psi) and 'dec0' (cos psi), proper motions 'pmra' (sin psi x dt) and
    'pmdec' (cos psi x dt) with dt in years about the channel's mean epoch, and
    parallax 'plx' (column f), and with astro_trend > 0 'acc_ra', 'acc_dec'
    (sin psi x dt^2, cos psi x dt^2) and with astro_trend = 2 'jerk_ra',
    'jerk_dec' (dt^3). Per RV channel: an 'offset'. If rv_trend > 0, shared
    trend columns are added.

    Channels of the same type share a companion's amplitudes (the search treats
    the astrometric and RV amplitudes of one companion as independent, which is
    conservative; the final Octofitter fit imposes the physical link).
    """

    def __init__(self, channels, rv_trend=0, astro_trend=0):
        """channels: list of data.Channel (astrometry and/or RV). rv_trend: 0, 1
        or 2, order of the polynomial trend shared by all RV rows. astro_trend:
        0, 1 (acceleration) or 2 (+ jerk) for the astrometric channels."""
        self.channels = list(channels)
        cat = lambda f: np.concatenate([f(c) for c in self.channels])
        self.n = sum(len(c) for c in self.channels)
        self.t = cat(lambda c: c.t)
        self.y = cat(lambda c: c.y)
        self.err = cat(lambda c: c.err)
        self.chan = np.concatenate([np.full(len(c), i) for i, c in enumerate(self.channels)]).astype(int)
        self.is_astro = cat(lambda c: np.full(len(c), c.kind == 'astro'))
        self.is_rv = ~self.is_astro
        zero = lambda c: np.zeros(len(c))
        self.sinpsi = cat(lambda c: c.sinpsi if c.kind == 'astro' else zero(c))
        self.cospsi = cat(lambda c: c.cospsi if c.kind == 'astro' else zero(c))
        self.fal = cat(lambda c: c.fal if c.kind == 'astro' else zero(c))
        self.has_astro = bool(self.is_astro.any())
        self.has_rv = bool(self.is_rv.any())
        self.tref = 0.5*(self.t.min() + self.t.max())
        self.span = self.t.max() - self.t.min()
        self.rv_trend = rv_trend if self.has_rv else 0
        self.astro_trend = astro_trend if self.has_astro else 0
        self.baseline, self.baseline_names = self._baseline()
        self.k_base = self.baseline.shape[1]
        self.k_comp_astro = 4 if self.has_astro else 0
        self.k_comp_rv = 2 if self.has_rv else 0

    def _baseline(self):
        """Baseline columns of all channels (see the class docstring). Returns
        (N x k_base matrix, list of column names)."""
        cols, names = [], []
        for i, c in enumerate(self.channels):
            m = (self.chan == i).astype(float)
            if c.kind == 'astro':
                dt = (self.t - c.t.mean())/365.25
                terms = [('ra0', self.sinpsi), ('dec0', self.cospsi),
                         ('pmra', self.sinpsi*dt), ('pmdec', self.cospsi*dt), ('plx', self.fal)]
                if self.astro_trend >= 1:
                    terms += [('acc_ra', self.sinpsi*dt**2), ('acc_dec', self.cospsi*dt**2)]
                if self.astro_trend >= 2:
                    terms += [('jerk_ra', self.sinpsi*dt**3), ('jerk_dec', self.cospsi*dt**3)]
                for nm, col in terms:
                    cols.append(m*col)
                    names.append(f'{c.name}:{nm}')
            else:
                cols.append(m)
                names.append(f'{c.name}:offset')
        if self.rv_trend:
            rv_t = self.is_rv*(self.t - self.tref)/self.span
            for k in range(1, self.rv_trend + 1):
                cols.append(self.is_rv*rv_t**k)
                names.append(f'rv_trend{k}')
        return np.column_stack(cols), names

    @property
    def trend(self):
        """Trend order of the model (0, 1 or 2)."""
        return max(self.rv_trend, self.astro_trend)

    def with_trend(self, order):
        """A copy of this stack (same data, including a `y` set in place) with
        trend order `order` for every data type present."""
        import copy
        s = copy.copy(self)
        s.rv_trend = order if self.has_rv else 0
        s.astro_trend = order if self.has_astro else 0
        s.baseline, s.baseline_names = s._baseline()
        s.k_base = s.baseline.shape[1]
        return s

    def n_params(self, n_comp):
        """Free parameters of a model with `n_comp` companions: baseline + per
        companion 3 nonlinear (P, e, M0) and the linear amplitudes (4
        astrometric, 2 RV); jitter not counted."""
        return self.k_base + n_comp*(3 + self.k_comp_astro + self.k_comp_rv)

    def sigma(self, jit):
        """Total error of every row, sqrt(err^2 + jitter^2); `jit` holds one
        jitter per channel (the channel's units)."""
        return np.sqrt(self.err**2 + np.asarray(jit)[self.chan]**2)

    def companion_columns(self, P, e, M0):
        """(N, k_comp) design columns of one companion with period P (d),
        eccentricity e and mean anomaly M0 (rad) at `tref`. Astrometric rows
        carry 4 columns, the coefficients of which are A, B, F, G (mas); RV rows
        carry 2, with coefficients K cos w and K sin w (m/s); a column is zero
        on rows of the other data type."""
        M = orbit.mean_anomaly(self.t, P, M0, self.tref)
        cols = []
        if self.has_astro:
            cols.append(orbit.astro_columns(M, e, self.sinpsi, self.cospsi)*self.is_astro[:, None])
        if self.has_rv:
            cols.append(orbit.rv_columns(M, e)*self.is_rv[:, None])
        return np.hstack(cols)

    def design(self, planets):
        """Design matrix: the baseline columns followed by the columns of each
        (P, e, M0) in `planets` (the order of `fit`'s coefficients)."""
        return np.hstack([self.baseline] + [self.companion_columns(*p) for p in planets])

    def fit(self, planets, jit):
        """Weighted least squares for all linear parameters with the jitter
        FIXED at `jit` (one value per channel; zeros = formal errors only).
        Returns coef (baseline coefficients, then each companion's amplitudes,
        as in `design`), residuals (data - model, native units) and chi2 = sum
        (resid/sigma)^2."""
        sig = self.sigma(jit)
        X = self.design(planets)
        coef = np.linalg.lstsq(X/sig[:, None], self.y/sig, rcond=1e-12)[0]
        r = self.y - X@coef
        return coef, r, float(np.sum((r/sig)**2))

    def lnL(self, resid, jit):
        """Gaussian log-likelihood -1/2 sum [r^2/s^2 + ln(2 pi s^2)] of
        residuals `resid` with s^2 = err^2 + jit^2."""
        s2 = self.sigma(jit)**2
        return float(-0.5*np.sum(resid**2/s2 + np.log(2*np.pi*s2)))

    def fit_jitter(self, resid, jit0=None, iters=1):
        """Maximum-likelihood jitter of each channel given fixed residuals (>=
        0), by a bounded 1-D minimisation per channel. An older, slower
        alternative to `jitter_fixed_point` (the argument `iters` is unused);
        the search uses `fit_free_jitter`."""
        from scipy.optimize import minimize_scalar
        jit = np.zeros(len(self.channels)) if jit0 is None else np.array(jit0, float)
        for i in range(len(self.channels)):
            m = self.chan == i
            r2, e2 = resid[m]**2, self.err[m]**2
            hi = 3*np.sqrt(np.mean(r2)) + 1e-9
            f = lambda s: 0.5*np.sum(r2/(e2 + s*s) + np.log(e2 + s*s))
            jit[i] = minimize_scalar(f, bounds=(0, hi), method='bounded',
                                     options={'xatol': 1e-6*hi}).x
        return jit

    def jitter_fixed_point(self, resid, s2, iters=6):
        """Maximum-likelihood jitter^2 per channel for fixed residuals, by the
        fixed-point iteration s2 <- sum w^2 (r^2 - err^2) / sum w^2 with w =
        1/(err^2 + s2), floored at 0."""
        e2, nch = self.err**2, len(self.channels)
        for _ in range(iters):
            w2 = 1/(e2 + s2[self.chan])**2
            s2 = np.maximum(np.bincount(self.chan, w2*(resid**2 - e2), nch)/np.bincount(self.chan, w2, nch), 0.0)
        return s2

    def fit_free_jitter(self, X, s2=None, iters=4):
        """Weighted least squares for the design matrix X (N x k) with each
        channel's jitter fitted jointly: `iters` rounds of a linear solve at the
        current jitter followed by the maximum-likelihood jitter update
        (`jitter_fixed_point`). `s2` = starting jitter^2 per channel (default 1%
        of the data variance; warm starts converge faster). This is the
        likelihood used by the periodogram, the refits and the simulated
        threshold. Returns coef, residuals, lnL and the fitted jitter^2 per
        channel."""
        nch = len(self.channels)
        s2 = np.full(nch, 0.01*np.var(self.y)) if s2 is None else np.asarray(s2, float)
        e2 = self.err**2
        for _ in range(iters):
            sig = np.sqrt(e2 + s2[self.chan])
            coef = np.linalg.lstsq(X/sig[:, None], self.y/sig, rcond=1e-12)[0]
            r = self.y - X@coef
            s2 = self.jitter_fixed_point(r, s2)
        v = e2 + s2[self.chan]
        return coef, r, float(-0.5*np.sum(r*r/v + np.log(2*np.pi*v))), s2


if __name__ == '__main__':
    import sys
    if '--test' in sys.argv:
        import data
        g = data.load_gaia_dr4(os.path.join(data.ROOT, 'data_conversion/GaiaDR4_data_Gaia4.csv'))
        s = Stack([g])
        assert s.baseline.shape == (s.n, 5) and s.companion_columns(100., 0.1, 0.).shape == (s.n, 4)
        print('model.py test passed (full tests: tests/test_blind_search.py)')
