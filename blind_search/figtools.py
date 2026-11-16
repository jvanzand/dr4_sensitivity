"""Figure helpers: legends that never cover plotted data.

`legend_clear(ax, **kw)` tries the standard in-axes legend positions in turn
and keeps the first whose box (in pixels, padded) contains none of the axes'
data: every marker, error bar and line vertex, with lines densified so a long
straight segment (e.g. a horizontal reference line) also counts. If every
position overlaps, it raises the top of the y range to make a free band; as a
last resort the legend goes outside the axes (right, or below) where it
overlaps no other axes.
`check_legends(fig)` reports any legend that still overlaps data; call it just
before saving.

Call both after all plotting and axis limits are final, and after every
panel of the figure exists (the check works in pixel coordinates of the
current layout, and an outside legend must avoid the other panels).

Test mode:  python -I blind_search/figtools.py --test
"""

import numpy as np

LOCS = ('upper left', 'upper right', 'lower left', 'lower right', 'upper center',
        'lower center', 'center left', 'center right')


def _densify(xy, step=3.0):
    """Points every `step` pixels along a polyline (display coordinates)."""
    out = [xy[:1]]
    for a, b in zip(xy[:-1], xy[1:]):
        n = int(np.hypot(*(b - a))//step)
        if n > 0:
            out.append(a + (b - a)*np.linspace(0, 1, n + 2)[1:, None])
        else:
            out.append(b[None])
    return np.concatenate(out)


def data_points(ax):
    """Display coordinates (N, 2) of everything drawn in `ax` with data:
    lines (densified), markers, error bars and filled regions."""
    pts = []
    for ln in ax.get_lines():
        xy = np.asarray(ln.get_xydata(), float)
        if len(xy) == 0 or not ln.get_visible():
            continue
        d = ln.get_transform().transform(xy)
        d = d[np.all(np.isfinite(d), axis=1)]
        if len(d):
            pts.append(_densify(d) if ln.get_linestyle() not in ('None', ' ', '') else d)
    for col in ax.collections:
        if not col.get_visible():
            continue
        offs = col.get_offsets()
        if len(offs):
            pts.append(col.get_offset_transform().transform(np.asarray(offs, float)))
        for p in col.get_paths():
            v = p.vertices
            if len(v):
                d = col.get_transform().transform(v)
                pts.append(_densify(d[np.all(np.isfinite(d), axis=1)]))
    for pa in ax.patches:
        v = pa.get_path().vertices
        if len(v):
            pts.append(_densify(pa.get_transform().transform(v)))
    if not pts:
        return np.zeros((0, 2))
    p = np.concatenate([q for q in pts if len(q)])
    box = ax.get_window_extent()                      # only what is inside the axes
    m = (p[:, 0] >= box.x0) & (p[:, 0] <= box.x1) & (p[:, 1] >= box.y0) & (p[:, 1] <= box.y1)
    return p[m]


def _n_inside(leg, pts, renderer, pad=4):
    bb = leg.get_window_extent(renderer).padded(pad)
    return int(np.sum((pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1)
                      & (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1)))


def _collides(leg, ax, renderer):
    """Whether a legend box overlaps another axes or leaves the figure."""
    bb = leg.get_window_extent(renderer)
    fb = ax.figure.bbox
    if bb.x0 < fb.x0 or bb.y0 < fb.y0 or bb.x1 > fb.x1 or bb.y1 > fb.y1:
        return True
    return any(o is not ax and o.get_visible() and bb.overlaps(o.get_window_extent())
               for o in ax.figure.axes)


def legend_clear(ax, locs=LOCS, grow=3, **kw):
    """Place the legend at the first position in `locs` where it covers no
    data. If none is free, raise the top of the y range by 15% of the range
    (up to `grow` times; not for fixed-aspect axes) and try the upper positions
    again. Last resort:
    outside the axes, right or below, wherever it hits no other axes."""
    fig = ax.figure
    fig.canvas.draw()                                  # final limits and layout
    renderer = fig.canvas.get_renderer()
    pts = data_points(ax)
    for loc in locs:
        leg = ax.legend(loc=loc, **kw)
        if _n_inside(leg, pts, renderer) == 0:
            return leg
    ## Equal-aspect axes re-adjust their limits, which would crop data: no growth.
    for _ in range(grow if ax.get_aspect() == 'auto' else 0):
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.15*(hi - lo))
        fig.canvas.draw()
        pts = data_points(ax)
        for loc in ('upper left', 'upper right', 'upper center'):
            leg = ax.legend(loc=loc, **kw)
            if _n_inside(leg, pts, renderer) == 0:
                return leg
    for anchor in ((1.01, 1), (0, -0.18)):
        leg = ax.legend(loc='upper left', bbox_to_anchor=anchor, **kw)
        fig.canvas.draw()
        if not _collides(leg, ax, renderer):
            return leg
    raise RuntimeError(f'no free place for the legend of {ax.get_title()!r}')


def check_legends(fig, raise_error=True):
    """Every legend in `fig` that covers data points of its own axes. Returns
    a list of (axes title, number of covered points); raises if any and
    `raise_error`."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bad = []
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None or not ax.axison:
            continue
        ax.apply_aspect()
        n = _n_inside(leg, data_points(ax), renderer, pad=0)
        if n:
            bad.append((ax.get_title() or ax.get_ylabel(), n))
    if bad and raise_error:
        raise RuntimeError(f'legend covers data: {bad}')
    return bad


if __name__ == '__main__':
    import sys
    if '--test' in sys.argv:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        x = np.linspace(0, 10, 200)
        ax.plot(x, 10 - x, 'o', label='falling points')     # data limits 0-10: the check
        ax.axhline(9.5, color='k', label='a line across the top')   # must autoscale first
        ax.legend(loc='upper right')
        assert check_legends(fig, raise_error=False), 'overlap not detected'
        legend_clear(ax)
        assert not check_legends(fig, raise_error=False), 'legend_clear left an overlap'
        print('figtools.py test passed')
