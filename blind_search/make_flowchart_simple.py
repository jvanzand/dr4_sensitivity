"""Draw the simple flowchart of the blind search (blind_search/search_flowchart.png).

Three boxes (Steps A, B, C), the stop/continue branches and the loop back to
Step A, in plain words for the README. The detailed version with the
mathematics and background is drawn by make_flowchart.py.

Usage: python -I blind_search/make_flowchart_simple.py [--test]
(--test writes /tmp/search_flowchart_simple_test.png)
"""

import os
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

W = 10.6                   # canvas width, inches
FS = 14.5                  # body font size
LINE = 0.325               # line height, inches
WRAP = 66                  # characters per line inside a wide box

STEP_A = ('Step A: compute the $\\Delta$BIC periodogram', [
    'Fit the current model (any detected planets + center-of-mass motion).',
    'For every period in a fine grid, fit the current model plus a candidate '
    'planet at that period.',
    'Calculate $\\Delta$BIC between the current model and the candidate '
    'model.',
])
STEP_B = ('Step B: set the detection threshold by simulation', [
    'Start with a modified version of the supplied data: the original epochs '
    'and errors, but all values set to 0 (a no-planet model).',
    'Add the signals of the current model (the detected planets and '
    'center-of-mass motion).',
    'Draw M = 50 noise-sampled data sets. Each contains all known signals '
    'plus white noise (the formal errors and the fitted jitter, added in '
    'quadrature).',
    'Run each noise-sampled data set through Step A and keep only its '
    'highest peak. Each peak is the most significant ($\\Delta$BIC) signal '
    'that emerges purely from noise.',
    'Fit an extreme-value (Gumbel) distribution to the peak $\\Delta$BIC '
    'values and calculate the 1% '
    'and 0.1% thresholds: the $\\Delta$BIC values that noise alone exceeds in '
    'only 1% or 0.1% of cases.',
])
STEP_C = ('Step C: decide', [
    'Compare the highest peak of the real periodogram (Step A) with the '
    'threshold (Step B; 0.1% by default).',
])
NO = ('No: stop', [
    'Report the number of detected companions.',
    'Provide starting parameters for an in-depth MCMC search.',
])
YES = ('Yes: add the planet', [
    'Add the new planet to the model.',
    'Refit all companions simultaneously, with jitter free.',
    'Return to Step A.',
])


def wrapped(items, width):
    """Bullet items -> list of lines (bullet on the first line, hanging indent)."""
    out = []
    for it in items:
        ls = textwrap.wrap(it, width=width, break_long_words=False)
        out.append([('• ' if i == 0 else '   ') + l for i, l in enumerate(ls)])
    return out


def box_height(content, width_chars):
    title, items = content
    n = sum(len(b) for b in wrapped(items, width_chars))
    return 0.75 + n*LINE + 0.12*len(items) + 0.1


def draw_box(ax, x0, ytop, w, h, content, width_chars, face, edge):
    from matplotlib.patches import FancyBboxPatch
    title, items = content
    ax.add_patch(FancyBboxPatch((x0, ytop - h), w, h, boxstyle='round,pad=0.02,rounding_size=0.12',
                                fc=face, ec=edge, lw=2.2))
    y = ytop - 0.45
    ax.text(x0 + 0.25, y, title, fontsize=FS + 1, fontweight='bold', va='center')
    y -= 0.45
    for b in wrapped(items, width_chars):
        for k, line in enumerate(b):
            ax.text(x0 + 0.25, y, line, fontsize=FS, va='center')
            y -= LINE
        y -= 0.12


def arrow(ax, p0, p1, text=None, color='#222222', rad=0.0):
    ax.annotate('', xy=p1, xytext=p0,
                arrowprops=dict(arrowstyle='-|>', lw=2.6, color=color, mutation_scale=26,
                                connectionstyle=f'arc3,rad={rad}'))


def main():
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    matplotlib.rcParams.update({'text.usetex': False, 'font.family': 'serif',
                                'mathtext.fontset': 'cm'})
    import matplotlib.pyplot as plt

    wide = 8.3
    narrow = 3.95
    cw_wide, cw_narrow = WRAP, 27
    hA, hB, hC = (box_height(c, cw_wide) for c in (STEP_A, STEP_B, STEP_C))
    hC += 0.85                 # room for the decision question
    hN, hY = box_height(NO, cw_narrow), box_height(YES, cw_narrow)
    gap = 0.75
    H = 1.0 + hA + hB + hC + max(hN, hY) + 3*gap + 0.4
    top = H - 1.0
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis('off')
    ax.text(W/2, H - 0.45, 'Blind search for companions: one pass of the loop', fontsize=FS + 4,
            fontweight='bold', ha='center', va='center')
    x0 = 0.5
    yA = top
    yB = yA - hA - gap
    yC = yB - hB - gap
    yD = yC - hC - gap
    draw_box(ax, x0, yA, wide, hA, STEP_A, cw_wide, '#fdf0cf', '#b88a00')
    draw_box(ax, x0, yB, wide, hB, STEP_B, cw_wide, '#fdf0cf', '#b88a00')
    draw_box(ax, x0, yC, wide, hC, STEP_C, cw_wide, '#ffd9d9', '#b00020')
    xm = x0 + wide/2
    arrow(ax, (xm, yA - hA - 0.03), (xm, yB + 0.03))
    ax.text(xm + 0.2, yA - hA - gap/2, 'periodogram', fontsize=FS - 1.5, va='center', color='#444444')
    arrow(ax, (xm, yB - hB - 0.03), (xm, yC + 0.03))
    ax.text(xm + 0.2, yB - hB - gap/2, 'threshold', fontsize=FS - 1.5, va='center', color='#444444')
    # decision label inside Step C
    ax.text(x0 + wide/2, yC - hC + 0.45, 'Is the highest peak above the threshold?', fontsize=FS + 1,
            fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#b00020', lw=2))
    # branches
    xn = x0 + 0.05
    xy = x0 + wide - narrow
    draw_box(ax, xn, yD, narrow, hN, NO, cw_narrow, '#d8f0d8', '#2e7d32')
    draw_box(ax, xy, yD, narrow, hY, YES, cw_narrow, '#fdf0cf', '#b88a00')
    arrow(ax, (xn + narrow/2, yC - hC - 0.03), (xn + narrow/2, yD + 0.03))
    ax.text(xn + narrow/2 - 0.2, yC - hC - gap/2, 'no', fontsize=FS + 1, fontweight='bold', ha='right', va='center')
    arrow(ax, (xy + narrow/2, yC - hC - 0.03), (xy + narrow/2, yD + 0.03))
    ax.text(xy + narrow/2 + 0.2, yC - hC - gap/2, 'yes', fontsize=FS + 1, fontweight='bold', ha='left', va='center')
    # loop back to Step A
    xr = x0 + wide + 0.55
    ax.plot([xy + narrow, xr], [yD - hY/2, yD - hY/2], color='#222222', lw=2.6, solid_capstyle='butt')
    ax.plot([xr, xr], [yD - hY/2, yA - hA/2], color='#222222', lw=2.6, solid_capstyle='butt')
    arrow(ax, (xr, yA - hA/2), (x0 + wide + 0.04, yA - hA/2))
    ax.text(xr + 0.12, (yD - hY/2 + yA - hA/2)/2, 'next iteration', fontsize=FS - 1.5, rotation=90,
            va='center', ha='left', color='#444444')
    out = '/tmp/search_flowchart_simple_test.png' if '--test' in sys.argv else os.path.join(HERE, 'search_flowchart.png')
    fig.savefig(out, dpi=110)
    print('wrote', out)


if __name__ == '__main__':
    main()
