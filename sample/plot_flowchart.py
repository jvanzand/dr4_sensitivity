"""Draft project flowchart: three lanes and four comparisons (PLAN.md, "Overview").

This is a reconstruction drawn from the Overview text in PLAN.md (the original sketch was
not saved in the repository); replace it with the original sketch if it differs.

Usage (project root, gaia_dr4_sensitivity environment):
    python -I sample/plot_flowchart.py                  # writes sample/project_flowchart.png
    python -I sample/plot_flowchart.py --out /tmp/t.png # local test
"""

import argparse
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

W, H = 2.9, 1.15            # box size in axis units
BOXES = {                   # key: (x, y, text, lane colour)
    'lit':   (1, 3, 'Published literature orbits\n(census, one solution per companion)', 'C7'),
    'cens':  (0, 2, 'Census of RV companions\n{stage1}\n{stage2}', 'C0'),   # counts filled from CSVs
    'rv':    (1, 2, 'RV fits\n(+ pre-DR4 astrometry, G23H)\nPLAN Step 3\n'
                    '(RV search run so far on 10 test systems)', 'C0'),
    'syn':   (2, 2, 'Synthetic DR4: our RV orbits\ninjected at our hosts\nPLAN Step 4', 'C0'),
    'synfit': (3, 2, 'Fits to synthetic DR4\n(predicted detection probability,\nposterior widths)', 'C0'),
    'dr4':   (1, 1, 'Real Gaia DR4 epoch astrometry\n(released 2 December 2026)', 'C1'),
    'dr4fit': (3, 1, 'DR4-only fits to real data\nPLAN Step 5', 'C1'),
    'epo':   (1, 0, 'Epochalypse simulated DR4-like\nepoch astrometry\n(Yahalomi+ 2026)', 'C2'),
    'epofit': (3, 0, 'DR4-only fits to the\nEpochalypse data', 'C2'),
    'truth': (2, -1.2, 'Epochalypse injected orbits\n(exact truth)', 'C2'),
}
LANES = [('cens', 'rv'), ('rv', 'syn'), ('syn', 'synfit'), ('dr4', 'dr4fit'), ('epo', 'epofit'),
         ('truth', 'epo')]   # truth -> epo: the injected orbits generate the data
COMPARISONS = [      # label, colour, from, to, arc, description
    ('A', 'C3', 'rv', 'dr4fit', 0.0, 'RV fits vs. real-DR4 fits: does DR4 alone recover the known '
                                      'companions? (main result)'),
    ('B', 'C4', 'rv', 'lit', 0.0, 'RV fits vs. literature fits: do we reproduce published orbits? '
                                  '(consistency of the RV fitting)'),
    ('C', 'C5', 'truth', 'epofit', -0.2, 'Epochalypse truths vs. our fits to its data: '
                                       'consistency of the DR4 fitting'),
    ('D', 'C6', 'synfit', 'dr4fit', 0.0, 'Synthetic-DR4 fits vs. real-DR4 fits: are the '
                                       'simulations realistic?'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(HERE, 'project_flowchart.png'))
    a = ap.parse_args()

    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    ## Stage-1 / stage-2 counts from the catalog files written by sample/catalog.py
    import pandas as pd
    cat = pd.read_csv(os.path.join(HERE, 'stage1_catalog.csv'))
    sysm = pd.read_csv(os.path.join(HERE, 'stage1_systems.csv'))
    fit = cat[cat.fit_set]
    stage1 = f'stage 1: {len(sysm)} systems, {len(cat)} companions'
    stage2 = (f'stage 2 (S/N >= 0.1, the fit set):\n{int(sysm.fit_set.sum())} systems, '
              f'{len(fit)} companions')
    BOXES['cens'] = BOXES['cens'][:2] + (BOXES['cens'][2].format(stage1=stage1, stage2=stage2),
                                         BOXES['cens'][3])

    fig, ax = plt.subplots(figsize=(24, 15))
    X = 3.5                                 # column spacing
    pos = {k: (v[0]*X, v[1]*1.9) for k, v in BOXES.items()}

    def edge(k, toward):
        """Point on the box edge closest to `toward`."""
        x, y = pos[k]
        tx, ty = pos[toward]
        if abs(ty - y) > abs(tx - x)*0.45:
            return (x, y + (H/2 if ty > y else -H/2))
        return (x + (W/2 if tx > x else -W/2), y)

    for k, (x, y, text, col) in BOXES.items():
        px, py = pos[k]
        ax.add_patch(FancyBboxPatch((px - W/2, py - H/2), W, H, boxstyle='round,pad=0.05',
                                    fc=col, ec=col, alpha=0.18, lw=2))
        ax.add_patch(FancyBboxPatch((px - W/2, py - H/2), W, H, boxstyle='round,pad=0.05',
                                    fc='none', ec=col, lw=2))
        ax.text(px, py, text, ha='center', va='center', fontsize=14)
    for a_, b_ in LANES:
        ax.add_patch(FancyArrowPatch(edge(a_, b_), edge(b_, a_), arrowstyle='-|>',
                                     mutation_scale=22, lw=2, color='0.25'))
    handles = []
    for lab, col, a_, b_, rad, desc in COMPARISONS:
        ax.add_patch(FancyArrowPatch(edge(a_, b_), edge(b_, a_), arrowstyle='<|-|>',
                                     mutation_scale=22, lw=3, color=col, ls='--',
                                     connectionstyle=f'arc3,rad={rad}'))
        (x0, y0), (x1, y1) = edge(a_, b_), edge(b_, a_)
        ax.text((x0 + x1)/2 + 0.45*rad*(y1 - y0), (y0 + y1)/2 - 0.45*rad*(x1 - x0), lab,
                fontsize=18, fontweight='bold', color=col, ha='center', va='center',
                bbox=dict(fc='white', ec=col, boxstyle='circle,pad=0.2'))
        handles.append(plt.Line2D([], [], color=col, ls='--', lw=3, label=f'{lab}: {desc}'))
    for k, t in (('cens', 'Lane 1: our pipeline'), ('dr4', 'Lane 2: real DR4'),
                 ('epo', 'Lane 3: Epochalypse')):
        pass
    ax.text(-0.9, 2*1.9 + 1.0, 'Lane 1: our pipeline', fontsize=15, fontweight='bold', color='C0')
    ax.text(X - W/2, 1.9 + 0.95, 'Lane 2: real DR4', fontsize=15, fontweight='bold', color='C1')
    ax.text(X - W/2, 0.95, 'Lane 3: Epochalypse', fontsize=15, fontweight='bold', color='C2')
    ax.legend(handles=handles, loc='lower left', fontsize=14, frameon=True,
              title='Dashed arrows: the four comparisons', title_fontsize=15,
              bbox_to_anchor=(0.0, 0.0))
    ax.set_xlim(-2.2, 3*X + 2.2)
    ax.set_ylim(-3.4, 3*1.9 + 1.2)
    ax.axis('off')
    ax.set_title('Project flowchart (draft, reconstructed from PLAN.md): three lanes feeding '
                 'four comparisons', fontsize=18)
    fig.savefig(a.out, dpi=300, bbox_inches='tight')
    print('wrote', a.out)


if __name__ == '__main__':
    main()
