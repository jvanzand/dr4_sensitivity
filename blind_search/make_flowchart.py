"""Draw the flowchart of the blind-search algorithm as a PNG
(blind_search/search_flowchart_detailed.png): the developer
version with the mathematics and background. The simple version for the README
is drawn by make_flowchart_simple.py.

Left column: the steps of one loop iteration (what each step does, its inputs
and outputs, the decision). Right column (grey): the background needed to
understand that step. Text is typeset with LaTeX (matplotlib `text.usetex`;
needs a TeX installation), so the symbols are the ones used in the notes. It
documents the algorithm as currently written (`search.run_search`, with the
simulated threshold). Edit the box texts below when the code changes. Usage:
python -I blind_search/make_flowchart.py [--test]   (--test writes
    /tmp/search_flowchart_test.png)
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

W = 20.2                                   # canvas width, inches (canvas units are inches)
XL0, XL1 = 0.9, 11.4                       # left (flowchart) column
XR0, XR1 = 12.55, 19.95                      # right (background) column
FS_L, FS_R = 15.0, 12.0                    # font sizes, points
TALL = re.compile(r'\\frac|\\dfrac|\\sum|\\underbrace|\\sqrt')

# ---------------------------------------------------------------------------------------------------------
# Content. Item kinds: title, io (inputs / outputs), bul (bullet; list of lines), eq (centred equation),
# dec (decision, boxed), tag (small grey code pointer), gap.
# ---------------------------------------------------------------------------------------------------------
INPUT_L = [
    ('title', r'\textbf{Start: data, trial periods, empty model}'),
    ('bul', [r'\textbf{Data}: $N$ transits (93 for Gaia-4), each with epoch $t_i$, along-scan position $y_i$ (mas),',
             r'formal error $\sigma_i$ (mas), scan angle $\psi_i$, along-scan parallax factor $f_i$.']),
    ('bul', [r'\textbf{Trial periods} $P_j$, $j=1,\dots,N_P$ ($N_P=564$ for Gaia-4): uniform in frequency, spacing',
             r'$\Delta f = 1/(2\pi T)$, from $20\,\mathrm{d}$ to $2T$, with $T$ the time span of the data.']),
    ('bul', [r'\textbf{Known companions} $\mathcal{K}=\varnothing$: the list of companions already in the model,',
             r'each stored as $(P,\,e,\,M_0)$.']),
    ('tag', r'code: \texttt{load\_gaia\_dr4} in \texttt{data.py}, \texttt{Stack} in \texttt{model.py}, '
            r'\texttt{period\_grid} in \texttt{periodogram.py}'),
]
INPUT_R = [
    ('title', r'\textbf{Background for the start: the model that is fitted}'),
    ('eq', r'$y_i \;=\; b_i \;+\; \sum_{k\in\mathcal{K}} \eta_k(t_i) \;+\; \epsilon_i,\qquad '
           r'\epsilon_i\sim\mathcal{N}\!\left(0,\;\sigma_i^2+s^2\right)$'),
    ('bul', [r'Gaia measures one coordinate per transit: the sky offset projected on',
             r'$(\sin\psi_i,\cos\psi_i)$ (east, north).']),
    ('bul', [r'$b_i=\Delta\alpha^*_0\sin\psi_i+\Delta\delta_0\cos\psi_i+\mu_{\alpha^*}\Delta t_i\sin\psi_i$',
             r'$\quad+\,\mu_{\delta}\Delta t_i\cos\psi_i+\varpi f_i$: the star\textquoteright s reference motion',
             r'(5 linear parameters: position, proper motion, parallax $\varpi$).']),
    ('bul', [r'$\eta_k(t)=(B_kX+G_kY)\sin\psi+(A_kX+F_kY)\cos\psi$ is the star\textquoteright s reflex orbit',
             r'around the barycentre, with $X=\cos E-e$, $Y=\sqrt{1-e^2}\,\sin E$, and $E$ from',
             r'$M=E-e\sin E=M_0+2\pi (t-t_\mathrm{ref})/P$ ($t_\mathrm{ref}$ = mid-baseline).']),
    ('bul', [r'$A,B,F,G$ (Thiele--Innes constants, mas) set the orbit\textquoteright s size and orientation.',
             r'Given $(P,e,M_0)$ \textbf{every other parameter is linear}.']),
    ('bul', [r'$s$ = \textbf{jitter}: extra white-noise error (mas) added in quadrature to every $\sigma_i$,',
             r'fitted from the data.']),
]

A_L = [
    ('title', r'\textbf{Step A \;|\; Compute the $\Delta$BIC periodogram}'),
    ('io', r'\textit{in}: data, $\mathcal{K}$, $\{P_j\}$ \qquad \textit{out}: $\Delta\mathrm{BIC}_j$ for $j=1,\dots,N_P$'),
    ('bul', [r'Fit the \textbf{current model} $\mathcal{M}_0$ (the $|\mathcal{K}|$-planet model) = reference motion + the companions in $\mathcal{K}$',
             r'(their $P,e,M_0$ fixed, their $A,B,F,G$ refitted):\; $\ln\mathcal{L}_0$.']),
    ('bul', [r'For every $P_j$ fit the \textbf{trial model} $\mathcal{M}_j$ (the $(|\mathcal{K}|+1)$-planet model) = $\mathcal{M}_0$ + one circular',
             r'companion ($e=0$) with free $A,B,F,G$:\; $\ln\mathcal{L}_j$.']),
    ('eq', r'$\Delta\mathrm{BIC}_j \;=\; 2\left(\ln\mathcal{L}_j-\ln\mathcal{L}_0\right) \;-\; 5\ln N$'),
    ('bul', [r'Jitter $s$ is \textbf{refitted in both models at every $P_j$}.']),
    ('tag', r'code: \texttt{delta\_bic} in \texttt{periodogram.py}; \texttt{fit\_free\_jitter} in \texttt{model.py}'),
]
A_R = [
    ('title', r'\textbf{Background for Step A: $\Delta$BIC, the penalty, the jitter fit}'),
    ('bul', [r'$\Delta$BIC (difference of Bayesian information criteria) is twice the gain in',
             r'log-likelihood minus $k\ln N$ for the $k$ new parameters: positive means the companion',
             r'is favoured even after paying for its parameters.']),
    ('bul', [r'$k=5$: the period and 4 linear amplitudes $A,B,F,G$. A circular template suffices',
             r'because free $A,B,F,G$ absorb the phase and the ellipse shape; eccentricity enters in Step E.']),
    ('bul', [r'\textbf{Fit for a given jitter}: weighted least squares, weights $1/(\sigma_i^2+s^2)$, residuals $r_i$,',
             r'$\ln\mathcal{L}=-\tfrac12\sum_i\left[\dfrac{r_i^2}{\sigma_i^2+s^2}+\ln 2\pi(\sigma_i^2+s^2)\right]$.']),
    ('bul', [r'\textbf{Jitter update} (maximum likelihood), alternated with the fit 4 times:',
             r'$s^2\leftarrow\max\left[0,\ \dfrac{\sum_i w_i^2\,(r_i^2-\sigma_i^2)}{\sum_i w_i^2}\right]$, $\;w_i=(\sigma_i^2+s^2)^{-1}$.']),
    ('bul', [r'\textbf{Why refit it everywhere}: a strong unmodelled signal inflates $s$ in the null',
             r'model, which bounds the spurious gain at wrong periods. With $s$ fixed at 0 that gain',
             r'swamped the eFAP threshold (an earlier bug).']),
    ('bul', [r'Noise only: the gain at one period is $\chi^2$ with 4 degrees of freedom (mean 4).']),
]

B_L = [
    ('title', r'\textbf{Step B \;|\; Set the detection threshold $T_\mathrm{sim}$ by simulation}'),
    ('io', r'\textit{in}: data, $\mathcal{K}$, fitted jitter $s$ \qquad \textit{out}: $T_\mathrm{sim}$'),
    ('bul', [r'Draw $M=50$ noise-sampled datasets: the fitted $|\mathcal{K}|$-planet model at the real epochs',
             r'plus noise $\mathcal{N}(0,\sigma_i^2+s^2)$, new random numbers for each.']),
    ('bul', [r'Run each one through Step A (known companions re-polished, jitter free) and keep only',
             r'its highest peak $h_m$, $m=1,\dots,M$. Fit a Gumbel distribution (location $\mu$, scale $\beta$) to $\{h_m\}$.']),
    ('eq', r'$T_\mathrm{sim} \;=\; \mu \;-\; \beta\,\ln\!\left[-\ln(1-\mathrm{FAP})\right],\qquad \mathrm{FAP}=10^{-3}$'),
    ('tag', r'code: \texttt{simulated\_threshold} in \texttt{simthreshold.py}'),
]
B_R = [
    ('title', r'\textbf{Background for Step B: where the threshold comes from}'),
    ('bul', [r'\textbf{Idea}: the false-alarm probability concerns the highest peak over all periods, so the',
             r'noise-only distribution of the highest peak is simulated directly.']),
    ('bul', [r'The noise-sampled datasets contain no unfound signal, so the real periodogram\textquoteright s body, and the',
             r'leakage of any real signal, never enters the threshold.']),
    ('bul', [r'\textbf{Gumbel (extreme-value) fit}: maxima of many roughly independent values follow it;',
             r'two parameters from $M=50$ maxima give a smooth quantile far in the tail, where a',
             r'direct quantile (the $1/\mathrm{FAP}$-th largest of $M$) does not exist.']),
    ('bul', [r'\textbf{Caveat}: the simulated noise is white (formal error plus fitted jitter);',
             r'correlated noise is not modelled.']),
    ('bul', [r'Unlike RVSearch\textquoteright s eFAP (tail of the periodogram\textquoteright s own log-histogram,',
             r'\texttt{efap.py}), $T_\mathrm{sim}$ does not rise when a strong signal is present.']),
]

D_L = [
    ('title', r'\textbf{Step C--D \;|\; Pick the peak and decide}'),
    ('io', r'\textit{in}: $\{\Delta\mathrm{BIC}_j\}$, $T_\mathrm{sim}$ \qquad \textit{out}: stop or continue'),
    ('bul', [r'Highest peak: $P_\star=\arg\max_j \Delta\mathrm{BIC}_j$.']),
    ('bul', [r'\textbf{yes}: add the companion at $P_\star$ (Step E).\;\; \textbf{no}: stop.']),
    ('dec', r'\textbf{Decision: is $\max_j \Delta\mathrm{BIC}_j > T_\mathrm{sim}$\,?}'),
]
D_R = [
    ('title', r'\textbf{Background for Step C--D: the decision}'),
    ('bul', [r'One companion is added per iteration, the one at the highest peak.']),
    ('bul', [r'A real search stops at the first non-detection (\texttt{run\_search} in \texttt{search.py}).']),
    ('bul', [r'In injection tests a detection counts as a match if its period is within 3\% of an',
             r'injected period.']),
]

STOP_L = [
    ('title', r'\textbf{Stop}'),
    ('bul', [r'Report $\mathcal{K}$:', r'$k-1$ companions.']),
    ('bul', [r'A real search ends here.']),
]
E_L = [
    ('title', r'\textbf{Step E \;|\; Refit the companions}'),
    ('io', r'\textit{in}: $\mathcal{K}$, $P_\star$ \quad \textit{out}: enlarged $\mathcal{K}$'),
    ('bul', [r'Maximise $\ln\mathcal{L}(\ln P,e,M_0)$ over the new',
             r'companion\textquoteright s three nonlinear parameters from',
             r'6 starting guesses $(e_0,M_{0,0})$, with $|P/P_\star-1|\le 0.1$',
             r'and $0\le e\le 0.7$.']),
    ('bul', [r'If $\mathcal{K}\neq\varnothing$: maximise again over all',
             r'$3(|\mathcal{K}|+1)$ nonlinear parameters jointly.']),
    ('bul', [r'$\mathcal{K}\leftarrow\mathcal{K}\cup\{(P,e,M_0)\}$; return to Step A.']),
    ('tag', r'code: \texttt{add\_and\_refit} in \texttt{fit.py}'),
]
E_R = [
    ('title', r'\textbf{Background for Step E: the refit}'),
    ('bul', [r'\textbf{Variable projection}: for fixed $(P,e,M_0)$ of all companions, the 5 baseline',
             r'parameters, every $A,B,F,G$ and $s$ are solved exactly as in Step A, so only',
             r'$3|\mathcal{K}|$ nonlinear numbers are searched.']),
    ('bul', [r'\textbf{Nelder--Mead}: a derivative-free simplex optimiser. The 6 starts guard against',
             r'local maxima: $(e_0,M_{0,0})=(0,0),(0.3,0),(0.3,2.1),(0.3,4.2),(0.6,1.0),(0.6,3.5)$.']),
    ('bul', [r'$e\le 0.7$ is the blind eccentricity cap (as in the RV search); $P$ stays within',
             r'$\pm10\%$ of $P_\star$ so the refit remains on the peak.']),
    ('bul', [r'The joint re-maximisation lets earlier companions adjust once the new one is present.']),
]

OUT_L = [
    ('title', r'\textbf{Hand-off (not written yet)}'),
    ('bul', [r'$\mathcal{K}$ and the fitted linear parameters become starting values for a',
             r'physical-orbit model in vanilla Octofitter v9.']),
]


def layout(items, fs):
    """Lines of one box: list of (x offset, baseline y offset below the box top,
    text, style) and height."""
    u = fs/13.0 if fs < 14 else fs/14.0
    y, out = 0.20, []
    for kind, val in items:
        if kind == 'gap':
            y += 0.12*u
        elif kind in ('title', 'io', 'dec', 'tag', 'eq'):
            dy = {'title': 0.44, 'io': 0.38, 'dec': 0.50, 'tag': 0.30, 'eq': 0.58}[kind]*u
            if kind == 'eq' and TALL.search(val):
                dy *= 1.0
            y += dy + (0.08*u if kind == 'tag' else 0.0)
            out.append((0.0 if kind != 'eq' else 'c', y, val, kind))
            y += 0.10*u if kind in ('title', 'io') else (0.22*u if kind == 'dec' else 0.0)
        elif kind == 'bul':
            y += 0.06*u
            for n, ln in enumerate(val):
                dy = (0.50 if TALL.search(ln) else 0.33)*u
                y += dy
                out.append((0.0 if n == 0 else 0.26, y, (r'$\bullet$\;' + ln) if n == 0 else ln, 'bul'))
    return out, y + 0.18


def render(out):
    import matplotlib
    matplotlib.rc_file(os.path.join(ROOT, 'matplotlibrc'))
    matplotlib.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
                                'text.latex.preamble': r'\usepackage{amsmath}\usepackage{amssymb}'})
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Rectangle

    rows = [  # (left items, right items, left x0, left x1, colour, dashed)
        ('input', INPUT_L, INPUT_R, XL0, XL1, '#dbe8fb'),
        ('A', A_L, A_R, XL0, XL1, '#fdf0cf'),
        ('B', B_L, B_R, XL0, XL1, '#fdf0cf'),
        ('D', D_L, D_R, XL0, XL1, '#ffd9d9'),
    ]
    GAP = 0.62
    # measure
    lay = {}
    for name, L, R, x0, x1, col in rows:
        lay[name] = (layout(L, FS_L), layout(R, FS_R))
    stopl, hstop = layout(STOP_L, FS_L)
    el, he = layout(E_L, FS_L)
    er, her = layout(E_R, FS_R)
    outl, hout = layout(OUT_L, FS_L)
    rowh = {n: max(lay[n][0][1], lay[n][1][1]) for n in lay}
    hE = max(he, her, hstop)
    X_STOP1, X_E0 = 3.7, 4.1
    top = 1.55
    H = top + sum(rowh.values()) + hE + hout + GAP*(len(rows) + 2) + 1.7
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)                                         # y downward, in inches
    ax.axis('off')

    def box(x0, x1, y0, h, col, lw=2.2, dashed=False, ec='k'):
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, h, boxstyle='round,pad=0,rounding_size=0.12', fc=col, ec=ec,
                                    lw=lw, ls='--' if dashed else '-', zorder=1))

    def put(x0, y0, lines, fs, x1):
        for xo, yo, txt, kind in lines:
            if xo == 'c':
                ax.text((x0 + x1)/2, y0 + yo, txt, fontsize=fs, ha='center', va='baseline', zorder=3)
            elif kind == 'tag':
                ax.text(x0 + 0.2, y0 + yo, txt, fontsize=fs - 2.5, ha='left', va='baseline', color='0.4', zorder=3)
            elif kind == 'dec':
                ax.text((x0 + x1)/2, y0 + yo, txt, fontsize=fs + 1, ha='center', va='baseline', zorder=3,
                        bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#b00020', lw=2))
            elif kind == 'title':
                ax.text(x0 + 0.2, y0 + yo, txt, fontsize=fs + 1.5, ha='left', va='baseline', zorder=3)
            else:
                ax.text(x0 + 0.2 + xo, y0 + yo, txt, fontsize=fs, ha='left', va='baseline', zorder=3)

    def arrow(x0, y0, x1, y1, lw=2.6, style='-|>'):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0), zorder=2,
                    arrowprops=dict(arrowstyle=style, lw=lw, color='0.15', mutation_scale=26, shrinkA=0, shrinkB=0))

    ax.text(W/2, 0.42, r'\textbf{Blind search for companions in Gaia DR4 epoch astrometry: one pass of the loop}',
            fontsize=21, ha='center', va='center')
    ax.text(W/2, 0.93, r'Left: what the algorithm does at each step (inputs, outputs, decisions). Grey, right: the background '
            r'needed to understand that step.  $N$ = number of transits, $N_P$ = number of trial periods.',
            fontsize=13, ha='center', va='center')

    y = top
    ys = {}
    xc = 6.1                                                  # x of the main flow arrows
    for name, L, R, x0, x1, col in rows:
        (ll, hl), (rl, hr) = lay[name]
        h = rowh[name]
        ys[name] = y
        box(x0, x1, y, hl, col)
        put(x0, y, ll, FS_L, x1)
        box(XR0, XR1, y, hr, '#f1f1f1', lw=1.6, ec='0.45')
        put(XR0, y, rl, FS_R, XR1)
        if name == 'input':
            lab = None
        bot = y + hl
        y += h + (GAP + 0.38 if name == 'input' else GAP)
        if name == 'input':
            arrow(xc, bot, xc, y)
        elif name == 'A':
            arrow(xc, bot, xc, y)
            ax.text(xc + 0.2, (bot + y)/2, r'$\{\Delta\mathrm{BIC}_j\}$, $j=1,\dots,N_P$', fontsize=13, va='center')
        elif name == 'B':
            arrow(xc, bot, xc, y)
            ax.text(xc + 0.2, (bot + y)/2, r'$T_\mathrm{sim}$', fontsize=13, va='center')
    yD = ys['D']
    yE = y
    # stop and step E side by side
    box(XL0, X_STOP1, yE, hstop, '#d8f0d4')
    put(XL0, yE, stopl, FS_L, X_STOP1)
    box(X_E0, XL1, yE, he, '#fdf0cf')
    put(X_E0, yE, el, FS_L, XL1)
    box(XR0, XR1, yE, her, '#f1f1f1', lw=1.6, ec='0.45')
    put(XR0, yE, er, FS_R, XR1)
    yDb = yD + lay['D'][0][1]
    arrow(3.2, yDb, 2.3, yE)
    ax.text(2.3, yDb + GAP/2 - 0.02, r'\textbf{no}', fontsize=15, ha='right', va='center')
    arrow(7.8, yDb, 7.8, yE)
    ax.text(8.0, yDb + GAP/2 - 0.02, r'\textbf{yes}', fontsize=15, ha='left', va='center')
    # loop frame and return arrow
    fy0, fy1 = ys['A'] - 0.38, yE + hE + 0.22
    ax.add_patch(Rectangle((0.15, fy0), XL1 + 0.95 - 0.15, fy1 - fy0, fc='none', ec='#7a3b00', lw=2.2, ls=':', zorder=0))
    ax.text(0.3, fy0 + 0.02, r'\textbf{LOOP over $k=1,2,3,\dots$}', fontsize=14, color='#7a3b00', va='bottom', ha='left',
            bbox=dict(fc='white', ec='none', pad=1))
    # return arrow: from the right edge of step E up the right margin to step A
    xr = XL1 + 0.3
    ax.plot([XL1, xr, xr], [yE + he - 0.5, yE + he - 0.5, ys['A'] + 0.45], color='0.15', lw=2.6, zorder=2)
    arrow(xr, ys['A'] + 0.45, XL1, ys['A'] + 0.45)
    ax.text(XL1 + 0.62, ys['A'] + 0.5, r'next iteration: $k\to k+1$', fontsize=12, rotation=90, ha='center', va='top',
            color='0.15')
    # output
    yO = fy1 + GAP
    arrow(2.3, yE + hstop, 2.3, yO)
    xo1 = XL1
    box(XL0, xo1, yO, hout, '#eeeeee', dashed=True)
    put(XL0, yO, outl, FS_L, xo1)
    ax.set_ylim(yO + hout + 0.3, 0)
    fig.set_size_inches(W, yO + hout + 0.3)
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main():
    out = '/tmp/search_flowchart_test.png' if '--test' in sys.argv else os.path.join(HERE, 'search_flowchart_detailed.png')
    render(out)
    print('wrote', out)


if __name__ == '__main__':
    main()
