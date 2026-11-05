import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/grieves_2017/'
L = open(R + 'src_1702.01784/grieves2017_arxiv.tex').read().split('\n')
i0 = [i for i, l in enumerate(L) if 'label{tab:compar}' in l][0]
j0 = [i for i, l in enumerate(L) if 'label{tab:starpar}' in l][0]
ms = {}
for l in L[j0:j0 + 25]:
    c = cells(l)
    if len(c) == 7 and '\\pm' in l:
        ms[re.sub(r'\$.*?\$', '', c[0]).split('(')[0].strip()] = val(c[4])
A, B = {}, {}
for l in L[i0:i0 + 45]:
    c = cells(l)
    if len(c) == 6 and 'MARVELS-' in c[1]:
        A[c[1]] = (c[0].split('(')[0].strip(), c[2:])
    if len(c) == 8 and c[0].startswith('MARVELS-'):
        B[c[0]] = c[1:]
rows = []
for d, (star, (msini, P, ecc, K)) in A.items():
    a, Tp, w, gam, dvdt, rms, chi = B[d]
    m = float(val(msini))
    note = '%s. MARVELS UF2D RVs (+ any follow-up as in paper), RVLIN fit, bootstrap errors. Tp JD-2450000 in paper (converted; JD not BJD). fit rms %s m/s; dv/dt %s m/s/d.' % (d, val(rms), val(dvdt))
    if m > 85:
        continue
    if m > 13 and m < 14:
        note += ' msini 13.1+/-1.6 near deuterium limit.'
    rows.append(dict(host=tyc(star), companion='b', P_d=val(P), e_P_d=err(P), K_ms=val(K), e=val(ecc), omega_deg=val(w),
        Tp_bjd=addbjd(Tp, 2450000), msini_mj=val(msini), a_au=val(a), mstar_msun=ms[star], fit_type='RV-only',
        table_ref='arXiv:1702.01784 Table 2 (Mstar Table 1)', notes=note + ' Candidate: paper notes some fits have weak phase coverage/aliases.'))
write(R + 'grieves_2017_orbits.csv', rows)
