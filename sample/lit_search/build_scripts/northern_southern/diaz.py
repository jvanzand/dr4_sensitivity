import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/diaz_2012/'
T = open(R + 'src_1111.1168/TableParameters_v3.tex').read().split('\n')
S = open(R + 'src_1111.1168/TableStellarParameters_v1.tex').read().split('\n')
def blocks(lines, keys):
    out = {}; names = None
    for l in lines:
        if l.strip().startswith('%'):
            continue
        c = cells(l)
        if len(c) == 5 and c[0] in ('', 'Parameters') and c[1].replace('\\tablefootmark{\\dagger}', '').startswith('HD'):
            names = [x.replace('\\tablefootmark{\\dagger}', '') for x in c[1:]]
            continue
        for k in keys:
            if names and c[0].startswith(k) and len(c) >= 5:
                for n, x in zip(names, c[1:5]):
                    out.setdefault(n, {})[k] = x
    return out
P = blocks(T, ['$P$', '$e$', '$\\omega$', '$K$', '$T_0$', '$a$', '$M_\\textrm{c} \\sin i', '$M_{c}$'])
M = blocks(S, ['$M_\\textrm{s}$'])
notes = {
 'HD22781': 'Likely in NASA Exoplanet Archive (msini 13.65). A second companion probably present (linear+quadratic drift fitted); parameters may change.',
 'HD156279': 'Massive planet (msini 9.71); likely in NASA Exoplanet Archive.',
 'HD14651': 'Hipparcos 3-sigma upper limit on Mc 0.28-0.40 Msun.',
 'HD30246': 'Poorly constrained orbit (K +1600/-560 m/s, e 0.84+/-0.08); msini +20.3/-8.2 MJ.',
 'HD92320': 'Hipparcos upper limit 0.26-0.34 Msun.',
 'HD137510': 'Fit combines SOPHIE with McDonald, TLS and Lick RVs (Endl+2004 discovery). Hipparcos constrains Mc to 0.019-0.057 Msun (99.7%); so a BD.',
 'HD16702': 'Hipparcos astrometric orbit: Mc = 0.35-0.45 Msun (99.7%): low-mass star, not BD.',
}
rows = []
for n in ['HD22781', 'HD156279', 'HD14651', 'HD30246', 'HD92320', 'HD137510', 'HD16702']:
    p = P[n]
    host = 'HD ' + n[2:]
    rows.append(dict(host=host, companion='b', P_d=val(p['$P$']), e_P_d=err(p['$P$']), K_ms=val(p['$K$']), e=val(p['$e$']),
        omega_deg=val(p['$\\omega$']), Tp_bjd=addbjd(p['$T_0$'], 2454000), msini_mj=val(p['$M_\\textrm{c} \\sin i']),
        a_au=val(p['$a$']), mstar_msun=val(M[n]['$M_\\textrm{s}$']), fit_type='RV-only',
        table_ref='arXiv:1111.1168 Table 5 (Mstar Table 3)',
        notes='SOPHIE RVs, single Keplerian; T0 = periastron (BJD-2454000 in paper, converted). ' + notes.get(n, '')))
write(R + 'diaz_2012_orbits.csv', rows)
