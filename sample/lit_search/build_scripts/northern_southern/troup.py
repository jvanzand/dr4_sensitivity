import sys; sys.path.insert(0, '/private/tmp/claude-501/-Users-judahvz-research-code-GitHub-latex-documents-companion-mass-vs-stellar/ceb6b3d7-e00e-4270-8185-44729956fa3f/scratchpad/scripts')
from lt import *
import numpy as np
from astropy.io import fits
R = '/Users/judahvz/research/code/GitHub/gaia_dr4_rv_sensitivity/sample/raw_lit/troup_2016/'
d = fits.open(R + 'src_1601.00688/goldOrbit-dr12.fits')[1].data
MJ = 1047.57  # MJup per Msun
rows = []
for r in d:
    msini = r['MSINI'] * MJ
    if not (5 <= msini <= 85):
        continue
    aid = r['APOGEE_ID'].strip()
    host = '2MASS J' + aid[2:]
    rows.append(dict(host=host, companion='b', P_d='%.6g' % r['PERIOD'], e_P_d='%.3g' % r['PERIOD_ERR'],
        K_ms='%.1f' % r['SEMIAMP'], e='%.4f' % r['ECC'], omega_deg='%.2f' % np.degrees(r['OMEGA']),
        Tp_bjd='%.5f' % r['TPERI'], msini_mj='%.2f' % msini, a_au='%.4f' % r['SEMIMAJ'], mstar_msun='%.3f' % r['MSTAR'],
        fit_type='RV-only', table_ref='arXiv:1601.00688 goldOrbit-dr12.fits (= Table 4, VizieR J/AJ/151/85)',
        notes='APOGEE-1 DR12 gold sample; class=%s; NRV=%d over %.0f d baseline; sigma_RV~%.0f m/s; e_K=%.3g m/s; e_e=%.3g; PUI=%.2f VUI=%.2f; slope=%.3f m/s/d. Tp is JD as given (time standard not stated). msini converted from Msun with 1047.57 MJ/Msun. Few-epoch, low-precision orbit: treat as candidate; period aliases possible.' % (
            r['STARTYPE'].strip(), r['NVISITS'], r['BASELINE'], r['SIGMA_V'], r['SEMIAMP_ERR'], r['ECC_ERR'], r['PUI'], r['VUI'], r['SLOPE'])))
write(R + 'troup_2016_orbits.csv', rows)
ms = np.array([float(x['msini_mj']) for x in rows]); a = np.array([float(x['a_au']) for x in rows])
print('msini 13-80:', ((ms >= 13) & (ms <= 80)).sum(), ' <13:', (ms < 13).sum(), ' a max', a.max())
from collections import Counter; print(Counter(x['notes'].split('class=')[1][:2] for x in rows))
