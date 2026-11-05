"""Build lee_2023_orbits.csv from the arXiv LaTeX source of Lee et al. 2023
(arXiv:2307.15897), Table 3 "Orbital solutions for HD 6860 b and HD 112300 b
probed by the MCMC simulation" (label tab3).

Source downloaded from https://arxiv.org/e-print/2307.15897 into src_2307.15897/.
Run from this directory:  python -I build.py
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp",
        "tp_system", "msini_mj", "a_au", "mstar_msun", "fit_type", "table_ref",
        "notes"]

tex = (HERE / "src_2307.15897" / "AA43725.tex").read_text()
i = tex.index("\\label{tab3}")
body = tex[tex.index("\\hline", i):tex.index("\\end{tabular}", i)]

VAL = r"\s*(-?[\d.]+)\s*(?:\$?\^\{\s*\+?\s*([\d.]+)\}_\{\s*-([\d.]+)\}\$?)?"
par = {}
for line in body.splitlines():
    line = line.split("%")[0]
    if "&" not in line:
        continue
    cells = [c.strip() for c in line.replace("\\\\", "").split("&")]
    vals = []
    for c in cells[1:]:
        m = re.match(VAL, c.replace("$", "$"))
        vals.append(m.groups() if m else None)
    par[cells[0]] = vals


def get(key, k):
    for name, v in par.items():
        if name.startswith(key):
            return v[k]
    raise KeyError(key)


# Stellar masses are stated in the text (Sect. 4), not in Table 3:
#  HD 6860:   "Assuming 2.49 M_sun (Dehaes et al. 2011) for the mass of HD 6860"
#  HD 112300: "Assuming a stellar mass of 1.4 +/- 0.3 M_sun"
mstar = {}
m = re.search(r"Assuming ([\d.]+) \$M_\{\\odot\}\$ \\citep\{2011A&A...533A.107D\} for the mass of HD~6860", tex)
mstar["HD 6860"] = m.group(1)
m = re.search(r"Assuming a stellar mass of ([\d.]+) \$\\pm\$ [\d.]+~\$M_\{\\odot\}\$", tex)
mstar["HD 112300"] = m.group(1)

extra = {
    "HD 6860": ("beta And (Mirach), M0 III, SIMBAD V=2.05: likely too bright for reliable Gaia astrometry. Paper attributes a "
                "separate ~1200 d signal in activity indicators to chromospheric activity. "),
    "HD 112300": ("delta Vir, SIMBAD V=3.38, M1 III semiregular variable (Hipparcos ~693 d photometric "
                  "period, distinct from the 467 d RV period); fit includes a linear RV "
                  "slope. "),
}
rows = []
for k, (host, comp) in enumerate([("HD 6860", "HD 6860 b"), ("HD 112300", "HD 112300 b")]):
    P, Pp, Pm = get("P (days)", k)
    rows.append(dict(
        host=host, companion=comp, P_d=P,
        e_P_d="%.2f" % ((float(Pp) + float(Pm)) / 2),
        K_ms=get("K (m", k)[0], e=get("$e$", k)[0], omega_deg=get("$\\omega$", k)[0],
        Tp=get("$T_{periastron}$", k)[0], tp_system="JD",
        msini_mj=get("$m$ sin", k)[0], a_au=get("$a$ (AU)", k)[0],
        mstar_msun=mstar[host], fit_type="RV-only",
        table_ref="arXiv:2307.15897 Table 3",
        notes=(extra[host] + "BOAO/BOES iodine RVs only (Exo-Striker MCMC); median values; "
               "e_P_d is the mean of the asymmetric +%s/-%s errors; Tp is 'T_periastron (JD)' "
               "as printed. M giant host with large RV jitter (Table 3), so an intrinsic "
               "origin cannot be fully excluded. mstar from text (Sect. 4)." % (Pp, Pm))))

with open(HERE / "lee_2023_orbits.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rows written")
