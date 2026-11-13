"""Build ren_2026_orbits.csv from the arXiv:2610.16662 LaTeX source
(Ren et al. 2026, AJ accepted): HD 203771 B.

Run from anywhere:  python -I build.py
Reads  src_2610.16662/main.tex, Table "allpriors" (two joint fits; no RV-only fit given).
Writes ren_2026_orbits.csv

The paper's 159-brown-dwarf sample is a compilation shown only in a figure (Gan+2025,
Barbato+2026, Unger+2023, Feng+2022, Bowler+2020); it is not tabulated, so only
HD 203771 B is recorded. We take the RV+Hip+GDR2+GDR3 column because it is the only one
that reports K; the RV+Gaia two-body column is summarized in notes.
Uncertainties (same column): the table prints +/- errors for K, M_c, a, e and P; asymmetric ones are averaged.
e_msini_mj = mean error on M_c * sin(i) (the printed true-mass error converted the same way as the mass; the
inclination error is not folded in).
"""
import csv
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "src_2610.16662" / "main.tex"
COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp", "tp_system",
        "msini_mj", "a_au", "mstar_msun", "fit_type", "table_ref", "notes",
        "e_K_ms", "e_msini_mj", "e_e", "e_a_au"]
NUM = r"[-+]?\d+(?:\.\d+)?"
TP_REF = 2457389.0  # "t_p (JD-2457389)" per the table row label


def parse(cell):
    c = cell.replace("$", "").replace("^{\\dag}", "").replace("\\dag", "").replace(" ", "")
    for pat, f in [
        (rf"^({NUM})\\pm({NUM})$", lambda m: (float(m[1]), float(m[2]))),
        (rf"^({NUM})\^\{{\+({NUM})\}}_\{{-({NUM})\}}$",
         lambda m: (float(m[1]), 0.5 * (float(m[2]) + float(m[3])))),
        (rf"^({NUM})_\{{-({NUM})\}}\^\{{\+({NUM})\}}$",
         lambda m: (float(m[1]), 0.5 * (float(m[2]) + float(m[3])))),
        (rf"^({NUM})$", lambda m: (float(m[1]), None)),
    ]:
        m = re.match(pat, c)
        if m:
            return f(m)
    raise ValueError(cell)


def main():
    tex = TEX.read_text()
    body = tex[tex.index("\\caption{Parameter priors and best-fits in two joint models"):
               tex.index("\\label{allpriors}")]
    rows = {}
    for line in body.splitlines():
        s = line.strip()
        if s.startswith("%") or s.count("&") < 3:
            continue
        cells = [c.strip() for c in s.rstrip("\\").split("&")]
        rows[cells[0]] = cells

    def get(label, col):  # col 2 = RV+Gaia two-body, col 3 = RV+Hip+GDR2+GDR3
        return parse(rows[label][col])

    P, eP = get("$P$ (days)", 3)
    tp = get("$t_{p}$ (JD-2457389)", 3)[0] + TP_REF
    Mc, eMc = get("$M_{c}$ ($M_{\\mathrm{J}}$)", 3)
    K, eK = get("$K$ (m~s$^{-1}$)", 3)
    a, ea = get("$a$ (AU)", 3)
    inc = get("$i$ (deg)", 3)[0]
    e, ee = get("$e$", 3)
    w_rad = get("$\\omega$ (rad)", 3)[0]
    mstar = get("$M_\\ast$ ($M_\\odot$)", 3)[0]
    msini = Mc * math.sin(math.radians(inc))

    # Two-body column for the notes.
    P2 = get("$P$ (days)", 2)[0]
    tp2 = get("$t_{p}$ (JD-2457389)", 2)[0] + TP_REF
    Mc2 = get("$M_{c}$ ($M_{\\mathrm{J}}$)", 2)[0]
    i2 = get("$i$ (deg)", 2)[0]
    e2 = get("$e$", 2)[0]
    w2 = get("$\\omega$ (rad)", 2)[0]

    notes = (
        "HIP 105740 = Gaia DR3 6827680034091744896; elements from the RV+Hip+GDR2+GDR3 joint fit "
        f"(HARPS RVs from HARPS-RVBank + Hipparcos IAD + Gaia DR2/DR3 catalog astrometry); true mass "
        f"Mc={Mc} MJ at i={inc} deg, msini_mj = Mc*sin(i) computed here; omega converted from "
        f"{w_rad} rad (paper uses the Gaia/Winn2022 ascending-node convention, consistent in the RV "
        f"model); Tp = {TP_REF:.0f} + table t_p (derived, not fitted, in this column), table "
        "labels it JD; e_P_d is the mean of asymmetric errors. "
        f"RV+Gaia DR3 two-body column: P={P2} d, tp=JD {tp2:.3f}, e={e2}, omega={w2} rad, "
        f"Mc={Mc2} MJ, i={i2} deg (no K reported). New fit; earlier RV+astrometry orbit in "
        "Feng+2022 (see feng_2022, stevenson_2023). The paper's 159-BD sample is not tabulated."
    )
    row = {
        "host": "HD 203771",
        "companion": "HD 203771 B",
        "P_d": f"{P:g}",
        "e_P_d": f"{eP:.3g}",
        "K_ms": f"{K:g}",
        "e": f"{e:g}",
        "omega_deg": f"{math.degrees(w_rad):.2f}",
        "Tp": f"{tp:.3f}",
        "tp_system": "JD",
        "msini_mj": f"{msini:.2f}",
        "a_au": f"{a:g}",
        "mstar_msun": f"{mstar:g}",
        "e_K_ms": f"{eK:.4g}" if eK else "",
        "e_msini_mj": f"{eMc * math.sin(math.radians(inc)):.4g}" if eMc else "",
        "e_e": f"{ee:.4g}" if ee else "",
        "e_a_au": f"{ea:.4g}" if ea else "",
        "fit_type": "RV+astrometry",
        "table_ref": "arXiv:2610.16662 Table 2 (allpriors), RV+Hip+GDR2+GDR3 column",
        "notes": notes,
    }
    with open(HERE / "ren_2026_orbits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerow(row)
    print(row)


if __name__ == "__main__":
    main()
