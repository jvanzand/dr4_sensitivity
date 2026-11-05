"""Build barbato_2026_orbits.csv from the arXiv:2606.29001 LaTeX source
(Barbato et al. 2026, GAPS LXXV, A&A accepted).

Run from anywhere:  python -I build.py
Reads  src_2606.29001/aa57585-25.tex
       table:rvfit          (HARPS-N RV-only PyORBIT Keplerian fits, bottom block)
       table:star-validated (Gaia DR3 IDs and SED stellar masses)
Writes barbato_2026_orbits.csv

The fits use mean longitude lambda_0, but the paper does not state its reference epoch,
so Tp is left blank.
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "src_2606.29001" / "aa57585-25.tex"

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp", "tp_system",
        "msini_mj", "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

NUM = r"[-+]?\d+(?:\.\d+)?"

# SIMBAD-style names for the column headers (the LaTeX uses thin spaces "\,").
EXTRA_NOTES = {
    "PM J13580+3141": (
        "Gaia-4 b of Stefansson+2025 (per the paper); fit combines HARPS-N with literature "
        "HPF, NEID and FIES RVs and includes a quadratic trend (adopted solution, dBIC=9.11); "
        "single-Keplerian alternative in text: P=574.46 d, K=309.97 m/s, e=0.336"
    ),
    "StKM 1-1494": (
        "companion also confirmed by GRAVITY imaging (Winterhalder+2024, per the paper); "
        "HARPS-N covers only 62% of the NSS period"
    ),
}


def parse_val(cell):
    """Return (value, err) where err is the symmetric or mean asymmetric 1-sigma error."""
    c = cell.strip().rstrip(".").replace("$", "").replace(" ", "")
    m = re.match(rf"^({NUM})\\pm({NUM})$", c)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.match(rf"^({NUM})_\{{-({NUM})\}}\^\{{\+?({NUM})\}}$", c)
    if m:
        return float(m.group(1)), 0.5 * (float(m.group(2)) + float(m.group(3)))
    m = re.match(rf"^({NUM})\^\{{\+({NUM})\}}_\{{-({NUM})\}}$", c)
    if m:
        return float(m.group(1)), 0.5 * (float(m.group(2)) + float(m.group(3)))
    raise ValueError(f"cannot parse {cell!r}")


def table_rows(tex, label):
    """Map row label -> list of cells for the sidewaystable with the given \\label."""
    start = tex.index(f"\\label{{{label}}}")
    body = tex[start:tex.index("\\end{tabular}", start)]
    rows = {}
    for line in body.splitlines():
        if "&" not in line:
            continue
        line = re.sub(r"\\\\(\[3pt\])?\s*$", "", line.strip())
        cells = [c.strip() for c in line.split("&")]
        rows[cells[0]] = cells[1:]
    return rows


def main():
    tex = TEX.read_text()
    fit = table_rows(tex, "table:rvfit")
    star = table_rows(tex, "table:star-validated")
    names = [n.replace("\\,", " ") for n in fit[""]]
    assert names == [n.replace("\\,", " ") for n in star[""]]

    def col(rows, key):
        hits = [k for k in rows if k.startswith(key)]
        assert len(hits) == 1, (key, hits)
        return rows[hits[0]]

    gaia = col(star, "\\gdr{3} ID")
    mstar = col(star, "$M_\\star$")
    P = col(fit, "$P_{\\rm b}$")
    K = col(fit, "$K_{\\rm b}$")
    E = col(fit, "$e_{\\rm b}$")
    W = col(fit, "$\\omega_{\\rm b}$")
    MS = col(fit, "$M_{\\rm b}\\sin{i}$")
    A = col(fit, "$a_{\\rm b}$")
    MI = col(fit, "$M_{b,i_{\\rm NSS}}$")
    NRV = col(fit, "$N_{\\rm RV}$")
    PNSS = col(fit, "$P_{\\rm NSS}$")

    out = []
    for i, name in enumerate(names):
        p, ep = parse_val(P[i])
        msini = parse_val(MS[i])[0]
        notes = (
            f"Gaia DR3 {gaia[i]}; RV-only PyORBIT fit to {NRV[i]} HARPS-N RVs; "
            f"Gaia DR3 NSS astrometric candidate (P_NSS={parse_val(PNSS[i])[0]} d); "
            f"mass using NSS inclination M_b,iNSS={parse_val(MI[i])[0]} MJ; "
            "Tp blank: paper fits mean longitude lambda_0 without stating its reference epoch; "
            "e_P_d is the mean of asymmetric errors where asymmetric; mstar from SED fit "
            "(table:star-validated)"
        )
        if name in EXTRA_NOTES:
            notes += "; " + EXTRA_NOTES[name]
        out.append({
            "host": name,
            "companion": f"{name} b",
            "P_d": f"{p:g}",
            "e_P_d": f"{ep:.3g}",
            "K_ms": f"{parse_val(K[i])[0]:g}",
            "e": f"{parse_val(E[i])[0]:g}",
            "omega_deg": f"{parse_val(W[i])[0]:g}",
            "Tp": "",
            "tp_system": "",
            "msini_mj": f"{msini:g}",
            "a_au": f"{parse_val(A[i])[0]:g}",
            "mstar_msun": f"{parse_val(mstar[i])[0]:g}",
            "fit_type": "RV-only",
            "table_ref": "arXiv:2606.29001 Table table:rvfit (RV fit block); mstar table:star-validated",
            "notes": notes,
        })

    with open(HERE / "barbato_2026_orbits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(out)
    nbd = sum(13 <= float(r["msini_mj"]) <= 80 for r in out)
    print(f"{len(out)} rows; {nbd} with 13 <= msini <= 80 MJ")


if __name__ == "__main__":
    main()
