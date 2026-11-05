"""Build cui_2026_orbits.csv from the arXiv:2603.01808 LaTeX source (Cui et al. 2026, PNAS).

Run from anywhere:  python -I build.py
Reads  src_2603.01808/BDdesert.tex  (longtable "Catalog of companions", label tab:sample)
Writes cui_2026_all.csv     (all 55 rows, with an in_range flag)
       cui_2026_orbits.csv  (rows with M sin i <= 80 MJ and (a <= 5 AU or P <= 4000 d))

The table gives true mass m (MJ), q, e, a (au), P (yr), inclination I (deg) and Omega (deg)
from joint RV + Hipparcos-Gaia fits (posteriors of Feng et al. 2022). No K, omega or Tp.
M sin i is computed here as m * sin(I).
"""
import csv
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEX = HERE / "src_2603.01808" / "BDdesert.tex"
DAYS_PER_YR = 365.25

COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp", "tp_system",
        "msini_mj", "a_au", "mstar_msun", "fit_type", "table_ref", "notes"]

NUM = r"[-+]?\d+(?:\.\d+)?(?:E[-+]?\d+)?"


def parse_cell(cell):
    """Return (value, err_plus, err_minus) from '$x^{+a}_{-b}$' or the q-column form."""
    c = cell.replace("\\mathrm", "").replace("$", "").replace(" ", "")
    m = re.match(rf"^({NUM})\^\{{?\+?({NUM})\}}?_\{{?-?({NUM})\}}?$", c)
    if not m:
        m2 = re.match(rf"^({NUM})_\{{?-?({NUM})\}}?\^\{{?\+?({NUM})\}}?$", c)
        if not m2:
            raise ValueError(f"cannot parse cell {cell!r}")
        return float(m2.group(1)), float(m2.group(3)), float(m2.group(2))
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


def host_name(raw):
    m = re.match(r"^(HD|HIP)(\d+)([A-Z]?)$", raw)
    if not m:
        raise ValueError(raw)
    return f"{m.group(1)} {m.group(2)}" + (f" {m.group(3)}" if m.group(3) else "")


def main():
    tex = TEX.read_text()
    block = tex.split("\\label{tab:sample}")[1].split("\\end{longtable}")[0]
    rows = []
    for line in block.splitlines():
        line = line.strip()
        if not re.match(r"^(HD|HIP)\d", line):
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        assert len(cells) == 10, line
        host = host_name(cells[0])
        m, mp, mm = parse_cell(cells[3])
        e, _, _ = parse_cell(cells[5])
        a, _, _ = parse_cell(cells[6])
        p_yr, pp, pm = parse_cell(cells[7])
        inc, _, _ = parse_cell(cells[8])
        msini = m * math.sin(math.radians(inc))
        p_d = p_yr * DAYS_PER_YR
        # The table rounds P errors to 0.01 yr; a printed 0.00 means < 0.005 yr, not zero.
        e_p_d = 0.5 * (pp + pm) * DAYS_PER_YR if (pp + pm) > 0 else None
        in_range = (msini <= 80.0) and (a <= 5.0 or p_d <= 4000.0)
        notes = (
            f"Cui+2026 Table S-catalog (tab:sample): true mass m={m:g} (+{mp:g}/-{mm:g}) MJ, "
            f"I={inc:g} deg, Omega={cells[9]}; msini_mj = m*sin(I) computed here; "
            f"P converted from {p_yr:g} yr (table gives 0.01 yr precision, ~3.7 d); "
            + ("e_P_d = mean of +/- errors; " if e_p_d is not None else
               "P error printed as 0.00 yr (< 0.005 yr), so e_P_d left blank; ")
            + "no K, omega or Tp in the paper; "
            "orbits are posteriors of the joint RV+Hipparcos-Gaia fits of Feng et al. 2022 "
            "(ApJS 262, 21), not new fits -- see feng_2022; stellar mass from TIC v8.2 per the "
            "paper but not tabulated (only q), so mstar left blank; RA/Dec="
            f"{cells[1]},{cells[2]}"
        )
        rows.append({
            "host": host,
            "companion": f"{host} (Cui+2026, P={p_yr:g} yr)",
            "P_d": f"{p_d:.1f}",
            "e_P_d": f"{e_p_d:.1f}" if e_p_d is not None else "",
            "K_ms": "",
            "e": f"{e:g}",
            "omega_deg": "",
            "Tp": "",
            "tp_system": "",
            "msini_mj": f"{msini:.2f}",
            "a_au": f"{a:g}",
            "mstar_msun": "",
            "fit_type": "RV+astrometry",
            "table_ref": "arXiv:2603.01808 SI Table 'Catalog of companions' (tab:sample)",
            "notes": notes,
            "_in_range": in_range,
        })
    assert len(rows) == 55, len(rows)

    with open(HERE / "cui_2026_all.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS + ["in_range"])
        w.writeheader()
        for r in rows:
            out = {k: r[k] for k in COLS}
            out["in_range"] = r["_in_range"]
            w.writerow(out)
    keep = [r for r in rows if r["_in_range"]]
    with open(HERE / "cui_2026_orbits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in keep:
            w.writerow({k: r[k] for k in COLS})
    nbd = sum(13 <= float(r["msini_mj"]) <= 80 for r in keep)
    print(f"{len(rows)} rows parsed; {len(keep)} in range; {nbd} with 13 <= msini <= 80 MJ")


if __name__ == "__main__":
    main()
