"""Build wittenmyer_2009_orbits.csv from the arXiv LaTeX sources.

Sources (downloaded from https://arxiv.org/e-print/<id> into src_<id>/):
  - arXiv:0903.0652, Wittenmyer et al. 2009, ApJS 182, 97 ("A Search for
    Multi-Planet Systems Using the Hobby-Eberly Telescope"): Table 3
    "Keplerian Orbital Solutions" (label planetparms) and Table 1 "Stellar
    Parameters" (label stellarparms).
  - arXiv:0812.4094, Wittenmyer et al. 2009, AJ 137, 3529 ("HD 91669b: A New
    Brown Dwarf Candidate from the McDonald Observatory Planet Search"):
    Tables 1 and 2.

Run from this directory:  python -I build.py
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLS = ["host", "companion", "P_d", "e_P_d", "K_ms", "e", "omega_deg", "Tp",
        "tp_system", "msini_mj", "a_au", "mstar_msun", "fit_type", "table_ref",
        "notes"]


def table_block(tex, label):
    """Return the \\startdata...\\enddata body of the deluxetable with label."""
    i = tex.index("\\label{%s}" % label)
    j = tex.index("\\enddata", i)
    k = tex.find("\\startdata", i, j)  # label may sit in the caption, before \startdata
    start = k + len("\\startdata") if k >= 0 else i + len(label) + 8
    return tex[start:j]


def strip_comments(s):
    return "\n".join(line.split("%")[0] if not line.lstrip().startswith("%")
                     else "" for line in s.splitlines())


def rows_of(body):
    body = strip_comments(body)
    out = []
    for r in body.split("\\\\"):
        r = " ".join(r.split())
        if r:
            out.append([c.strip() for c in r.split("&")])
    return out


def num(s):
    s = re.sub(r"\\tablenotemark\{.*?\}", "", s).strip()
    return s


# ---------------------------------------------------------------- ApJS 182, 97
tex = (HERE / "src_0903.0652" / "ms.tex").read_text()

mstar = {}
for r in rows_of(table_block(tex, "stellarparms")):
    star = r[0]
    m = re.match(r"([\d.]+)", num(r[3]))
    mstar[star] = m.group(1)

rows = []
notes_by_planet = {
    "HD 74156 b": "Table 3 note a: results from two-planet fit. ",
    "47 UMa b": "Table 3 note b: results from one-planet fit (47 UMa c not fit). ",
    "HD 130322 b": "Table 3 note c: excludes Udry et al. (2000) data. ",
    "HD 128311 c": "T0=JD 2456987 lies after the data span (as printed). ",
}
for r in rows_of(table_block(tex, "planetparms")):
    name = num(r[0])
    host = name.rsplit(" ", 1)[0]
    # columns: P, eP, T0, eT0, e, ee, w, ew, K, eK, msini, emsini, a, ea, chi2, rms
    P, eP, T0, eT0, e, ee, w, ew, K, eK, ms, ems, a, ea = r[1:15]
    tp = "%.4f" % (float(T0) + 2400000.0)
    tp = tp.rstrip("0").rstrip(".") if "." in tp else tp
    notes = notes_by_planet.get(name, "")
    notes += ("Tp from Table 3 'T_0' column (JD-2400000; +2400000 added); the paper "
              "does not define T_0 explicitly, but it matches the periastron epoch "
              "adopted for these orbits elsewhere (e.g., exoplanets.org/Kiefer+2021 "
              "HD 136118 b); sigma_T0=%s d. GaussFit Keplerian fit to HET + McDonald "
              "2.7m + published RVs; mstar from Table 1." % eT0)
    if host == "HD 178911B":
        host_s = "HD 178911 B"
    else:
        host_s = host
    rows.append(dict(host=host_s, companion=name, P_d=P, e_P_d=eP, K_ms=K, e=e,
                     omega_deg=w, Tp=tp, tp_system="JD", msini_mj=ms, a_au=a,
                     mstar_msun=mstar[host], fit_type="RV-only",
                     table_ref="arXiv:0903.0652 Table 3 (Keplerian Orbital Solutions)",
                     notes=notes))

# ---------------------------------------------------------------- AJ 137, 3529
tex2 = (HERE / "src_0812.4094" / "ms.tex").read_text()
body = table_block(tex2, "tbl-2")
par = {}
for r in rows_of(body):
    par[r[0]] = r[1]


def first(s):
    return re.match(r"\s*(-?[\d.]+)", s).group(1)


def err(s):
    return re.search(r"\\pm\$?\s*([\d.]+)", s).group(1)


m1 = first(re.search(r"Mass & (.*?)\\\\", tex2).group(1))
T0 = float(first(par["$T_0$"])) + 2400000.0
rows.append(dict(host="HD 91669", companion="HD 91669 b",
                 P_d=first(par["Period"]), e_P_d=err(par["Period"]),
                 K_ms=first(par["$K$"]), e=first(par["$e$"]),
                 omega_deg=first(par["$\\omega$"]), Tp="%.1f" % T0, tp_system="JD",
                 msini_mj=first(par["M sin $i$"]), a_au=first(par["$a$"]),
                 mstar_msun=m1, fit_type="RV-only",
                 table_ref="arXiv:0812.4094 Table 2 (Wittenmyer+2009, AJ 137, 3529)",
                 notes=("Tp from Table 2 'T_0' (JD-2400000; +2400000 added; sigma 0.6 d); "
                        "not explicitly defined as periastron, but Unger+2023's "
                        "independent Tp (2457279) differs by 8.00 periods, consistent "
                        "with periastron. 18 McDonald 2.7m RVs only, 4.3 yr. "
                        "mstar 0.914 (+0.018/-0.087) from Table 1.")))

with open(HERE / "wittenmyer_2009_orbits.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rows written")
