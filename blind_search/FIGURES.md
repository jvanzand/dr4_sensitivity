# Figures

Only figures that belong to core usage are kept in git. Testing and comparison figures
are written to the gitignored `results/` folders and can be redrawn from saved files.

## Figure code

| Code | What it draws | Output |
|---|---|---|
| `plots.plot_step_sim`, called by `search.run_search(out_dir=..., label=...)` | one search step: the real ΔBIC periodogram with example noise-sampled periodograms and the simulated threshold; the histogram of the noise-sampled datasets' highest peaks with the Gumbel fit; the residuals before and after adding the detected companion | `<out_dir>/<label>_step<k>.png` |
| `plots.redraw_step` | the same figure, redrawn from a step saved by `search.run_search(save_dir=...)` (no search) | any path |
| `figtools.legend_clear`, `figtools.check_legends` | legends placed where they cover no data, and the check run before every save | (helpers) |
| `injection/injection_map.py plot` (`plot_map`) | the injection/recovery map: every injection with smoothed completeness contours | `injection/figures/<results-file-stem>.png` |
| `validation/synthetic_dr4.py` (`process_figure`) | how a benchmark dataset is made: RVs and reference fit, ground truth, zero-planet model, injected data | `validation/results/figures/synthetic_<system>_seed<k>.png` |
| `validation/benchmark.py run` / `redraw` | the step figures of every benchmark dataset | `validation/results/<system>_seed<k>/` |
| `docs/make_flowchart_simple.py` | the simple algorithm flowchart used in the README | `docs/search_flowchart.png` |
| `docs/make_flowchart.py` | the detailed developer flowchart with the mathematics (needs TeX) | `docs/search_flowchart_detailed.png` |

## Figure files in git

- `docs/gaia4_real_step1.png`, `gaia4_real_step2.png`: the step figure on the real
  Gaia-4 data (step 1 finds Gaia-4 b; step 2 finds nothing).
- `docs/search_flowchart.png` (simple), `docs/search_flowchart_detailed.png` (developer).
- `../injection/figures/injection_map_M50.png`, `injection_map_M100.png`,
  `injection_map_efap.png`: the 2000-injection maps (simulated threshold with 50 and 100
  noise-sampled datasets; eFAP threshold, same injections).

Removed code and figures are kept in the backup archives next to the repository
(`gaia_dr4_sensitivity_blind_search_backup_2026-11-13.tar.gz`,
`gaia_dr4_sensitivity_backup_2026-11-16.tar.gz`).
