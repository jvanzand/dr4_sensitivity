# Figures in blind_search (2026-11-13)

Only the figures that belong to core usage are kept. Testing, validation and
method-comparison figures were removed on 2026-11-13 (Judah's decision); new
testing figures will be added later and pruned again when the time is right.

## Core figure code

| Code | What it draws | Output |
|---|---|---|
| `plots.plot_step_sim`, called by `search.run_search(out_dir=..., label=...)` | one search step: the real Delta-BIC periodogram with example noise-sampled periodograms and the simulated threshold; the histogram of the noise-sampled datasets' highest peaks with the Gumbel fit; the residuals before and after adding the detected companion | `<out_dir>/<label>_step<k>.png` |
| `injection_map.py plot` (`plot_map`) | the injection/recovery map: every injection (recovered or not) with smoothed completeness contours, and the smoothed completeness image | `figures_search/<results-file-stem>.png` |
| `make_flowchart.py` | the algorithm flowchart (needs TeX) | `search_flowchart.png` |

## Figure files kept

- `figures_search/injection_map_M50.png`, `injection_map_M100.png`,
  `injection_map_efap.png`: the 2000-injection maps (simulated threshold with
  50 and 100 noise-sampled datasets; eFAP threshold, same injections).
- `figures_search/gaia4_real_step1.png`, `gaia4_real_step2.png`: the step
  figure on the real Gaia-4 data (step 1 finds Gaia-4 b; step 2 finds nothing).
- `search_flowchart.png`.

## Removed (copies of the code and figures are in
`../../gaia_dr4_sensitivity_blind_search_backup_2026-11-13.tar.gz`)

Code: `make_figures.py`, `make_sim_figures.py`, `make_method_diagram.py`,
`make_slide_figures.py`, `make_threshold_comparison_figure.py`,
`run_injection_tests.py`, `run_threshold_comparison.py`,
`injection_map_compare.py`, `injection_map_compare_methods.py`,
`plots.plot_step` (the eFAP step figure) and `plots.plot_injection`; the
figure in `reliability_tests.py`. Figure folders `figures/`, `figures_compare/`,
`figures_sim/`, `figures_slides/` and the other files of `figures_search/`.
The result files in `results/` (JSON and JSON-lines) were not touched.
`search.run_search` now writes step figures only for the simulated threshold,
and the options that existed only for the removed scripts (the fixed-floor
detection rule and the arguments `floor`, `truth_info`, `title_prefix`,
`resolution_tol`) are gone. The eFAP detection rule remains (`method='efap'`).
