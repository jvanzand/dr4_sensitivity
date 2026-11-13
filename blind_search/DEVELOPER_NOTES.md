# blind_search: developer notes

Details for the people who work on the code (Judah and Claude). The short, user-facing
description is in [README.md](README.md). Figure code and what was removed:
[FIGURES.md](FIGURES.md). Decisions and results, with dates: `../PLAN.md`, `../QUESTIONS.md`.

A standalone Python tool (no Octofitter, no MCMC) that decides, from a star's
Gaia DR4 epoch astrometry (and optionally RVs), how many companions the data
support and at which periods. It was built to answer: *what detection threshold
decides when to add a planet, so that we neither miss real signals nor accept weak
 ones?* DR4 is released on 2026-12-02; until then the tool runs on the
12-source pre-release epoch data in `data_conversion/` (Gaia-4, HD 114762) and
on simulated data with the same epochs and errors.

## Terminology (Judah, 2026-11-13; use these terms in all text)

- **Real data**: the measured transits of a star.
- **Model+noise dataset**: a dataset built as a known model plus noise, for
  example the zero-planet model (no signal) plus white noise at the formal
  errors, or a model with injected planets plus noise. It stands in for real
  data in tests. The reliability test T1 uses N = 45 of them; the injection map
  makes 2000 (one per injected planet).
- **Noise-sampled dataset**: one of the M datasets drawn from the model fitted
  to a given dataset (real or model+noise) plus noise, so any signal in it
  arises from noise only. The simulated threshold of that dataset is built from
  its M = 50 noise-sampled datasets. (Earlier notes, figures and result files
  call them "fake" or "simulated" datasets, "fakes".)
- **M** counts noise-sampled datasets per threshold; **N** counts model+noise
  datasets in a test. They are different numbers: a test with N model+noise
  datasets draws N x M noise-sampled datasets in total.
- **Simulated threshold**: the detection threshold built from the
  noise-sampled datasets (method `'sim'`); **eFAP**: RVSearch's empirical
  threshold from a dataset's own periodogram (method `'efap'`), which uses no
  noise-sampled datasets.

## How the search works

Each search step has a model with k - 1 companions already found (step 1: none)
and asks whether a k-th companion is warranted.

1. **Periodogram.** For every trial period P (20 d to twice the baseline; grid
   uniform in frequency, spacing 1/(2 pi T) for baseline T) fit the current
   model plus one *circular* companion. Because every amplitude is linear for
   fixed (P, e, M0), this is a small least-squares problem per period. Jitter
   is refitted at every period. The statistic is
   Delta-BIC(P) = Delta chi^2(P) - k_new ln N (`periodogram.py`).
2. **Threshold.** Draw M noise-sampled datasets from the current k - 1 companion model
   (its prediction at the real epochs plus Gaussian noise of standard deviation
   sqrt(err^2 + s^2), s the fitted jitter). Run each through the identical
   pipeline (re-polish the known companions, compute the periodogram) and keep
   only its highest peak. Fit a Gumbel (extreme-value) distribution to the M
   maxima; the threshold is its (1 - FAP) quantile (`simthreshold.py`). The
   real data's own periodogram is never used, so an undetected real signal
   cannot raise the threshold. (The alternative, RVSearch's empirical FAP in
   `efap.py`, fits the periodogram's own tail; a strong signal's leakage lifts
   that fit, so it misses companions.)
3. **Decision.** If the highest real peak exceeds the threshold, add a
   companion there and refit all companions' (P, e, M0) with free jitter
   (`fit.py`), then go to step k + 1. Otherwise stop (`search.py`).

## Modules

| file | what it does |
|---|---|
| `orbit.py` | Kepler's equation and the linear columns of an orbit (Thiele-Innes, RV) |
| `data.py` | data channels (Gaia DR4 transits, RVs) and `select` |
| `model.py` | the stacked weighted linear model (baseline, companion columns, jitter) |
| `periodogram.py` | the Delta-BIC periodogram and period grid |
| `efap.py` | RVSearch's empirical false-alarm threshold (kept for comparison) |
| `simthreshold.py` | the simulated threshold (noise-sampled datasets, Gumbel fit) |
| `fit.py` | refit of (P, e, M0) for all companions; add a companion at a peak |
| `search.py` | `run_search`: the iterative search (simulated threshold; eFAP optional) |
| `injection.py` | physical companions (masses, orbit angles) to linear amplitudes; Gaia-4 b |
| `plots.py` | `plot_step_sim`: the per-step figure written by `run_search(out_dir=...)` |
| `injection_map.py` | injection/recovery map on the Gaia-4 epochs (`run`, `plot`, `--method sim\|efap`) |
| `reliability_tests.py` | the reliability tests below |
| `make_flowchart_simple.py` | draws `search_flowchart.png`, the simple flowchart shown in the README |
| `make_flowchart.py` | draws `search_flowchart_detailed.png`, the detailed flowchart with the mathematics (needs TeX) |

Outputs: `results/` (JSON and JSON-lines results), `figures_search/` (the maps
and example step figures), `search_flowchart.png`, `search_flowchart_detailed.png`. [FIGURES.md](FIGURES.md)
lists the figure code and what was removed.

## Running

Use the project environment (`environment.yml`, Python 3.14) from the repo
root, with `python -I`:

```bash
python -I tests/test_blind_search.py --quick      # unit and integration tests (~1.5 min)
python -I blind_search/<module>.py --test         # every script has a quick test mode
python -I blind_search/reliability_tests.py       # the reliability tests (~8 min)
python -I blind_search/injection_map.py run --n 2000 --M 50 --out blind_search/results/injection_map_M50.jsonl
```

Compute (12 cores, 11 workers): one injection with M = 50 noise-sampled datasets takes about 33 s
of one core, M = 100 about 63 s, the eFAP rule about 1 s; a 2000-injection map
takes about 95 min, 3.3 h and 4 min respectively.

## What has been shown

- On 2000 injected planets (0.5-20 M_J, 0.5-5 AU, circular orbits, Gaia-4
  epochs and formal errors, FAP 1e-3), the simulated threshold detects 1351 and
  recovers 1092 (period and mass within 25%); eFAP detects 813 and recovers 725
  of the same injections. The 50% detection SNR is 6.4 against 16.3. M = 50 and
  M = 100 noise-sampled datasets give the same map within noise.
- Reliability tests: see the next section.

## Reliability tests

The suite (`python -I blind_search/reliability_tests.py`, 492 s on 11 workers,
12 of 12 checks passed on 2026-11-13) asks two questions: does the search rule
out noise, and does it find signals? In all of it the threshold is the
simulated one with M = 50 noise-sampled datasets and FAP 1e-3 unless stated.

| test | what is checked | result |
|---|---|---|
| T1 noise only | N = 45 model+noise datasets (zero-planet model plus white noise at the Gaia-4 epochs and errors); each gets its own threshold from M = 50 noise-sampled datasets; its highest peak is compared with it | 2 of 45 above the FAP 0.1 threshold (4.5 expected, 95% interval 0.5-15%); 0 of 45 above FAP 0.01 (0.45 expected; a weak check); p-values of the highest peaks uniform (KS statistic 0.17, p = 0.13) |
| T2a signals | detection versus SNR from the stored 2000-injection map | detected: 0% below SNR 2, 0.4% at SNR 2-4, 10% at 4-6, 65% at 6-8, 98% at 8-12, 100% above 12 (1020 of 1020); recovered (period and mass within 25%) 84-95% above SNR 12 |
| T2b regression | the first 5 injections of the map rerun with the current code | identical peaks, thresholds, detections and recoveries |
| T3 iterative search | 3 noise draws each of: no companion; Gaia-4 b; Gaia-4 b + a 200 d companion; a faint companion (SNR 3) | nothing found; b found then stop; both found then stop; nothing found; no spurious detection in any draw |
| T4 real Gaia-4 | the pre-release data | Gaia-4 b found (step 1: peak 181.3 against threshold 14.3), nothing at step 2 (peak -13.6 against 16.9) |

Reading the results honestly: with N = 45 the p-value test is the real check of
T1, and it shows no miscalibration, but it is weak (p = 0.13); an earlier run
with N = 100 (the 45 were its first 45, same seeds) gave 5 false alarms at FAP
0.1 (10 expected) and KS p = 0.61. Both runs fall below the expected rate, a
hint that the threshold is slightly conservative, not a significant one. FAP
1e-3 cannot be tested at this N. Among detected
signals above SNR 12, about 11% fail the 25% period-and-mass recovery window;
the reason has not been investigated. The detection curve is measured for
circular injections, one companion, and white noise only.

## Limitations (known, deliberate for now)

- The noise-sampled datasets carry white noise only; correlated noise would raise false alarms
  (the earlier correlated-noise experiment is described in `QUESTIONS.md`).
- Injections are circular; the refit of a detected peak is circular in the map.
  The periodogram template is circular; eccentricity enters only in the refit.
- The map has one companion per star; multi-companion behaviour is tested only
  in the reliability tests (T3).
- A companion absorbed by a poor first fit of an earlier one can go undetected
  (`QUESTIONS.md`, "Accepted limitation 2026-11-12").
- The BIC penalty counts 4 astrometric amplitudes plus the period; the 4 vs 5
  convention is unsettled. Whether the injected orbit orientation matches
  Octofitter's angle conventions has not been verified (test "V1"); sizes,
  periods and shapes do not depend on it.

Decisions, results and the glossary with dates live in `../PLAN.md` and
`../QUESTIONS.md`.
