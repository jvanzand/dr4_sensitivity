# validation: benchmarks and reliability tests

Checks that the blind search behaves as intended, on data whose answer is known.

## Benchmarks

Synthetic Gaia DR4 data of known RV systems. For each system, the planets' RV orbits
(period, eccentricity, phase, semi-amplitude) are drawn from a fit to the real RVs; the
inclination and node are drawn at random. Gaia's scanning law and error model
(`gaiamock`) give realistic epochs and noise. The datasets are frozen in `benchmarks/`,
so every version of the search runs on identical data.

![How a benchmark dataset is made](benchmark_example.png)

| system | what it tests |
|---|---|
| π Men | strong long-period planet; the two inner planets must not be claimed |
| HD 111232 | moderate signal |
| HD 168443 | strong outer planet; a weak inner planet as a stress test |
| HD 89707 | very high eccentricity (e = 0.9) |
| TOI-1736 | weak signal (optional) |

The list and the search settings are in `benchmarks.toml`.

```bash
python -I validation/benchmark.py generate    # make datasets listed but missing (never overwrites)
python -I validation/benchmark.py run         # search them; products in validation/results/ (not in git)
python -I validation/benchmark.py redraw      # redraw every figure from saved products, no search
python -I validation/benchmark.py compare     # table against the stored reference results
```

## Reliability tests

```bash
python -I validation/reliability_tests.py     # 12 checks, about 8 minutes
```

Noise alone exceeds the threshold at the stated rate; injected planets are found as
expected; the real Gaia-4 planet is found and nothing more.

| file | contents |
|---|---|
| `synthetic_dr4.py` | makes a benchmark dataset and its figure |
| `benchmark.py`, `benchmarks.toml` | runs, redraws and compares the benchmarks |
| `reliability_tests.py` | the reliability tests |
| `rv_hd191939.py` | the RV-only search of HD 191939 |
