# Rules for duplicate RVs within a system

One exposure can reach us several times: from two databases (Jump and DACE), from one
database's own reduction and a publication's table, or from two reductions of the same
spectra. A fit must see each exposure once. These rules decide which copies to drop.
They are standard for every system; `merge_rvs.py` applies them, and the decisions they
encode were made on 2026-11-04.

## What can and can't be a copy

- **Only within an instrument family.** RVs from different instruments are never copies.
  A family groups every name and upgrade era of one spectrograph across sources (Jump's
  `hires_j`, `hires_k` and DACE's `HIRES` are all `hires`; DACE's `HARPS03`, `HARPS15`
  and Jump's `harps03` are all `harps`). The mapping is an explicit table in
  `merge_rvs.py`. An instrument missing from it becomes its own family, with a warning
  so that it gets added.
- **Never within one dataset.** Repeated exposures in one dataset (back-to-back shots)
  are real measurements and are all kept.
- **Exact repeats.** One exception: a row repeated in one dataset with the same time
  (to 1 s) and the same RV is a data-entry duplicate, and is dropped.

## Which copy is kept

Datasets in a family are ranked, and the higher-ranked copy is kept:
1. **Source:** Jump's CPS pipeline, then DACE, then Jump's imported datasets.
2. **Within a source:** a reduction of the spectra before a publication's table; for
   Jump's imports, the dataset Jump marks preferred first.
3. **Size:** the dataset with more RVs of the star, since one longer consistent series
   needs fewer offsets.
4. **Name,** alphabetically, so the order is reproducible.

## How a copy is recognized

Evidence, strongest first:
- **E1, same exposure ID.** Two RVs come from the same raw frame: the same timestamp in
  their raw file names (normalized to 1 s; e.g., DACE's `HARPN.2019-11-17T19:49:17.761`
  and Jump's imported `r.HARPN.2019-11-17T19-49-17.761_CCF_A.fits.gz`).
- **E2, same time.** Otherwise, two RVs are paired by time:
  - One to one. All candidate pairs between a lower-ranked dataset and the
    higher-ranked RVs are sorted by time difference and assigned nearest first, so no RV
    is paired twice by the same dataset.
  - Within a window that depends on the time stamps: **2 min** when both datasets have
    exposure-level times (a pipeline's times from the raw frames), **15 min** when either
    is a published or uploaded series. Published times can be defined differently from
    the reduction's: Butler+2017's HIRES times differ from Jump's by up to 9 min for the
    same exposures.
  - A pair with more than one candidate within the window is flagged **ambiguous** (a
    crowded night); it is still paired nearest first.

## Checking the pairs

True copies of one dataset by another differ by a constant zero-point offset. For each
pair of datasets:
- **Offset and scatter:** the median RV difference, and its robust scatter
  (1.4826 × median absolute deviation), compared with the expected scatter (the median
  of the two RVs' combined errors).
- **Rejected pairs:** a pair whose RV difference is more than 5 × max(scatter, expected)
  from the offset is not a copy. Both RVs are kept and the pair is flagged for review.
  This needs at least 3 pairs; with fewer, time-matched pairs are accepted but flagged
  **unvalidated** (exposure-ID pairs need no check).
- **Inconsistent datasets:** if the scatter exceeds 3 × the expected scatter, the
  dataset pair is flagged for review.
- **Remainders:** RVs of a lower-ranked dataset that copy nothing stay, as that dataset.
  If fewer than 5 remain, the dataset is flagged as a **small remainder**, with the
  offset measured from its pairs, so it could be tied to the higher-ranked dataset
  instead of getting its own offset.

## Records

- `rv_duplicates.csv`: each dropped RV, the evidence (E1, E2, or exact repeat), the RV
  it copies, their time and RV differences, the number of candidates, and flags.
- `rv_duplicate_pairs.csv`: per pair of datasets, the number of pairs, the fraction of
  the lower-ranked dataset they cover, offset, scatter, expected scatter, rejected pairs,
  and flags.
