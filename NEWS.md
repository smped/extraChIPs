# Changes in version 0.99.0 (2022-03-28)

- Submitted to Bioconductor

# Changes in version 1.0.1 (2022-06-01)

- Bugfix so `as_tibble()` respects original column names

# Changes in version 1.0.2 (2022-07-22)

- Added `collapseTranscripts = "auto"` as the default for `plotHFGC()`
- `getProfileData()` now returns log2 transformed data by default 
- `getProfileData()` now uses `bplapply()` internally
- Bugfixes for `plotPie()`, `distinctMC()` `colToRanges()` and `stitchRanges()`
