# Changes in version 0.99.0 (2022-03-28)

- Submitted to Bioconductor

# Changes in version 1.1.1 (2022-06-01)

- Bugfix so `as_tibble()` respects original column names

# Changes in version 1.1.2 (2022-07-22)

- Added `collapseTranscripts = "auto"` as the default for `plotHFGC()`
- `getProfileData()` now returns log2 transformed data by default 
- `getProfileData()` now uses `bplapply()` internally
- Bugfixes for `plotPie()`, `distinctMC()` `colToRanges()` and `stitchRanges()`

# Changes in version 1.1.3 (2022-08-01)

- Added `plotOverlaps()` for generation of Venn Diagrams and ComplexUpset plots

# Changes in version 1.1.5 (2022-10-11)

- Added `makeConsensus()` and updated vignette

# Changes in version 1.2.4 (2023-01-31)

- Added `plotSplitDonut()`
- Significant bugfix & rewrite for `plotAssayDensities()`, `plotAssayPCA()` and `plotAssayRle()`
