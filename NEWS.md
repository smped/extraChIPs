# extraChIPs 1.17.2

## Bug Fixes

- Improved memory consumption for `dualFilter()`
- Deprecated `keep.totals` & set default bin size to be 10kb


# extraChIPs 1.17.1

## Bug Fixes

- Better handled label background colours for plotSplitDonut
- Updated `fitAssayDiff()` for compatability with edgeR v4.0.0

# extraChIPs 1.12.1

## Bug Fixes

- Patched for ggplot2 v4.0.0
- Switched to `SimpleUpset` for Upset plots

# extraChIPs 1.11.2

## Changes

- Changed handling of arguments in plotting functions, setting NULL as the primary default value

# extraChIPs 1.11.0

Bioconductor 3.20 release

# extraChIPs 1.9.6

## New Features

- Added `centrePeaks()` to recentre peaks using any files with coverage

# extraChIPs 1.7.7

## Improvements

- Added `merge_within` to `makeConsensus()` for better handling when `method = "coverage"`

# extraChIPs 1.7.6

## New Features

- Added the DESeq2 Wald statistic to options for `fitAssayDiff()`

# extraChIPs 1.5.14

## Improvements

- Added `p_mu0` to output of `fitAssayDiff()`
- Added `respectLevels` and filtering to `plotProfileHeatmap()`

# extraChIPs 1.5.13

## Improvements

- Added handling of unquoted column names to most plotting functions
- Added passing of specific columns to `dualFilter()`
- Added `drop` argument to `addDiffStatus()`

# extraChIPs 1.5.12

## New Features

- Added `plotGrlCol()`

# extraChIPs 1.5.11

## New Features

- Added `defineSeqinfo()`

# extraChIPs 1.5.10

## Improvements

- Added bed format to `importPeaks()`

# extraChIPs 1.5.8

## Improvements

- Added control of side-axis label position for `plotProfileHeatmap()`
- Added option to return merged key-value ranges for `mergeByHMP()`

# extraChIPs 1.5.7

## Improvements

- Matched DiffBind and csaw settings for `fitAssayDiff()`
- Added `min_win` to all merging functions
- Added `n_max` to `getProfileData()`

# extraChIPs 1.5.6

## New Features

- Added `mapGrlCols()`

# extraChIPs 1.5.5

## New Features

- Added `plotPairwise()` and `addDiffStatus()`

# extraChIPs 1.4.2

## New Features

- Added fixed-width vignette and edited sliding window vignette
- Added `se` and `peaks` as example data for man pages and vignettes
- Added `defineRegions()`

## Improvements

- Enabled `plotHFGC()` without Ideogram tracks
- Enabled use of offsets for normalisation in `fitAssayDiff()`

# extraChIPs 1.3.9

## Improvements

- Changed labelling strategy for `plotPie()`

# extraChIPs 1.3.8

## Improvements

- Added coverage option to `makeConsensus()`

# extraChIPs 1.3.7

## New Features

- Added `plotAssayHeatmap()`
- Added `fitAssayDiff()` and coercion of `TopTags` objects

# extraChIPs 1.3.6

## New Features

- Added `mergeByHMP()` for merging overlapping windows using the harmonic mean p-value

## Improvements

- Expanded arguments for `plotSplitDonut()`

# extraChIPs 1.3.5

## New Features

- Added `plotSplitDonut()`

## Bug Fixes

- Fixed bug in `plotAssayDensities()` and `plotAssayRle()` and enabled plotting by group

# extraChIPs 1.3.4

## New Features

- Added `mergeBySig()`

# extraChIPs 1.1.5

## New Features

- Added `makeConsensus()` and updated vignette

# extraChIPs 1.1.3

## New Features

- Added `plotOverlaps()` for generation of Venn Diagrams and ComplexUpset plots

# extraChIPs 1.1.2

## Improvements

- Added `collapseTranscripts = "auto"` as default for `plotHFGC()`
- `getProfileData()` now returns log2 transformed data by default
- `getProfileData()` now uses `bplapply()` internally

## Bug Fixes

- Fixed bugs in `plotPie()`, `distinctMC()`, `colToRanges()`, and `stitchRanges()`

# extraChIPs 1.1.1

## Bug Fixes

- Fixed `as_tibble()` to respect original column names

# extraChIPs 0.99.0

## Major Changes

- Submitted to Bioconductor
