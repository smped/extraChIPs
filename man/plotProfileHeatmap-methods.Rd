% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotProfileHeatmap.R
\name{plotProfileHeatmap}
\alias{plotProfileHeatmap}
\alias{plotProfileHeatmap,GenomicRangesList-method}
\alias{plotProfileHeatmap,GenomicRanges-method}
\title{Draw a coverage Profile Heatmap}
\usage{
plotProfileHeatmap(object, ...)

\S4method{plotProfileHeatmap}{GenomicRangesList}(
  object,
  profileCol = "profile_data",
  xValue = "bp",
  fillValue = "score",
  facetX = NULL,
  facetY = NULL,
  colour = facetY,
  linetype = NULL,
  summariseBy = c("mean", "median", "min", "max", "none"),
  xLab = xValue,
  yLab = NULL,
  fillLab = fillValue,
  labelFunX = waiver(),
  relHeight = 0.3,
  sortFilter = NULL,
  maxDist = 100,
  ...
)

\S4method{plotProfileHeatmap}{GenomicRanges}(
  object,
  profileCol = "profile_data",
  xValue = "bp",
  fillValue = "score",
  facetX = NULL,
  facetY = NULL,
  colour = facetY,
  linetype = NULL,
  summariseBy = c("mean", "median", "min", "max", "none"),
  xLab = xValue,
  yLab = NULL,
  fillLab = fillValue,
  labelFunX = waiver(),
  relHeight = 0.3,
  summaryLabelSide = "left",
  respectLevels = FALSE,
  sortFilter = NULL,
  maxDist = 100,
  ...
)
}
\arguments{
\item{object}{A GRanges or GRangesList object}

\item{...}{Passed to \link[ggplot2]{facet_grid} internally. Can be utilised
for switching panel strips or passing a labeller function}

\item{profileCol}{Column name specifying where to find the profile DataFrames}

\item{xValue, fillValue}{Columns within the profile DataFrames for heatmaps}

\item{facetX, facetY}{Columns used for faceting across the x- or y-axis
respectively}

\item{colour}{Column used for colouring lines in the summary panel. Defaults
to any column used for facetY}

\item{linetype}{Column used for linetypes in the summary panel}

\item{summariseBy}{Function for creating the summary plot in the top panel.
If set to 'none', no summary plot will be drawn. Otherwise the top panel will
contain a line-plot representing this summary value for each x-axis bin}

\item{xLab, yLab, fillLab}{Labels for plotting aesthetics. Can be overwritten
using labs() on any returned object}

\item{labelFunX}{Function for formatting x-axis labels}

\item{relHeight}{The relative height of the top summary panel.
Represents the fraction of the plotting area taken up by the summary panel.}

\item{sortFilter}{If calling on a GRangesList, a method for subsetting the
original object (e.g. 1:2). If calling on a GRanges object should be and
expression able to be parsed as a filtering expression using
\link[rlang]{eval_tidy}. This is applied when sorting the range order down
the heatmap such that ranges can be sorted by one or specific samples, or
all. Ranges will always be sorted such that those with the strongest signal
are at the top of the plot}

\item{maxDist}{Maximum distance from the centre to find the strongest signal
when arranging the ranges}

\item{summaryLabelSide}{Side to place y-axis for the summary plot in the top
panel}

\item{respectLevels}{logical(1) If FALSE, facets along the y-axis will be
arranged in descending order of signal, otherwise any original factor levels
will be retained}
}
\value{
A \code{ggplot2} object, able to be customised using standard \code{ggplot2} syntax
}
\description{
Plot a coverage Profile Heatmap across multiple ranges
}
\details{
Convenience function for plotting coverage heatmaps across a common set of
ranges, shared between one or more samples. These are most commonly the
coverage values from merged samples within a treatment group. THe input
data structure is based on that obtained from \link{getProfileData}, and can
be provided either as a GRanges object (generally for one sample) or as a
GRangesList.

A 'profile DataFrame' here refers to a data.frame (or tibble, or DataFrame)
with a coverage value in one column that corresponds to a genomic bin of a
fixed size denoted in another, as generated by \link{getProfileData}.
Given that multiple ranges are most likely to be drawn, each profile data
frame must be the same size in terms of the number of bins, each of which
represent a fixed number of nucleotides. At a minimum this is a two column
data frame although getProfileData will provide three columns for each
specified genomic region.

If using a GRangesList, each list element will be drawn as a separate panel
by default. These panels will appear in the same order as the list elements
of the GRangesList, although this can easily be overwritten by passing a
column name to the facetX argument. The default approach will add the
original element names as the column "name" which can be seen in the $data
element of any resultant ggplot object produced by this function.
}
\examples{
\donttest{
library(rtracklayer)
fl <- system.file(
"extdata", "bigwig", c("ex1.bw", "ex2.bw"), package = "extraChIPs"
)
bwfl <- BigWigFileList(fl)
names(bwfl) <- c("ex1", "ex2")

gr <- GRanges(
  c(
    "chr10:103880281-103880460", "chr10:103892581-103892760",
    "chr10:103877281-103877460"
   )
)
pd <- getProfileData(bwfl, gr)
plotProfileHeatmap(pd, "profile_data") +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  labs(fill = "Coverage")
}

}
