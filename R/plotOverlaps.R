
plotOverlaps <- function(data, col = hcl.colors(3, "Zissou")) {

  ## Collect areas & totals, then separate
  areas <- .calculateAreas(data)
  totals <- areas[grepl("total", names(areas))]
  names(totals) <- gsub("_total", "", names(totals))
  areas <- areas[!grepl("total", names(areas))]

  ## Define Colours & Labels
  col <- rep(col, length(data))[seq_along(data)]
  nm <- names(data)
  if (length(nm) != length(data)) nm <- c("X", "Y", "Z")[seq_along(data)]

  ## Below is taken largely from Biovenn::draw.venn. By Tim Hulsen
  ## https://cran.r-project.org/web/packages/BioVenn/index.html

  ## Amplification
  total <- sum(unlist(areas))
  amp <- 1e5 / total
  grps <- names(areas)
  amps <- lapply(
    grps,
    function(x) {
      i <- grepl(x, names(areas))
      sum(unlist(areas[i])) * amp
    }
  )
  names(amps) <- grps

  single_grps <- intersect(c("x", "y", "z"), grps)
  pairwise_grps <- intersect(c("xy", "xz", "yz"), grps)
  three_way <- "xyz" %in% grps

  ## Radius calculations
  r <- lapply(amps[single_grps], function(x) .sqr(x/pi))
  names(r) <- single_grps

  ## Distances
  d <- lapply(
    pairwise_grps,
    function(x) {
      i <- gsub("^([xyz])([xyz])", "\\1", x)
      j <- gsub("^([xyz])([xyz])", "\\2", x)
      if (areas[[x]] > 0) {
        .pairwiseDistance(r[[i]], r[[j]], areas[[x]])
      }
    }
  )
  names(d) <- pairwise_grps
  d <- d[vapply(d, length, integer(1)) > 0]

  ## Dist for horizontal plots
  d <- sapply(
    names(d),
    function(x) {
      d1 <- d[[x]]
      d2 <- sum(unlist(d)) - d1
      if (d1 > d2) d1 <- d2
      d1
    },
    simplify = FALSE
  )

  browser()

  ## Angle calculations

  ## PPU

  ## Centres

  ## Intersection Points

  ## Special Cases with missing overlaps


  ## Get circle parameters

  ## Draw circles

  ## Add group labels

  ## Add counts

}
ex <- list(
  x = letters[1:10],
  y = letters[6:15],
  z = letters[10:26]
)
plotOverlaps(ex)

#' @importFrom methods is
.calculateAreas <- function(data) {

  ## Limit this to 3way interactions. Anything biggeris better as an UpSet plot
  stopifnot(is(data, "list"))
  stopifnot(length(data) <= 3)
  stopifnot(all(vapply(data, is.character, logical(1))))
  n <- length(data)
  areas <- lapply(data, length)
  names(areas) <- c("x_total", "y_total", "z_total")[seq_along(data)]
  if (n == 2) {
    areas$xy <- length(intersect(data[[1]], data[[2]]))
    areas$x <- areas$x_total - areas$xy
    areas$y <- areas$y_total - areas$xy
    i <- c("x_total", "y_total", "x", "y", "xy")
    areas <- areas[i]
  }
  if (n == 3) {
    areas$xyz <- length(intersect(intersect(data[[1]], data[[2]]), data[[3]]))
    areas$xy <- length(intersect(data[[1]], data[[2]])) - areas$xyz
    areas$xz <- length(intersect(data[[1]], data[[3]])) - areas$xyz
    areas$yz <- length(intersect(data[[2]], data[[3]])) - areas$xyz
    areas$x <- areas$x_total - areas$xy - areas$xz - areas$xyz
    areas$y <- areas$y_total - areas$xy - areas$yz - areas$xyz
    areas$z <- areas$z_total - areas$xz - areas$yz - areas$xyz
    i <- c("x_total", "y_total", "z_total", "x", "y", "z", "xy", "xz", "yz", "xyz")
    areas <- areas[i]
  }

  areas

}
ex <- list(
  x = letters[1:10],
  y = letters[6:15],
  z = letters[10:26]
)
.calculateAreas(ex)


.pairwiseDistance <- function(r1, r2, ol, .err = 1e-3) {
  .getVals <- function(r1, r2, d) {
    s1 <- .sq(r1) * .arccos((.sq(d) + .sq(r1) - .sq(r2))/(2 * d * r1))
    s2 <- .sq(r2) * .arccos((.sq(d) + .sq(r2) - .sq(r1))/(2 * d * r2))
    s3 <- 0.5 * .sqr(
      round(
        (d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (-d + r1 + r2), 5
      )
    )
    s1 + s2 - s3
  }
  d <- r1 + r2
  thresh <- ol * (1 - .err)
  while(
    thresh > .getVals(r1, r2, d)
  ) {
    d <- d - min(r1, r2) * .err;
  }
  d
}



#' Utility functions
.sq <- function(x) x^2
.sqr <- function(x) sqrt(round(abs(x)))
.arccos <- function(x, .digits = 5) {
  acos(
    min(
      max(-1, round(x, .digits)),
      1
    )
  )
}
