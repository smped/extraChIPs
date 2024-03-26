#' Try finding viable ways for defining NFRs within histone peaks.
#' The idea behind this is taken directly from HisTrader [https://www.biorxiv.org/content/10.1101/2020.03.12.989228v1.full]
#' using a slow & fast moving averages & taking the difference.
#' Effectively, this is a PPO/MACD-H from stock trading
#'
#' The basic steps need to be:
#' 1. Merge peaks within a certain distance.
#'     + Merging based on a typical NFR may be the wisest strategy
#' 2. Split the merged peaks using GPos
#' 3. Define the two MAs for each GPos element
#'     + Start with 100/25, or looks at HisTrader defaults
#' 4. Calculate the difference as PPO-H
#' 5. Look for regions where the PPOH drops below a certain threshold
#'     + Might need to set a minimum width for ranges to be retained
#'     + Also might need to set a minimum coverage
#' 6. Scale to minimum width ranges? Maybe...
#'
#' Notably, the required inputs are a set of peaks, and a coverage file.
#' The coverage file can be FE or SMPR from macs2 bdcmp or callpeak
#'

## Load some peaks
library(extraChIPs)
library(rtracklayer)
library(ggplot2)
library(tidyverse)
sq <- defineSeqinfo("GRCh37")
peaks <- importPeaks(
    "~/DRMCRL/GRAVI_testing/output/macs2/H3K27ac/H3K27ac_consensus_peaks.bed.gz",
    type = "bed", seqinfo = sq
) |>
    unlist() |>
    GenomicRanges::reduce(min.gapwidth = 300)
bw <- "~/DRMCRL/GRAVI_testing/output/macs2/H3K27ac/H3K27ac_E2DHT_merged_treat_pileup.bw" |>
    BigWigFile()

gr <- peaks[4]
cov <- import.bw(bw, which = gr)
scores_df <- tibble(
    score = Rle(cov$score, width(cov)) |> as.numeric(),
    start = start(GPos(cov))
)
a <- 150 # Nucleosomes are 147
b <- 25
scores_df %>%
    mutate(
        slow = c(
            rep(NA, a - 1), zoo::rollmean(score, a)
        ),
        fast = c(
            rep(NA, b - 1), zoo::rollmean(score, b)
        ),
        diff = fast / slow - 1
    ) %>%
    dplyr::filter(!is.na(diff)) %>%
    mutate(
        region = ifelse(diff < 0, "NFR", "NOR"),
        change = region != lag(region, 1),
        id = c(0, cumsum(change[-1]))
    ) %>%
    dplyr::mutate(
        region = ifelse(dplyr::n() < 100, "NOR", region), .by = id
    ) %>%
    dplyr::mutate(
        region = case_when(
            id == max(id) ~ "NOR",
            id == min(id) ~ "NOR",
            TRUE ~ region
        )
    ) %>%
    ggplot(aes(start, score)) +
    geom_line(linewidth = 1) +
    geom_line(aes(y = slow), colour = "blue") +
    geom_line(aes(y = fast), colour = "red") +
    geom_vline(aes(xintercept = start), data = . %>% dplyr::filter(region == "NFR"), alpha = 0.2) +
    theme_bw()

%>%
    pivot_longer(cols = all_of(c("score", "diff"))) %>%
    ggplot(aes(start, value)) +
    geom_line(linewidth = 1) +
    geom_line(
        aes(y = fast),
        data = . %>% dplyr::filter(name == "score"), colour = "red", alpha = 0.5
    ) +
    geom_line(
        aes(y = slow),
        data = . %>% dplyr::filter(name == "score"), colour = "blue", alpha = 0.58
    ) +
    geom_hline(
        aes(yintercept = y),
        data = tibble(name = "diff", y = 0)
    ) +
    geom_vline(
        aes(xintercept = start),
        data = . %>% dplyr::filter(centre), colour = "grey70", linetype = 2
    ) +
    geom_rect(
        aes(xmin = start - 73.5, xmax = start + 73.5, ymin = -Inf, ymax = Inf),
        data = . %>% dplyr::filter(centre), fill = "grey70", alpha = 0.4
    ) +
    facet_grid(rows = vars(name), scales = "free_y", ) +
    theme_bw()

