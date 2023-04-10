#' Setup the new data for the differential binding vignette
library(tidyverse)
library(rtracklayer)
library(extraChIPs)

h3k <- read_rds("../GRAVI_MDA-MB-453_AR_GATA3_H3K27ac/output/differential_h3k27ac/H3K27ac/H3K27ac_Veh_DHT-filtered_windows.rds")
chr10 <- subset(h3k, seqnames == "chr10")
rowData(chr10) <- rowData(chr10)["overlaps_ref"]
assay(chr10, "logCPM") <- NULL
assay(chr10, "qsmooth") <- NULL
colnames(chr10) <- colData(chr10)$label
colData(chr10)$bam.files <- file.path("../data/bam", paste0(colnames(chr10), ".bam"))
colData(chr10)$input <- rep("../data/bam/input.bam", 6)
colData(chr10)$design <- NULL
write_rds(chr10[1:20000,], here::here("data/chr10.rds"), compress = "gz")
