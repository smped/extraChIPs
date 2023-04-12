# Creation of data objects

## HiC Data

This was artificially generated using:

```r
library(GenomicInteractions)
ex_hic <- GInteractions(
  GRanges("chr10:103880000-103882500"), 
  GRanges("chr10:103892500-103895000")
)
ex_hic$distance <- calculateDistances(ex_hic)
```

## GRanges Objects

```r
library(rtracklayer)
library(plyranges)
library(extraChIPs)
gr <- GRanges("chr10:103860000-103920000")
gtf <- import.gff("path/to/file", type = "gene")
ids <- subsetByOverlaps(gtf, gr)$gene_id
ex_genes <- subset(gtf, gene_id %in% ids) %>% 
  select(gene = gene_id, symbol = gene_name) %>% 
  reduceMC()
```

```r
gtf <- import.gff("path/to/file", type = "exon")
ex_trans <- gtf %>%
  subset(gene_id %in% ids) %>% 
  as_tibble(rangeAsChar = FALSE) %>%
  group_by(transcript_id) %>%
  mutate(exon = paste0(transcript_id, "_", seq_along(transcript_id))) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = sq) %>%
  mutate(feature = as.character(type)) %>%
  select(type, gene = gene_id, exon, transcript = transcript_id, symbol = gene_name) %>%
  sort() %>%
  setNames(.$transcript)
```

```r
gtf <- import.gff("path/to/file", type = "transcript")
ex_prom <- gtf %>% 
  promoters(upstream = 1500, downstream = 500) %>% 
  reduce()
```

## Bam Files

These represent ER binding data from the HCI-005 PDX taken from 
[Hickey _et al_ (2015) The androgen receptor is a tumor suppressor in estrogen receptor-positive breast cancer _Nat Med_](https://doi.org/10.1038/s41591-020-01168-7),
along with the Input sample. Regions were extracted using 

```bash
samtools view -h filename.bam "chr10:103860000-103920000" > fileout.bam
samtools index fileout.bam
```

## BigWig Files

BigWig files were created by taking bedGraph files as output by `macs2 callpeak`,
converting to BigWig using the [GRAVI](https://github.com/smped/GRAVI) workflow,
loading into R then exporting using

```r
fl <- BigWigFile(file.path("path_to_complete_bigwig"))
bw <- import.bw(fl, which = gr)
export.bw(bw, "path/to/output.bw")
```

## Bed Files

This was taken as the merged sample `macs2 callpeak` output in `narrowPeak`
format and restricted to the above range

```r
peaks <- importPeaks("path/to/source.bed")
export.bed(subsetByOverlaps(peaks[[1]], gr), "path/to/output.bed.gz")
```

