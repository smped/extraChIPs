## Cytogenetic bands for GRCh37
url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
f <- tempfile()
download.file(url, f)
ln <- readLines(gzfile(f))
ln <- ln[grepl("chr[0-9XY]+\\t", ln)]
df <- do.call(rbind, strsplit(ln, split = "\t", fixed = TRUE))
df <- as.data.frame(df)
names(df) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
df$chromStart <- as.integer(df$chromStart)
df$chromEnd <- as.integer(df$chromEnd)
grch38.cytobands <- df
usethis::use_data(grch37.cytobands)

## Cytogenetic bands for GRCh38
url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
f <- tempfile()
download.file(url, f)
ln <- readLines(gzfile(f))
ln <- ln[grepl("chr[0-9XY]+\\t", ln)]
df <- do.call(rbind, strsplit(ln, split = "\t", fixed = TRUE))
df <- as.data.frame(df)
names(df) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
df$chromStart <- as.integer(df$chromStart)
df$chromEnd <- as.integer(df$chromEnd)
grch38.cytobands <- df
usethis::use_data(grch38.cytobands)
