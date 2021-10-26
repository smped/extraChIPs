# We need a combination of ranges which do & don't map
genes <- GRanges(c("chr1:1-10:*", "chr1:21-30:-", "chr1:31-40:+"))
genes$gene_id <- paste0("gene", seq_along(genes))
## Add a promoter for each gene
prom <- promoters(genes, upstream = 1, downstream = 1)
## Add an enhancer in gene2 & to the right of gene3
enh <- GRanges(c("chr1:25", "chr1:50"))
enh$which <- c("gene2", "gene1 (gi)")
## Map the last enhancer to gene1's promoter
gi <- InteractionSet::GInteractions(prom[1], enh[2])
mcols(gi) <- c()
## Now define key ranges
gr <- GRanges(paste0("chr1:", seq(0, 50, by = 10)))
gr$via_prom <- CharacterList(
  list("gene1", NULL, NULL, c("gene2", "gene3"), NULL, "gene1")
  )
gr$via_enh_gap0 <- CharacterList(
  list(NULL, NULL, NULL, NULL, NULL, NULL)
)
gr$via_enh_gap5 <- CharacterList(
  list(NULL, NULL, "gene2", "gene2", NULL, NULL)
)
gr$via_gi <- CharacterList(
  list("gene1", NULL, NULL, NULL, NULL, "gene1")
)

