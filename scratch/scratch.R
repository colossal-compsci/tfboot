suppressPackageStartupMessages({
  library(tidyverse)
  library(plyranges)
  library(BSgenome.Cfamiliaris.UCSC.canFam3)
  library(TxDb.Cfamiliaris.UCSC.canFam3.ncbiRefSeq)
  library(motifbreakR)
})

bs <- BSgenome.Cfamiliaris.UCSC.canFam3
bs

tx <- TxDb.Cfamiliaris.UCSC.canFam3.ncbiRefSeq
tx

# Extract the transcription start sites (TSS) by GENE.
# Note if you want to do this for every TRANSCRIPT, do somethign like this first:
# GenomicFeatures::transcriptsBy(tx, by="gene")
tss <-
  suppressMessages(genes(tx)) %>%
  filter(seqnames=="chr38") %>%
  arrange(start) %>%
  resize(width=1, fix="start")
tss

# Use plyranges to get some kb upstream
width=5000
upstream <-
  tss %>%
  anchor_5p() %>%
  stretch(width-1)

set.seed(1); upstream <- upstream %>% sample(20)

vcf_file <- system.file("extdata", "rosie38.vcf.gz", package="tfboot", mustWork = TRUE)

snps <- motifbreakR::snps.from.file(vcf_file, format="vcf", search.genome=bs)
upstreamsnps <- snps %>% join_overlap_intersect(upstream)
mb1 <- motifbreakR(snpList = upstreamsnps, pwmList = MotifDb, threshold = 0.9)
mb1

upstream |> split()

