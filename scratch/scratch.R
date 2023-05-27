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


# statsig -----------------------------------------------------------------

library(dplyr)
mbres <- tfboot::vignettedata$mbres
all <- tfboot::vignettedata$all
summarize_mb <- function(x) {
  x |>
    dplyr::summarize(nsnps=dplyr::n(),
                     nstrong=sum(effect=="strong"),
                     alleleDiffAbsMean=mean(abs(alleleDiff)),
                     alleleDiffAbsSum=sum(abs(alleleDiff)),
                     alleleEffectSizeAbsMean=mean(abs(alleleEffectSize)),
                     alleleEffectSizeAbsSum=sum(abs(alleleEffectSize)),
    )
}
summarize_mb(mbres)
ngenes <- length(unique(mbres$gene_id))
allgenes <- unique(all$gene_id)
boots=100
bootwide <-
  tibble::tibble(boot=1:1000, ngenes=ngenes) |>
  dplyr::mutate(genes = purrr::map(ngenes, ~tibble::tibble(gene_id=sample(allgenes, size=., replace=TRUE)))) |>
  dplyr::mutate(mbmbres = purrr::map(genes,  ~dplyr::inner_join(all, ., by="gene_id", relationship='many-to-many'))) |>
  dplyr::mutate(mbsum = purrr::map(mbmbres, summarize_mb)) |>
  dplyr::mutate(genes = purrr::map_chr(genes, function(x) paste(x[[1]], collapse=";"))) |>
  tidyr::unnest_wider(col=mbsum)
bootwide
bootnest <-
  bootwide |>
  dplyr::select(-ngenes, -genes, -mbres) |>
  tidyr::gather("metric", "bootdist", -boot) |>
  group_by(metric) |>
  summarize(bootdist=list(sort(c(bootdist))))
bootnest
mbres
summarize_mb(mbres) |>
  tidyr::gather("metric", "stat") |>
  dplyr::inner_join(bootnest, by="metric") |>
  mutate(max=purrr::map_dbl(bootdist, max)) |>
  mutate(p=purrr::map2_dbl(stat, bootdist, function(x, y) sum(x<y)/length(y)))
