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
mball <- tfboot::vignettedata$mball
attributes(mbres)$mb_to_tibble <- TRUE
attributes(mball)$mb_to_tibble <- TRUE

mbsmry <- mb_summarize(mbres)
mbboot <- mb_bootstrap(mball, ngenes=5)

if (is.null(attributes(mbsmry)$mb_summarize)) warning("mbsmry is usually the result of `mb_summarize()`")
if (is.null(attributes(mbboot$bootwide)$mb_bootstrap)) warning("mbboot is usually the result of `mb_summarize()`")
if (is.null(attributes(mbboot$bootdist)$mb_bootstrap)) warning("mbboot is usually the result of `mb_summarize()`")
if (length(unique(mbboot$bootwide$ngenes))!=1L) stop("You have different numbers of genes in some of your bootstrap resamples.")
if (unique(mbboot$bootwide$ngenes)!=mbsmry$ngenes) stop("You have different number of genes in mbsmry and mbboot.")
nboots <- unique(purrr::map_int(mbboot$bootdist$bootdist, length))
if (length(nboots)!=1L) stop("Something went wrong. You have different numbers of bootstraps for different genes.")
bootstats <-
  mbsmry |>
  tidyr::gather("metric", "stat") |>
  dplyr::inner_join(mbboot$bootdist, by="metric") |>
  dplyr::mutate(bootmax=purrr::map_dbl(bootdist, max)) |>
  dplyr::mutate(p=purrr::map2_dbl(stat, bootdist, function(x, y) sum(x<y)/length(y)))
bootstats |>
  tidyr::unnest(bootdist) |>
  ggplot::ggplot(ggplot::aes(bootdist)) +
  ggplot::geom_histogram(bins=30) +
  ggplot::facet_wrap(~metric, scales="free") +
  ggplot::geom_vline(data=bootstats, ggplot::aes(xintercept=stat), col="red") +
  ggplot::theme_bw()






mbboot |>
  dplyr::select(-boot, -genes, -ngenes) |>
  tidyr::gather("metric", "bootdist") |>
  dplyr::group_by(.data$metric) |>
  summarize(bootdist=list(sort(c(bootdist))))
bootnest
mbres
summarize_mb(mbres) |>
  tidyr::gather("metric", "stat") |>
  dplyr::inner_join(bootnest, by="metric") |>
  mutate(max=purrr::map_dbl(bootdist, max)) |>
  mutate(p=purrr::map2_dbl(stat, bootdist, function(x, y) sum(x<y)/length(y)))


mbtibble <- mball
boots=10
key_col="gene_id"
mb_bootstrap <- function(mbtibble, ngenes, boots=100, key_col="gene_id") {
  names(mbtibble)[names(mbtibble)==key_col] <- ".key"
  allgenes <- unique(mbtibble$.key)
  bootwide <-
    tibble::tibble(boot=1:boots, ngenes=ngenes) |>
    dplyr::mutate(genes = purrr::map(ngenes, ~tibble::tibble(.key=sample(allgenes, size=., replace=TRUE)))) |>
    dplyr::mutate(mbres = purrr::map(genes,  ~dplyr::inner_join(., mbtibble, by=".key", relationship='many-to-many'))) |>
    dplyr::mutate(mbsum = purrr::map(mbres, mb_summarize)) |>
    dplyr::mutate(genes = purrr::map_chr(genes, function(x) paste(x[[1]], collapse=";"))) |>
    tidyr::unnest_wider(col=mbsum)
  return(bootwide)
}

mb_bootstrap(mbtibble, ngenes=5, boots=100)

x <- tibble::tibble(boot=1:boots, ngenes=ngenes) |>
  dplyr::mutate(genes = purrr::map(ngenes, ~tibble::tibble(.key=sample(allgenes, size=., replace=TRUE)))) |>
  dplyr::mutate(mbres = purrr::map(genes,  ~dplyr::inner_join(., mbtibble, by=".key", relationship='many-to-many')))
x |>
  mutate(mbsum=map(mbres, mb_summarize)) |>
  dplyr::mutate(mbsum = purrr::map(mbres, mb_summarize)) |>
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
