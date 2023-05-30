#' Read in SNPs from a VCF
#'
#' Helper function to read in data with [motifbreakR::snps.from.file].
#'
#' Note that the VCF **must** be filtered to only contain variant sites (i.e.,
#' no `0/0`), or only homozygous alt sites if you choose (`0/1` or `1/1`). This
#' can be accomplished with bcftools:
#'
#' `# Filter to any variant sites:` \cr
#' `bcftools view -i 'GT="alt"' ...` \cr
#' `# Filter to homozygous alt sites:` \cr
#' `bcftools view -i 'GT="AA"' ...`
#'
#' @param file File path to a VCF file. See details.
#' @param bsgenome An object of class BSgenome for the species you are interrogating; see [BSgenome::available.genomes] for a list of species.
#'
#' @return A GRanges object containing `SNP_id`, `REF`, and `ALT` columns.
#' @export
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Cfamiliaris.UCSC.canFam3)
#' vcf_file <- system.file("extdata", "rosie38.vcf.gz", package="tfboot", mustWork = TRUE)
#' snps <- read_vcf(vcf_file, BSgenome.Cfamiliaris.UCSC.canFam3)
#' snps
#' }
read_vcf <- function(file, bsgenome) {
  stopifnot(file.exists(file))
  stopifnot(inherits(bsgenome, "BSgenome"))
  motifbreakR::snps.from.file(file, format = "vcf", search.genome = bsgenome)
}

#' Get upstream intervals
#'
#' Get upstream intervals from a GRanges object.
#'
#' @param gr A GenomicRanges object.
#' @param width Width of the upstream interval. Default 5000.
#'
#' @return A GRanges object with intervals of width specified by `width` upstream of the start position of the `gr` input GRanges object.
#' @export
#'
#' @examples
#' gr <- data.frame(seqnames = "chr1",
#'                  strand = c("+", "-"),
#'                  start = c(42001, 67001),
#'                  width = c(1000, 1000),
#'                  gene = c("A", "B")) |>
#'   plyranges::as_granges()
#' gr
#' get_upstream(gr)
get_upstream <- function(gr, width=5000L) {
  stopifnot(inherits(gr, "GRanges"))
  upstream <-
    gr |>
    plyranges::mutate(ranges_original=GenomicRanges::ranges(gr)) |>
    GenomicRanges::resize(width=1) |>
    plyranges::anchor_3p() |>
    plyranges::stretch(width-1) |>
    plyranges::shift_left(1)
  # Set an attribute to check later if needed
  # attr(upstream, "upstream") <- TRUE
  return(upstream)
}

#' Get SNPs in upstream regions
#'
#' Get SNPs upstream of gene regions. Input SNPs read in with [read_vcf] and an appropriate TxDb object.
#'
#' @param snps SNPs read in with [read_vcf].
#' @param txdb A [GenomicFeatures::TxDb-class] object with transcript annotations for the organism of interest. Must match the organism specified in the `bsgenome` argument of [read_vcf].
#' @param level Currently only `"genes"` and `"transcripts"` are supported, which run [GenomicFeatures::genes] and [GenomicFeatures::transcripts], respectively.
#' @param ... Additional arguments passed to [get_upstream] (e.g., `width=`).
#'
#' @return A GRanges object containing SNPs in regions upstream of intervals in the specified TxDb. See [get_upstream].
#' @export
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Cfamiliaris.UCSC.canFam3)
#' library(TxDb.Cfamiliaris.UCSC.canFam3.ncbiRefSeq)
#' bs <- BSgenome.Cfamiliaris.UCSC.canFam3
#' tx <- TxDb.Cfamiliaris.UCSC.canFam3.ncbiRefSeq
#' vcf_file <- system.file("extdata", "rosie38.vcf.gz", package="tfboot", mustWork = TRUE)
#' snps <- read_vcf(vcf_file, BSgenome.Cfamiliaris.UCSC.canFam3)
#' snps
#' upstreamsnps <- get_upstream_snps(snps=snps, txdb=tx, level="genes")
#' upstreamsnps
#' }
get_upstream_snps <- function(snps, txdb, level="genes", ...) {
  # Check types
  stopifnot(inherits(snps, "GRanges"))
  stopifnot(inherits(txdb, "TxDb"))
  # Get features from the txdb
  if (level=="genes") {
    features <- suppressMessages(GenomicFeatures::genes(txdb))
  } else if (level=="transcripts") {
    features <- suppressMessages(GenomicFeatures::transcripts(txdb))
  } else {
    stop("level must be either 'genes' or 'transcripts'")
  }
  # Get upstream intervals of those features
  upstream <- get_upstream(features, ...)
  # Intersect with SNPs to get SNPs in those upstream regions
  upstreamsnps <- plyranges::join_overlap_intersect(snps, upstream)
  # Set an attribute to check later if needed
  # attr(upstream, "upstream") <- TRUE
  return(upstreamsnps)
}


#' Split GRanges by gene
#'
#' Splits a GRanges object into a GRangesList by a column (typically `gene_id`).
#' This function is deprecated and generally has no good use case. Originally
#' written to split up a GRanges object into a list to iterate over using
#' `furrr::future_map()`, but deprecated in favor of using built-in
#' parallelization in motifbreakR.
#'
#' @param gr A GRanges object returned by [get_upstream_snps].
#' @param key_col The name of the column in `gr` to split by (default `gene_id`).
#'
#' @return A list of genomic ranges split by `split_col`.
#' @export
#'
#' @examples
#' \dontrun{
#' gr <- data.frame(seqnames=rep(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'            start=1:10,
#'            width=1,
#'            gene_id = rep(c("gene1", "gene2", "gene3", "gene4"), c(4, 2, 1, 3))) |>
#'   plyranges::as_granges()
#' gr
#' split_gr_by_id(gr, key_col="gene_id")
#' }
split_gr_by_id <- function(gr, key_col="gene_id") {
  .Deprecated("motifbreakR(..., BPPARAM = bpparam())")
  stopifnot(inherits(gr, "GRanges"))
  stopifnot(key_col %in% names(GenomicRanges::mcols(gr)))
  grsplit <- split(gr, GenomicRanges::mcols(gr)[[key_col]])
  return(grsplit)
}

#' motifbreakR results to tibble
#'
#' Make a compact tibble with only select columns from motifbreakR results GRanges objects.
#'
#' It's a good idea to run motifbreakR _once_ on the background set of _all_ genes and save this as an RData (or .rds) file, and read in when you need them.
#'
#' @param mb motifbreakR results GRanges object.
#' @param key_col The name of the column used to key the txdb. Default `gene_id`. May be `transcript_id` or otherwise if you use a different value of `level` in [get_upstream_snps].
#'
#' @return A tibble containing the key column (usually `gene_id`), and a select number of other columns needed for downstream statistical analysis.
#' @export
#'
#' @examples
mb_to_tibble <- function(mb, key_col="gene_id") {
  stopifnot(inherits(mb, "GRanges"))
  standardcols <- c("SNP_id",
                    "geneSymbol",
                    "pctRef",
                    "pctAlt",
                    "scoreRef",
                    "scoreAlt",
                    "effect",
                    "alleleDiff",
                    "alleleEffectSize")
  mycols <- c(key_col, standardcols)
  stopifnot(all(mycols %in% colnames(GenomicRanges::mcols(mb))))
  tib <- tibble::as_tibble(GenomicRanges::mcols(mb)[,mycols])
  names(tib)[names(tib)=="geneSymbol"] <- "tf"
  tib <- tib |> dplyr::arrange(.data$gene_id, .data$tf)
  attr(tib, "mb_to_tibble") <- TRUE
  return(tib)
}

#' Summarize motifbreakR results
#'
#' Summarizes motifbreakR results as a tibble from [mb_to_tibble]. See details.
#'
#' Summarizes motifbreakR results. Returns a tibble with columns indicating:
#' 1. `ngenes`: The number of genes in the SNP set.
#' 1. `nsnps`: The number of SNPs total.
#' 1. `nstrong`: The number of SNPs with a "strong" effect.
#' 1. `alleleDiffAbsMean` The mean of the absolute values of the `alleleDiff` scores.
#' 1. `alleleDiffAbsSum` The sum of the absolute values of the `alleleDiff` scores.
#' 1. `alleleEffectSizeAbsMean` The mean of the absolute values of the `alleleEffectSize` scores.
#' 1. `alleleEffectSizeAbsSum` The sum of the absolute values of the `alleleEffectSize` scores.
#'
#' @param mbtibble motifbreakR results summarized with [mb_to_tibble].
#' @param key_col The name of the column used to key the txdb. Default `gene_id`. May be `transcript_id` or otherwise if you use a different value of `level` in [get_upstream_snps].
#'
#' @return A tibble. See description.
#' @export
#'
#' @examples
#' data(vignettedata)
#' vignettedata$mbres
#' mb_summarize(vignettedata$mbres)
mb_summarize <- function(mbtibble, key_col="gene_id") {
  # if (is.null(attributes(mbtibble)$mb_to_tibble)) warning("This function is usually run after mb_to_tibble()")
  mbsum <-
    mbtibble |>
    dplyr::summarize(ngenes=dplyr::n_distinct(.data[[key_col]]),
                     nsnps=dplyr::n(),
                     nstrong=sum(.data$effect=="strong"),
                     alleleDiffAbsMean=mean(abs(.data$alleleDiff)),
                     alleleDiffAbsSum=sum(abs(.data$alleleDiff)),
                     alleleEffectSizeAbsMean=mean(abs(.data$alleleEffectSize)),
                     alleleEffectSizeAbsSum=sum(abs(.data$alleleEffectSize)),
    )
  attr(mbsum, "mb_summarize") <- TRUE
  return(mbsum)
}

#' Bootstrap motifbreakR results
#'
#' Bootstrap motifbreakR results. Takes motifbreakR results as a tibble (from
#' [mb_to_tibble]), draws `boots` random samples of `ngenes` genes, and returns
#' as a list (1) a wide tibble with results for each bootstrap, and (2) another
#' tibble with the distribution of each metric in a listcol. See examples.
#'
#' Typically, you want to run this function on a full set of all genes to create
#' an empirical null distribution. You should run motifbreakR _once_ on all
#' genes, save this as an RData object, and read it in during bootstrap
#' resampling. See vignettes.
#'
#' @param mbtibble A tibble of motifbreakR results from [mb_to_tibble].
#' @param ngenes The number of genes to sample with each bootstrap resample.
#' @param boots The number of bootstrap resamples.
#' @param key_col The name of the column used to key the txdb. Default `gene_id`. May be `transcript_id` or otherwise if you use a different value of `level` in [get_upstream_snps].
#'
#' @return A list of two tibbles. See description and examples.
#' @export
#'
#' @examples
#' data(vignettedata)
#' vignettedata$mball
#' mb_bootstrap(vignettedata$mball, ngenes=5, boots=5)
mb_bootstrap <- function(mbtibble, ngenes, boots=100, key_col="gene_id") {
  names(mbtibble)[names(mbtibble)==key_col] <- ".key"
  allgenes <- unique(mbtibble$.key)
  bootwide <-
    tibble::tibble(boot=1:boots, ngenes=ngenes) |>
    dplyr::mutate(genes = purrr::map(ngenes, ~tibble::tibble(.key=sample(allgenes, size=., replace=FALSE)))) |>
    dplyr::mutate(mbres = purrr::map(.data$genes,  ~dplyr::inner_join(., mbtibble, by=".key", relationship='many-to-many'))) |>
    dplyr::mutate(mbsum = purrr::map(.data$mbres, mb_summarize, key_col=".key")) |>
    dplyr::mutate(genes = purrr::map_chr(.data$genes, function(x) paste(x[[1]], collapse=";"))) |>
    # dplyr::select(.data$boot, .data$genes, .data$mbsum) |>
    dplyr::select("boot", "genes", "mbsum") |>
    tidyr::unnest_wider(col="mbsum")
  attr(bootwide, "mb_bootstrap") <- TRUE
  bootdist <-
    bootwide |>
    # dplyr::select(-.data$boot, -.data$genes, -.data$ngenes) |>
    dplyr::select(-"boot", -"genes", -"ngenes") |>
    tidyr::gather("metric", "bootdist") |>
    dplyr::group_by(.data$metric) |>
    dplyr::summarize(bootdist=list(sort(c(bootdist))))
  attr(bootdist, "mb_bootstrap") <- TRUE
  return(tibble::lst(bootwide, bootdist))
}

#' Bootstrap statistics
#'
#' Get statistics (p-values) on your gene set's motifbreakR results compared to the bootstrapped empirical null distribution.
#'
#' @param mbsmry Results from running [mb_summarize] on motifbreakR results from your gene set of interest.
#' @param mbboot Results from running [mb_bootstrap] on a full background of all genes. Typically this should be performed once, with results read in from a file.
#'
#' @return A tibble with each metric from your bootstrap resampling (see [mb_summarize]), and p-value comparing the actual value of your gene set (`stat`) against the empirical null distribution (`bootdist`).
#' @export
#'
#' @examples
#' data(vignettedata)
#' mbres <- vignettedata$mbres
#' mball <- vignettedata$mball
#' mbsmry <- mb_summarize(mbres)
#' mbsmry
#' set.seed(42)
#' mbboot <- mb_bootstrap(mball, ngenes=5, boots = 100)
#' mbboot
#' mb_bootstats(mbsmry, mbboot)
mb_bootstats <- function(mbsmry, mbboot) {
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
    dplyr::mutate(bootmax=purrr::map_dbl(.data$bootdist, max)) |>
    dplyr::mutate(p=purrr::map2_dbl(.data$stat, .data$bootdist, function(x, y) sum(x<y)/length(y)))
  attr(bootstats, "mb_bootstats") <- TRUE
  return(bootstats)
}


#' Plot bootstrap distributions
#'
#' Plot bootstrap distributions of motifbreakR results with your critical value highlighted by a vertical red line.
#'
#' @param bootstats Output from running [mb_bootstats] on your results and a bootstrap resampling.
#'
#' @return A ggplot2 plot object. See description and examples.
#' @export
#'
#' @examples
#' data(vignettedata)
#' mbres <- vignettedata$mbres
#' mball <- vignettedata$mball
#' mbsmry <- mb_summarize(mbres)
#' set.seed(42)
#' mbboot <- mb_bootstrap(mball, ngenes=5, boots = 250)
#' bootstats <- mb_bootstats(mbsmry, mbboot)
#' plot_bootstats(bootstats)
plot_bootstats <- function(bootstats) {
  if (is.null(attributes(bootstats)$mb_bootstats)) stop("Call this function on the results of `mb_bootstats()`")
  bootstats$metric <- sprintf("%s\n(p=%.3f)", bootstats$metric, bootstats$p)
  bootstats |>
    tidyr::unnest(.data$bootdist) |>
    ggplot2::ggplot(ggplot2::aes(.data$bootdist)) +
    ggplot2::geom_histogram(bins=30) +
    ggplot2::facet_wrap(~metric, scales="free") +
    ggplot2::geom_vline(data=bootstats, ggplot2::aes(xintercept=.data$stat), col="red") +
    ggplot2::theme_bw()
}
