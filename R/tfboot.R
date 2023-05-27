#' Read in SNPs from a VCF
#'
#' Helper function to read in data with [motifbreakR::snps.from.file]
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
#' @return FIXME
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
#' @return FIXME
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

#' Get upstream SNPs
#'
#' Get SNPs upstream of gene regions. Input SNPs read in with [read_vcf] and an appropriate TxDb object.
#'
#' @param snps FIXME
#' @param txdb FIXME
#' @param level FIXME
#' @param ... Additional arguments passed to [get_upstream].
#'
#' @return FIXME
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
  class(tib) <- append(class(tib), "mbtibble")
  return(tib)
}

#' Summarize motifbreakR results
#'
#' Summarizes motifbreakR results as a tibble from [mb_to_tibble]. See details.
#'
#' Summarizes motifbreakR results. Returns a tibble with columns indicating:
#' 1. `nsnps`: The number of SNPs total
#' 1. `nstrong`: The number of SNPs with a "strong" effect
#' 1. `alleleDiffAbsMean` The mean of the absolute values of the `alleleDiff` scores
#' 1. `alleleDiffAbsSum` The sum of the absolute values of the `alleleDiff` scores
#' 1. `alleleEffectSizeAbsMean` The mean of the absolute values of the `alleleEffectSize` scores
#' 1. `alleleEffectSizeAbsSum` The sum of the absolute values of the `alleleEffectSize` scores
#'
#' @param mbtibble motifbreakR results summarized with [mb_to_tibble].
#'
#' @return A tibble. See description.
#' @export
#'
#' @examples
summarize_mb <- function(mbtibble) {
  if (!inherits(mbtibble, "mbtibble")) warning("This function is usually called on the results of mb_to_tibble().")
  mbtibble |>
    dplyr::summarize(nsnps=dplyr::n(),
                     nstrong=sum(.data$effect=="strong"),
                     alleleDiffAbsMean=mean(abs(.data$alleleDiff)),
                     alleleDiffAbsSum=sum(abs(.data$alleleDiff)),
                     alleleEffectSizeAbsMean=mean(abs(.data$alleleEffectSize)),
                     alleleEffectSizeAbsSum=sum(abs(.data$alleleEffectSize)),
    )
}
