#' Title
#'
#' @param gr A GenomicRanges object
#' @param width Width of the upstream interval
#'
#' @return fixme
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
#' get_upstream_interval(gr)
get_upstream_interval <- function(gr, width=5000L) {
  gr |>
    plyranges::mutate(ranges_original=GenomicRanges::ranges(gr)) |>
    GenomicRanges::resize(width=1) |>
    plyranges::anchor_3p() |>
    plyranges::stretch(width-1) |>
    plyranges::shift_left(1)
}
