FUN <- function(x, ...) {
  require(BSgenome.Cfamiliaris.UCSC.canFam3)
  require(MotifDb)
  motifbreakR::motifbreakR(x, ...)
}
furrr::future_map(grl, function(x) FUN(x, pwmList=motifs, threshold=.9))
library(furrr)
furrr::future_map(grl, function(x) motifbreakR(x, pwmList=motifs, threshold=.9))


mb_pois_list <- parallel_motifbreakR(poi_list, pwmList=motifs, threshold=0.9)
parallel_motifbreakR <- function(grl, cpus=NULL, ...) {
  stopifnot(inherits(grl, "GRangesList"))
  grllist <- as.list(grl)
  maxcpus <- parallel::detectCores()-1
  if (is.null(cpus) || cpus>maxcpus) cpus <- maxcpus
  message(sprintf("Parallelizing using %s CPUs", cpus))
  future::plan(future::multisession, workers=cpus)
  furrr::future_map(grl, function(x) motifbreakR::motifbreakR(x, ...))
}

parallel_motifbreakR <- function(grl, cpus=NULL, ...) {
  stopifnot(inherits(grl, "list"))
  stopifnot(inherits(grl[[1]], "GRanges"))
  genome.package <- unique(unlist(lapply(grl, function(x) x@genome.package)))
  stopifnot(length(genome.package)==1L)
  maxcpus <- parallel::detectCores()-1
  if (is.null(cpus) || cpus>maxcpus) cpus <- maxcpus
  message(sprintf("Parallelizing using %s CPUs", cpus))
  f <- function(x, ...) {
    suppressMessages(loadNamespace(genome.package))
    suppressMessages(loadNamespace("MotifDb"))
    # suppressMessages(require(genome.package, character.only = TRUE))
    # suppressMessages(require(MotifDb))
    motifbreakR::motifbreakR(x, ...)
  }
  future::plan(future::multisession, workers=cpus)
  furrr::future_map(grl, function(x) f(x, ...))
}
mb_pois_list <- parallel_motifbreakR(poi_list, pwmList=motifs, threshold=0.9)


library(furrr)
plan(multisession, workers=5)
furrr::future_map(grl, function(x) motifbreakR::motifbreakR(x, pwmList=motifs, threshold=.9))

