# plyranges vignette ------------------------------------------------------

suppressPackageStartupMessages(library(plyranges))
set.seed(100)
df <- data.frame(start=c(2:-1, 13:15),
                 width=c(0:3, 2:0))
df
rng <- df %>% as_iranges()
rng

grng <- df %>%
  transform(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
            strand = sample(c("+", "-"), 7, replace = TRUE),
            gc = runif(7)) %>%
  as_granges()
grng

rng <- as_iranges(data.frame(start=c(1, 2, 3), end=c(5, 2, 8)))
rng
grng <- as_granges(data.frame(start=c(1, 2, 3), end=c(5, 2, 8),
                              seqnames = "seq1",
                              strand = c("+", "*", "-")))
grng
mutate(rng, width = 10)
mutate(anchor_end(rng), width = 10)
grng
mutate(anchor_3p(grng), width = 10) # leave negative strand fixed
mutate(anchor_5p(grng), width = 10) # leave positive strand fixed

grng <- data.frame(seqnames = sample(c("chr1", "chr2"), 7, replace = TRUE),
                   strand = sample(c("+", "-"), 7, replace = TRUE),
                   gc = runif(7),
                   start = 1:7,
                   width = 10) %>%
  as_granges()
grng
grng_by_strand <- grng %>%
  group_by(strand)
grng_by_strand

grng %>% filter(gc < 0.3)

ir0 <- data.frame(start = c(5,10, 15,20), width = 5) %>%
  as_iranges()
ir1 <- data.frame(start = 2:6, width = 3:7) %>%
  as_iranges()
ir0
ir1

ir1 <- ir1 %>%
  mutate(gc = runif(length(.)))
ir1
ir0 %>%
  group_by_overlaps(ir1) %>%
  summarise(gc = mean(gc))

query <- data.frame(seqnames = "chr1",
                    strand = c("+", "-"),
                    start = c(1, 9),
                    end =  c(7, 10),
                    key.a = letters[1:2]) %>%
  as_granges()
query

subject <- data.frame(seqnames = "chr1",
                      strand = c("-", "+"),
                      start = c(2, 6),
                      end = c(4, 8),
                      key.b = LETTERS[1:2]) %>%
  as_granges()
subject


# nullranges --------------------------------------------------------------

library(nullranges)
library(GenomicRanges)
set.seed(1)
gr <- GRanges("chr1", IRanges(0:4 * 10 + 1, width=5), seqlengths=c(chr1=50))
gr
bootRanges(gr, R=1, blockLength=5, type="permute")
library(plyranges)
gr <- GRanges("chr1", IRanges(0:99 * 10 + 1, width=5), seqlengths=c(chr1=1000))
gr <- gr %>% mutate(x=rnorm(n()))

# Simple stat test
gr %>%
  sample(10, replace=TRUE) %>%
  summarize(sum(x))


# motifbreakr -------------------------------------------------------------

library(motifbreakR)
library(plyranges)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
pca.snps.file <- system.file("extdata", "pca.enhancer.snps", package = "motifbreakR")
pca.snps <- as.character(read.table(pca.snps.file)[,1])
variants <- snps.from.rsid(rsid = pca.snps,
                           dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19)
variants
motifbreakr.results <- motifbreakR(snpList = variants, pwmList = MotifDb, threshold = 0.9)
plotMB(results = motifbreakr.results, rsid = "rs7837328", effect = "strong")
library(BSgenome)
library(BSgenome.Cfamiliaris.UCSC.canFam3)
bsg <- BSgenome.Cfamiliaris.UCSC.canFam3
snps <- snps.from.file("inst/extdata/rosie.vcf.gz",
                       format="vcf",
                       search.genome = bsg)
snps
seqinfo(bsg)
bsg %>%
  seqinfo() %>%
  as("GRanges")
subseq(bsg$chr1, start=10225, end=10225)
getSeq(bsg, snps[1:10])


