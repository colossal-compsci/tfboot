library(BSgenome.Ggallus.UCSC.galGal6)
library(TxDb.Ggallus.UCSC.galGal6.refGene)
library(org.Gg.eg.db)
library(VariantAnnotation)

makesnp <- function(x) {
  if (x=="N") return("N")
  nt <- c("A", "C", "G", "T")
  stopifnot(x %in% nt)
  nt <- setdiff(nt, x)
  sample(nt, 1)
}

# Set the BSGenome
bsgenome <- BSgenome.Ggallus.UCSC.galGal6

# Get the length of each chromosome from the BSGenome object
chromosome_lengths <- seqlengths(bsgenome)

# Define how many SNPs you want to generate per chromosome
num_snps <- 100  # Change as needed

# Initialize an empty GRanges object to store your simulated SNPs
simulated_snps <- GRanges()

chr_name <- "chr33"
# Loop through each chromosome
# for (chr_name in names(chromosome_lengths)) {
for (chr_name in "chr33") {
  # Generate random positions for the SNPs
  snp_positions <- sort(sample(1:chromosome_lengths[chr_name], num_snps))

  # Generate GRanges object for this chromosome
  chr_snps <- GRanges(seqnames = chr_name,
                      ranges = IRanges(start = snp_positions, width = 1),
                      strand = "*")

  chr_snps$REF <- getSeq(bsgenome,
                         names = rep(chr_name, length(snp_positions)),
                         start = snp_positions,
                         end = snp_positions)

  chr_snps$ALT <-
    chr_snps$REF |>
    unlist() |>
    as.character() |>
    strsplit(split="") |>
    unlist() |>
    sapply(makesnp) |>
    setNames(rep(chr_name, length(snp_positions))) |>
    DNAStringSet()

  # Add the simulated SNPs for this chromosome to the total
  simulated_snps <- c(simulated_snps, chr_snps)
}

vcf <- VCF(rowRanges = simulated_snps)
ref(vcf) <- simulated_snps$REF
alt(vcf) <- simulated_snps$ALT |> as.character() |> VariantAnnotation:::.toDNAStringSetList()

writeVcf(vcf, here::here("inst/extdata/galGal6-chr38.vcf"))
system(paste("head -n10", here::here("inst/extdata/galGal6-chr38.vcf")))
