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

# Loop through each chromosome
# for (chr_name in names(chromosome_lengths)) {
for (chr_name in "chr33") {
  # Generate random positions for the SNPs
  snp_positions <- sort(sample(1:chromosome_lengths[chr_name], num_snps))

  # Generate random nucleotide substitutions
  snp_substitutions <- sample(c("A", "T", "C", "G"), num_snps, replace = TRUE)

  # Generate GRanges object for this chromosome
  chr_snps <- GRanges(seqnames = chr_name,
                      ranges = IRanges(start = snp_positions, width = 1),
                      strand = "*")

  chr_snps$REF <- getSeq(bsgenome, names = rep(chr_name, length(snp_positions)), start = snp_positions, end = snp_positions)

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

writeVcf(vcf, here::here("inst/extdata/galGal6-chr38.vcf.gz"))



# Define the new sample name
new_sample <- "sample1"

# Create a new DataFrame for the new sample
new_sample_df <- DataFrame(sampleID = new_sample)

# Add the new sample to colData(vcf)
colData(vcf) <- new_sample_df

# Check if the new sample has been added
colData(vcf)


gt <- matrix(sample(c("0/1", "1/1"), size=length(vcf), replace=TRUE))
geno(vcf)$GT <- gt
