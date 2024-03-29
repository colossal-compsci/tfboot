library(BSgenome.Ggallus.UCSC.galGal6)
library(TxDb.Ggallus.UCSC.galGal6.refGene)
library(VariantAnnotation)

set.seed(2023-06-19)

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
num_snps <- 20000  # Change as needed

# Initialize an empty GRanges object to store your simulated SNPs
simulated_snps <- GRanges()

# Loop through each chromosome
for (chr_name in names(chromosome_lengths)) {

  # skip everything except chr33
  if (chr_name!="chr33") next

  # Generate random positions for the SNPs
  snp_positions <- sort(sample(10000:(chromosome_lengths[chr_name]-10000), num_snps))

  # Generate GRanges object for this chromosome
  chr_snps <- GRanges(seqnames = chr_name,
                      ranges = IRanges(start = snp_positions, width = 1),
                      strand = "*")

  # Add the simulated SNPs for this chromosome to the total
  simulated_snps <- c(simulated_snps, chr_snps)
}

# Make the VCF
vcf <- VCF(rowRanges = simulated_snps,
           colData=DataFrame(Samples=1L, row.names = "sample1"),
           exptData=list(header=VCFHeader()))

# Set ref alleles (DNAStringSet)
ref(vcf) <- getSeq(bsgenome, simulated_snps)

# Set ref alleles (DNAStringSetList)
alt(vcf) <-
  ref(vcf) |>
  as.character() |>
  sapply(makesnp) |>
  unname() |>
  VariantAnnotation:::.toDNAStringSetList()

# Set the genotypes
assays(vcf, withDimnames=FALSE) <- list(GT=matrix(rep("1/1", length(vcf)), ncol=1, dimnames = list(NULL, "sample1")))


# Write the VCF
vcffile <- here::here("inst/extdata/galGal6-chr33.vcf")
writeVcf(vcf, vcffile)
system(paste0("head -n10 ", vcffile))

# Make it spec compliant
system(paste0("gsed -i '1 i##fileformat=VCFv4.3' ", vcffile))
system(paste0("bgzip -f ", vcffile))
system(paste0("tabix -f ", vcffile, ".gz"))
system(paste0("bcftools view ", vcffile, ".gz | head -n10"))

