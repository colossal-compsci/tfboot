# Get only SNPs on chr38 which are 0/1 or 1/1, and fix the reference allele

# Rosie
tabix ~/Downloads/RosieMorrill-PapillonDog.vcf.gz
bcftools view -t chr38 -i 'GT="alt"' ~/Downloads/RosieMorrill-PapillonDog.vcf.gz \
  | bcftools norm --check-ref s --fasta-ref ~/genomes/canFam3.fasta \
  | bcftools sort -Oz -o rosie38.vcf.gz \
  && tabix -f rosie38.vcf.gz \
  && bcftools index -s rosie38.vcf.gz

# Todd
tabix ~/Downloads/TodMorrill-PapillonDog.vcf.gz
bcftools view -t chr38 -i 'GT="alt"' ~/Downloads/TodMorrill-PapillonDog.vcf.gz \
  | bcftools norm --check-ref s --fasta-ref ~/genomes/canFam3.fasta \
  | bcftools sort -Oz -o todd38.vcf.gz \
  && tabix -f todd38.vcf.gz
