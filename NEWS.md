# tfboot 1.0.0

- Public release. Documentation built with pkgdown at https://colossal-compsci.github.io/tfboot/.
- Preprint on bioRxiv: Turner, S.D., et al. (2023). tfboot: Bootstrapping and statistical analysis for transcription factor binding site-disrupting variants in gene sets. bioRxiv 2023.07.14.549004. DOI: [10.1101/2023.07.14.549004](https://doi.org/10.1101/2023.07.14.549004).

# tfboot 0.2.1

- Use simulated chicken in all `@examples` where previously used "Rosie" requiring a nonstandard package to run vignette examples.
- Rewrite `read_vcf()` function to no longer use `motifbreakR::snps.from.file()` and use `VariantAnnotation::readVcf()` instead, followed by some post-processing. 
- Move motifbreakR to `Suggests` avoiding the hard dependency, since motifbreakR isn't technically a requirement to install and use tfboot. 
- Added a changelog / NEWS.md.

# tfboot 0.2.0

Use chicken in vignette.

# tfboot 0.1.0

Initial version.
