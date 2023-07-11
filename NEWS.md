# tfboot 0.2.1

- Use simulated chicken in all `@examples` where previously used "Rosie" requiring a nonstandard package to run vignette examples.
- Rewrite `read_vcf()` function to no longer use `motifbreakR::snps.from.file()` and use `VariantAnnotation::readVcf()` instead, followed by some post-processing. 
- Move motifbreakR to `Suggests` avoiding the hard dependency, since motifbreakR isn't technically a requirement to install and use tfboot. 
- Added a changelog / NEWS.md.

# tfboot 0.2.0

Use chicken in vignette.

# tfboot 0.1.0

Initial version.
