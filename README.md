
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tfboot

<!-- badges: start -->
<!-- badges: end -->

The goal of tfboot is to facilitate statistical analysis of SNPs
disrupting transcription factor binding sites (TFBS) using bootstrapping
resampling to create empirical null distributions.

## Installation

You can install the development version of tfboot from
[GitHub](https://github.com/) with the code below. You’ll need to
replace `gh_pat` with your [personal access
token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token).

``` r
# First install devtools if you don't have it already
# install.packages("devtools")

# Set your personal access token
gh_pat="12345-abcde-CHANGEME-67890-xyz"

# Install directly from GitHub:
devtools::install_github(repo = "https://github.com/colossal-compsci/tfboot", 
                         auth_token = gh_pat, 
                         upgrade = FALSE, 
                         build_vignettes = TRUE, 
                         force = TRUE, 
                         dependencies = c("Imports", "Suggests"))
```

## Documentation

See the full introductory vignette with:

``` r
vignette("intro", package="tfboot")
```

## Example

Let’s use an example from the vignette. First, let’s load some pre-baked
motifbreakR results. `mbres` is a set of motifbreakR results run on SNPs
in the 5kb promoter region of a random selection of 5 genes. `mball` is
the precomputed set of motifbreakR results run on SNPs in the promoter
region of *all* genes.

``` r
library(tfboot)
```

``` r
mbres <- vignettedata$mbres
mbres
#> # A tibble: 449 × 10
#>    gene_id  SNP_id       tf    pctRef pctAlt scoreRef scoreAlt effect alleleDiff
#>    <chr>    <chr>        <chr>  <dbl>  <dbl>    <dbl>    <dbl> <chr>       <dbl>
#>  1 ARHGAP30 38:21390808… ACE2   0.754  0.957     3.73     4.71 strong      0.980
#>  2 ARHGAP30 38:21394731… ASCL1  0.882  0.774     4.68     4.13 weak       -0.558
#>  3 ARHGAP30 38:21394731… ASCL1  0.857  0.726     5.75     4.91 strong     -0.846
#>  4 ARHGAP30 38:21390808… AT5G…  0.861  0.692     4.29     3.48 strong     -0.808
#>  5 ARHGAP30 38:21394731… Atoh1  0.997  0.848     6.19     5.29 strong     -0.901
#>  6 ARHGAP30 38:21390808… BZIP…  0.765  0.860     4.98     5.57 weak        0.591
#>  7 ARHGAP30 38:21394731… Bhlh…  0.746  0.853     4.98     5.68 weak        0.696
#>  8 ARHGAP30 38:21394731… ERF0…  0.741  0.927     3.47     4.28 strong      0.811
#>  9 ARHGAP30 38:21394731… ERF13  0.723  0.903     3.99     4.95 strong      0.966
#> 10 ARHGAP30 38:21394731… ERF7   0.738  0.942     2.72     3.40 weak        0.678
#> # ℹ 439 more rows
#> # ℹ 1 more variable: alleleEffectSize <dbl>
mball <- vignettedata$mball
mball
#> # A tibble: 34,053 × 10
#>    gene_id SNP_id        tf    pctRef pctAlt scoreRef scoreAlt effect alleleDiff
#>    <chr>   <chr>         <chr>  <dbl>  <dbl>    <dbl>    <dbl> <chr>       <dbl>
#>  1 ACKR1   38:22769951:… ABF3   0.920  0.733     4.88     3.94 strong     -0.947
#>  2 ACKR1   38:22769951:… ARG80  0.910  0.727     3.93     3.14 strong     -0.788
#>  3 ACKR1   38:22769951:… CAMT…  0.916  0.703     4.09     3.17 strong     -0.921
#>  4 ACKR1   38:22769951:… CAMT…  0.940  0.741     4.27     3.38 strong     -0.893
#>  5 ACKR1   38:22769951:… HES1   0.863  0.682     4.90     3.90 strong     -0.997
#>  6 ACKR1   38:22769951:… NHP10  0.881  0.716     5.41     4.43 strong     -0.981
#>  7 ACKR1   38:22769951:… OJ10…  0.851  0.701     5.57     4.59 strong     -0.974
#>  8 ACKR1   38:22769951:… SUT1   0.888  0.639     3.61     2.62 strong     -0.993
#>  9 ACKR1   38:22769951:… ZNF8…  0.840  0.896     7.97     8.47 weak        0.498
#> 10 ACKR1   38:22769951:… nuc-1  0.856  0.697     5.42     4.42 strong     -0.998
#> # ℹ 34,043 more rows
#> # ℹ 1 more variable: alleleEffectSize <dbl>
```

Summarize motifbreakR results on our 5 genes of interest. This shows us
the actual values for the number of SNPs in the upstream regions of
these genes, and summary statistics on the allele differences, effect
sizes, etc. See the vignette and `?mb_summarize` for details.

``` r
mbsmry <- mb_summarize(mbres)
mbsmry
#> # A tibble: 1 × 7
#>   ngenes nsnps nstrong alleleDiffAbsMean alleleDiffAbsSum alleleEffectSizeAbsM…¹
#>    <int> <int>   <int>             <dbl>            <dbl>                  <dbl>
#> 1      5   449     363             0.837             376.                  0.157
#> # ℹ abbreviated name: ¹​alleleEffectSizeAbsMean
#> # ℹ 1 more variable: alleleEffectSizeAbsSum <dbl>
```

Bootstrap resample motifbreakR results for all genes. Resample sets of 5
genes 250 times.

``` r
set.seed(42)
mbboot <- mb_bootstrap(mball, ngenes=5, boots = 250)
mbboot$bootwide
#> # A tibble: 250 × 9
#>     boot genes           ngenes nsnps nstrong alleleDiffAbsMean alleleDiffAbsSum
#>    <int> <chr>            <int> <int>   <int>             <dbl>            <dbl>
#>  1     1 HSD17B7;SLAMF9…      5  1463    1100             0.808            1182.
#>  2     2 SLAMF7;LOC1110…      5   631     419             0.774             489.
#>  3     3 GPATCH2;CDK18;…      5   342     258             0.810             277.
#>  4     4 LOC111094614;L…      5   789     539             0.785             620.
#>  5     5 LOC111094522;L…      5   608     435             0.800             486.
#>  6     6 LOC111094613;S…      5   568     381             0.779             443.
#>  7     7 PIGM;LOC111094…      5   991     706             0.799             792.
#>  8     8 LEMD1;SH2D1B;F…      5   472     336             0.790             373.
#>  9     9 FCER1G;LOC1110…      5   685     490             0.802             549.
#> 10    10 EIF2D;LOC10215…      5   808     486             0.756             611.
#> # ℹ 240 more rows
#> # ℹ 2 more variables: alleleEffectSizeAbsMean <dbl>,
#> #   alleleEffectSizeAbsSum <dbl>
mbboot$bootdist
#> # A tibble: 6 × 2
#>   metric                  bootdist   
#>   <chr>                   <list>     
#> 1 alleleDiffAbsMean       <dbl [250]>
#> 2 alleleDiffAbsSum        <dbl [250]>
#> 3 alleleEffectSizeAbsMean <dbl [250]>
#> 4 alleleEffectSizeAbsSum  <dbl [250]>
#> 5 nsnps                   <dbl [250]>
#> 6 nstrong                 <dbl [250]>
```

Compare the values from our five genes of interest to the empirical null
distribution from bootstrap resampling.

``` r
bootstats <- mb_bootstats(mbsmry, mbboot)
bootstats
#> # A tibble: 6 × 5
#>   metric                     stat bootdist     bootmax     p
#>   <chr>                     <dbl> <list>         <dbl> <dbl>
#> 1 nsnps                   449     <dbl [250]> 1646     0.82 
#> 2 nstrong                 363     <dbl [250]> 1253     0.732
#> 3 alleleDiffAbsMean         0.837 <dbl [250]>    0.845 0.02 
#> 4 alleleDiffAbsSum        376.    <dbl [250]> 1347.    0.796
#> 5 alleleEffectSizeAbsMean   0.157 <dbl [250]>    0.173 0.556
#> 6 alleleEffectSizeAbsSum   70.7   <dbl [250]>  264.    0.816
```

Visualize the results:

``` r
plot_bootstats(bootstats)
```

<img src="man/figures/README-plot_bootstats-1.png" width="100%" />

See `vignette("intro", package="tfboot")` for more.
