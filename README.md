
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
#>    gene_id  SNP_id          tf        pctRef pctAlt scoreRef scoreAlt effect alleleDiff alleleEffectSize
#>    <chr>    <chr>           <chr>      <dbl>  <dbl>    <dbl>    <dbl> <chr>       <dbl>            <dbl>
#>  1 ARHGAP30 38:21390808:T:C ACE2       0.754  0.957     3.73     4.71 strong      0.980           0.199 
#>  2 ARHGAP30 38:21394731:T:C ASCL1      0.882  0.774     4.68     4.13 weak       -0.558          -0.105 
#>  3 ARHGAP30 38:21394731:T:C ASCL1      0.857  0.726     5.75     4.91 strong     -0.846          -0.127 
#>  4 ARHGAP30 38:21390808:T:C AT5G04390  0.861  0.692     4.29     3.48 strong     -0.808          -0.163 
#>  5 ARHGAP30 38:21394731:T:C Atoh1      0.997  0.848     6.19     5.29 strong     -0.901          -0.145 
#>  6 ARHGAP30 38:21390808:T:C BZIP52     0.765  0.860     4.98     5.57 weak        0.591           0.0918
#>  7 ARHGAP30 38:21394731:T:C Bhlha15    0.746  0.853     4.98     5.68 weak        0.696           0.105 
#>  8 ARHGAP30 38:21394731:T:C ERF096     0.741  0.927     3.47     4.28 strong      0.811           0.176 
#>  9 ARHGAP30 38:21394731:T:C ERF13      0.723  0.903     3.99     4.95 strong      0.966           0.177 
#> 10 ARHGAP30 38:21394731:T:C ERF7       0.738  0.942     2.72     3.40 weak        0.678           0.189 
#> # ℹ 439 more rows
mball <- vignettedata$mball
mball
#> # A tibble: 34,053 × 10
#>    gene_id SNP_id          tf           pctRef pctAlt scoreRef scoreAlt effect alleleDiff alleleEffectSize
#>    <chr>   <chr>           <chr>         <dbl>  <dbl>    <dbl>    <dbl> <chr>       <dbl>            <dbl>
#>  1 ACKR1   38:22769951:C:T ABF3          0.920  0.733     4.88     3.94 strong     -0.947          -0.179 
#>  2 ACKR1   38:22769951:C:T ARG80         0.910  0.727     3.93     3.14 strong     -0.788          -0.182 
#>  3 ACKR1   38:22769951:C:T CAMTA2        0.916  0.703     4.09     3.17 strong     -0.921          -0.207 
#>  4 ACKR1   38:22769951:C:T CAMTA3        0.940  0.741     4.27     3.38 strong     -0.893          -0.197 
#>  5 ACKR1   38:22769951:C:T HES1          0.863  0.682     4.90     3.90 strong     -0.997          -0.176 
#>  6 ACKR1   38:22769951:C:T NHP10         0.881  0.716     5.41     4.43 strong     -0.981          -0.160 
#>  7 ACKR1   38:22769951:C:T OJ1058_F05.8  0.851  0.701     5.57     4.59 strong     -0.974          -0.149 
#>  8 ACKR1   38:22769951:C:T SUT1          0.888  0.639     3.61     2.62 strong     -0.993          -0.245 
#>  9 ACKR1   38:22769951:C:T ZNF816        0.840  0.896     7.97     8.47 weak        0.498           0.0530
#> 10 ACKR1   38:22769951:C:T nuc-1         0.856  0.697     5.42     4.42 strong     -0.998          -0.158 
#> # ℹ 34,043 more rows
```

Summarize motifbreakR results on our 5 genes of interest. This shows us
the actual values for the number of SNPs in the upstream regions of
these genes, and summary statistics on the allele differences, effect
sizes, etc. See the vignette and `?mb_summarize` for details.

``` r
mbsmry <- mb_summarize(mbres)
mbsmry
#> # A tibble: 1 × 7
#>   ngenes nsnps nstrong alleleDiffAbsMean alleleDiffAbsSum alleleEffectSizeAbsMean alleleEffectSizeAbsSum
#>    <int> <int>   <int>             <dbl>            <dbl>                   <dbl>                  <dbl>
#> 1      5   449     363             0.837             376.                   0.157                   70.7
```

Bootstrap resample motifbreakR results for all genes. Resample sets of 5
genes 250 times.

``` r
set.seed(42)
mbboot <- mb_bootstrap(mball, ngenes=5, boots = 250)
mbboot$bootwide
#> # A tibble: 250 × 9
#>     boot genes                                                    ngenes nsnps nstrong alleleDiffAbsMean alleleDiffAbsSum alleleEffectSizeAbsMean alleleEffectSizeAbsSum
#>    <int> <chr>                                                     <int> <int>   <int>             <dbl>            <dbl>                   <dbl>                  <dbl>
#>  1     1 HSD17B7;SLAMF9;LOC102151674;LOC111094585;LOC102155451         5  1463    1100             0.808            1182.                   0.162                  236. 
#>  2     2 SLAMF7;LOC111094577;LOC111094537;HSD17B7;LOC111094548         5   631     419             0.774             489.                   0.159                  100. 
#>  3     3 GPATCH2;CDK18;LOC102153439;LOC111094501;LOC111094466          5   342     258             0.810             277.                   0.157                   53.6
#>  4     4 LOC111094614;LOC111094516;CD1E;LOC111094588;LOC111094527      5   789     539             0.785             620.                   0.160                  126. 
#>  5     5 LOC111094522;LOC111094553;FAM177B;LOC111094466;CNTN2          5   608     435             0.800             486.                   0.151                   91.8
#>  6     6 LOC111094613;SLAMF9;RXRG;LOC111094514;APOA2                   5   568     381             0.779             443.                   0.147                   83.4
#>  7     7 PIGM;LOC111094605;LOC111094471;LOC111094509;AIDA              5   991     706             0.799             792.                   0.159                  158. 
#>  8     8 LEMD1;SH2D1B;FCER1A;CDK18;LOC111094593                        5   472     336             0.790             373.                   0.159                   74.9
#>  9     9 FCER1G;LOC111094573;LOC111094582;LEMD1;LOC111094560           5   685     490             0.802             549.                   0.152                  104. 
#> 10    10 EIF2D;LOC102152982;PPP1R15B;LOC111094577;LOC111094514         5   808     486             0.756             611.                   0.148                  120. 
#> # ℹ 240 more rows
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
