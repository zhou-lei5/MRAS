
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRAS

<!-- badges: start -->
<!-- badges: end -->

MRAS is designed to identify crucial RNA-binding proteins (RBPs)
responsible for splicing variations in diverse scenarios, including
cancer vs. normal, primary vs. recurrence, and more, not just in bulk
data but also in single-cell data.

## Installation

You can install the development version of MRAS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zhou-lei5/MRAS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MRAS)
#> Loading required package: data.table
#> Loading required package: fgsea
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
## "expr" is RBP expression matrix.
data("expr")
expr[1:5,1:5]
#>        SRR5962198 SRR5962199  SRR5962200  SRR5962201 SRR5962202
#> A1CF    0.1005414  0.0000000  0.03379508  0.04432587  1.1422434
#> ANKHD1  0.8390107  0.7166435  0.77311709  0.77295568  2.2460320
#> CELF1  19.4532827 18.6930904 11.78120909 16.90546963 21.4689576
#> CELF2  29.7642191 25.5567231 11.75549788 26.64966474 16.2681018
#> CELF3   0.0000000  0.0840603  0.06387838  0.10054003  0.3084332
```

``` r
## "psi" is events psi matrix."expr" and "psi" should have same column names.
data("psi")
psi[1:5,1:5]
#>                                                                                               SRR5962198
#> ENSG00000147403.18_RPL10_chrX_+_154400463_154400626_154399802_154399941_154400701_154400811        0.984
#> ENSG00000196924.19_FLNA_chrX_-_154357250_154357274_154354824_154355072_154357433_154357623         0.407
#> ENSG00000101972.20_STAG2_chrX_+_124090853_124090964_124090574_124090764_124094017_124094144        0.213
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678600_103685170_103685260      0.402
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678651_103685170_103685260      0.402
#>                                                                                               SRR5962199
#> ENSG00000147403.18_RPL10_chrX_+_154400463_154400626_154399802_154399941_154400701_154400811        0.969
#> ENSG00000196924.19_FLNA_chrX_-_154357250_154357274_154354824_154355072_154357433_154357623         0.399
#> ENSG00000101972.20_STAG2_chrX_+_124090853_124090964_124090574_124090764_124094017_124094144        0.213
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678600_103685170_103685260      0.402
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678651_103685170_103685260      0.000
#>                                                                                               SRR5962200
#> ENSG00000147403.18_RPL10_chrX_+_154400463_154400626_154399802_154399941_154400701_154400811        0.969
#> ENSG00000196924.19_FLNA_chrX_-_154357250_154357274_154354824_154355072_154357433_154357623         0.322
#> ENSG00000101972.20_STAG2_chrX_+_124090853_124090964_124090574_124090764_124094017_124094144        0.255
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678600_103685170_103685260      0.251
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678651_103685170_103685260      0.192
#>                                                                                               SRR5962201
#> ENSG00000147403.18_RPL10_chrX_+_154400463_154400626_154399802_154399941_154400701_154400811        0.978
#> ENSG00000196924.19_FLNA_chrX_-_154357250_154357274_154354824_154355072_154357433_154357623         0.463
#> ENSG00000101972.20_STAG2_chrX_+_124090853_124090964_124090574_124090764_124094017_124094144        0.172
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678600_103685170_103685260      0.173
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678651_103685170_103685260      0.322
#>                                                                                               SRR5962202
#> ENSG00000147403.18_RPL10_chrX_+_154400463_154400626_154399802_154399941_154400701_154400811        0.961
#> ENSG00000196924.19_FLNA_chrX_-_154357250_154357274_154354824_154355072_154357433_154357623         0.514
#> ENSG00000101972.20_STAG2_chrX_+_124090853_124090964_124090574_124090764_124094017_124094144        0.237
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678600_103685170_103685260      0.198
#> ENSG00000123562.18_MORF4L2_chrX_-_103684680_103684729_103678498_103678651_103685170_103685260      0.242
```

``` r
## Users can utilize the MRAS function for a streamlined analysis, or execute individual steps separately if they prefer to have more control over specific aspects of the analysis.
a<-MRAS(
  expr,
  psi,
  rbp_interested = "ESRP1",
  m = 22,
  n = 22,
  DS_pvalue = 0.05,
  DS_dPSI = 0.1,
  method =  "spearman",
  num1 = 0.2,
  num2 = 10,
  dpsi_network_threshold = 0.1,
  Regulate_threshold = 0.5,
  threads = 6,
  path_use="./tests/",
  result_type = "Top10"
)
#> Step1:Performing differential splicing analysis...
#> Step2:Constructing RBP-Event regulatory relationship network...
#> Step3:Preparing for enrichment analysis...
#> Step4:Performing enrichment analysis...
#> Finish!
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

# `{r pressure, echo = FALSE} # plot(pressure) #`

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
