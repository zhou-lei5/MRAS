
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

This is a basic example which shows you how to use MRAS:

``` r
library(MRAS)
#> Loading required package: data.table
#> Loading required package: fgsea
## basic example code
```

Test data is included in MRAS and you can import it using data(), which
contains the RBP expression matrix as well as the event PSI matrix.

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

The easiest way to use it: directly use the function `MRAS()`. The
specific parameters are detailed in `??MRAS` or `help(MRAS)`.

``` r
## Users can utilize the MRAS function for a streamlined analysis, or execute individual steps separately if they prefer to have more control over specific aspects of the analysis.
result<-MRAS(
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

After running `MRAS()`, there are three ways to display the results. In
addition to setting the form directly in the parameters, users can also
obtain other result display forms through the functions `get_Top10()`,
`get_tab_all()`, and `get_tab_simple()`. This allows users to access
additional result display formats without having to rerun `MRAS()`.

``` r
result_Top10<-get_Top10(path_use = "./tests/")
result_tab_simple<-get_tab_simple(path_use = "./tests/")
result_tab_all<-get_tab_all(path_use = "./tests/")
head(result_tab_simple[1:5,])
#>     RBP     logFC   score1  m1 score1_nor nes1_size   nes1_es nes1_nes
#> 1 ESRP1 1.6128311 32.14689 106  0.8416978       193 0.8403875 1.726525
#> 2 RBM11 1.7321882 29.95712  80  0.7843613       148 0.8169598 1.676204
#> 3 RBMS3 1.6080483 32.37771  99  0.8477417       187 0.8369700 1.710362
#> 4 RBMS1 0.6918446 12.59383  89  0.3297241       133 0.8669439 1.761774
#> 5 FKBP4 1.1147372 19.50317  88  0.5106369       143 0.8463503 1.738259
#>        nes1_p nes2_size   nes2_es nes2_nes      nes2_p overlap total_size
#> 1 0.000999001       146 0.9149168 1.636088 0.000999001     106       1857
#> 2 0.000999001       146 0.8429471 2.799824 0.000999001      80       1857
#> 3 0.000999001       146 0.9027223 1.685116 0.000999001      99       1857
#> 4 0.000999001       146 0.9045543 2.384146 0.000999001      89       1857
#> 5 0.000999001       146 0.9023890 1.984417 0.000999001      88       1857
#>         OR         pval    score3
#> 1 49.20177 3.606096e-84 116.98127
#> 2 29.12594 4.097512e-58 107.21450
#> 3 38.62165 1.166684e-74  94.36527
#> 4 58.80437 2.355515e-78  81.44105
#> 5 45.34175 2.331438e-72  79.86525
```

Of course, you can also distribute the MRAS analysis, which will help
you better understand the principles of MRAS. Finally, if you have any
more questions, you can submit them in Github and we will do our best to
answer them.(<https://github.com/zhou-lei5/MRAS>)
