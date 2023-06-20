
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
## basic example code
```

Test data is included in MRAS and you can import it using data(), which
contains the RBP expression matrix as well as the event PSI matrix.

``` r
## "hcc_expr" is RBP expression matrix.
data("hcc_expr")
hcc_expr[1:5,1:5]
#>        SRR3182261 SRR3129836 SRR3129837 SRR3129838 SRR3129839
#> A1CF   21.7575656 116.042384  47.071361  37.644303 54.8133944
#> ANKHD1  0.5555836   1.244645   0.727923   1.429632  0.8184703
#> CELF1  12.5222101  33.131067  29.052948  25.443336 26.5795453
#> CELF2   1.1010763   1.974816   3.688240   1.537169  1.4341667
#> CNOT4   2.6800668   4.430709   5.144975   3.858476  3.3844203
```

``` r
## "hcc_psi" is events psi matrix."hcc_expr" and "hcc_psi" should have same column names.
data("hcc_psi")
hcc_psi[1:5,1:3]
#>                                                                                SRR3182261 SRR3129836 SRR3129837
#> FAM3A_ES_chrX_-_153740638_153740735_153740181_153740204_153741147_153741260         0.118      0.150      0.389
#> FAM3A_ES_chrX_-_153740638_153740755_153740201_153740204_153741147_153741260         0.084      0.072      0.202
#> TAFAZZIN_ES_chrX_+_153642438_153642527_153641819_153641904_153647882_153647962      0.310      0.330      0.447
#> RPL10_ES_chrX_+_153628805_153628967_153628144_153628282_153629043_153629152         0.994      0.997      0.992
#> SSR4_ES_chrX_+_153061383_153061581_153060131_153060209_153061889_153062007          0.049      0.013      0.010
```

The easiest way to use it: directly use the function `MRAS()`. The
specific parameters are detailed in `??MRAS` or `help(MRAS)`.

``` r
## Users can utilize the MRAS function for a streamlined analysis, or execute individual steps separately if they prefer to have more control over specific aspects of the analysis.
result1<-MRAS(
  expr = hcc_expr,
  psi = hcc_psi,
  rbp_interested = "SF3B4",
  m = 50,
  n = 50,
  DS_pvalue = 0.05,
  DS_dPSI = 0.1,
  method =  "spearman",
  num1 = 0.15, num2 = 0.15,
  BS = NULL, net_type = "1",
  dpsi_network_threshold = 0.1,
  Regulate_threshold = 0.5,
  threads = 6,
  path_use="./tests/",
  group=T,
  result_type = "Top10"
)
#> Step1:Performing differential splicing analysis...
#> Step2:Preparing data...
#> Step3:Constructing RBP-Event regulatory relationship network...
#> Joining with `by = join_by(rbp, BS)`
#> 
#> Step4:Performing enrichment analysis...
#> Finish!
result1
#>      rbp_interested rank RBP1    rank1              RBP2  rank2      RBP3      rank3              RBP4   rank4              RBP5  
#> [1,] "SF3B4"        "1"  "SF3B4" "24.6550702640961" "PKM" "22.88593" "IGF2BP2" "22.1797574937597" "RRP9" "20.9632303515802" "BOP1"
#>      rank5              RBP6   rank6             RBP7    rank7              RBP8    rank8              RBP9    rank9             
#> [1,] "20.7894630961655" "XPO5" "18.900063049553" "NELFE" "18.2031305015691" "RBM42" "18.0584531954273" "SNRPA" "18.0043547945609"
#>      RBP10  rank10           
#> [1,] "RALY" "17.358981454952"
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
#>       RBP    logFC   score1  m1 score1_nor nes1_size   nes1_es nes1_nes      nes1_p nes2_size   nes2_es nes2_nes      nes2_p overlap
#> 1   SF3B4 2.407401 163.3349 447  0.7065393       486 0.9013645 1.810043 0.000999001       632 0.9271163 1.363220 0.000999001     447
#> 2     PKM 3.611442 231.1632 417  1.0000000       488 0.8812575 1.762094 0.000999001       632 0.9080232 1.320125 0.000999001     417
#> 3 IGF2BP2 3.386165 161.8666 284  0.7001862       297 0.8887830 1.772951 0.000999001       632 0.9253928 1.263373 0.000999001     284
#> 4    RRP9 2.166866 150.7078 451  0.6519076       500 0.8958257 1.798930 0.000999001       632 0.9269796 1.326146 0.000999001     451
#> 5    BOP1 2.686065 181.3979 434  0.7846890       498 0.8916058 1.784817 0.000999001       632 0.9222807 1.346828 0.000999001     434
#>   total_size        OR          pval   score3
#> 1       4259 200.00000  0.000000e+00 24.65507
#> 2       4259  96.79408  0.000000e+00 22.88593
#> 3       4259 200.00000 3.645079e-243 22.17976
#> 4       4259 181.69070  0.000000e+00 20.96323
#> 5       4259 121.47310  0.000000e+00 20.78946
```

If using sc-RNA seq data, some methods should change.

``` r
data("sc_brca_expr")
data("sc_brca_psi")
result2<-MRAS(
  expr = sc_brca_expr,
  psi = sc_brca_psi,
  rbp_interested = "ESRP1",
  m = 198,
  n = 317,
  DS_pvalue = 0.05,
  DS_dPSI = 0.1,
  method =  "spearman",
  num1 = 0.1, num2 = 0.1,
  BS = NULL, net_type = "1",
  dpsi_network_threshold = 0.1,
  Regulate_threshold = 0.5,
  threads = 6,
  path_use="./tests/",
  group=T,
  sc=T,
  result_type = "Top10"
)
#> Step1:Performing differential splicing analysis...
#> Step2:Preparing data...
#> Step3:Constructing RBP-Event regulatory relationship network...
#> Joining with `by = join_by(rbp, BS)`
#> 
#> Step4:Performing enrichment analysis...
#> Finish!
result2
#>      rbp_interested rank RBP1    rank1      RBP2    rank2              RBP3       rank3              RBP4    rank4             
#> [1,] "ESRP1"        "1"  "ESRP1" "16.80279" "RBM47" "3.75512809116142" "APOBEC3C" "1.77086642393804" "MBNL1" "1.39991934961294"
#>      RBP5      rank5              RBP6    rank6               RBP7    rank7             RBP8     rank8               RBP9  
#> [1,] "HNRNPH2" "1.27600655740444" "RBM28" "0.774343871900044" "DDX24" "0.7702568970529" "PABPC1" "0.536357017226172" "SND1"
#>      rank9               RBP10   rank10             
#> [1,] "0.439591932068094" "CELF2" "0.439122487144725"
```

Of course, you can also distribute the MRAS analysis, which will help
you better understand the principles of MRAS. Finally, if you have any
more questions, you can submit them in Github and we will do our best to
answer them (<https://github.com/zhou-lei5/MRAS>).
