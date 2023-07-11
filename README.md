
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRAS

<!-- badges: start -->
<!-- badges: end -->

MRAS is designed to identify crucial RNA-binding proteins (RBPs)
responsible for splicing variations in diverse scenarios, including
cancer vs. normal, primary vs. recurrence, and more, not just in bulk
data but also in single-cell data.

## The Overview of MRAS

<figure>
<img src="png/Fig1_0707.png" data-margin="10px"
alt="The Overview of MRAS" />
<figcaption aria-hidden="true">The Overview of MRAS</figcaption>
</figure>

## Installation and Library

You can install the development version of MRAS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zhou-lei5/MRAS")
library(MRAS)
```

## Usage and Examples

There are some basic descriptions which shows you how to use MRAS:

Input type 1: Direct input of the set of alternative splicing events.
Firstly, you need to download the pre-constructed regulation network
generated by MRAS.In general, you will get two matrix:
“rbp_event_deal_all_total” and “rbp_event_deal_all”.

``` r
#The basic code of input type 1.
MRAS(input_type = "1",
     expr,psi,
     rbp_interested,
     m,n,
     DS_pvalue,DS_dPSI,
     rbp_event_deal_all_total,rbp_event_deal_all,
     result_type,threads,path_use)
```

Input type 2: Reconstruct the network using inferred relationships.
Firstly, you need to download the inferred relationships generated by
MRAS. In general, yo will get one matrix: “rbp_net_mat_group”.

``` r
#The basic code of input type 2.
MRAS(input_type = "2",
     expr,psi,
     rbp_interested,
     m,n,
     DS_pvalue,DS_dPSI,
     rbp_net_mat_group,group,
     result_type,threads,path_use)
```

Input type 3: Construct a network using the user’s own data. Here, we
have prepared some test data for users to better understand the usage of
MRAS. Test data is included in MRAS and you can import it using data(),
which contains the RBP expression matrix as well as the event PSI
matrix.

Use MRAS in bulk rna-seq data:

``` r
## "hcc_expr" is RBP expression matrix.
library(MRAS)
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

<!-- # ```{r MRAS_BULK} -->

``` r
## Users can utilize the MRAS function for a streamlined analysis, or execute individual steps separately if they prefer to have more control over specific aspects of the analysis.
result_bulk<-MRAS(input_type = "3",
  expr = hcc_expr,
  psi = hcc_psi,
  rbp_interested = "SF3B4",
  m = 50, n = 50,
  DS_pvalue = 0.05, DS_dPSI = 0.1,
  method =  "spearman",
  num1 = 0.15, num2 = 0.15,
  BS = NULL, smooth = F,
  dpsi_network_threshold = 0.1, Regulate_threshold = 0.5,
  group = T,
  result_type = "Top10", threads = 6, path_use = "./tests/"
)
result_bulk
```

After running `MRAS()`, there are three ways to display the results. In
addition to setting the form directly in the parameters, users can also
obtain other result display forms through the functions `get_Top10()`,
`get_tab_all()`, and `get_tab_simple()`. This allows users to access
additional result display formats without having to rerun `MRAS()`.

<!-- ```{r MRAS_result} -->

``` r
result_Top10<-get_Top10(path_use = "./tests/")
result_tab_simple<-get_tab_simple(path_use = "./tests/")
result_tab_all<-get_tab_all(path_use = "./tests/")
head(result_tab_simple[1:5,])
```

Use MRAS in bulk rna-seq data:

<!-- # ```{r MRAS_sc} -->

``` r
data("sc_brca_expr")
data("sc_brca_psi")
result_sc<-MRAS(input_type = "3",
  expr = sc_brca_expr,
  psi = sc_brca_psi,
  rbp_interested = "ESRP1",
  m = 198, n = 317,
  DS_pvalue = 0.05, DS_dPSI = 0.1,
  method =  "spearman",
  num1 = 0.1, num2 = 0.1,
  BS = NULL, smooth = F,
  dpsi_network_threshold = 0.1, Regulate_threshold = 0.5,
  group = T, sc = T,
  result_type = "Top10", threads = 6, path_use = "./tests/"
)
result_sc
```

## Tools: AS Events ID converter

MRAS provides an ID converter specifically designed for splice events.
This converter facilitates the matching of splice event coordinates
obtained from different software, allowing seamless integration with the
pre-constructed regulatory network generated by MRAS. This functionality
simplifies the process of mapping splice events to the existing
regulatory network, increasing the usability and versatility of MRAS.
<span style="display: block; margin-top: 20px; margin-bottom: 10px;">
<img src="png/ID_format.png" data-margin="10px"
alt="AS Events ID format" /> </span>

MRAS provides the following features for ID conversion of splice
events: 1. “id_find: This function allows the user to input the output
path of commonly used splicing event identification software such as
rMATS, SUPPA, and JUM. MRAS will directly return the PSI matrix or
canonical splice event ID associated with the input data. This allows
for seamless integration into the MRAS pre-built regulatory network.
2.”id_normalization”: This function guides the user step-by-step through
the input of the corresponding column coordinates, allowing for the
standardized output of splicing event IDs. The process ensures
consistency and compatibility in the representation of splicing events.
3. “id_change: This function converts splicing event IDs recognized by
two different splicing event identification software. By default, a
mismatch coordinate difference of 2 is used to account for potential
differences in coordinate systems between the software. This allows
users to bridge the gap between different software outputs and harmonize
the representation of splicing events. You can get help by `??MRAS::FUN`
or `help("FUN")`.

## Help

Of course, you can also distribute the MRAS analysis, which will help
you better understand the principles of MRAS. Finally, if you have any
more questions, you can submit them in Github and we will do our best to
answer them (<https://github.com/zhou-lei5/MRAS>).
