## code to prepare `DATASET` dataset goes here
sample<-read.table("./inst/extdata/sample.txt",header = F)
sample<-do.call(c,sample)

brca_expr<-read.table("./inst/extdata/expr.txt",header = T)

brca_psi<-read.table("./inst/extdata/psi.txt",header = T)

brca_expr<-brca_expr[,sample]
brca_psi<-brca_psi[,sample]

load("./inst/extdata/hcc_expr.RData")
load("./inst/extdata/hcc_psi.RData")
usethis::use_data(hcc_psi,hcc_expr, overwrite = TRUE)

load("./inst/extdata/sc_brca_expr.RData")
load("./inst/extdata/sc_brca_psi.RData")
usethis::use_data(sc_brca_expr,sc_brca_psi, overwrite = TRUE)


# load("./inst/extdata/TCGA_BRCA_expr.RData")
# load("./inst/extdata/TCGA_BRCA_psi.RData")
# load("./inst/extdata/TCGA_BRCA_net.RData")

usethis::use_data(brca_expr,brca_psi, overwrite = TRUE)

# usethis::use_data(TCGA_BRCA_expr, overwrite = TRUE)
# usethis::use_data(TCGA_BRCA_psi, overwrite = TRUE)
# usethis::use_data(TCGA_BRCA_net, overwrite = TRUE)

load("./inst/extdata/string_net.RData")

usethis::use_data(string_net,overwrite = TRUE,internal = TRUE)

