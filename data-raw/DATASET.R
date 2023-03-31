## code to prepare `DATASET` dataset goes here
sample<-read.table("./inst/extdata/sample.txt",header = F)
sample<-do.call(c,sample)

expr<-read.table("./inst/extdata/expr.txt",header = T)

psi<-read.table("./inst/extdata/psi.txt",header = T)

expr<-expr[,sample]
psi<-psi[,sample]

usethis::use_data(expr,psi, overwrite = TRUE)
