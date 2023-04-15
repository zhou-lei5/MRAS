#' RBP expression matrix pre-preocessing.
#'
#' @param expr RBP expression matrix
#' @param num,cutoff Threshold.
#'
#' @return RBP expression matrix after filtering out the low-expression RBPs.
#' @export
#'

data_preprocess_expr<-function(expr,num,cutoff){
  f1<-filter_expr(expr = expr,n = num,cutoff = cutoff)
  expr<-expr[f1,]
  return(expr)
}



filter_expr<-function(expr,n,cutoff=1){
  results<-matrix(TRUE,nrow = nrow(expr),ncol = 1)
  for (i in 1:nrow(expr)) {
    length0<-length(which(as.numeric(expr[i,])<=cutoff))
    results[i,1]<-ifelse((length0)>n,FALSE,TRUE)
  }
  return(results[,1])
}
