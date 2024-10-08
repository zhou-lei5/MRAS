#' This function provides the means of processing the expression matrix.
#'
#' @param expr RBP expression matrix.
#' @param m Number of the first condition.
#' @param n Number of the second condition.
#' @param RBP_cutoff The cutoff of RBP expression between two conditions.
#'
#' @return Preparing for enrichment analysis.
#' @export
#'

get_RBP_use<-function(expr,m,n,RBP_cutoff=0.05){
  RBP_use<-as.data.frame(expr)
  RBP_use$logFC<-unlist(apply(expr,1,function(x){
    return(abs(log2(mean(as.numeric(x[1:m]))/mean(as.numeric(x[(m+1):(m+n)])))))
  }))
  P<-unlist(apply(expr,1,function(x){
    p<-t.test(as.numeric(x[1:m]),as.numeric(x[(m+1):(m+n)]))
    return(p$p.value)
  }))
  if_P<-(P<RBP_cutoff)
  RBP_use<-RBP_use[if_P,]
  RBP_use<-as.data.frame(RBP_use)
  RBP_use$logFC[is.infinite(RBP_use$logFC)]<-0
  RBP_use$logFC[is.nan(RBP_use$logFC)]<-0
  return(RBP_use)
}
