#' This This function provides methods for finding correlations and constructing them into conforming matrices.
#'
#' @param x,y numeric data.frame or matrix. x and y must have the same dimension.
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#'
#' @return The returned correlation matrix is spoofed with a total of 4 columns, the first column is the row name of the input x, the second is the row name of the input y, the third column is the magnitude of the correlation, and the fourth column is the significance level p-value. If x and y do not have row names, the returned value is the row to which they belong.
#' @export
#' @import checkmate
#' @importFrom  stats cor pt
#'
#' @examples
#' x <- matrix(sample(1:100,100),10,10)
#' y <- matrix(sample(1:100,100),10,10)
#' rownames(x)=rownames(y)=paste0("row",1:10)
#' colnames(x)=colnames(y)=paste0("col",1:10)
#' cor_spearman(x,y)
cor_spearman<-function(x,y,method = "spearman"){
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(y)
  rr<-cor(t(x),t(y),method = method)
  df<-ncol(x)-2
  pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
  pp[which(pp>1)]<-2-pp[which(pp>1)]
  corr_mat<-cor_mat2l(rr,pp)
  return(corr_mat)
}


cor_mat2l<-function(cormat,pmat){
  if ((is.null(rownames(cormat)))|(is.null(colnames(cormat)))){
    data.frame(row = c(1:nrow(cormat))[row(cormat)],
               col = c(1:ncol(cormat))[col(cormat)],
               cor = c(unlist(cormat)),p = c(unlist(pmat)))
  }else{
    data.frame(row = rownames(cormat)[row(cormat)],
             col = colnames(cormat)[col(cormat)],
             cor = c(unlist(cormat)),p = c(unlist(pmat)))
  }
}
