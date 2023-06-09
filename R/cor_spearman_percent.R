#' This function provides the method to find the correlation based on the samples with high and low expressions in each row, taking a certain proportion and constructing it into a matrix that meets the requirements.
#'
#' @param x,y numeric data.frame or matrix. x and y must have the same dimension.
#' @param num1 Denotes the proportion/number of samples to take the expression high or low and to construct the regulatory relationship, respectively.
#'             Input requirements: if the input is less than 1, the input is considered as the sample proportion, and greater than 1, the input is considered as the number of samples.
#'             When choosing the input sample proportion, it is required to be no greater than 0.5.
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#' @param threads If you want to use multiple threads, you can set threads, where the default parameter is 2.
#'
#' @return The returned correlation matrix is spoofed with a total of 4 columns,
#'         the first column is the row name of the input x,
#'         the second is the row name of the input y,
#'         the third column is the magnitude of the correlation,
#'         the fourth column is the significance level p-value.
#' @export
#' @import doParallel
#' @import foreach
#' @import parallel
#' @importFrom  stats cor
#'

#'
#' @examples
#' x <- matrix(sample(1:10000,10000),100,100)
#' y <- matrix(sample(1:10000,10000),100,100)
#' rownames(x)=rownames(y)=paste0("row",1:100)
#' colnames(x)=colnames(y)=paste0("col",1:100)
#' cor_spearman(x,y)
cor_spearman_percent<-function(x,y,method = "spearman",num1 = num1,threads = 2){
  if (num1<1) num1<-ceiling(ncol(x)*num1)
  use_mat<-t(apply(x, 1, function(q){
    use<-c(order(q,decreasing=TRUE)[1:num1],
           order(q,decreasing=FALSE)[1:num1])
    return(use)
  }))
  ii=NULL
  cll<-parallel::makeCluster(threads)
  doParallel::registerDoParallel(cll)
  # clusterEvalQ(cll, .libPaths("D:/R"))
  re<-foreach::foreach (ii = 1:nrow(x), .combine = "rbind") %dopar% {
    t1<-matrix(as.numeric(unlist(x[ii,use_mat[ii,]])),nrow = 1)
    t2<-y[,use_mat[ii,]]
    rr<-cor(t(t1),t(t2),method = method)
    df<-ncol(t1)-2
    rownames(rr)<-rownames(x)[ii]
    pp<-2*pt(rr*sqrt(df)/sqrt(1-rr*rr),df,lower.tail = F)
    pp[which(pp>1)]<-2-pp[which(pp>1)]
    corr_mat<-data.frame(row = rownames(rr)[row(rr)],
                         col = colnames(rr)[col(rr)],
                         cor = c(unlist(rr)),p = c(unlist(pp))
    )
    return(corr_mat)
  }
  parallel::stopCluster(cll)
  return(re)
}
