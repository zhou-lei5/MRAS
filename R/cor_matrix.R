#' This function provides methods for finding correlations and constructing them into conforming matrices.
#'
#' @param RBP_matrix RBP expression matrix
#' @param Events_matrix Events splicing matrix
#' @param method "pearson" or "spearman"
#'
#' @return The result include correlation coefficient and significance level matrix.
#' @export
#'
#' @importFrom stats cor.test complete.cases
#' @examples
#' a<-matrix(rnorm(25),nrow=5,ncol=5)
#' rownames(a)<-paste("gene",1:5,sep="")
#' b<-matrix(rnorm(25),nrow=5,ncol=5)
#' rownames(b)<-paste("event",1:5,sep="")
#' c<-cor_matrix(a,b)
#' head(c[1:5,])
cor_matrix<-function(RBP_matrix,Events_matrix,method="pearson"){
  RBP_Events<-matrix(NA,nrow = nrow(RBP_matrix)*nrow(Events_matrix),ncol = 5)
  colnames(RBP_Events)<-c("rbp_ID","Events","type","corr","p")
  for (p in 1:nrow(RBP_matrix)) {
    for (q in 1:nrow(Events_matrix)) {
      rank_number<-(p-1)*nrow(Events_matrix)+q
      cor_result<-cor_method(as.numeric(RBP_matrix[p,]),as.numeric(Events_matrix[q,]),method)
      RBP_Events[rank_number,1]<-as.character(rownames(RBP_matrix)[p])
      RBP_Events[rank_number,2]<-as.character(rownames(Events_matrix)[q])
      RBP_Events[rank_number,4]<-cor_result[["corr"]]
      RBP_Events[rank_number,5]<-cor_result[["pvalue"]]
    }
  }
  RBP_Events[,3]<-ifelse(RBP_Events[,4]>0,"up","down")
  RBP_Events<-RBP_Events[complete.cases(RBP_Events),]
  return(RBP_Events)
}


#' This function is used to provide a method for finding correlation.
#'
#' @param x a numeric vector
#' @param y a numeric vector
#' @param method "pearson" or "spearman"
#'
#' @return The result include correlation coefficient and significance level.
#' @export
#'
#' @examples
#'
#' a<-rnorm(10)
#' b<-rnorm(10)
#' cor_method(a,b,"pearson")
cor_method<-function(x,y,method="pearson"){
  a<-stats::cor.test(x,y,method=method)
  return(list(corr=a$estimate,pvalue=a$p.value))
}

