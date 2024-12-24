#' This function provides methods to find the differential splicing events and construct them into a conforming matrix.
#'
#' @param psi Events splicing matrix（column should be ordered by conditions）
#' @param m The number of condition A samples
#' @param n The number of condition B samples
#'
#' @return differential splicing matrix
#' @export
#'
DS_matrix<-function(psi,m,n){
  Events_DS<-matrix(NA,nrow = nrow(psi),ncol = 5)
  Events_DS<-t(apply(psi,1,function(x){
    a<-DS(x,m,n)
    bb<-c(rownames(x),a[["pvalue"]],a[["FDR"]],a[["dPSI"]],a[["abs_dPSI"]])
    return(bb)
  }))
  Events_DS<-Events_DS[complete.cases(Events_DS),]
  Events_DS[,2]<-stats::p.adjust(as.numeric(Events_DS[,1]),method = "BH")
  Events_DS<-cbind(rownames(Events_DS),Events_DS)
  colnames(Events_DS)<-c("Events_ID","pvalue","FDR","dPSI","abs_dPSI")
  return(Events_DS)
}
