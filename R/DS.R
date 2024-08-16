#' This function provides a method to find the differential splicing events.
#'
#' @param x a numeric vector about your events splicing level.And x should be ordered by conditions
#' @param m number of the first condition.
#' @param n number of the second condition.
#'
#' @return The result include p-value, FDR, dPSI and abs_dPSI.
#' @export
#' @importFrom  stats t.test var
#'
#'
#' @examples
#' a<-abs(rnorm(10))
#' DS(a,5,5)
#'
#'
DS<-function(x,m,n){
  x1<-as.numeric(x[1:m])
  x2<-as.numeric(x[(m+1):(m+n)])
  va1<-var(x1)
  va2<-var(x2)
  if (va1==0) x1[1]<-x1[1]+0.001
  if (va2==0) x2[1]<-x2[1]+0.001

  a<-t.test(x1,x2)
  dPSI<-mean(as.numeric(x2))-mean(as.numeric(x1))
  FDR<-stats::p.adjust(a$p.value,method = "BH")
  return(list(pvalue=a$p.value,FDR=FDR,dPSI=dPSI,abs_dPSI=abs(dPSI)))
}


