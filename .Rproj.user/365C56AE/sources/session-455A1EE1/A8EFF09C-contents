#' This function provides a way to get significantly differentially spliced events from the differentially spliced event matrix.
#'
#' @param Events_DS The difference matrix of the clipping events, obtained by the MRAS function DS_matrix().
#' @param DS_pvalue Significance level of differential splicing analysis.
#' @param DS_dPSI Threshold for differential splicing.
#'
#' @return Events with significant differences clipped by threshold filtering in Events_DS.
#' @export
#'

get_Event_DS_sig<-function(Events_DS,DS_pvalue=0.05,DS_dPSI=0.1){
  Event_DS_sig<-Events_DS[which((as.numeric(Events_DS[,2])<0.05)&(as.numeric(Events_DS[,5])>0.1)),]
  colnames(Event_DS_sig)<-c("Events","pvalue","FDR","dPSI","abs_dPSI")
  return(Event_DS_sig)
}
