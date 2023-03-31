#' This function provides a way to obtain the matrix of all RBP-Events regulatory relationships.
#'
#' @param path_use The path to the file used to store the output.
#'
#' @return Get rbp_event_deal_all file and save RData.
#'         rbp_event_deal_all_total contains four columns indicating the event ID, the degree of moderation by RBP, the RBP, the degree of difference of splicing events during inference.
#' @export
#'

get_rbp_event_deal_all_total<-function(path_use){
  if (file.exists(paste0(path_use,"deal/rbp_event_deal_all.txt"))){
    rbp_event_deal_all_total<-read.table(paste0(path_use,"deal/rbp_event_deal_all_total.txt"))
    colnames(rbp_event_deal_all_total)<-c("ID","score","RBP","dPSI")
    save(rbp_event_deal_all_total,file = paste0(path_use,"deal/rbp_event_deal_all_total.RData"))
    return(rbp_event_deal_all_total)
  }else{
    stop("Please run this step after running MRAS_network()")
  }
}
