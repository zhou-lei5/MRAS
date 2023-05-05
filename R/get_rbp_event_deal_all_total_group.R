#' This function provides a way to obtain the matrix of all RBP-Events regulatory relationships.
#'
#' @param path_use The path to the file used to store the output.
#'
#' @return Get rbp_event_deal_all file and save RData.
#'         rbp_event_deal_all_total contains four columns indicating the event ID, the degree of moderation by RBP, the RBP, the degree of difference of splicing events during inference.
#' @export
#'

get_rbp_event_deal_all_total_group<-function(path_use){
  if (file.exists(paste0(path_use,"deal/rbp_event_deal_all_group.txt"))){
    rbp_event_deal_all_total_group<-read.table(paste0(path_use,"deal/rbp_event_deal_all_total_group.txt"))
    colnames(rbp_event_deal_all_total_group)<-c("ID","score","RBP","dPSI")
    save(rbp_event_deal_all_total_group,file = paste0(path_use,"deal/rbp_event_deal_all_total_group.RData"))
    return(rbp_event_deal_all_total_group)
  }else{
    stop("Please run this step after running MRAS_network()")
  }
}
