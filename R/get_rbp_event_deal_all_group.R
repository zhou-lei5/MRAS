#' This function provides a way to obtain the matrix of all RBP-Events regulatory relationships with P > 0.5.
#'
#' @param path_use The path to the file used to store the output.
#'
#' @return Get rbp_event_deal_all file and save RData.
#'         rbp_event_deal_all contains three columns indicating the event ID, the degree to which the RBP regulates it, the RBP
#' @export
#'

get_rbp_event_deal_all_group<-function(path_use){
  if (file.exists(paste0(path_use,"deal/rbp_event_deal_all_group.txt"))){
    rbp_event_deal_all_group<-read.table(paste0(path_use,"deal/rbp_event_deal_all_group.txt"))
    colnames(rbp_event_deal_all_group)<-c("ID","score","RBP")
    save(rbp_event_deal_all_group,file = paste0(path_use,"deal/rbp_event_deal_all_group.RData"))
    return(rbp_event_deal_all_group)
  }else{
    stop("Please run this step after running MRAS_network()")
  }
}
