#' This function provides a way to obtain the matrix of all RBP-Events regulatory relationships with P > 0.5.
#'
#' @param path_useful The path to the file used to store the output.
#'
#' @return Get rbp_event_deal_all file and save RData.
#'         rbp_event_deal_all contains three columns indicating the event ID, the degree to which the RBP regulates it, the RBP
#' @export
#'

get_rbp_event_deal_all<-function(path_useful){
  if (file.exists(paste0(path_useful,"/rbp_event_deal_all.txt"))){
      rbp_event_deal_all<-read.delim(paste0(path_useful,"/rbp_event_deal_all.txt"),header = F)
      colnames(rbp_event_deal_all)<-c("ID","score","RBP")
      save(rbp_event_deal_all,file = paste0(path_useful,"/rbp_event_deal_all.RData"))
      return(rbp_event_deal_all)
  }else{
    stop("Please run this step after running MRAS_network()")
  }
}
