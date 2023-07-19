#' Looking for interacting RBPs.
#'
#' @param rbp_interested The name of the RBP interested.
#'
#' @return interacting RBPs table.
#' @export
#'

get_group<-function(rbp_interested){
  tab<-string_net[which(rownames(string_net)==rbp_interested),which(string_net[which(rownames(string_net)==rbp_interested),]!=0)]
  if (!(nrow(tab)>0)){
    stop("This RBP may have no interacting RBP.")
  }else{
    return(tab)
  }
}
#' Searching for splicing events co-regulated by interacting RBPs.
#'
#' @param rbp_interested The name of the RBP interested.
#' @param rbp_event_deal_all the network built by MRAS.
#'
#' @return splicing events co-regulated by interacting RBPs.
#' @export
#'

get_group_events<-function(rbp_interested,rbp_event_deal_all){
  tab<-string_net[which(rownames(string_net)==rbp_interested),which(string_net[which(rownames(string_net)==rbp_interested),]!=0)]
  if (!(nrow(tab)>0)){
    stop("This RBP may have no interacting RBP.")
  }else{
    a1<-rbp_event_deal_all[which(rbp_event_deal_all$RBP==rbp_interested),]

    a2<-rbp_event_deal_all[which(rbp_event_deal_all$RBP %in% c(rbp_interested,colnames(tab))),]

    a3<-a2[which(a3$ID %in% a1$ID),]

    if (nrow(a3)>nrow(a1)){
      return(a3)
    }else{
      stop("This RBP may not have other RBPs that co-regulate splicing events.")
    }
  }
}
