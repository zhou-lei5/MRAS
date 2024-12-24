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
#' @param rbp_event_deal_all The network built by MRAS.
#' @param Event_DS_sig The differetial splicing events got by MRAS.
#'
#' @return splicing events co-regulated by interacting RBPs.
#' @export
#'

get_group_events<-function(rbp_interested,rbp_event_deal_all,Event_DS_sig){
  tab<-string_net[which(rownames(string_net)==rbp_interested),which(string_net[which(rownames(string_net)==rbp_interested),]!=0)]
  if (!(nrow(tab)>0)){
    stop("This RBP may have no interacting RBP.")
  }else{
    a1<-rbp_event_deal_all[which(rbp_event_deal_all$RBP==rbp_interested),]

    a2<-rbp_event_deal_all[which(rbp_event_deal_all$RBP %in% c(rbp_interested,colnames(tab))),]

    a3<-a2[which(a2$ID %in% a1$ID),]
    a3_DE<-a3[which(a3$ID %in% Event_DS_sig[,1]),]

    if (length(unique(a3_DE$RBP)) == 1){
      stop("This RBP may not have other RBPs that co-regulate splicing events.")
    }else{
      rbp_all<-c(rbp_interested,colnames(tab))
      tab_t<-cbind(RBP=rownames(t(tab)),t(tab))
      num<-as.data.frame(table(a3_DE$RBP))
      colnames(num)<-c("RBP","number")
      num<-merge(num,tab_t,all = T)
      rownames(num)<-num[,1]
      num[rbp_interested,3]<-1
      num[is.na(num)]<-0
      num$percent<-num$number/num[rbp_interested,2]
      # num$inter_stenghth<-ifelse(num$RBP == rbp_interested,1,
      #                            ifelse(num$RBP %in% colnames(tab),tab[1,setdiff(unique(a3_DE$RBP),rbp_interested)],0))

      num<-num[order(num$number,decreasing = T),]
      group_event<-list(summary = num,
                        events  = a3_DE)
      return(group_event)
    }
  }
}
