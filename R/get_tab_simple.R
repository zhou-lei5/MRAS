#' This function provides a way to get the results in the form of "tab_simple".
#'
#' @param path_useful The path to the file used to store the output.
#'
#' @return Return results in "tab_simple" format
#' @export
#'
get_tab_simple<-function(path_useful){
  dat<-read.table(paste0(path_useful,"/result_tab_simple.txt"),header = T)
  tmp<-dat[,c("RBP","score1_nor","nes1_nes","nes2_nes","OR","score3")]
  colnames(tmp)<-c("RBP","D","NES1","NES2","odds","MRAS_Score")
  return(tmp)
}


#' This function provides a way to get the results in the form of "tab_simple".
#'
#' @param expr The path to the file used to store the output.
#' @param m Number of the first condition.
#' @param n Number of the second condition.
#' @param result_tab_simple The output of MRAS.
#' @param RBP_cutoff The cutoff of RBP expression between two conditions.
#'
#' @return Return results with direction in "tab_simple" format
#' @export
#'
get_direaction<-function(expr,m,n,result_tab_simple,RBP_cutoff=0.05){
  RBP_use<-as.data.frame(expr)
  RBP_use$logFC<-unlist(apply(expr,1,function(x){
    return(log2(mean(as.numeric(x[1:m]))/mean(as.numeric(x[(m+1):(m+n)]))))
  }))
  P<-unlist(apply(expr,1,function(x){
    p<-t.test(as.numeric(x[1:m]),as.numeric(x[(m+1):(m+n)]))
    return(p$p.value)
  }))
  if_P<-(P<RBP_cutoff)
  RBP_use<-RBP_use[if_P,]
  RBP_use<-as.data.frame(RBP_use)
  RBP_use$logFC[is.infinite(RBP_use$logFC)]<-0
  RBP_use$logFC[is.nan(RBP_use$logFC)]<-0
  result_tab_simple<-cbind(result_tab_simple,logFC_new=0)
  result_tab_simple<-cbind(result_tab_simple,score3_new=0)

  for (i in 1:nrow(result_tab_simple)) {
    result_tab_simple$logFC_new[i]<-RBP_use[result_tab_simple$RBP[i],"logFC"]
    result_tab_simple$score3_new[i]<-ifelse(result_tab_simple$logFC_new[i]>0,result_tab_simple$score3[i],-result_tab_simple$score3[i])
  }
  # result_tab_simple<-result_tab_simple[order(result_tab_simple$score3_new,decreasing = T),]
  result_tab_simple$direction[i]<-ifelse(result_tab_simple$logFC_new[i]>0,"Upregulation","Downregulation")
  return(result_tab_simple)
}
