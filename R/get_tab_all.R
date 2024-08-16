#' This function provides a way to get the results in the form of "tab_all".
#'
#' @param path_useful The path to the file used to store the output.
#'
#' @return Return results in "tab_all" format
#' @export
#'
get_tab_all<-function(path_useful){
  rbp_uni_mat_2=NULL
  load(paste0(path_useful,"/result_tab_all.RData"))
  tmp<-rbp_uni_mat_2[,c("RBP","Events","score1_nor","nes1_nes","nes2_nes","OR","score3")]
  colnames(tmp)<-c("RBP","Targets","D","NES1","NES2","odds","MRAS_Score")
  return(tmp)
}
