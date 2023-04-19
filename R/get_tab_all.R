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
  return(rbp_uni_mat_2)
}
