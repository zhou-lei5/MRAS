#' This function provides a way to get the results in the form of "tab_simple".
#'
#' @param path_useful The path to the file used to store the output.
#'
#' @return Return results in "tab_simple" format
#' @export
#'
get_tab_simple<-function(path_useful){
  dat<-read.table(paste0(path_useful,"/result_tab_simple.txt"),header = T)
  return(dat)
}
