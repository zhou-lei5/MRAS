#' This function provides a way to get the results in the form of "tab_simple".
#'
#' @param path_use The path to the file used to store the output.
#'
#' @return Return results in "tab_simple" format
#' @export
#'
get_tab_simple<-function(path_use){
  dat<-read.table(paste0(path_use,"/result_tab_simple.txt"))
  return(dat)
}
