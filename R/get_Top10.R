#' This function provides a way to get the results in the form of "Top10".
#'
#' @param path_useful The path to the file used to store the output.
#'
#' @return Return results in "top10" format
#' @export
#'
get_Top10<-function(path_useful){
  dat<-read.table(paste0(path_useful,"/result_top10.txt"),header = T)
  return(dat)
}
