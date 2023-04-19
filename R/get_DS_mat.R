#' This function provides a way to differential splicing matrix.
#'
#' @param path_useful The path to the file used to store the output.
#'
#' @return Get DS_mat file.
#'
#' @export

get_DS_mat<-function(path_useful){
  DS_mat<-read.table(paste0(path_useful,"DS_mat.txt"),header = T)
  return(DS_mat)
}
