#' This is a small function for pasting objects together.
#'
#' @param name_list a vector you want to paste
#' @param type which collapse you want to use.If not, then use ","
#'
#' @return a character
#' @export
#'
#' @examples
#' a<-c("you","me","he","she","it")
#' tj(a,"_")
#' tj(a)
tj<-function(name_list,type = ","){
  name_list<-name_list[complete.cases(name_list)]
  name<-paste0(name_list,collapse = type)
  return(name)
}
