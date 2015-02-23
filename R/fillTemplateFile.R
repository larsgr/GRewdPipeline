#' Fill out template file
#'
#' @param x A named character vector where the names are the placeholders and values is the replacements
#' @param templateFile Path of template file
#' @value Character vector of filled out template
#' 
#' @details The template file must contain placholder names of the format ${name}

fillTemplateFile <- function(x, templateFile ){
  txt <- readLines(con=templateFile)
  
  for( placeholder in names(x) ){
    txt <- gsub(paste("${",placeholder,"}",sep=""),x[placeholder],txt,fixed=T)
  }
  
  return(txt)
}