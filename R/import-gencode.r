#' Load a gencode file. 
#'  
#' @param path Path to the gencode file. 
#' @return A data frame containing the rows from the file that match the 
#' specified feature and contain the given attribute. 
#' @examples
#' \dontrun{
#' import_gencode(path = 'path/to/file.tsv')
#' }
#' @export

import_gencode <- function(path){
  if (missing(path)) stop("path not defined")
}
