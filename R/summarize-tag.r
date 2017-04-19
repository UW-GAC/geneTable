#' Summarize the features in a gencode-derived tibble.
#'
#' summarize_tag returns a tibble containing the number of features included
#' in a tibble produced by genetable::import_gencode() (or similar) that are
#' annotated with the tag provided in the argument. Legal tags are listed at
#' https://www.gencodegenes.org/gencode_tags.html
#'
#'  
#' @param gtf A tibble containing gene definitions in the formate made by
#'   genetable::import_gencode()
#' @param tag The tag to summarize features by - a string.
#' @return A table containing the cross tabulation of features by tag.
#' @examples
#' \dontrun{
#' summarize_tag(gtf = gtf, tag = "basic")
#' }
#' @export

summarize_tag <- function(gtf, tag = "basic"){
  # check arguments
  if (missing(gtf)) stop("gtf not defined")
  return(table(gtf$feature, gtf$tag == tag))
}
