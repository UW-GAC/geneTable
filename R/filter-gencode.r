#' Filter the gencode tibble by feature and tag.
#'
#' This function filters a gencode tibble (such as returned from
#' import_gencode() to select the given feature type, taged with the given tag.
#'
#' @param gtf The gencode tibble
#' @param featurearg The feature type of interest. Must be one of {gene,
#'   transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine}.
#'   Default = "transcript"
#' @param tagarg The attribute of interest. Legal tags are listed at
#' https://www.gencodegenes.org/gencode_tags.html Default = "basic"
#' @return A tibble containing the rows matching the feature and tag of interest
#' @examples
#' \dontrun{
#' filter_gencode(gtf = gtf, feature = "transcript", tag = "basic")
#' }
#' @importFrom dplyr filter
#' @importFrom tibble is.tibble
#' @export

filter_gencode <- function(gtf,
                           featurearg = "transcript",
                           tagarg){
  # check arguments
  if (missing(gtf)) stop("gtf not defined")

  if (!(featurearg %in%
        c("gene",
          "transcript",
          "exon",
          "CDS",
          "UTR",
          "start_codon",
          "stop_codon",
          "Selenocysteine"))) {
    msg <- paste("feature must be one of gene, transcript, exon, CDS, UTR, ",
                "start_codon, stop_codon, or Selenocysteine", sep = "")
    stop(msg)
  }
  if (missing(tagarg)) {
                       return(filter(gtf, feature == featurearg)) #nolint
  } else {
         return(filter(gtf, feature == featurearg, tag == tagarg)) #nolint
  }
}
