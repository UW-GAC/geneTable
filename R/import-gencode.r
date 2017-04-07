#' Get the features with attributes of interest from a gencode .gtf file 
#'
#' The gtf file format is described at 
#' https://www.gencodegenes.org/data_format.html
#'  
#' @param path Path to the gencode file. 
#' @param feature The feature type of interest. Must be one of {gene, 
#'   transcript, exon, CDS, UTR, start_codon, stop_codon, Selenocysteine}. 
#'   Default = "transcript"
#' @param attribute The attribute of interest (from the 9th column of the gtf 
#'   file). Default = "tag basic"
#' @return A data frame containing the rows from the file that match the 
#'   specified feature and contain the given attribute. 
#' @examples
#' \dontrun{
#' import_gencode(path = 'path/to/file.tsv')
#' }
#' @export

import_gencode <- function(path,
                           feature = "transcript",
                           attribute = "tag basic"){
  # check arguments
  if (missing(path)) stop("path not defined")

  if (!(feature %in%
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

  if (!is.character(attribute)) stop("attribute must be a string")

  # read the file to make a tibble for tidyverse work.
  gtf <- readr::read_tsv(path,
                         comment = "#",
                         col_names = c("seqname",
                                       "source",
                                       "feature",
                                       "start",
                                       "end",
                                       "score",
                                       "strand",
                                       "frame",
                                       "attribute")
                                       )
}
