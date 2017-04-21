#' Get the features with attributes of interest from a gencode .gtf file 
#'
#' The gtf file format is described at 
#' https://www.gencodegenes.org/data_format.html
#'  
#' @param path Path to the gencode file.
#' @param featuretag Tag to include in tag variable (defined at
#'  https://www.gencodegenes.org/gencode_tags.html) Default = "basic"
#' @return A tidied tibble
#' @examples
#' \dontrun{
#' import_gencode(path = 'path/to/file.tsv')
#' }
#' @importFrom magrittr %>%
#' @export

import_gencode <- function(path, featuretag = "basic"){
  # check arguments
  if (missing(path)) stop("path not defined")

  # read the file to make a tibble for tidyverse work.
  gtf <- suppressMessages(readr::read_tsv(path,
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
  )

  # define regular expression to extract fields of interest
  rx <- paste('gene_id "(.+?)"; ', # required gene_id field in gtf
              'transcript_id "(.+?)"; ', # required in gtf
              'gene_type "(.+?)"; ', # required in gtf
              ".*; ", # gene_status not needed
              'gene_name "(.+?)"; ', # required in gtf
              'transcript_type "(.+?)"; ', # required in gtf
              'transcript_status "(.+?)"; ', # required in gtf
              'transcript_name "(.+?)"; ', # required in gtf
              ".*; ", # some optional fields begin
              'tag "(', featuretag, ')"', # match the tag from args
              sep = "") # join regex string without spaces

  # parse the attribute column of gtf with regular expression
  parsed_gtf <- gtf %>%
  tidyr::extract(attribute,
          c("gene_id",
            "transcript_id",
            "gene_type",
            "gene_name",
            "transcript_type",
            "transcript_status",
            "transcript_name",
            "tag"),
          rx, remove = FALSE) %>%
  dplyr::select(-attribute)

  return(parsed_gtf)
}
