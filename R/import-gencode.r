#' Get the features with attributes of interest from a gencode .gtf file 
#'
#' The gtf file format is described at 
#' https://www.gencodegenes.org/data_format.html
#'  
#' @param path Path to the gencode file. 
#' @return A tidied tibble
#' @examples
#' \dontrun{
#' import_gencode(path = 'path/to/file.tsv')
#' }
#' @importFrom magrittr %>%
#' @export

import_gencode <- function(path){
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

  # extract fields of interest
  rx <- paste('gene_id "(.+?)"; transcript_id "(.+?)"; gene_type "(.+?)";',
              '.*; gene_name "(.+?)"; transcript_type "(.+?)";',
              'transcript_status "(.+?)";',
              'transcript_name "(.+?)";( .*; tag "(.+?)";)*', sep = " ")

  parsed_gtf <- gtf %>%
  tidyr::extract(attribute,
          c("gene_id",
            "transcript_id",
            "gene_type",
            "gene_name",
            "transcript_type",
            "transcript_status",
            "transcript_name",
            "extra",
            "tag"),
          rx, remove = FALSE) %>%
  dplyr::select(-attribute, -extra)

  return(parsed_gtf)
}
