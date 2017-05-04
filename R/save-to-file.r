#' Save the tibble defining feature ranges to a tab-delimited text file.
#'
#' Save the tibble of genomic feature ranges to a tab-delimited file with
#' appropriate header columns.
#'
#' The output file includes the following fields as comments in the header:
#'   Date: $DATE
#'   genetable Version: $VERSION
#'   Notes: $NOTES
#'   Data Dictionary:
#'     variable type description
#'     chr character chromosome in genome build GRCh37/hg19
#'     start integer 5' start position of the gene in genome build GRCh37/hg19
#'     end integer 3' end position of the gene in genome build GRCh37/hg19
#'     strand character DNA strand. values {+,-}
#'     gene_id character gene identifier with ENSG prefix
#'     gene_name character Gene name
#'     type character biotype of the feature
#'
#' @param feature_bounds a tibble such as produced by define_boundaries()
#' @param file_name the target file name for writing (default 
#'   "feature_bounds_DATE.tsv")
#' @param notes a character string of notes to include in the file header.
#' @examples
#' \dontrun{
#' notes <- paste("Features selected by filtering gtf for feature =",
#'                '"transcript" and tag = "basic"', sep = " ")
#' save_to_file(gene_bounds, file_name = "gene_bounds_20170504", 
#'              notes = notes)
#' }
#' @importFrom lubridate today
#' @importFrom readr write_tsv
#' @export

save_to_file <- function(feature_bounds, file_name, notes){
  if (missing(file_name)) {
    file_name <- paste("feature_bounds_",
                       format(today(), "%Y%m%d"),
                       ".tsv",
                       sep = "")
  }

  if (missing(notes)) {
    notes <- ""
  }

  header <- paste("# Date: ", today(), "\n",
                  "# genetable Version: ",
                  packageDescription("genetable")$Version, "\n",
                  "# Notes:", notes, "\n#\n",
                  "# Data Dictionary:\n",
                  "#   variables type description\n",
                  "#   chr character chromosome in genome build GRCh37/hg19\n",
                  "#   strand character DNA strand. values {+,-}\n",
                  "#   gene_id character gene identifier with ENSG prefix\n",
                  "#   gene_name character Gene name\n",
                  "#   agg_start integer 5' start position of the feature in",
                  " genome build GRCh37/hg19\n",
                  "#   agg_end integer 3' end position of the gene in genome",
                  " build GRCh37/hg19\n",
                  "#   source comma seprated list annotation sources of",
                  " transcripts in an aggregation unit\n",
                  "#   type character biotype of the feature\n",
                  "#   merge_count integer number of features joined in",
                  " aggregation unit\n",
                  "#   agg_size integer length of aggregation unit\n",
                  sep = "")

  fileconn <- file(file_name)
  writeLines(header, fileconn)
  close(fileconn)

  write_tsv(feature_bounds,
            file_name,
            append = TRUE,
            col_names = TRUE)

  return(invisible(feature_bounds))
}
