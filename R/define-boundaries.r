#' Define the genomic region spanning all the feature types of interest.
#'
#' This function processes a filtered gencode tibble (such as returned from 
#' filter_gencode()) to define the genomic region that spans all the features 
#' of interest.
#'
#' NOTE: for now, this works for grouping on gene_id. I may generalize later
#' if we want to group on other features.
#'
#' @param filtered_gtf The filtered gencode tibble
#' @param feature The feature type of interest. Must be a field of the
#'  filtered_gtf tibble (e.g., feature %in% names(filtered_gtf).
#'  Default = "gene_id".
#' @return A tibble containing the regions of interest. Variables in this 
#'  tibble include:
#'  - chr - chromosome in genome build GRCh37/hg19
#'  - strand - DNA strand. values {+,-}
#'  - gene_id - gene identifier with ENSG prefix
#'  - gene_name - Gene name
#'  - agg_start - lower bound of genic unit
#'  - agg_end - upper bound of genic unit
#'  - source - comma seprated lists values of annotation sources of transcripts
#'    in a genic unit
#'  - transcript_type - comma seprated values of biotypes of transcripts in a
#'    genic unit
#'  - merge_count - number of transcripts merged in the genic unit
#'  - agg_size - length of the genic unit derived as a difference of agg_end 
#'    agg_start

#' @examples
#' \dontrun{
#' define_boundaries(filtered_gtf, feature = "gene_id")
#' }
#' @importFrom tibble is.tibble
##' @import dplyr
#' @importFrom purrr map_chr
#' @importFrom tidyr nest
#' @importFrom stringr str_replace_all
#' @export

define_boundaries <- function(filtered_gtf,
                           grouping = "gene_id"){
  # check arguments
  if (missing(filtered_gtf)) stop("filtered gtf not defined")

  # select fields of interest from gtf and group by feature
  grouped <- filtered_gtf %>%
    select(chr = seqname,
           strand,
           gene_id,
           gene_name,
           start,
           end,
           source,
           transcript_type
           )

  grouped <- group_by_(grouped, "chr", grouping) # apparently pipes break env

  # reformat the source and transcript_type variables from tidy to comma-
  # separated list
  reformatted <- .parse_sources(grouped)

  # add agg_start, agg_end, merge_count, agg_size variables
  expanded <- .add_ranges(reformatted)

  # remove duplicate rows from the expanded gtf tibble
  return(
         distinct(expanded, .keep_all = TRUE)
         )

}

#
#------------------------------------------------------------------------------
#
#' Parse source and transcript_types columns to merge rows for transcripts with
#' same transcript_id. In such cases, unique sources and transcript_types
#' should be merged to comma-separated lists in the source and transcript_type
#' variables.
#'
#' @param grouped_gtf a filtered gtf field grouped by feature
#' @return a tibble

.parse_sources <- function(grouped_gtf){
  listed_source <- grouped_gtf %>%
      nest(source, .key = source) %>%
      mutate(source = map_chr(.$source,
                              function(x){
                                stringr::str_c(unique(x),
                                               collapse = ", ")
                              }
                              )
      ) %>%
      mutate(source = str_replace_all(.$source, '(c\\(\")|"|,|\\)', "")) %>%
      mutate(source = str_replace_all(.$source, " ", ", "))

    listed_transcript_type <- grouped_gtf %>%
      nest(transcript_type, .key = transcript_type) %>%
      mutate(transcript_type = map_chr(.$transcript_type,
                                       function(x) {
                                         stringr::str_c(unique(x),
                                                        collapse = ", ")
                                       }
                                       )
      ) %>%
      mutate(transcript_type = str_replace_all(
                                               .$transcript_type,
                                               '(c\\(\")|"|,|\\)', ""
                                               )
      ) %>%
      mutate(transcript_type = str_replace_all(.$transcript_type, " ", ", "))

    # add the nested columns to grouped
    with_listed <- grouped_gtf %>%
      inner_join(listed_source,
                 by = c("chr", "gene_id")) %>%
      inner_join(listed_transcript_type, by = c("chr", "gene_id")) %>%
      select(-source.x, -transcript_type.x)

    return(with_listed)
}

#
#------------------------------------------------------------------------------
#
#' Get feature ranges and add agg_start, agg_end, and merge_count variables
#' @param reformatted_gtf a filtered gtf field grouped by feature with comma-
#'   separated lists for the source and transcript_type fields
#' @return a tibble 

.add_ranges <- function(reformatted_gtf){
  summary_columns <- reformatted_gtf %>%
    summarize(agg_start = min(start),
              agg_end = max(end),
              merge_count = n())
  # add summary columns to final table
  combined <- reformatted_gtf %>% inner_join(summary_columns)

  # add agg_size derived variable
  all <- combined %>%
    mutate(agg_size = agg_end - agg_start) %>%
    select(chr,
           strand,
           gene_id,
           gene_name,
           agg_start,
           agg_end,
           source = source.y,
           transcript_type = transcript_type.y,
           merge_count,
           agg_size)

    return(all)
}
