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
#' @import dplyr
#' @export

define_boundaries <- function(filtered_gtf,
                           feature = "gene_id"){
  # check arguments
  if (missing(filtered_gtf)) stop("filtered gtf not defined")

  by_feature <- select(filtered_gtf, chr = seqname, strand, gene_id,
                       gene_name, start, end, source, transcript_type) %>%
  group_by(gene_id, chr, strand, gene_name, source, transcript_type) %>%
  summarize(agg_start = min(start), agg_end = max(end),
            merge_count = n()) %>%
  mutate(agg_size = agg_end - agg_start)

  return(select(by_feature, chr, strand, gene_id, gene_name, agg_start,
                agg_end, source, transcript_type, merge_count, agg_size))

}
