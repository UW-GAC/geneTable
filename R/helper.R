#' hack to pass devtools::check()
#' see: https://stackoverflow.com/questions/9439256/
#' @noRd
  utils::globalVariables(c("start", "end", "agg_end", "agg_start", "chr",
                           "strand", "gene_id", "gene_name", "source.y",
                           "transcript_type.y", "merge_count", "agg_size",
                           "optional", ".", "transcript_type", "source.x",
                           "feature", "tag", "packageDescription", "attribute",
                           "tags", "ccdsids", "onts", "seqname",
                           "transcript_type.x"))