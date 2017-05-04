#' Import and tidy a gencode .gtf file 
#'
#' The gtf file format is described at 
#' https://www.gencodegenes.org/data_format.html. This function imports a gtf
#' and returns a tidied tibble, in which each column is a variable and each
#' row is an observation
#'
#' @param path Path to the gencode file.
#' @return A tidied tibble
#' @examples
#' \dontrun{
#' import_gencode(path = 'path/to/file.tsv')
#' }
#' @importFrom tidyr extract
#' @importFrom readr read_tsv
#' @importFrom stringr str_trim
#' @importFrom stringr str_replace_all
#' @importFrom tidyr unnest
#' @export

import_gencode <- function(path, featuretag = "basic"){
  # check arguments
  if (missing(path)) stop("path not defined")

  # read the file to make a tibble for tidyverse work.
  raw_gtf <- .read_gtf(path)

  # pick off the required fields from the 9th column.
  intermediate_gtf <- .pull_required(raw_gtf)

  # Next, peel off the remaining optional fields.
  messy_gtf <- .parse_optional(intermediate_gtf)

  # Finally, tidy the messy tibble (that has fields with 0 to many entries)
  tidied_gtf <- .tidy_lists(messy_gtf)

  return(tidied_gtf)
}

#
#------------------------------------------------------------------------------
#
#' Read GENCODE gtf file to tibble. Assign fields based on tab separation
#' with no further parsing.
#'
#' @param path PATH to gencode file.
#' @return A tibble

.read_gtf <- function(path){
  suppressMessages(read_tsv(path,
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
}

#
#------------------------------------------------------------------------------
#
#' Preliminary parsing for the attribute column of raw gtf tibble: use regular
#' expression to pick off the easier required fields, leave rest in "optional"
#'
#' @param raw_gtf A tibble from .read_gtf().
#' @return A tibble with attribute column partially parsed and replaced with
#' the "optional" column with the remaining unparsed fields.

.pull_required <- function(raw_gtf){
  # define regular expression to use with tidyr::extract to get required fields
  # from attribute column, and split off optional for further tidying
  rx <- paste('gene_id "(.+?)"; ', # required gene_id field in gtf
                'transcript_id "(.+?)"; ', # required in gtf
                'gene_type "(.+?)";', # required in gtf
                ".*", # gene_status removed from newer releases
                'gene_name "(.+?)"; ', # required in gtf
                'transcript_type "(.+?)";', # required in gtf
                ".*", # transcript_status removed from newer releases
                'transcript_name "(.+?)"; ', # required in gtf
                "(?:exon_number (.+?); )?", # exon_number (not in all lines)
                '(?:exon_id "(.+?)"; )?', # exon_id (not in all lines)
                "level ([123]); ?", # required in gtf
                "(.*)$", # additional optional fields
                sep = "") # join regex string without spaces

  # parse the attribute column of gtf with regular expression - pick off
  # the easier required fields, leave rest in "optional"
  intermediate_gtf <- raw_gtf %>% extract(attribute,
                c("gene_id",
                  "transcript_id",
                  "gene_type",
                  "gene_name",
                  "transcript_type",
                  "transcript_name",
                  "exon_number",
                  "exon_id",
                  "level",
                  "optional"
                  ),
                rx, remove = FALSE
                ) %>% select(-attribute)
  return(intermediate_gtf)
}

#
#------------------------------------------------------------------------------
#
#' Intermediate parsing - parse optional fields from "optional" column of 
#' intermediate_gtf (which originated from the 9th column of the GENCODE gtf
#' file).
#'
#' @param intermediate_gtf A tibble from .pull_required().
#' @return A tibble with columns for optional fields from .gtf file, including
#' variables that may include multiple values.

.parse_optional <- function(intermediate_gtf){
# TODO: refactor this - lots of repeated typing
  messy_gtf <-
  intermediate_gtf %>% extract(optional,
                               "tags",
                               '((?:tag "(?:.+?)"; ?)+)',
                               remove = FALSE
                               ) %>%
  extract(optional, "ccdsids",
          '((?:ccsid "(?:.+?)"; ?)+)',
          remove = FALSE
  ) %>%
  extract(optional, "havana_gene",
          'havana_gene "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "havana_transcript",
          'havana_transcript "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "protein_id",
          'protein_id "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "onts",
          '((?:ont "(?:.+?)"; ?)+)',
          remove = FALSE
  ) %>%
  extract(optional, "transcript_support_level",
          'transcript_support_level "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "remap_status",
          'remap_status "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "remap_original_id",
          'remap_original_id "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "remap_original_location",
          'remap_original_location "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "remap_num_mappings",
          'remap_num_mappings "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "remap_target_status",
          'remap_target_status "(.+?)"',
          remove = FALSE
  ) %>%
  extract(optional, "remap_substituted_missing_target",
          'remap_substituted_missing_target "(.+?)"',
          remove = TRUE
  )

  return(messy_gtf)
}

#
#------------------------------------------------------------------------------
#
#' Tidy fields with 0 to many entries: tags, ccdsids, onts.
#'
#' @param messy_gtf A tibble from .parse_optional().
#' @return A tidied tibble with no list-columns.

.tidy_lists <- function(messy_gtf){
  trimmed_gtf <- messy_gtf %>%
  # remove leading and trailing whitespace from multi-value fields
  mutate(tags = str_trim(tags, side = "both"),
         ccdsids = str_trim(ccdsids, side = "both"),
         onts = str_trim(onts, side = "both")) %>%
  # remove quotes, semicolons, and key strings from multi-value fields
  mutate(tags = str_replace_all(tags,
                                c('"' = "", # remove "
                                  ";" = "", # remove ;
                                  "tag " = "") # remove 'tag 's
                                ),
         ccdsids = str_replace_all(ccdsids,
                                   c('"' = "",
                                     ";" = "",
                                     "ccdsid " = "")
                                   ),
         onts = str_replace_all(onts,
                                c('"' = "",
                                  ";" = "",
                                  "ont " = "")
                                )) %>%
  # split the strings by single spaces to make list-columns
  mutate(tags = strsplit(tags, " "),
         ccdsids = strsplit(ccdsids, " "),
         onts = strsplit(onts, " "))

  # finally, unnest the list-columns
  # (have to do individually b/c different #s in each list)
  unnested_gtf <-
  trimmed_gtf %>%
  unnest(ccdsids, .drop = FALSE) %>%
  unnest(onts, .drop = FALSE) %>%
  unnest(tags, .drop = FALSE) %>%
  rename(tag = tags, ccdsid = ccdsids, ont = onts)

  return(unnested_gtf)
}
