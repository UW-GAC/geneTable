% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save-to-file.r
\name{save_to_file}
\alias{save_to_file}
\title{Save the tibble defining feature ranges to a tab-delimited text file.}
\usage{
save_to_file(feature_bounds, file_name, notes)
}
\arguments{
\item{feature_bounds}{a tibble such as produced by define_boundaries()}

\item{file_name}{the target file name for writing (default 
"feature_bounds_DATE.tsv")}

\item{notes}{a character string of notes to include in the file header.}
}
\description{
Save the tibble of genomic feature ranges to a tab-delimited file with
appropriate header columns.
}
\details{
The output file includes the following fields as comments in the header:
  Date: $DATE
  genetable Version: $VERSION
  Notes: $NOTES
  Data Dictionary:
    variable type description
    chr character chromosome in genome build GRCh37/hg19
    start integer 5' start position of the gene in genome build GRCh37/hg19
    end integer 3' end position of the gene in genome build GRCh37/hg19
    strand character DNA strand. values {+,-}
    gene_id character gene identifier with ENSG prefix
    gene_name character Gene name
    type character biotype of the feature
}
\examples{
\dontrun{
notes <- paste("Features selected by filtering gtf for feature =",
               '"transcript" and tag = "basic"', sep = " ")
save_to_file(gene_bounds, file_name = "gene_bounds_20170504", 
             notes = notes)
}
}
