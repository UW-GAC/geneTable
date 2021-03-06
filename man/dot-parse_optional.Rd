% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import-gencode.r
\name{.parse_optional}
\alias{.parse_optional}
\title{Intermediate parsing - parse optional fields from "optional" column of
intermediate_gtf (which originated from the 9th column of the GENCODE gtf
file).}
\usage{
.parse_optional(intermediate_gtf)
}
\arguments{
\item{intermediate_gtf}{A tibble from .pull_required().}
}
\value{
A tibble with columns for optional fields from .gtf file, including
variables that may include multiple values.
}
\description{
Intermediate parsing - parse optional fields from "optional" column of
intermediate_gtf (which originated from the 9th column of the GENCODE gtf
file).
}
