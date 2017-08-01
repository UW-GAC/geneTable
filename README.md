# genetable

The goal of genetable is to provide some utility functions for making a database table of gene boundary definitions. For now, it is focused on GeneCode gene definitions, as defined in the GENCODE .gtf files, such as availble at https://www.gencodegenes.org/releases/

## Installation

You can install genetable from github with:

```R
# install.packages("devtools")
devtools::install_github("UW-GAC/genetable")
```

## Example

This example assumes you've downloaded the gencode gtf file to the working directory.

```R
# install genetable package
devtools::install_github("UW-GAC/genetable")

# load libraries for convenience
library(genetable)
library(tidyverse)

path <- "gencode.v19.annotation.gtf.gz"

# import the gtf file to a tidy data frame (a tibble)
# this is slow
gtf <- import_gencode(path)

# look at the tibble
glimpse(gtf)

# summarize the number of features by tag.
summarize_tag(gtf, tag = "basic")

# filter gtf file to return transcript features tagged basic
basic_transcripts <- filter_gencode(gtf, featurearg = "transcript", tagarg = "basic")

# or filter for features == "gene"
genes <- filter_gencode(gtf, featurearg = "gene")

# define the boundaries of the feature of interest
# this is slow
#gene_bounds <- define_boundaries(basic_transcripts, "gene_id")
gene_bounds <- define_boundaries(genes, "gene_id")

# can check the resulting tibble for sanity
glimpse(gene_bounds)

# save to file
note <- 'This file includes starting and ending ranges for feature = "gene" in the gtf file.'
save_to_file(gene_bounds, notes = note) # will automatically make file called feature_bounds_DATE.tsv
```

## Additional Information

We're documenting our process of creating this package. You may be interested in reading information at [https://bheavner.github.io/new_package.html](https://bheavner.github.io/new_package.html).
