context("test_filter_gencode - unit tests")

test_that("filter_gencode returns an error if the file path is not defined", {
            expect_error(filter_gencode(), "gtf not defined")
})

test_that("filter-gencode returns an error if disallowed feature is given", {
            gtf <- import_gencode("gencode_test.tsv")
            err <- paste("feature must be one of gene, transcript, exon, ",
                         "CDS, UTR, start_codon, stop_codon, or ",
                         "Selenocysteine", sep = "")
            expect_error(filter_gencode(gtf, feature = "foo"), err)
})

test_that("filter-gencode returns a tibble", {
            gtf <- import_gencode("gencode_test.tsv")
            expect_true(is.tibble(filter_gencode(gtf)))
})

test_that("filter-gencode returns the right size tibble for known input", {
            gtf <- import_gencode("gencode_test.tsv")
            expect_equal(dim(dplyr::filter(gtf,
                                    feature == "exon",
                                    tag == "basic")),
                         c(5, 16))
})
