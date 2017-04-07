context("test_import_gencode - unit tests")

test_that("import-gencode returns an error if the file path is not defined", {
            expect_error(import_gencode(), "path not defined")
})

test_that("import-gencode returns an error if disallowed feature is given", {
            err <- paste("feature must be one of gene, transcript, exon, ",
                         "CDS, UTR, start_codon, stop_codon, or ",
                         "Selenocysteine", sep = "")
            expect_error(import_gencode(path = "path/to/file.tsv",
                                        feature = "foo"), err)
})

test_that("import-gencode returns an error if the attribute is not a string", {
            expect_error(import_gencode(path = "pat/to/file.tsv",
                                        attribute = 7),
                         "attribute must be a string")
})
