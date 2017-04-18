context("test_import_gencode - unit tests")

test_that("import-gencode returns an error if the file path is not defined", {
            expect_error(import_gencode(), "path not defined")
})

test_that("import-gencode returns a tibble", {
            expect_true(tibble::is.tibble(import_gencode("gencode_test.tsv")))
})

test_that("import-gencode returns 16 columns", {
            expect_equal(ncol(import_gencode("gencode_test.tsv")), 16)
})
