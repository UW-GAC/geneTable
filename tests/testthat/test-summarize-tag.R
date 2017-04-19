context("test_summarize_tag - unit tests")

test_that("summarize-tag returns an error if the gtf tibble is not defined", {
            expect_error(summarize_tag(), "gtf not defined")
})

test_that("summarize-tag returns a table", {
            gtf <- import_gencode("gencode_test.tsv")
            expect_true(is.table(summarize_tag(gtf)))
})
