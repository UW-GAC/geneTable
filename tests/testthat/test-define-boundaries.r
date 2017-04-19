context("test_define_boundaries - unit tests")

test_that(
  "define-boundaries returns an error if the file path is not defined", {
            expect_error(define_boundaries(), "filtered gtf not defined")
})

test_that("define-boundaries returns a tibble", {
            gtf <- import_gencode("gencode_test.tsv")
            filtered <- filter(gtf, feature == "exon", tag == "basic")
            expect_true(is.tibble(define_boundaries(filtered)))
})

#test_that("filter-gencode returns the right size tibble for known input", {
 #           gtf <- import_gencode("gencode_test.tsv")
 #           expect_equal(dim(dplyr::filter(gtf,
 #                                   feature == "exon",
 #                                   tag == "basic")),
 #                        c(5, 16))
#})
