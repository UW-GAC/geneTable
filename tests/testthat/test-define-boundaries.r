context("test_define_boundaries - unit tests")

test_that(
  "define-boundaries returns an error if the file path is not defined", {
            expect_error(define_boundaries(), "filtered gtf not defined")
})

test_that("define-boundaries returns a tibble", {
            load("filtered_gtf.RData")
            expect_true(is.tibble(define_boundaries(filtered)))
})

test_that("define_boundaries returns the right size tibble for known input", {
            load("filtered_gtf.RData")
            expect_equal(dim(define_boundaries(filtered)), c(41, 10))
})
