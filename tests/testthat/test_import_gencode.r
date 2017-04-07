context("test_import_gencode - unit tests")

test_that("import-gencode returns an error if the file path is not defined", {
              expect_error(import_gencode(), "path not defined")
              })
