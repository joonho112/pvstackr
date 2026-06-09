test_that("package scaffold loads", {
  expect_true(requireNamespace("pvstackr", quietly = TRUE))
})
