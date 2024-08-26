test_that("hello_ctrex outputs the correct message.", {
  testthat::expect_output(hello_ctrex(), "Welcome to the ctrex package.")
})
