context("Test of input check functions")

## Check (0, 1) ################################################################
test_that("Interval check gives appropriate errors", {
  expect_error(check.0.1(c(0.5, 1), "errmess"), "errmess")
  expect_error(check.0.1(c(0.5, 0), "errmess"), "errmess")
  expect_error(check.0.1(c(0.5, -1), "errmess"), "errmess")
  check.0.1(c(0.1, 0.3, 0.99), "error")
})

## Check positive integer ######################################################
test_that("Positive integer check gives appropriate errors", {
  expect_error(check.pos.int(c(0, 1, 2), "errmess"), "errmess")
  expect_error(check.pos.int(c(-1, 1, 2), "errmess"), "errmess")
  expect_error(check.pos.int(c(1, 2, 3.7), "errmess"), "errmess")
  check.pos.int(c(1, 2, 10))
})