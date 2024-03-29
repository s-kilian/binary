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

## Check consistency of delta with method ######################################
test_that("Out of range delta gives error", {
  # Risk difference
  ## error
  for (delta in c(-2, -1, 1, 1.1)) {
    for (better in c("high", "low")) {
      expect_error(
        check.delta.method.better(
          delta = delta,
          method = "RD",
          better = better
        )
      )
    }
  }
  ## no error
  for (delta in c(-0.1, -0.9)) {
    expect_silent(
      check.delta.method.better(
        delta = delta,
        method = "RD",
        better = "high"
      )
    )
  }
  
  # Risk ratio
  ## error
  for (delta in c(-1, 0)) {
    for (better in c("high", "low")) {
      expect_error(
        check.delta.method.better(
          delta = delta,
          method = "RR",
          better = better
        )
      )
    }
  }
  ## no error
  for (delta in c(0.1, 0.9)) {
    expect_silent(
      check.delta.method.better(
        delta = delta,
        method = "RR",
        better = "high"
      )
    )
  }
  
  # Odds ratio
  ## error
  for (delta in c(-1, 0)) {
    for (better in c("high", "low")) {
      expect_error(
        check.delta.method.better(
          delta = delta,
          method = "OR",
          better = better
        )
      )
    }
  }
  ## no error
  for (delta in c(0.1, 0.9)) {
    expect_silent(
      check.delta.method.better(
        delta = delta,
        method = "OR",
        better = "high"
      )
    )
  }
})
test_that("Unfitting delta gives warning",{
  # warnings
  data.frame(
    method = c("RD", "RD", "RD", "RD", "RR", "RR", "RR", "RR", "OR", "OR", "OR", "OR"),
    better = c("high", "high", "low", "low", "high", "high", "low", "low", "high", "high", "low", "low"),
    delta = c(0.1, 0.9, -0.2, -0.7, 1.1, 100, 0.8, 0.001, 1.1, 100, 0.8, 0.001)
  )->
    df
  for (i in 1:nrow(df)) {
    expect_warning(
      check.delta.method.better(
        delta = df$delta[i],
        method = df$method[i],
        better = df$better[i]
      )
    )
  }
  # no warnings
  data.frame(
    method = c("RD", "RD", "RD", "RD", "RR", "RR", "RR", "RR", "OR", "OR", "OR", "OR"),
    better = c("high", "high", "low", "low", "high", "high", "low", "low", "high", "high", "low", "low"),
    delta = c(-0.1, -0.9, 0.2, 0.7, 0.8, 0.001, 1.1, 100, 0.8, 0.001, 1.1, 100)
  )->
    df
  for (i in 1:nrow(df)) {
    expect_silent(
      check.delta.method.better(
        delta = df$delta[i],
        method = df$method[i],
        better = df$better[i]
      )
    )
  }
})

## Check behavior of delta.null ###############################################
test_that("Check handling of delta = NULL",{
  # warnings and default delta setting
  data.frame(
    method = c("RD", "RR", "OR"),
    default.delta = c(0, 1, 1)
  ) ->
    df
  for(i in 1:nrow(df)){
    expect_warning(
      check.delta.null(
        method = df$method[i],
        delta = NULL
      )
    )
    expect_equal(
      check.delta.null(
        method = df$method[i],
        delta = NULL
      ),
      df$default.delta[i]
    )
  }
  
  # no changes
  for (method in c("RD", "RR", "OR")) {
    delta <- rnorm(1)
    expect_equal(
      check.delta.null(
        method = method,
        delta = delta
      ),
      delta
    )
  }
})
