context("Tests for test statistic construction and default test statistics.")

test_that("new_exbin_ts works",{
  new_exbin_ts() -> ts.1
  new_exbin_ts(
    ts_fun = function(x_E, x_C, n_E, n_C) x_E/n_E-x_C/n_C,
    ts_fun_args = list(n_E = 10, n_C = 12),
    quant_fun = qnorm,
    quant_fun_args = list(mean = 0.1)
  ) -> ts.2
  for (ts in list(ts.1, ts.2)) {
    expect_equal(length(ts), 4)
    expect_s3_class(ts, "exbin_ts")
  }
})

test_that("validate_exbin_ts work",{
  # ts_fun not a function
  new_exbin_ts(
    ts_fun = "function"
  ) ->
    ts.1
  # ts_fun a function with false argument names
  new_exbin_ts(
    ts_fun = function(x, y) x-y
  ) ->
    ts.2
  # ts_fun_args not a list
  new_exbin_ts(
    ts_fun = function(x_E, x_C) 0,
    ts_fun_args = c(n_E = 12)
  ) ->
    ts.3
  # quant_fun a function with false argument names
  new_exbin_ts(
    ts_fun = function(x_E, x_C) 0,
    ts_fun_args = list(n_E = 12),
    quant_fun = function(q) q
  ) ->
    ts.4
  # quant_fun_args not a list
  new_exbin_ts(
    ts_fun = function(x_E, x_C) 0,
    ts_fun_args = list(n_E = 12),
    quant_fun = function(p) p,
    quant_fun_args = c(q = 13)
  ) ->
    ts.5
  # correct specification
  new_exbin_ts(
    ts_fun = function(x_E, x_C) x_E*0,
    ts_fun_args = list(n_E = 12),
    quant_fun = function(p) p,
    quant_fun_args = list(q = 13)
  ) ->
    ts.6
  
  for (ts in list(ts.1, ts.2, ts.3, ts.4, ts.5)) {
    expect_error(validate_exbin_ts(ts))
  }
  validate_exbin_ts(ts.6)
})

test_that("exbin_ts work",{
  expect_error(exbin_ts())
  exbin_ts(
    ts_fun = function(x_E, x_C, n_E, n_C) x_E/n_E-x_C/n_C
  ) ->
    ts.1
  exbin_ts(
    ts_fun = function(x_E, x_C, n_E, n_C) x_E/n_E-x_C/n_C,
    ts_fun_args = list(n_E = 12, n_C = 12, useless = "bla"),
    quant_fun = function(p) p,
    quant_fun_args = list(useless = "blabla")
  ) ->
    ts.2
  for (ts in list(ts.1, ts.2)) {
    expect_equal(length(ts), 4)
    expect_s3_class(ts, "exbin_ts")
  }
})

test_that("calc_ts.exbin_ts and calc_quant.exbin_ts work",{
  exbin_ts(
    ts_fun = function(x_E, x_C, n_E, n_C) x_E/n_E-x_C/n_C,
    ts_fun_args = list(n_E = 10),
    quant_fun = function(p, fact) fact*p,
    quant_fun_args = list(fact = 2, false.fact = 7)
  ) ->
    ts.1
  
  # calculate test statistic
  calc_ts.exbin_ts(
    ts = ts.1,
    x_E = 3:6,
    x_C = 1:4,
    n_C = 10
  ) ->
    ts.1.ts
  expect_equal(ts.1.ts, c(0.2, 0.2, 0.2, 0.2))
  
  # calculate quantile
  calc_quant.exbin_ts(
    ts = ts.1,
    p = 0.07,
    fact = 3
  ) ->
    ts.1.quant
  expect_equal(ts.1.quant, 0.21)
  
  # use predefined factor
  calc_quant.exbin_ts(
    ts = ts.1,
    p = 0.07
  ) ->
    ts.1.quant
  expect_equal(ts.1.quant, 0.14)
  
  # argument n_C not given
  expect_error(
    calc_ts.exbin_ts(
      ts = ts.1,
      x_E = 3:6,
      x_C = 1:4
    )
  )
  
  # arguments have unequal length
  expect_error(
    calc_ts.exbin_ts(
      ts = ts.1,
      x_E = 3:6,
      x_C = 1:3,
      n_C = 10
    )
  )
  
  # arguments have wrong type
  expect_error(
    calc_ts.exbin_ts(
      ts = ts.1,
      x_E = letters[1:3],
      x_C = 1:3,
      n_C = 10
    )
  )
  
  # argument p not given
  expect_error(
    calc_quant.exbin_ts(
      ts = ts.1,
      fact = 0.3
    )
  )
  
  # argument has wrong type
  expect_error(
    calc_quant.exbin_ts(
      ts = ts.1,
      p = "a"
    )
  )
  
  # argument is in wrong range
  expect_error(
    calc_quant.exbin_ts(
      ts = ts.1,
      p = 1.2
    )
  )
})

# test_that("default test statistics work",{
# 
# })
