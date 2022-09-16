context("Tests for samplesize_appr and samplesize_exact.")

test_that("samplesize_exact works.",{
  expect_silent(samplesize_exact(
    p_EA = 0.7,
    p_CA = 0.4,
    delta = 0,
    alpha = 0.025,
    beta = 0.2,
    r = 1,
    size_acc = 3,
    eff_meas = "RD",
    test_stat = NULL,
    better = "high"
  ))
})
