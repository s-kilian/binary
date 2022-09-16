context("Tests for RD and RR teststat, critval, power, appr. samplesize, and exact sample size")

test_that("Functions work", {
  data.frame(
    eff_meas = c("RD", "RD", "RR", "RR", "RR"),
    better = c("high", "low", "high", "low", "high"),
    n_E = rep(15, 5),
    n_C = rep(23, 5),
    delta = c(-0.1, 0.8, 0.9, 4, 1/6),
    alpha = rep(0.025, 5),
    power = rep(0.8, 5),
    size_acc = rep(2, 5),
    p_EA = c(rep(0.5, 4), 0.3),
    p_CA = c(rep(0.3, 4), 0.5),
    r = c(rep(1.2, 4), 0.8)
  ) ->
    df
  
  for(i in 1:5){
    ts <- default_ts(eff_meas = df$eff_meas[i], better = df$better[i])
    critval(
      alpha = df$alpha[i],
      n_E = df$n_E[i],
      n_C = df$n_C[i],
      delta = df$delta[i],
      eff_meas = df$eff_meas[i],
      size_acc = df$size_acc[i],
      test_stat = ts
    )
    samplesize_appr(
      p_EA = df$p_EA[i],
      p_CA = df$p_CA[i],
      delta = df$delta[i],
      eff_meas = df$eff_meas[i],
      better = df$better[i],
      alpha = df$alpha[i],
      beta = 1-df$power[i],
      r = df$r[i]
    )
    samplesize_exact(
      p_EA = df$p_EA[i],
      p_CA = df$p_CA[i],
      delta = df$delta[i],
      eff_meas = df$eff_meas[i],
      better = df$better[i],
      alpha = df$alpha[i],
      beta = 1-df$power[i],
      r = df$r[i],
      size_acc = df$size_acc[i]
    )
  }
})
