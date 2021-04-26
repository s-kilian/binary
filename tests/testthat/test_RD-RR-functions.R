context("Tests for RD and RR teststat, critval, power, appr. samplesize, and exact sample size")

test_that("Functions work", {
  data.frame(
    method = c("RD", "RD", "RR", "RR", "RR"),
    better = c("high", "low", "high", "low", "high"),
    n_E = rep(15, 5),
    n_C = rep(23, 5),
    delta = c(-0.1, 0.8, 0.9, 4, 1/6),
    alpha = rep(0.025, 5),
    power = rep(0.8, 5),
    size_acc = rep(3, 5),
    p_EA = c(rep(0.5, 4), 0.3),
    p_CA = c(rep(0.3, 4), 0.5),
    r = c(rep(2, 4), 1/2)
  ) ->
    df
  
  for(i in 1:3){
    teststat(
      df = expand.grid(
        x_E = 0:df$n_E[i],
        x_C = 0:df$n_C[i]
      ),
      n_E = df$n_E[i],
      n_C = df$n_C[i],
      delta = df$delta[i],
      method = df$method[i],
      better = df$better[i]
    )
    critval(
      alpha = df$alpha[i],
      n_E = df$n_E[i],
      n_C = df$n_C[i],
      delta = df$delta[i],
      method = df$method[i],
      better = df$better[i],
      size_acc = df$size_acc[i]
    )
    samplesize_appr(
      p_EA = df$p_EA[i],
      p_CA = df$p_CA[i],
      delta = df$delta[i],
      method = df$method[i],
      better = df$better[i],
      alpha = df$alpha[i],
      beta = 1-df$power[i],
      r = df$r[i]
    )
    samplesize_exact(
      p_EA = df$p_EA[i],
      p_CA = df$p_CA[i],
      delta = df$delta[i],
      method = df$method[i],
      better = df$better[i],
      alpha = df$alpha[i],
      beta = 1-df$power[i],
      r = df$r[i],
      size_acc = df$size_acc[i]
    )
  }
})
