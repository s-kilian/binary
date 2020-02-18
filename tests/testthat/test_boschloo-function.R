context("Test results from literature")

library(tidyverse)

## Power                    ####################################################
# Superiority

test_that("Boschloo-Power can be reproduced", {
# Raised nominal level and exact maximal size: According to [Boschloo, 1970],
# the raised nominal level for sample sizes 10 and 15 and true level alpha = 0.05
# can be chosen as approximately 0.09
#Raise.level(0.05, 10, 15, 4)

# Reject probability: According to [Boschloo, 1970], the power (reject prob.)
# for sample size 15 in both groups and true rates p1 = 0.3 (0.6, 0.7, 0.8),
# p0 = 0.1(0.1, 0.1, 0.2, 0.2) and alpha = 0.01 (0.05) is:
data.frame(
  n_E = rep(15, 8),
  n_C = rep(15, 8),
  p1 = rep(c(0.3, 0.6, 0.7, 0.8), 2),
  p0 = rep(c(0.1, 0.1, 0.2, 0.2), 2),
  alpha = rep(c(0.01, 0.05), each = 4),
  power = c(
    0.1493,
    0.7394,
    0.6796,
    0.8723,
    0.3531,
    0.9189,
    0.8962,
    0.9744
  )
) ->
  df
# Compute power with function power_boschloo():
df %>%
  mutate(
    power2 = as.vector(
      sapply(
      1:length(n_E),
      function(i) power_boschloo(
        df = expand.grid(
          x_C = 0:n_C[i],
          x_E = 0:n_E[i]
        ) %>%
          teststat_boschloo(n_C = n_C[i], n_E = n_E[i]) %>%
          mutate(
            reject = cond_p <= critval_boschloo(alpha = alpha[i], n_C = n_C[i], n_E = n_E[i], size_acc = 3)["nom_alpha_mid"]
          ),
        n_C = n_C[i],
        n_E = n_E[i],
        p_CA = p0[i],
        p_EA = p1[i]
      )
    )
    )
  ) -> df

expect_equal(df$power, df$power2, tolerance = .001, scale = 1)

})



test_that("Approximate sample size can be reproduced", {
  
# Approximate sample size for chi-quare test: According to [Kieser, 2018], the
# approximate sample size for alpha = 0.025, power = 0.8, r = 1 (2) and true rates
# p0 = 0.1 (0.2, 0.3, 0.4), p1 = p0 + 0.2 is given by:
data.frame(
  p1 = rep(c(0.3, 0.4, 0.5, 0.6), 2),
  p0 = rep(c(0.1, 0.2, 0.3, 0.4), 2),
  power = rep(0.8, 8),
  r = rep(c(1, 2), each = 4),
  alpha = rep(0.025, 8),
  n.appr = c(
    124,
    164,
    186,
    194,
    147,
    189,
    213,
    219
  )
) ->
  df
# Compute approximate sample size with function samplesize_normal_appr()
df %>%
  mutate(
    n.appr.2 = as.vector(
      sapply(
      1:length(p1),
      function(i) samplesize_normal_appr(alpha = alpha[i], beta = 1-power[i], r = r[i], p_CA = p0[i], p_EA = p1[i]) %>%
        unlist() %>%
        sum()
    )
    )
  ) -> df

expect_equal(df$n.appr, df$n.appr.2, tolerance = .5, scale = 1)

})



test_that("Exact sample size can be reproduced", {
  
# Exact sample size for Fisher-Boschloo test: According to [Kieser, 2018], the
# exact sample size for alpha = 0.025, power = 0.8, r = 1 (2) and true rates
# p0 = 0.1 (0.2, 0.3, 0.4), p1 = p0 + 0.2 is given by:
data.frame(
  p1 = rep(c(0.3, 0.4, 0.5, 0.6), 2),
  p0 = rep(c(0.1, 0.2, 0.3, 0.4), 2),
  power = rep(0.8, 8),
  r = rep(c(1, 2), each = 4),
  alpha = rep(0.025, 8),
  n.ex = c(
    126,
    168,
    190,
    204,
    147,
    186,
    213,
    219
  )
) ->
  df
# Compute exact sample size with function samplesize_exact_boschloo()
df %>%
  mutate(
    n.ex.2 = as.vector(
      sapply(
      1:length(p1),
      function(i) samplesize_exact_boschloo(alpha = alpha[i], beta = 1-power[i], r = r[i], p_CA = p0[i], p_EA = p1[i], size_acc = 3) %>%
        (function(x) x[["n_C"]]+x[["n_E"]])
    ))
  ) -> df


expect_equal(df$n.ex, df$n.ex.2, tolerance = .5, scale = 1)

})




test_that("Second exact sample size can be reproduced", {


# Exact sample size for Fisher-Boschloo test: According to [Kieser, 2018], the
# exact sample size for alpha = 0.025, power = 0.8, r = 1 (2) and true rates
# p0 = 0.1 (0.2, 0.3, 0.4), p1 = p0 + 0.2 is given by:
data.frame(
  p1 = rep(c(0.3, 0.4, 0.5, 0.6), 2),
  p0 = rep(c(0.1, 0.2, 0.3, 0.4), 2),
  power = rep(0.8, 8),
  r = rep(c(1, 2), each = 4),
  alpha = rep(0.025, 8),
  n.ex = c(
    126,
    168,
    190,
    204,
    147,
    186,
    213,
    219
  )
) ->
  df
# Compute exact sample size with function samplesize_exact_boschloo()
df %>%
  mutate(
    n.ex.2 = as.vector(
      sapply(
      1:length(p1),
      function(i) samplesize_exact_boschloo(alpha = alpha[i], beta = 1-power[i], r = r[i], p_CA = p0[i], p_EA = p1[i], size_acc = 3) %>%
        (function(x) x[["n_C"]]+x[["n_E"]])
    ))
  ) -> df

expect_equal(df$n.ex, df$n.ex.2, tolerance = .5, scale = 1)

})




test_that("Non-inferiority example", {

# Non-inferiority
# Raised nominal level: According to [Wellek, 2010], the raised nominal level
# for delta = 1/3 (1/2), n_C = n_E = 10 (25, 50), alpha = 0.05 is (Table 6.22, p.178):
data.frame(
  n_C = rep(c(10, 25, 50), 2),
  n_E = rep(c(10, 25, 50), 2),
  alpha = c(0.05, 0.05084, 0.05, 0.05, 0.05065, 0.05),
  delta = rep(c(2/3, 1/2), each = 3),
  nom.alpha = c(
    0.08445,
    0.09250,
    0.07716,
    0.10343,
    0.07445,
    0.06543
  ),
  size = c(
    0.04401,
    0.05083,
    0.04991,
    0.04429,
    0.05064,
    0.04999
  ),
  nom.alpha.2 = NA,
  size.2 = NA
) ->
  df
# Calculate raised nominal level with function Raise.level.NI():
for (i in 1:nrow(df)) {
  suppressWarnings(critval_boschloo_NI(
    alpha = df$alpha[i],
    n_C = df$n_C[i],
    n_E = df$n_E[i],
    gamma = df$delta[i],
    size_acc = 3
  )) ->
    result
  df$nom.alpha.2[i] <- result[["nom_alpha_mid"]]
  df$size.2[i] <- result[["size"]]
}


expect_equal(df$nom.alpha, df$nom.alpha.2, tolerance = .5, scale = 1)


expect_equal(df$size, df$size.2, tolerance = .5, scale = 1)


})
