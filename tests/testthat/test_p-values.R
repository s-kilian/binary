context("Test calculation of p-values")

n_E <- 10
n_C <- 10
delta <- -0.1
method <- "RD"
better <- "high"
size_acc <- 3
alpha <- 0.05

expand.grid(
  x_E = 0:n_E,
  x_C = 0:n_C
) %>%
  teststat(
    n_E = n_E,
    n_C = n_C,
    delta = delta,
    method = method,
    better = better
  ) %>%
  mutate(
    reject = stat >= qnorm(1-alpha)
  ) ->
  df
find_max_prob(
  df = df,
  n_E = n_E,
  n_C = n_C,
  method = method,
  delta = delta,
  calc_method = "grid search",
  size_acc = size_acc
)
x_E. <- 7
x_C. <- 3
p_value(
  x_E. = x_E.,
  x_C. = x_C.,
  n_E = n_E,
  n_C = n_C,
  better = better,
  method = method,
  delta = delta,
  calc_method = "uniroot",
  size_acc = size_acc
)

# Test confidence region
conf_region(
  x_E. = x_E.,
  x_C. = x_C.,
  n_E = n_E,
  n_C = n_C,
  alpha = 0.05,
  method = method,
  delta_acc = 2,
  size_acc = size_acc
)

