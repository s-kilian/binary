context("Test calculation of p-values")

p_value(
  7, 4, 100, 110, "RD", -0.1, size_acc = 3, better = "high"
)$p_max
critval(alpha = 0.05, 110, 100, method = "RD", delta = -0.1, size_acc = 3, better = "high")

n_E <- 100
n_C <- 100
