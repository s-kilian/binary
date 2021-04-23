context("Test calculation of p-values")

p_value(
  3, 4, 10, 11, "OR", 0.9, size_acc = 2, better = "high"
)$p_max
