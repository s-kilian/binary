

# FUnction to create Table for Reference Values

data.frame(
  method = c(rep("RD", 10), rep("RR", 10)),
  better = c(rep("high",5), rep("low",5), rep("high", 5), rep("low", 5)),
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

