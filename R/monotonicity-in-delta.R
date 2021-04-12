##..............................................................................
##               Project: Binary exact NI tests
##               Purpose: Evaluate monotonicity of p-value in NI margin
##                 Input: 
##                Output: 
##      Date of creation: 2021-03-26
##   Date of last update: 2021-03-26
##                Author: Samuel Kilian
##..............................................................................

# Röhmel 2005 showed that the p-value of Chan's exact risk difference test
# is not decreasing in delta.

source("R/Non_Inferiority_RD_RR.r")

# Values from Röhmel's example. H0: high values of x_E are better
n_E <- 170
n_C <- 248
x_E. <- 76
x_C. <- 130
delta.vec <- seq(0.02, 0.028, by = 0.001)

# Example doesn't work. Observed risk difference is 0.077 and not able to
# produce a small p-value with a NI-margin of 0.02
x_E./n_E - x_C./n_C

# Alternative NI margin
delta.vec <- seq(0.15, 0.16, by = 0.001)


p.max.vec <- c()
for (delta in delta.vec) {
  p_value(
    x_E. = x_E.,
    x_C. = x_C.,
    n_E = n_E,
    n_C = n_C,
    method = "RD",
    delta = -delta,
    size_acc = 2,
    better = "low"
  ) ->
    res
  
  p.max.vec <- c(p.max.vec, max(res))
}

plot(delta.vec, p.max.vec)

p_value(
  x_E. = x_E.,
  x_C. = x_C.,
  n_E = n_E,
  n_C = n_C,
  method = "RD",
  delta = -0.02,
  size_acc = 2,
  better = "low"
) %>%
  max()
