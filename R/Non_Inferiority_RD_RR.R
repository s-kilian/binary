##..............................................................................
##               Project: Package binary (working title)
##               Purpose: Provide template to implement exact test and sample
##                        size calculation for arbitrary test statistic
##                 Input: None
##                Output: None
##      Date of creation: 2019-07-03
##   Date of last update: 2019-07-04
##                Author: Samuel Kilian
##..............................................................................

## Functions ###################################################################

#' Calculate RD test statistic
#' 
#' \code{test_RD} returns the value of the Farrington-Manning test statistic
#' for non-inferiority of the risk difference between two proportions.
#' 
#' If higher values of $x_E$ favour the alternative hypothesis, we are interested
#' in testing the null hypothesis
#' $$H_0: p_E - p_C \le \delta ,$$
#' where the NI-margin is usually non-positive: $\delta \le 0$.
#' The test statistic for this hypothesis is
#' $$T_{\RD, \delta}(x_E, x_C) = \frac{\hat p_E - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},$$
#' where $\tilde p_C = \tilde p_C(x_E, x_C)$ is the MLE of $p_C$ and
#' $\tilde p_E = \tilde p_C + \delta$ is the MLE of $p_E$ under $p_E - p_C = \delta$.
#' High values of $T_{\RD, \delta}$ favour the alternative hypothesis.
#' 
#' @param x_E Vector of number of events in experimental group.
#' @param x_C Vector of number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of x_E favour the alternative 
#' hypothesis and "low" vice versa.
#' @return Vector of values of the RD test statistic.
#' 
#' @export
#' 
#' @examples
#' test_RD(3, 4, 10, 10, 0.2, "high")
test_RD <- function(x_E, x_C, n_E, n_C, delta, better = c("high", "low")){
  p_E <- x_E / n_E
  p_C <- x_C / n_C
  
  theta <- n_C / n_E
  p1hat <- x_E / n_E
  p2hat <- x_C / n_C
  a <- 1 + theta
  b <- -(1 + theta + p1hat + theta * p2hat + delta * (theta + 2))
  c <- delta^2 + delta * (2*p1hat + theta + 1) + p1hat + theta * p2hat
  d <- - p1hat * delta * (1 + delta)
  # Define the parameters for solving the equation
  v <- b^3/(3*a)^3 - b*c/(6*a^2) + d/(2*a)
  u <- sign(v) * sqrt(b^2/(3*a)^2 - c/(3*a))
  w <- 1/3 * (pi + acos(ifelse(u == 0, 0, round(v/u^3, 10))))   # das round wurde eingebaut, weil wenn bei v/u^3 etwas minimal gr??er als 1 rauskommt es zu einer Fehlermeldung kommt! 
  # Define the solution
  p_E0 <- 2*u*cos(w) - b/(3*a)
  p_C0 <- round(p_E0 - delta, 10)                  # round eingebaut aus gleichem Grund wie oben.
  
  denom <- ifelse(p_E - p_C - delta == 0, 1, sqrt(p_E0*(1-p_E0)/n_E + p_C0*(1-p_C0)/n_C))
  num <- p_E - p_C - delta
  if (better == "high"){
    return <-  num/denom
  }
  if (better == "low"){
    return <- -num/denom
  }
  return(return)
  
}


#' Calculate RR test statistic
#' 
#' \code{test_RR} returns the value of the Farrington-Manning test statistic
#' for non-inferiority of the risk ratio between two proportions.
#' 
#' If higher values of $x_E$ favour the alternative hypothesis, we are interested
#' in testing the null hypothesis
#' $$H_0: p_E / p_C \le \delta ,$$
#' where the NI-margin is usually smaller than 1: $\delta < 1$.
#' The test statistic for this hypothesis is
#' $$T_{\RD, \delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},$$
#' where $\tilde p_C = \tilde p_C(x_E, x_C)$ is the MLE of $p_C$ and
#' $\tilde p_E = \tilde p_C + \delta$ is the MLE of $p_E$ under $p_E / p_C = \delta$.
#' High values of $T_{\RD, \delta}$ favour the alternative hypothesis.
#' 
#' @param x_E Vector of number of events in experimental group.
#' @param x_C Vector of number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of x_E favour the alternative 
#' hypothesis and "low" vice versa.
#' @return Vector of values of the RD test statistic.
#' 
#' @export
#' 
#' @examples
#' test_RD(3, 4, 10, 10, 0.2, "high")
test_RR <- function(x_E, x_C, n_E, n_C, delta, better){
  p_E <- x_E / n_E
  p_C <- x_C / n_C
  
  theta <- n_C / n_E
  p1hat <- x_E / n_E
  p2hat <- x_C / n_C
  a <- 1 + theta
  b <- -(delta*(1 + theta*p2hat) + theta + p1hat)
  c <- delta * (p1hat + theta * p2hat)
  # Define the solution
  p_E0 <- (-b - sqrt(round(b^2 - 4*a*c,10)))/(2*a)
  p_C0 <- round(p_E0 / delta, 10)
  
  denom <- ifelse(p_E - delta * p_C == 0, 1, sqrt(round(p_E0*(1-p_E0)/n_E + p_C0*(1-p_C0)*delta^2/n_C, 10)))
  num <- p_E - delta * p_C
  if (better == "high"){
    return <- num/denom
  }
  if (better == "low"){
    return <- -num/denom
  }
    return(return)
}




