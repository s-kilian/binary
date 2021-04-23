##..............................................................................
##               Project: Binary exact NI tests
##               Purpose: Provide wrapper functions
##                 Input: None
##                Output: Wrapper functions
##      Date of creation: 2021-04-23
##   Date of last update: 
##                Author: Samuel Kilian
##..............................................................................

#' Calculate test statistic
#' 
#' \code{teststat} takes a data frame with variables \code{x_E} and \code{x_C}
#' and adds the variable \code{stat} with the value of the test statistic
#' specified by \code{n_E}, \code{n_C}, \code{method}, \code{delta} and \code{better}.
#' 
#' If higher values of $x_E$ favour the alternative hypothesis (\code{better = "high"}), we are interested
#' in testing the null hypothesis
#' $$H_0: e(p_E, p_C) \le \delta ,$$
#' where $e$ is one of the effect measures risk difference (\code{method = "RD"}),
#' risk ratio (\code{method = "RR"}), or odds ratio (\code{method = "OR"}).
#' The test statistic for risk difference is
#' $$T_{\RD, \delta}(x_E, x_C) = \frac{\hat p_E - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},$$
#' where $\tilde p_C = \tilde p_C(x_E, x_C)$ is the MLE of $p_C$ and
#' $\tilde p_E = \tilde p_C + \delta$ is the MLE of $p_E$ under $p_E - p_C = \delta$.
#' High values of $T_{\RD, \delta}$ favour the alternative hypothesis.
#' The test statistic for risk ratio is
#' $$T_{\RR, \delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},$$
#' where $\tilde p_C = \tilde p_C(x_E, x_C)$ is the MLE of $p_C$ and
#' $\tilde p_E = \tilde p_C + \delta$ is the MLE of $p_E$ under $p_E / p_C = \delta$.
#' High values of $T_{\RR, \delta}$ favour the alternative hypothesis.
#' The test statistic for Odds Ratio
#' $$ T_{\OR, \delta} = \sum\limits_{k = X_E}^\infty  f_{\ncHg(X_E+X_C, n_E, n_C, \delta)}(k) $$
#' is based on Fisher's non-central hypergeometric distribution with density
#' $$ f_{\ncHg(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, $$
#' where $A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}$.
#' The density is zero if $k < \max(0, s-n_C)$ or $k > \min(n_E, s)$.
#' Small values of $T_{\OR, \delta}$ favour the alternative hypothesis.
#' 
#'  
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param method Specifies the effect measure/test statistic. One of "RD", "RR", or "OR".
#' @param better "high" if higher values of x_E favour the alternative 
#' hypothesis and "low" vice versa.
#' @return
#' A list with the two elements \code{p_max} and \code{p_vec}.
#' \code{p_max} is the maximum p-value and most likely servers as "the one" p-value.
#' \code{p_vec} is a named vector. The names indicate the true proportion pairs
#' $(p_E, p_C)$ with $e(p_E, p_C) = \delta$ that underly the calculation of
#' the p-values. It can be used for plotting the p-value versus the true proportions.
#' 
#' @examples
#' n_E <- 10
#' n_C <- 11
#' df <- expand.grid(
#'   x_E = 0:n_E,
#'   x_C = 0:n_C
#' )
#' teststat(
#'   df = df,
#'   n_E = n_E,
#'   n_C = n_C,
#'   method = "RD",
#'   delta = -0.1,
#'   better = "high"
#' )
teststat <- function(df, n_E, n_C, delta, method, better){
  if (method == "RR") {
    return <- df %>%
      dplyr::mutate(
        stat = -test_RR(x_E, x_C, n_E, n_C, delta, better)
      ) 
  }
  
  if (method == "RD") {
    return = df %>%
      dplyr::mutate(
        stat = -test_RD(x_E, x_C, n_E, n_C, -delta, better)
      )
  }
  
  if (method == "OR") {
    return <- df %>%
      dplyr::mutate(s = x_C+x_E) %>%
      dplyr::group_by(s) %>%
      dplyr::mutate(
        stat = BiasedUrn::pFNCHypergeo(x_C, n_C, n_E, s[1], 1/delta)
      ) %>%
      ungroup()
  }
  
  return(return)
}

# function to compute p_E from p_C and NI-margin delta s.t. effect(p_E, p_C) = delta
p_C.to.p_E <- function(p_C, method, delta){
  
  if (method == "RR") {
    p_E <- p_C * delta
  }
  if (method == "RD") {
    p_E <- p_C + delta
  }
  if (method == "OR") {
    p_E <- 1/(1+(1-p_C)/(delta*p_C))
  }
  
  return(p_E)
}

#' Calculate p-value(s)
#' 
#' \code{p_value} returns the vector of p-values (dependent on the true $H_0$ proportions)
#' and its maximum of the test for non-inferiority of two proportions specified by \code{method}.
#' 
#' If higher values of $x_E$ favour the alternative hypothesis (\code{better = "high"}), we are interested
#' in testing the null hypothesis
#' $$H_0: e(p_E, p_C) \le \delta ,$$
#' where $e$ is one of the effect measures risk difference (\code{method = "RD"}),
#' risk ratio (\code{method = "RR"}), or odds ratio (\code{method = "OR"}).
#' The test statistic for risk difference is
#' $$T_{\RD, \delta}(x_E, x_C) = \frac{\hat p_E - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},$$
#' where $\tilde p_C = \tilde p_C(x_E, x_C)$ is the MLE of $p_C$ and
#' $\tilde p_E = \tilde p_C + \delta$ is the MLE of $p_E$ under $p_E - p_C = \delta$.
#' High values of $T_{\RD, \delta}$ favour the alternative hypothesis.
#' The test statistic for risk ratio is
#' $$T_{\RR, \delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},$$
#' where $\tilde p_C = \tilde p_C(x_E, x_C)$ is the MLE of $p_C$ and
#' $\tilde p_E = \tilde p_C + \delta$ is the MLE of $p_E$ under $p_E / p_C = \delta$.
#' High values of $T_{\RR, \delta}$ favour the alternative hypothesis.
#' The test statistic for Odds Ratio
#' $$ T_{\OR, \delta} = \sum\limits_{k = X_E}^\infty  f_{\ncHg(X_E+X_C, n_E, n_C, \delta)}(k) $$
#' is based on Fisher's non-central hypergeometric distribution with density
#' $$ f_{\ncHg(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, $$
#' where $A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}$.
#' The density is zero if $k < \max(0, s-n_C)$ or $k > \min(n_E, s)$.
#' Small values of $T_{\OR, \delta}$ favour the alternative hypothesis.
#' 
#'  
#' @param x_E Number of events in experimental group.
#' @param x_C Number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of x_E favour the alternative 
#' hypothesis and "low" vice versa.
#' @return
#' A list with the two elements \code{p_max} and \code{p_vec}.
#' \code{p_max} is the maximum p-value and most likely servers as "the one" p-value.
#' \code{p_vec} is a named vector. The names indicate the true proportion pairs
#' $(p_E, p_C)$ with $e(p_E, p_C) = \delta$ that underly the calculation of
#' the p-values. It can be used for plotting the p-value versus the true proportions.
#' 
#' @examples
#' p_value(
#'   x_E. = 3,
#'   x_C. = 4,
#'   n_E = 10,
#'   n_C = 10,
#'   method = "RD",
#'   delta = -0.1,
#'   size_acc = 3,
#'   better = "high"
#' )
p_value <- function(
  x_E.,
  x_C.,
  n_E,
  n_C,
  method,
  delta = NULL,
  size_acc = 3,
  better = c("high", "low")
){
  # Check if input is correctly specified
  check.pos.int(
    size_acc,
    "size_acc has to be positive integer."
  )
  check.pos.int(
    c(x_E.+1, x_C.+1, n_E, n_C, n_E-x_E.+1, n_C-x_C.+1),
    "n_E, n_C have to be positive integers and x_E. in {0, ..., n_E}, x_C. in {0, ..., n_C}."
  )
  
  if (
    any(
      sapply(
        list(x_E., x_C., n_E, n_C, gamma, size_acc, better),
        length
      ) != 1
    )
  ) {
    stop("Input values have to be single values.")
  }
  
  # Set delta to default if not specified and throw warning
  delta <- check.delta.null(
    method = method,
    delta = delta
  )
  
  # Check if delta is in correct range
  check.delta.method.better(
    delta = delta, 
    method = method,
    better = better
  )
  
  # Define grid for p_C
  p_C <- seq(10^-size_acc, 1-10^-size_acc, by = 10^-size_acc)
  
  # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
  # the null hypothesis
  p_E <- p_C.to.p_E(p_C, method, delta)
  
  p_C <- p_C[p_E >= 0 & p_E <= 1]
  p_E <- p_E[p_E >= 0 & p_E <= 1]
  
  df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
  
  df %>%
    teststat(n_E, n_C, delta, method, better) %>%
    dplyr::mutate(
      reject = stat <= stat[x_E == x_E. & x_C == x_C.]
    ) %>%
    power(n_C, n_E, p_C, p_E) ->      # function power actually computes the rejection probability which in this case is the p-value
    p.values
  
  # return vector of p-values. The "one" p-value would be max(p.values).
  return(
    list(
      p_max = max(p.values),
      p_vec = p.values
    )
  )
}