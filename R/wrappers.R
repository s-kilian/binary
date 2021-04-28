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
#' $$ T_{\OR, \delta} = 1-(1 - F_{\ncHg(X_E+X_C, n_E, n_C, \delta)}(x_E-1)) $$
#' is based on Fisher's non-central hypergeometric distribution with density
#' $$ f_{\ncHg(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, $$
#' where $A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}$.
#' The density is zero if $k < \max(0, s-n_C)$ or $k > \min(n_E, s)$.
#' High values of $T_{\OR, \delta}$ favour the alternative hypothesis (due to "1-...").
#' 
#' @param df data frame with variables x_E and x_C
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
        stat = test_RR(x_E, x_C, n_E, n_C, delta, better)
      ) 
  }
  
  if (method == "RD") {
    return = df %>%
      dplyr::mutate(
        stat = test_RD(x_E, x_C, n_E, n_C, delta, better)
      )
  }
  
  if (method == "OR") {
    return <- df %>%
      dplyr::mutate(s = x_C+x_E) %>%
      dplyr::group_by(s) %>%
      dplyr::mutate(
        stat = BiasedUrn::pFNCHypergeo(x_E-1, n_E, n_C, s[1], delta)
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

# function to compute d p_E/d p_C from p_C and NI-margin delta
d.p_E.p_C <- function(p_C, method, delta){
  
  if (method == "RR") {
    d.p_E.p_C <- rep(delta, length(p_C))
  }
  if (method == "RD") {
    d.p_E.p_C <- rep(1, length(p_C))
  }
  if (method == "OR") {
    d.p_E.p_C <- delta/(delta*p_C + 1 - p_C)^2
  }
  
  return(d.p_E.p_C)
}

# Calculate derivative after p_E of probability of a specific region under H0
prob.derivative <- function(
  p_C.vec,      # vector of true values of p_C
  x_E,          # vector of values belonging to region
  x_C,          # vector of values belonging to region
  n_E,
  n_C,
  method,
  delta
){
  if(!(length(x_E) == length(x_C))){
    stop("x_E and x_C must have same length.")
  }
  
  p_E.vec <- p_C.to.p_E(p_C = p_C.vec, method = method, delta = delta)
  d.p_E.p_C <- d.p_E.p_C(p_C = p_C.vec, method = method, delta = delta)
  
  result <- c()
  for (i in 1:length(p_C.vec)) {
    result <- c(result,
                sum(
                  dbinom(x_E, n_E, p_E.vec[i])*dbinom(x_C, n_C, p_C.vec[i])*
                    (x_E/p_E.vec[i]*d.p_E.p_C[i] - (n_E-x_E)/(1-p_E.vec[i])*d.p_E.p_C[i] + x_C/p_C.vec[i] - (n_C-x_C)/(1-p_C.vec[i]))
                )
    )
  }
  return(result)
}

find_max_prob_uniroot <- function(
  df,             # data frame with variables x_E, x_C and reject
  n_E,
  n_C,
  method,
  delta
){
  # # for testing
  # method = "RR"
  # n_E <- 100
  # n_C <- 100
  # delta <- 1.1
  # alpha <- 0.05
  # expand.grid(
  #   x_E = 0:n_E,
  #   x_C = 0:n_C
  # ) %>%
  #   teststat(
  #     n_E = n_E,
  #     n_C = n_C,
  #     delta = delta,
  #     method = method,
  #     better = "high"
  #   ) %>%
  #   mutate( reject = stat >= qnorm(1-alpha)) ->
  #   df
  df %>%
    filter(reject) ->
    df.
  
  if(method == "RD") interval.p_C <- c(max(0, -delta), min(1, 1-delta))
  if(method == "RR") interval.p_C <- c(0, min(1, 1/delta))
  if(method == "OR") interval.p_C <- c(0, 1)
  
  rootSolve::uniroot.all(
    f = prob.derivative,
    interval = interval.p_C,
    n = 100,
    maxiter = 10^4,
    x_E = df.$x_E,
    x_C = df.$x_C,
    n_E = n_E,
    n_C = n_C,
    method = method,
    delta = delta
  ) ->
    roots
  
  p_CA <- c(interval.p_C, roots)
  p_EA <- p_C.to.p_E(p_C = p_CA, method = method, delta = delta)
  
  return(
    power(
      df = df,
      n_C = n_C,
      n_E = n_C,
      p_CA = p_CA,
      p_EA = p_EA
    )
  )
}

power <- function(df, n_C, n_E, p_CA, p_EA){
  # Take data frame df with variable x_C and x_E representing all possible
  # response pairs for group sizes n_C and n_E and variable reject indicating
  # whether coordinates belong to rejection region.
  # Compute exact prob. of rejection region for all pairs (p_CA, p_EA).
  
  if (
    n_C+1 != df %>% dplyr::pull(x_C) %>% unique() %>% length() |
    n_E+1 != df %>% dplyr::pull(x_E) %>% unique() %>% length()
  ) {
    stop("Values of x_C and x_E have to fit n_C and n_E.")
  }
  
  if (
    length(p_CA) != length(p_EA) |
    !all(p_CA >= 0 & p_CA <= 1 & p_EA >= 0 & p_EA <= 1)
  ) {
    stop("p_CA and p_EA must have same length and values in [0, 1].")
  }
  
  
  # compute uncond. size for every p
  df %>%
    dplyr::filter(reject) ->
    df.reject
  
  sapply(
    1:length(p_CA),
    function(i) {
      sum(stats::dbinom(df.reject$x_C, n_C, p_CA[i])*stats::dbinom(df.reject$x_E, n_E, p_EA[i]))
    }
  ) ->
    result
  names(result) <- paste(p_CA, p_EA, sep = ", ")
  return(result)
}

find_max_prob <- function(
  df,             # data frame with variables x_E, x_C and reject
  n_E,
  n_C,
  method,
  delta,
  calc_method = c("uniroot", "grid search"),     # method to find maximum
  size_acc = 3
){
  if(calc_method == "uniroot"){
    result <- max(find_max_prob_uniroot(
      df = df,
      n_E = n_E,
      n_C = n_C,
      method = method,
      delta = delta
    ))
  }
  if(calc_method == "grid search"){
    # Define grid for p_C
    p_C <- seq(10^-size_acc, 1-10^-size_acc, by = 10^-size_acc)
    
    # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
    # the null hypothesis
    p_E <- p_C.to.p_E(p_C, method, delta)
    
    p_C <- p_C[p_E >= 0 & p_E <= 1]
    p_E <- p_E[p_E >= 0 & p_E <= 1]
    
    result <- max(power(
      df = df,
      n_E = n_E,
      n_C = n_C,
      p_EA = p_E,
      p_CA = p_C
    ))
  }
  return(result)
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
      reject = stat >= stat[x_E == x_E. & x_C == x_C.]
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

conf_int <- function(
  x_E.,
  x_C.,
  n_E,
  n_C,
  alpha = 0.05,
  method,
  delta = NULL,
  size_acc = 3,
  better = c("high", "low")
){
  # # for testing
  # x_E. = 10
  # x_C. = 20
  # n_E <- 100
  # n_C <- 100
  # alpha <- 0.05
  # method <- "RD"
  # delta <- 
  # work in progress
}

effect <- function(p_E, p_C, method){
  if(method == "RD") effect <- p_E-p_C
  if(method == "RR") effect <- p_E/p_C
  if(method == "OR") effect <- (p_E/(1-p_E))/(p_C/(1-p_C))
  return(effect)
}

samplesize_appr <- function(p_EA, p_CA, delta, alpha, beta, r, method, better){
  check.effect.delta.better(
    p_E = p_EA,
    p_C = p_CA,
    method = method,
    delta = delta,
    better = better
  )
  
  if(method == "RD"){
    # if low values of p_EA favor the alternative, flip probabilites and delta
    if(better == "low"){
      p_EA <- 1-p_EA
      p_CA <- 1-p_CA
      delta <- -delta
    }
    theta <- 1/r
    
    a <- 1 + theta
    b <- -(1 + theta + p_EA + theta * p_CA + delta * (theta + 2))
    c <- delta^2 + delta * (2*p_EA + theta + 1) + p_EA + theta * p_CA
    d <- -p_EA * delta * (1 + delta)
    
    # Define the parameters for solving the equation
    v <- b^3/(3*a)^3 - b*c/(6*a^2) + d/(2*a)
    u <- sign(v) * sqrt(b^2/(3*a)^2 - c/(3*a))
    w <- 1/3 * (pi + acos(v/u^3))
    
    # Define the solution
    p_E0 <- 2*u*cos(w) - b/(3*a)
    p_C0 <- p_E0 - delta
    
    z_alpha <- stats::qnorm(1 - alpha, mean = 0, sd = 1)
    z_beta <- stats::qnorm(1 - beta, mean = 0, sd = 1)
    
    # Calculating the sample size
    n <- (1+r) / r * (z_alpha * sqrt(r * p_C0 * (1 - p_C0) + p_E0 * (1-p_E0))
                      + z_beta * sqrt(r * p_CA * (1 - p_CA) + p_EA * (1-p_EA)))^2 / 
      (p_EA - p_CA - delta)^2
    
    n_C <- ceiling(1/(1+r) * n)
    n_E <- ceiling(r/(1+r) * n)
    n_total <- n_E + n_C
  }
  
  
  if(method == "RR"){
    # if low values of p_EA favor the alternative, swap groups and flip delta and r
    if(better == "low"){
      p_CA. <- p_EA
      p_EA <- p_CA
      p_CA <- p_CA.
      delta <- 1/delta
      r <- 1/r
    }
    theta <- 1/r
    
    a <- 1 + theta
    b <- -(delta*(1 + theta*p_CA) + theta + p_EA)
    c <- delta * (p_EA + theta * p_CA)
    
    # Define the solution
    p_E0 <- (-b - sqrt(round(b^2 - 4*a*c,10)))/(2*a)
    p_C0 <- p_E0 / delta
    
    z_alpha <- stats::qnorm(1 - alpha, mean = 0, sd = 1)
    z_beta <- stats::qnorm(1 - beta, mean = 0, sd = 1)
    
    # Calculating the sample size
    n <- (1+r) / r * (z_alpha * sqrt(r * delta^2 * p_C0 * (1 - p_C0) + p_E0 * (1 - p_E0))
                      + z_beta * sqrt(r * delta^2 * p_CA * (1 - p_CA) + p_EA * (1 - p_EA)))^2 / 
      (p_EA - delta * p_CA)^2
    
    n_C <- ceiling(1/(1+r) * n)
    n_E <- ceiling(r/(1+r) * n)
    if(better == "low"){
      n_E. <- n_C
      n_C <- n_E
      n_E <- n_E.
    }
    n_total<- n_E + n_C
    
  }
  
  if(method == "OR"){
    if(better == "low"){
      p_CA. <- p_EA
      p_EA <- p_CA
      p_CA <- p_CA.
      delta <- 1/delta
      r <- 1/r
    }
    samplesize_Wang(
      p_EA = p_EA,
      p_CA = p_CA,
      gamma = delta,
      alpha = alpha,
      beta = beta,
      r = r
    ) ->
      res
    
    n_C <- res$n_C
    n_E <- res$n_E
    
    if(better == "low"){
      n_E. <- n_C
      n_C <- n_E
      n_E <- n_E.
    }
    n_total<- n_E + n_C
  }
  
  return(data.frame(n_total = n_total, n_E = n_E, n_C = n_C))
}
