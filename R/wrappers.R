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
#' \loadmathjax
#' \code{teststat} takes a data frame with variables \code{x_E} and \code{x_C}
#' and adds the variable \code{stat} with the value of the test statistic
#' specified by \code{n_E}, \code{n_C}, \code{method}, \code{delta} and \code{better}.
#' 
#' If higher values of \mjseqn{x_E} favour the alternative hypothesis (\code{better = "high"}), we are interested
#' in testing the null hypothesis
#' \mjsdeqn{H_0: e(p_E, p_C) \le \delta ,}
#' where \mjseqn{e} is one of the effect measures risk difference (\code{method = "RD"}),
#' risk ratio (\code{method = "RR"}), or odds ratio (\code{method = "OR"}).
#' The test statistic for risk difference is
#' \mjsdeqn{T_{\mbox{RD}, \delta}(x_E, x_C) = \frac{\hat{p_E} - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
#' where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
#' \mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E - p_C = \delta}.
#' High values of \mjseqn{T_{\mbox{RD}, \delta}} favour the alternative hypothesis.
#' The test statistic for risk ratio is
#' \mjsdeqn{T_{\mbox{RR}, \delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
#' where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
#' \mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E / p_C = \delta}.
#' High values of \mjseqn{T_{\mbox{RR}, \delta}} favour the alternative hypothesis.
#' The test statistic for Odds Ratio
#' \mjsdeqn{ T_{\mbox{OR}, \delta} = 1-(1 - F_{\mbox{ncHg}(X_E+X_C, n_E, n_C, \delta)}(x_E-1)) }
#' is based on Fisher's non-central hypergeometric distribution with density
#' \mjsdeqn{ f_{\mbox{ncHg}(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, }
#' where \mjseqn{A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}}.
#' The density is zero if \mjseqn{k < \max(0, s-n_C)} or \mjseqn{k > \min(n_E, s)}.
#' High values of \mjseqn{T_{\mbox{OR}, \delta}} favour the alternative hypothesis (due to "1-...").
#' 
#' @param df data frame with variables \mjseqn{x_E} and \mjseqn{x_C}
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param method Specifies the effect measure/test statistic. One of "RD", "RR", or "OR".
#' @param better "high" if higher values of x_E favour the alternative 
#' hypothesis and "low" vice versa.
#' @return
#' A list with the two elements \code{p_max} and \code{p_vec}.
#' \code{p_max} is the maximum p-value and most likely serves as "the one" p-value.
#' \code{p_vec} is a named vector. The names indicate the true proportion pairs
#' \mjseqn{(p_E, p_C)} with \mjseqn{e(p_E, p_C) = \delta} that underly the calculation of
#' the p-values. It can be used for plotting the p-value versus the true proportions.
#' 
#' @export
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
    if(better == "high"){
      return <- df %>%
        dplyr::mutate(s = x_C+x_E) %>%
        dplyr::group_by(s) %>%
        dplyr::mutate(
          stat = BiasedUrn::pFNCHypergeo(x_E-1, n_E, n_C, s[1], delta)
        ) %>%
        dplyr::ungroup()
    }
    if(better == "low"){
      return <- df %>%
        dplyr::mutate(s = x_C+x_E) %>%
        dplyr::group_by(s) %>%
        dplyr::mutate(
          stat = 1 - BiasedUrn::pFNCHypergeo(x_E-1, n_E, n_C, s[1], delta)
        ) %>%
        dplyr::ungroup()
    }
  }
  
  return(return)
}

# Calculate approximate quantile function of test statistic under H0
appr_teststat_quantile <- function(
  p,
  method,
  better
){
  check.0.1(
    values = p,
    message = "p has to be in interval (0, 1)."
  )
  
  if(method == "RD") result <- stats::qnorm(p)
  if(method == "RR") result <- stats::qnorm(p)
  if(method == "OR") result <- p
  
  return(result)
}

# function to create grid for p_C
p_C.grid <- function(
  method,
  delta,
  acc
){
  if(method == "RD") grid <- seq(max(0, -delta), min(1, 1-delta), length.out = 10^acc+1)
  if(method == "RR") grid <- seq(0, min(1, 1/delta), length.out = 10^acc+1)
  if(method == "OR") grid <- seq(0, 1, length.out = 10^acc+1)
  
  return(grid)
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

#' Calculate derivative of rejection probability
#' 
#' Calculate derivative after \mjseqn{p_E} of probability of a specific region under \mjseqn{H_0}.
#' 
#' 
#' @param p_C.vec data frame with variables \mjseqn{x_E} and \mjseqn{x_C}
#' @param x_E vector of values belonging to region
#' @param x_C vector of values belonging to region
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param method Specifies the effect measure/test statistic. One of "RD", "RR", or "OR".
#'
#' @return
#' Vector of values of the derivative
#' 
#' @export
#' 
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
    dplyr::filter(reject) ->
    df.
  
  # Interval of possible values for p_C
  interval.p_C <- p_C.grid(
    method = method,
    delta = delta,
    acc = 0
  )
  
  # find roots of derivative of exact probability of rejection region
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
      n_E = n_E,
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
    result <- find_max_prob_uniroot(
      df = df,
      n_E = n_E,
      n_C = n_C,
      method = method,
      delta = delta
    )
  }
  if(calc_method == "grid search"){
    # Define grid for p_C
    p_C <- p_C.grid(method = method, delta = delta, acc = size_acc)
    
    # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
    # the null hypothesis
    p_E <- p_C.to.p_E(p_C, method, delta)
    
    result <- power(
      df = df,
      n_E = n_E,
      n_C = n_C,
      p_EA = p_E,
      p_CA = p_C
    )
  }
  return(
    list(
      p_max = max(result),
      p_vec = result
    )
  )
}

#' Calculate p-value(s)
#' 
#' \code{p_value} returns the vector of p-values (dependent on the true \mjseqn{H_0} proportions)
#' and its maximum of the test for non-inferiority of two proportions specified by \code{method}.
#' 
#' If higher values of \mjseqn{x_E} favour the alternative hypothesis (\code{better = "high"}), we are interested
#' in testing the null hypothesis
#' \mjsdeqn{H_0: e(p_E, p_C) \le \delta ,}
#' where \mjseqn{e} is one of the effect measures risk difference (\code{method = "RD"}),
#' risk ratio (\code{method = "RR"}), or odds ratio (\code{method = "OR"}).
#' The test statistic for risk difference is
#' \mjsdeqn{T_{\mbox{RD}, \delta}(x_E, x_C) = \frac{\hat{p_E} - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
#' where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
#' \mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E - p_C = \delta}.
#' High values of \mjseqn{T_{\mbox{RD}, \delta}} favour the alternative hypothesis.
#' The test statistic for risk ratio is
#' \mjsdeqn{T_{\mbox{RR}, \delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
#' where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
#' \mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E / p_C = \delta}.
#' High values of \mjseqn{T_{\mbox{RR}, \delta}} favour the alternative hypothesis.
#' The test statistic for Odds Ratio
#' \mjsdeqn{ T_{\mbox{OR}, \delta} = 1-(1 - F_{\mbox{ncHg}(X_E+X_C, n_E, n_C, \delta)}(x_E-1)) }
#' is based on Fisher's non-central hypergeometric distribution with density
#' \mjsdeqn{ f_{\mbox{ncHg}(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, }
#' where \mjseqn{A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}}.
#' The density is zero if \mjseqn{k < \max(0, s-n_C)} or \mjseqn{k > \min(n_E, s)}.
#' High values of \mjseqn{T_{\mbox{OR}, \delta}} favour the alternative hypothesis (due to "1-...").
#' 
#'  
#' @param x_E. Number of events in experimental group.
#' @param x_C. Number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of x_E favour the alternative 
#' hypothesis and "low" vice versa.
#' @param method Specifies the effect measure/test statistic. One of "RD", "RR", or "OR".
#' @param size_acc Accuracy of grid
#' @param calc_method "grid search" or "uniroot"
#' 
#' 
#' @return
#' A list with the two elements \code{p_max} and \code{p_vec}.
#' \code{p_max} is the maximum p-value and most likely servers as "the one" p-value.
#' \code{p_vec} is a named vector. The names indicate the true proportion pairs
#' \mjseqn{(p_E, p_C)} with \mjseqn{e(p_E, p_C) = \delta} that underly the calculation of
#' the p-values. It can be used for plotting the p-value versus the true proportions.
#' 
#' @export
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
  better = c("high", "low"),
  calc_method = c("uniroot", "grid search")
){
  calc_method <- match.arg(calc_method)
  
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
  
  df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
  
  df %>%
    teststat(
      n_E = n_E,
      n_C = n_C,
      delta = delta,
      method = method,
      better = better
    ) %>%
    dplyr::mutate(
      reject = stat >= stat[x_E == x_E. & x_C == x_C.]
    ) %>%
    find_max_prob(
      df = .,
      n_E = n_E,
      n_C = n_C,
      method = method,
      delta = delta,
      calc_method = calc_method,
      size_acc = size_acc
    ) ->
    result
  
  return(result)
}

# calculate intervals where bounds of approximate confidence interval lie inside
appr_confint_bound_intervals <- function(
  x_E.,
  x_C.,
  n_E,
  n_C,
  alpha = 0.05,
  method
){
  if(method == "RD"){
    interval_lb <- c(
      max(
        x_E./n_E - x_C./n_C - appr_teststat_quantile(
          p = 1-alpha,
          method = "RD",
          better = "high"
        ) *
          0.5*sqrt(1/n_E + 1/n_C),
        -1
      ),
      min(
        x_E./n_E - x_C./n_C,
        1
      )
    )
    interval_ub <- c(
      max(
        -1,
        x_E./n_E - x_C./n_C
      ),
      min(
        1,
        x_E./n_E - x_C./n_C + appr_teststat_quantile(
          p = 1-alpha,
          method = "RD",
          better = "high"
        ) *
          0.5*sqrt(1/n_E + 1/n_C)
      )
    )
  } 
  if(method == "RR"){
    # has to be refined
    interval_lb <- c(0, 10^1)
    interval_ub <- c(0, 10^1)
    # interval_lb <- c(
    #   max(
    #     0,
    #     x_E./n_E / (x_C./n_C) - appr_teststat_quantile(
    #       p = 1-alpha,
    #       method = "RR",
    #       better = "high"
    #     ) *
    #       0.5*sqrt(1/n_E + 1/n_C) / (x_C./n_C)
    #   ),
    #   x_E./n_E / (x_C./n_C)
    # )
    # interval_ub <- c(
    #   max(
    #     0,
    #     x_E./n_E / (x_C./n_C)
    #   ),
    #   x_E./n_E / (x_C./n_C) + appr_teststat_quantile(
    #     p = 1-alpha,
    #     method = "RR",
    #     better = "high"
    #   ) *
    #     0.5*sqrt(1/n_E + 1/n_C) / (x_C./n_C)
    # )
  } 
  if(method == "OR"){
    # has to be refined
    interval_lb <- c(0, 10^1)
    interval_ub <- c(0, 10^1)
  } 
  
  return(
    list(
      interval_lb = interval_lb,
      interval_ub = interval_ub
    )
  )
}

# approximate confidence interval
conf_region_appr <- function(
  x_E.,
  x_C.,
  n_E,
  n_C,
  alpha = 0.05,
  method,
  size_acc = 3,
  acc = 3
){
  # define function to find root
  f <- function(delta.vec, better){
    result <- c()
    for(delta in delta.vec){
      data.frame(
        x_E = x_E.,
        x_C = x_C.
      ) %>%
        teststat(
          n_E = n_E,
          n_C = n_C,
          delta = delta,
          method = method,
          better = better
        ) %>%
        `[[`("stat") ->
        u
      result <- c(
        result,
        u - appr_teststat_quantile(
          p = 1-alpha,
          method = method,
          better = better
        )
      )
    }
    
      return(result)
  }
  
  bound_intervals <- appr_confint_bound_intervals(
    x_E = x_E.,
    x_C = x_C.,
    n_E = n_E,
    n_C = n_C,
    alpha = alpha,
    method = method
  )
  
  # lower and upper bound of approximate conf. int.
  lb <- rootSolve::uniroot.all(
    f = f,
    interval = bound_intervals$interval_lb,
    better = "high"
  ) %>% min()
  ub <- rootSolve::uniroot.all(
    f = f,
    interval = bound_intervals$interval_ub,
    better = "low"
  ) %>% max()
  
  return(
    c(lb, ub)
  )
}

conf_region <- function(
  x_E.,
  x_C.,
  n_E,
  n_C,
  alpha = 0.05,
  method,
  size_acc = 3,
  delta_acc = 2
){
  # # for testing
  # x_E. = 10
  # x_C. = 20
  # n_E <- 100
  # n_C <- 100
  # alpha <- 0.05
  # method <- "RD"
  # delta_acc <- 2

  cr_appr_outer <- conf_region_appr(
    x_E. = x_E.,
    x_C. = x_C.,
    n_E = n_E,
    n_C = n_C,
    alpha = alpha/2,
    method = method
  )
  cr_appr_inner <- conf_region_appr(
    x_E. = x_E.,
    x_C. = x_C.,
    n_E = n_E,
    n_C = n_C,
    alpha = alpha*2,
    method = method
  )
  delta_vec_low <- c(
    seq(cr_appr_outer[1], cr_appr_inner[1], by = 10^-delta_acc)
  )
  delta_vec_high <- c(
    seq(cr_appr_inner[2], cr_appr_outer[2], by = 10^-delta_acc)
  )
    
  delta_in_region <- c()
  suppressWarnings(
    for (i in 1:length(delta_vec_low)) {
      p_val <- p_value(
        x_E. = x_E.,
        x_C. = x_C.,
        n_E = n_E,
        n_C = n_C,
        better = "high",
        method = method,
        delta = delta_vec_low[i],
        calc_method = "uniroot",
        size_acc = size_acc
      )
      delta_in_region <- c(delta_in_region, p_val$p_max > alpha)
    }
  )
  suppressWarnings(
    for (i in 1:length(delta_vec_high)) {
      p_val <- p_value(
        x_E. = x_E.,
        x_C. = x_C.,
        n_E = n_E,
        n_C = n_C,
        better = "low",
        method = method,
        delta = delta_vec_high[i],
        calc_method = "uniroot",
        size_acc = size_acc
      )
      delta_in_region <- c(delta_in_region, p_val$p_max > alpha)
    }
  )
  
  delta_vec <- c(delta_vec_low, delta_vec_high)
  
  (1:length(delta_in_region))[delta_in_region] %>% 
    diff() %>%
    max() ->
    max.diff
  
  connected <- max.diff <= 1
  
  largest.hole <- NA
  if(!connected) largest.hole <- max.diff*10^-delta_acc
  
  return(
    list(
      delta_vec = delta_vec,
      delta_in_region = delta_in_region,
      lb = min(delta_vec[delta_in_region]),
      ub = max(delta_vec[delta_in_region]),
      connected = connected,
      largest.hole = largest.hole,
      ci_appr = conf_region_appr(
        x_E. = x_E.,
        x_C. = x_C.,
        n_E = n_E,
        n_C = n_C,
        alpha = alpha,
        method = method
      )
    )
  )
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
  
  # Check whether alpha lies in the interval (0, 0.5]
  check.alpha(alpha = alpha)
  
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
    if(v != 0) u <- sign(v) * sqrt(b^2/(3*a)^2 - c/(3*a)) else u <- 1
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

# function to compute the critical value of a specific test statistic
critval <- function(alpha, n_C, n_E, method, delta, size_acc = 4, better, start_value = NULL){
  # method defines the used test with corresponding statistic and, size_acc defines
  # the accuracy of the grid used for the nuisance parameter p_C
  
  # Create data frame of all test statistics ordered by test statistic
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat(n_E, n_C, delta, method, better) %>%
    dplyr::arrange(desc(stat)) ->
    df.stat
  
  # Extract stat, x_C and x_E as vector
  stat <- df.stat$stat
  x_C <- df.stat$x_C
  x_E <- df.stat$x_E
  
  # Find starting value for the search of critical value. Take the
  # quantile of the approximate distribution of stat of no value is provided.
  if(
    is.null(start_value)
    ){
    start_value <- appr_teststat_quantile(
      p = 1-alpha,
      method = method,
      better = better
    ) 
  }
  
  
  # Find row number of df.stat corresponding to starting value
  # <- row of df.stat where stat is maximal with stat <= start_value
  # Special case with very small sample sizes can lead to stat > start_value 
  # for all rows. Then set i <- 1
  i <- max(sum(stat<start_value), 1)
  i.max <- length(stat)
  
  # Define rough grid for p_C
  acc <- 1
  p_C <- p_C.grid(method = method, delta = delta, acc = acc)
  
  # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
  # the null hypothesis
  p_E <- p_C.to.p_E(p_C, method, delta)
  
  
  # Calculate exact size for pair (p_C[i], p_E[i])
  pr <- sapply(1:length(p_C), function(j)
    sum(stats::dbinom(x_C[1:i], n_C, p_C[j]) * stats::dbinom(x_E[1:i], n_E, p_E[j])) )
  
  # Increase index if maximal size is too low
  while (max(pr) <= alpha & i < i.max) {
    i <- i+1
    pr <- pr + stats::dbinom(x_C[i], n_C, p_C) * stats::dbinom(x_E[i], n_E, p_E)
  }
  
  
  # Decrease index if maximal size is too high and iterate grid accuracy
  for (acc in 1:size_acc) {
    
    # Define grid for p_C
    p_C <- p_C.grid(method = method, delta = delta, acc = acc)
    
    # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
    # the null hypothesis
    p_E <- p_C.to.p_E(p_C, method, delta)
    
    pr <- sapply(1:length(p_C), function(j)
      sum(stats::dbinom(x_C[1:i], n_C, p_C[j]) * stats::dbinom(x_E[1:i], n_E, p_E[j])) )
    
    # Decrease index if maximal size is too high
    while (max(pr) > alpha & i >= 1) {
      # Compute new sizes
      pr <- pr - stats::dbinom(x_C[i], n_C, p_C) * stats::dbinom(x_E[i], n_E, p_E)
      i <- i-1
    }
    
  }
  if(i == 0){
    result <- list(
      crit.val.lb = Inf,
      crit.val.mid = Inf,
      crit.val.ub = stat[1],
      max.size = 0
    )
  } else {
    # Decrease index further as long as rows have the same test statistic value
    while (stat[i+1] == stat[i] & i >= 1) {
      pr <- pr - stats::dbinom(x_C[i], n_C, p_C) * stats::dbinom(x_E[i], n_E, p_E)
      i <- i-1
    }
    
    # Critical value can now be chosen between stat[i+1] and stat[i]
    crit.val.mid <- (stat[i+1] + stat[i])/2
    
    result <- list(
      crit.val.lb = stat[i],
      crit.val.mid = crit.val.mid,
      crit.val.ub = stat[i+1],
      max.size = max(pr)
    )
  }
  
  
  return(result)
}

# Function to compute exact sample size
samplesize_exact <- function(p_EA, p_CA, delta, alpha, beta, r, size_acc = 3, method, better){
  # Calculate exact sample size for method "X and specified
  # level alpha, beta, allocation ratio r = n_E/n_C, true rates p_CA, p_EA and
  # NI-margin delta
  # Accuracy of calculating the critical value can be specified by size_acc.
  # Output: Sample sizes per group (n_C, n_E), critical value and exact power.
  
  # Check whether p_EA, p_CA lie in the alternative hypothesis
  check.effect.delta.better(
    p_E = p_EA,
    p_C = p_CA,
    method = method,
    delta = delta,
    better = better
  )
  
  # Check whether alpha lies in the interval (0, 0.5]
  check.alpha(alpha = alpha)
  
  # Estimate sample size with approximate formula
  n_appr <- samplesize_appr(
    p_EA = p_EA,
    p_CA = p_CA,
    delta = delta,
    alpha = alpha,
    beta = beta,
    r = r,
    method = method,
    better = better
  )
  
  # Use estimates as starting values
  n_C <- n_appr[["n_C"]]
  n_E <- n_appr[["n_E"]]
  
  # Initiate data frame
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat(n_C = n_C, n_E = n_E, delta = delta, method = method, better) ->
    df
  
  # Calculate critical value
  crit.val <- critval(alpha = alpha, n_C = n_C, n_E = n_E, delta = delta, size_acc = size_acc, method = method, better = better)["crit.val.mid"]
  
  # Calculate exact power
  df %>%
    dplyr::mutate(reject = stat >= crit.val) %>%
    power(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
    exact_power
  
  # Decrease sample size if power is too high
  if(exact_power > 1-beta){
    while(exact_power > 1-beta){
      # Store power and nominal level of last iteration
      last_power <- exact_power
      last_crit.val <- crit.val
      
      # Decrease sample size by minimal amount possible with allocation ratio r
      if (r >= 1) {
        n_C <- n_C - 1
        n_E <- ceiling(r*n_C)
      } else {
        n_E <- n_E - 1
        n_C <- ceiling(1/r*n_E)
      }
      
      # Initiate data frame
      expand.grid(
        x_C = 0:n_C,
        x_E = 0:n_E
      ) %>%
        teststat(n_C = n_C, n_E = n_E, delta = delta, method = method, better = better) ->
        df
      
      # Calculate raised nominal level
      crit.val <- critval(alpha = alpha, n_C = n_C, n_E = n_E, delta = delta, size_acc = size_acc, method = method, better = better, start_value = crit.val)["crit.val.mid"]
      
      # Calculate exact power
      df %>%
        dplyr::mutate(reject = stat >= crit.val) %>%
        power(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
        exact_power
    }
    # Go one step back
    if (r >= 1) {
      n_C <- n_C + 1
      n_E <- ceiling(r*n_C)
    } else {
      n_E <- n_E + 1
      n_C <- ceiling(1/r*n_E)
    }
    exact_power <- last_power
    crit.val <- last_crit.val
  }
  
  # If power is too low: increase sample size until power is achieved
  while (exact_power < 1-beta) {
    if (r >= 1) {
      n_C <- n_C + 1
      n_E <- ceiling(r*n_C)
    } else {
      n_E <- n_E + 1
      n_C <- ceiling(1/r*n_E)
    }
    
    # Initiate data frame
    expand.grid(
      x_C = 0:n_C,
      x_E = 0:n_E
    ) %>%
      teststat(n_C = n_C, n_E = n_E, delta = delta, method = method, better = better) ->
      df
    
    # Calculate raised nominal level
    crit.val <- critval(alpha = alpha, n_C = n_C, n_E = n_E, delta = delta, size_acc = size_acc, method = method, better = better, start_value = crit.val)["crit.val.mid"]
    
    # Calculate exact power
    df %>%
      dplyr::mutate(reject = stat >= crit.val) %>%
      power(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
      exact_power
  }
  
  return(
    list(
      n_C = n_C,
      n_E = n_E,
      crit.val = crit.val,
      exact_power = exact_power
    )
  )
}
