##..............................................................................
##               Project: Package binary (working title)
##               Purpose: Provide test statistics for default tests
##      Date of creation: 2019-07-03
##                Author: Samuel Kilian
##..............................................................................

# Define class exbin_ts ####
# Constructor
new_exbin_ts <- function(
  ts_fun = function(x_E, x_C){},
  ts_fun_args = list(),
  quant_fun = function(p){},
  quant_fun_args = list()
){
  structure(
    list(
      ts_fun = ts_fun,
      ts_fun_args = ts_fun_args,
      quant_fun = quant_fun,
      quant_fun_args = quant_fun_args
    ),
    class = "exbin_ts"
  )
}

# Validator
validate_exbin_ts <- function(ts){
  # Check of basic properties
  if(!is.function(ts$ts_fun) | !is.function(ts$quant_fun)) stop("ts_fun and quant_fun must be functions.")
  if(!is.list(ts$ts_fun_args) | !is.list(ts$quant_fun_args)) stop("ts_fun_args and quant_fun_args must be lists.")
  if(!all(c("x_E", "x_C") %in% formalArgs(ts$ts_fun))) stop("ts_fun must have arguments x_E and x_C.")
  if(!(c("p") %in% formalArgs(ts$quant_fun))) stop("quant_fun must have argument p.")
  
  # Minimal example for ts_fun
  calc_ts.exbin_ts(
    ts = ts,
    x_E = c(0, 0, 1, 1),
    x_C = c(0, 1, 0, 1),
    n_E = 2,
    n_C = 2,
    delta = 0,
    better = "high"
  ) ->
    ts_fun_res
  if(!is.numeric(ts_fun_res) | length(ts_fun_res) != 4) stop("ts_fun didn't produce a numeric vector of the correct length when applied to a minimal example with n_E = 2, n_C = 2, delta = 0, better = \"high\".")
  
  # Minimal example for quant_fun
  calc_quant.exbin_ts(
    ts = ts,
    p = 0.05,
    n_E = 2,
    n_C = 2,
    delta = 0,
    better = "high"
  ) ->
    quant_fun_res
  if(!(is.null(quant_fun_res) | (is.numeric(quant_fun_res) & length(quant_fun_res) == 1))) stop("quant_fun didn't produce a numeric vector of the correct length (or NULL) when applied to a minimal example with n_E = 2, n_C = 2, delta = 0, better = \"high\".")
  
  ts
}

# Helper
exbin_ts <- function(ts_fun, ts_fun_args = list(), quant_fun = function(p){}, quant_fun_args = list()){
    ts <- new_exbin_ts(
    ts_fun = ts_fun,
    ts_fun_args = ts_fun_args,
    quant_fun = quant_fun,
    quant_fun_args = quant_fun_args
  )
  
  validate_exbin_ts(ts)
}

# Methods
calc_ts.exbin_ts <- function(ts, x_E, x_C, ...){
  # Check if x_E and x_C are vectors of non-negative integers of same length
  check.pos.int(
    values = c(x_E, x_C) + 1,
    message = "x_E and x_C have to be vectors of non-negative integers."
  )
  if(length(x_E) != length(x_C)) stop("x_E and x_C must have equal length.")
  
  arguments <- list(...)
  do.call(
    what = ts$ts_fun,
    args = c(
      list(
        x_E = x_E,
        x_C = x_C
      ),
      arguments[intersect(names(arguments), formalArgs(ts$ts_fun))],   # use given arguments that are taken by ts_fun
      ts$ts_fun_args[names(ts$ts_fun_args) %in% formalArgs(ts$ts_fun) & !(names(ts$ts_fun_args) %in% names(arguments))]
    )
  )
}
calc_quant.exbin_ts <- function(ts, p, ...){
  # Check if p is a vector of values between 0 and 1
  check.0.1(
    values = p,
    message = "p has to be a vector of values between 0 and 1."
  )
  
  arguments <- list(...)
  do.call(
    what = ts$quant_fun,
    args = c(
      list(
        p = p
      ),
      arguments[intersect(names(arguments), formalArgs(ts$quant_fun))],
      ts$quant_fun_args[names(ts$quant_fun_args) %in% formalArgs(ts$quant_fun) & !(names(ts$quant_fun_args) %in% names(arguments))]
    )
  )
}

update.exbin_ts <- function(ts, ...){
  arguments <- list(...)
  ts_fun_args.to.update <- intersect(formalArgs(ts$ts_fun), names(arguments))
  ts$ts_fun_args[ts_fun_args.to.update] <- arguments[ts_fun_args.to.update]
  quant_fun_args.to.update <- intersect(formalArgs(ts$quant_fun), names(arguments))
  ts$quant_fun_args <- arguments[quant_fun_args.to.update]
  ts
}

#' Calculate Farrington-Manning RD test statistic
#' 
#' \code{test_stat_FM_RD} returns the value of the Farrington-Manning test statistic
#' for non-inferiority of the risk difference between two proportions.
#' 
#' If higher values of \mjseqn{x_E} favor the alternative hypothesis, we are interested
#' in testing the null hypothesis
#' \mjsdeqn{H_0: p_E - p_C \le \delta ,}
#' where the NI-margin is usually non-positive: \mjseqn{\delta \le 0}.
#' The test statistic for this hypothesis is
#' \mjsdeqn{T_{\delta}(x_E, x_C) = \frac{\hat p_E - \hat p_C - \delta}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
#' where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
#' \mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E} under \mjseqn{p_E - p_C = \delta}.
#' High values of \mjseqn{T_{\delta}} favor the alternative hypothesis.
#' 
#' @param x_E Vector of number of events in experimental group.
#' @param x_C Vector of number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of \mjseqn{x_E} favor the alternative 
#' hypothesis and "low" vice versa.
#' @return Vector of values of the RD test statistic.
#' 
#' @export
#' 
#' @examples
#' test_stat_FM_RD(3, 4, 10, 10, 0.2, "high")
test_stat_FM_RD <- function(x_E, x_C, n_E, n_C, delta, better){
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
  p_E0 <- round(2*u*cos(w) - b/(3*a), 10)
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


#' Calculate Farrington-Manning RR test statistic
#' 
#' \code{test_stat_FM_RR} returns the value of the Farrington-Manning test statistic
#' for non-inferiority of the risk ratio between two proportions.
#' 
#' If higher values of \mjseqn{x_E} favor the alternative hypothesis, we are interested
#' in testing the null hypothesis
#' \mjsdeqn{H_0: p_E / p_C \le \delta ,}
#' where the NI-margin is usually smaller than 1: \mjseqn{\delta < 1}.
#' The test statistic for this hypothesis is
#' \mjsdeqn{T_{\delta}(x_E, x_C) = \frac{\hat p_E - \delta \cdot \hat p_C}{\sqrt{\frac{\tilde p_E(1 - \tilde p_E)}{n_E} + \delta^2\frac{\tilde p_C(1 - \tilde p_C)}{n_C}}},}
#' where \mjseqn{\tilde p_C = \tilde p_C(x_E, x_C)} is the MLE of \mjseqn{p_C} and
#' \mjseqn{\tilde p_E = \tilde p_C + \delta} is the MLE of \mjseqn{p_E$ under \mjseqn{p_E / p_C = \delta}.
#' High values of \mjseqn{T_{\delta}} favor the alternative hypothesis.
#' 
#' @param x_E Vector of number of events in experimental group.
#' @param x_C Vector of number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of \mjseqn{x_E} favor the alternative 
#' hypothesis and "low" vice versa.
#' @return Vector of values of the RD test statistic.
#' 
#' @export
#' 
#' @examples
#' test_stat_FM_RR(3, 4, 10, 10, 0.2, "high")
test_stat_FM_RR <- function(x_E, x_C, n_E, n_C, delta, better){
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


#' Calculate Boschloo's test statistic
#' 
#' \code{test_stat_Boschloo_OR} returns the value of Boschloo's test statistic
#' for non-inferiority of the odds ratio between two proportions.
#' 
#' If higher values of \mjseqn{x_E} favor the alternative hypothesis, we are interested
#' in testing the null hypothesis
#' \mjsdeqn{H_0: p_E/(1-p_E) / (p_C/(1-p_C)) \le \delta ,}
#' where the NI-margin is usually smaller than 1: \mjseqn{\delta < 1}.
#' The test statistic
#' \mjsdeqn{ T_{\delta} = 1-(1 - F_{\mbox{ncHg}(X_E+X_C, n_E, n_C, \delta)}(x_E-1)) }
#' is based on Fisher's non-central hypergeometric distribution with density
#' \mjsdeqn{ f_{\mbox{ncHg}(s, n_E, n_C, \delta)}(k) = \frac{\binom{n_E}{k}\cdot \binom{n_C}{s-k}\cdot \delta^k}{\sum\limits_{l \in A_{s, n_E, n_C}} \binom{n_E}{l}\cdot \binom{n_C}{s-l}\cdot \delta^l}, }
#' where \mjseqn{A_{s, n_E, n_C} = \{\max(0, s-n_C), \dots, \min(n_E, s)\}}.
#' The density is zero if \mjseqn{k < \max(0, s-n_C)} or \mjseqn{k > \min(n_E, s)}.
#' High values of \mjseqn{T_{\delta}} favor the alternative hypothesis.
#' 
#' @param x_E Vector of number of events in experimental group.
#' @param x_C Vector of number of events in control group.
#' @param n_E Sample size in experimental group.
#' @param n_C Sample size in control group.
#' @param delta Non-inferiority margin.
#' @param better "high" if higher values of \mjseqn{x_E} favor the alternative 
#' hypothesis and "low" vice versa.
#' @return Vector of values of the RD test statistic.
#' 
#' @export
#' 
#' @examples
#' test_stat_Boschloo_OR(3, 4, 10, 10, 0.2, "high")
test_stat_Boschloo_OR <- function(x_E, x_C, n_E, n_C, delta, better){
  if(better == "high"){
    return <- sapply(
      1:length(x_E),
      function(i) BiasedUrn::pFNCHypergeo(x_E[i]-1, n_E, n_C, x_E[i]+x_C[i], delta)
    )
  }
  if(better == "low"){
    return <- sapply(
      1:length(x_E),
      function(i) 1 - BiasedUrn::pFNCHypergeo(x_E[i]-1, n_E, n_C, x_E[i]+x_C[i], delta)
    )
  }
  return(return)
}

