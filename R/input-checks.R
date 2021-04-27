##..............................................................................
##               Project: Fisher-Boschloo test
##               Purpose: Provide functions to check input values.
##                 Input: none
##                Output: Functions:
##                        check.0.1: Value lies in open interval (0, 1)?
##                        check.pos.int: Value is positive integer?
##      Date of creation: 2020-02-18
##   Date of last update: 2020-02-18
##                Author: Samuel Kilian
##..............................................................................


## Checks       ################################################################
check.0.1 <- function(
  values,           # vector of values that shall lie in the interval (0, 1)
  message           # error message to display if one of the values doesn't lie
  # in (0,1)
){
  # Check whether values all lie in interval (0, 1) and output error message if not
  if (any(values <= 0 | values >= 1)) {
    stop(message)
  }
}

check.pos.int <- function(
  values,           # vector of values that shall be positive integers
  message           # error message to display if one of the values is not a
  # positive integer
){
  # Check whether values all are positive integers and output error message if not
  if (any(values %% 1 != 0 | values <= 0)) {
    stop(message)
  }
}

check.delta.method.better <- function(
  delta,
  method,
  better
){
  # Check whether delta fits the selected method
  if(method == "RD"){
    if(abs(delta) >= 1) stop("The NI-margin for risk difference has to be in (-1, 1).")
    if(better == "high" & delta > 0) warning("If better == \"high\", the NI-margin is usually not greater than 0.")
    if(better == "low" & delta < 0) warning("If better == \"low\", the NI-margin is usually not smaller than 0.")
  }
  # Check whether delta fits specification of better
  if(method %in% c("RR", "OR")){
    if(delta <= 0) stop("The NI-margin for risk ratio and odds ratio has to be positive.")
    if(better == "high" & delta > 1) warning("If better == \"high\", the NI-margin is usually not greater than 1.")
    if(better == "low" & delta < 1) warning("If better == \"low\", the NI-margin is usually not smaller than 1.")
  }
}

check.delta.null <- function(
  method,
  delta
){
  if(is.null(delta)){
    if(method == "RD"){
      warning("delta not specified and set to 0 (superiority test regarding risk difference).")
      delta <- 0
    }
    if(method %in% c("RR", "OR")){
      warning("delta not specified and set to 1 (superiority test regarding risk ratio or odds ratio).")
      delta <- 1
    }
  }
  return(delta)
}

check.effect.delta.better <- function(
  p_E,
  p_C,
  method,
  delta,
  better
){
  # check if specified prob. lie in region of alternative hyp.
  effect <- effect(p_E = p_E, p_C = p_C, method = method)
  if(better == "high") fits.alt <- effect > delta
  if(better == "low") fits.alt <- effect < delta
  if(!fits.alt) stop(paste0("delta = ", delta, " is ", better, "er than assumed effect ", method, " = ", effect))
}

