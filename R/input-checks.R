##..............................................................................
##               Project: binary
##               Purpose: Provide functions to check input values.
##                 Input: none
##                Output: Functions:
##                        check.0.1: Value lies in open interval (0, 1)?
##                        check.pos.int: Value is positive integer?
##                        check.delta.eff_meas.better: Delta fits eff_meas and better?
##                        check.delta.null: Default delta value
##                        check.effect.delta.better: Effect fits delta and better?
##      Date of creation: 2020-02-18
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

check.alpha <- function(
  alpha           # alpha that shall lie in the interval (0, 0.5]
){
  # Check whether alpha lie in interval (0, 0.5] and output error message if not
  if (any(alpha <= 0 | alpha > 0.5)) {
    stop("alpha has to be in interval (0, 0.5]")
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

check.delta.eff_meas.better <- function(
  delta,
  eff_meas,
  better
){
  # Check whether delta fits the selected eff_meas
  if(eff_meas == "RD"){
    if(abs(delta) >= 1) stop("The NI-margin for risk difference has to be in (-1, 1).")
    if(better == "high" & delta > 0) warning("If better == \"high\", the NI-margin is usually not greater than 0.")
    if(better == "low" & delta < 0) warning("If better == \"low\", the NI-margin is usually not smaller than 0.")
  }
  # Check whether delta fits specification of better
  if(eff_meas %in% c("RR", "OR")){
    if(delta <= 0) stop("The NI-margin for risk ratio and odds ratio has to be positive.")
    if(better == "high" & delta > 1) warning("If better == \"high\", the NI-margin is usually not greater than 1.")
    if(better == "low" & delta < 1) warning("If better == \"low\", the NI-margin is usually not smaller than 1.")
  }
}

check.delta.null <- function(
  eff_meas,
  delta
){
  if(is.null(delta)){
    if(eff_meas == "RD"){
      warning("delta not specified and set to 0 (superiority test regarding risk difference).")
      delta <- 0
    }
    if(eff_meas %in% c("RR", "OR")){
      warning("delta not specified and set to 1 (superiority test regarding risk ratio or odds ratio).")
      delta <- 1
    }
  }
  return(delta)
}

check.effect.delta.better <- function(
  p_E,
  p_C,
  eff_meas,
  delta,
  better
){
  # check if specified prob. lie in region of alternative hyp.
  effect <- effect(p_E = p_E, p_C = p_C, eff_meas = eff_meas)
  if(better == "high") fits.alt <- effect > delta
  if(better == "low") fits.alt <- effect < delta
  if(!fits.alt) stop(paste0("delta = ", delta, " is ", better, "er than assumed effect ", eff_meas, " = ", effect))
}

