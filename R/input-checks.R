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