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


library(tidyverse)

# function to calculate test statistic for Risk Difference
test_RD <- function(x_E, x_C, n_E, n_C, delta, better){
  p_E <- x_E / n_E
  p_C <- x_C / n_C
  
  theta <- n_C / n_E
  p1hat <- x_E / n_E
  p2hat <- x_C / n_C
  a <- 1 + theta
  b <- -(1 + theta + p1hat + theta * p2hat + (-delta) * (theta + 2))
  c <- (-delta)^2 + (-delta) * (2*p1hat + theta + 1) + p1hat + theta * p2hat
  d <- - p1hat * (-delta) * (1 + (-delta))
  # Define the parameters for solving the equation
  v <- b^3/(3*a)^3 - b*c/(6*a^2) + d/(2*a)
  u <- sign(v) * sqrt(b^2/(3*a)^2 - c/(3*a))
  w <- 1/3 * (pi + acos(ifelse(u == 0, 0, round(v/u^3, 10))))   # das round wurde eingebaut, weil wenn bei v/u^3 etwas minimal gr??er als 1 rauskommt es zu einer Fehlermeldung kommt! 
  # Define the solution
  p_E0 <- 2*u*cos(w) - b/(3*a)
  p_C0 <- p_E0 - (-delta)
  
  if (better == "high"){
   return <-  ifelse(p_E - p_C + delta == 0, 0, -(p_E - p_C + delta) /sqrt(p_E0*(1-p_E0)/n_E + p_C0*(1-p_C0)/n_C))
  }
  if (better == "low"){
    return <- -ifelse(p_E - p_C + delta == 0, 0, -(p_E - p_C + delta) /sqrt(p_E0*(1-p_E0)/n_E + p_C0*(1-p_C0)/n_C))
  }
  return(return)
  
}


# function to calculate test statistic for Relative Risk
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
  p_C0 <- p_E0 / delta
  
  if (better == "high"){
    return <- ifelse(p_E - delta * p_C == 0, 0, -(p_E - delta * p_C) /sqrt(p_E0*(1-p_E0)/n_E + p_C0*(1-p_C0)*delta^2/n_C))
  }
  if (better == "low"){
    return <- -ifelse(p_E - delta * p_C == 0, 0, -(p_E - delta * p_C) /sqrt(p_E0*(1-p_E0)/n_E + p_C0*(1-p_C0)*delta^2/n_C))
  }
    return(return)
}



# function to compute test statistic for every response pair
teststat <- function(df, n_E, n_C, delta, method, better){
  # df is a data frame with all possible response pairs (x_E, x_C), i.e.
  # df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
  # n_E, n_C are the sample sizes in the groups
  # delta is the NI-margin.
  # method defines the specific test (statistic), i.e. "Risk ratio".
  # small values of stat favor the alternative.
  if (method == "RR") {
    return = df %>%
      mutate(
        stat = test_RR(x_E, x_C, n_E, n_C, delta, better)
      ) 
  }
   
  if (method == "RD") {
    return = df %>%
      mutate(
        stat = test_RD(x_E, x_C, n_E, n_C, delta, better)
      )
  }
  
  return(return)
}


# function to compute p_E from p_C and NI-margin delta s.t. effect(p_E, p_C) = delta
p_C.to.p_E <- function(p_C, method, delta){
  
  if (method == "RR") {
    p_E <- p_C * delta
  }
  if (method == "RD") {
    p_E <- p_C - delta
  }
  
  return(p_E)
}


# function to compute the critical value of a specific test statistic
critval <- function(alpha, n_C, n_E, method, delta, size_acc = 4, better){
  # method defines the used test with corresponding statistic and, size_acc defines
  # the accuracy of the grid used for the nuisance parameter p_C
  
  # Create data frame of all test statistics ordered by test statistic
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat(n_E, n_C, delta, method, better) %>%
    arrange(stat) ->
    df.stat
  
  # Extract stat, x_C and x_E as vector
  stat <- df.stat$stat
  x_C <- df.stat$x_C
  x_E <- df.stat$x_E
  
  # Find starting value for the search of critical value. E.g. take the
  # quantile of the approximate distribution of stat
  start_value <- -qnorm(1-alpha) 
  
  # Find row number of df.stat corresponding to starting value
  # <- row of df.stat where stat is maximal with stat <= start_value
  # Special case with very small sample sizes can lead to stat > start_value 
  # for all rows. Then set i <- 1
  i <- sum(stat<start_value)
  
  # Define rough grid for p_C
  acc <- 1
  p_C <- seq(10^-acc, 1-10^-acc, by = 10^-acc)
  
  # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
  # the null hypothesis
  p_E <- p_C.to.p_E(p_C, method, delta)
  
  p_C <- p_C[p_E >= 0 & p_E <= 1]
  p_E <- p_E[p_E >= 0 & p_E <= 1]
  
  # Calculate exact size for every pair (p_C, p_E)
  sapply(
    1:length(p_C),
    function(j) dbinom(x_C[1:i], n_C, p_C[j])*dbinom(x_E[1:i], n_E, p_E[j])
  ) ->
    size.vec
  
  # Increase index if maximal size is too low
  while (max(apply(size.vec, 2, sum)) <= alpha) {
    i <- i+1
    
    # Compute new sizes
    size.vec <- rbind(size.vec, dbinom(x_C[i], n_C, p_C)*dbinom(x_E[i], n_E, p_E))
  }
  
  # Decrease index if maximal size is too high and iterate grid accuracy
  for (acc in 1:size_acc) {
    
    # Define grid for p_C
    p_C <- seq(10^-acc, 1-10^-acc, by = 10^-acc)
    
    # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
    # the null hypothesis
    p_E <- p_C.to.p_E(p_C, method, delta)
    
    p_C <- p_C[p_E >= 0 & p_E <= 1]
    p_E <- p_E[p_E >= 0 & p_E <= 1]
    
    sapply(
      1:length(p_C),
      function(j) dbinom(x_C[1:i], n_C, p_C[j])*dbinom(x_E[1:i], n_E, p_E[j])
    ) ->
      size.vec
    
    # Decrease index if maximal size is too high
    while (max(apply(size.vec, 2, sum)) > alpha & i >= 1) {
      # Compute new sizes
      size.vec <- size.vec[-i,]
      i <- i-1
    }
  }
 
  # Decrease index further as long as rows have the same test statistic value
  while (stat[i+1] == stat[i] & i >= 1) {
    size.vec <- size.vec - dbinom(x_C[i], n_C, p_C)*dbinom(x_E[i], n_E, p_E)
    i <- i-1
  }
  
  # Critical value can now be chosen between stat[i+1] and stat[i]
  crit.val.mid <- (stat[i+1] + stat[i])/2
  
  # Return range of critical values and maximal size
  return(
    list(
      crit.val.lb = stat[i],
      crit.val.mid = crit.val.mid,
      crit.val.ub = stat[i+1],
      max.size = max(apply(size.vec, 2, sum))
    )
  )
}


power <- function(df, n_C, n_E, p_CA, p_EA, better){
  # Take data frame df with variable x_C and x_E representing all possible
  # response pairs for group sizes n_C and n_E and variable reject indicating
  # whether coordinates belong to rejection region.
  # Compute exact prob. of rejection region for all pairs (p_CA, p_EA).
  
  if (
    n_C+1 != df %>% pull(x_C) %>% unique() %>% length() |
    n_E+1 != df %>% pull(x_E) %>% unique() %>% length()
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
  sapply(
    1:length(p_CA),
    function(i) {
      df %>%
        filter(reject) %>%
        mutate(prob = dbinom(x_C, n_C, p_CA[i])*dbinom(x_E, n_E, p_EA[i])) %>%
        pull(prob) %>%
        sum()
    }
  ) ->
    result
  names(result) <- paste(p_CA, p_EA, sep = ", ")
  return(result)
}





# Function to compute exact sample size
samplesize_exact <- function(p_EA, p_CA, delta, alpha, beta, r, size_acc = 3, method, better){
  # Calculate exact sample size for method "X and specified
  # level alpha, beta, allocation ratio r = n_E/n_C, true rates p_CA, p_EA and
  # OR-NI-margin delta
  # Accuracy of calculating the critical value can be specified by size_acc.
  # Output: Sample sizes per group (n_C, n_E), critical value and exact power.
  
  # Check wheter p_EA, p_CA lie in the alternative hypothesis
  #if (
  #  p_C.to.p_E(p_CA, method, delta) >= p_EA
  #) {
  #  stop("p_EA, p_CA has to belong to alternative hypothesis.")
  #}
  
  # Estimate sample size with approximate formula
  n_appr <- samplesize_appr(
    p_EA = p_EA,
    p_CA = p_CA,
    delta = delta,
    alpha = alpha,
    beta = beta,
    r = r,
    method = method
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
    mutate(reject = stat <= crit.val) %>%
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
      crit.val <- critval(alpha = alpha, n_C = n_C, n_E = n_E, delta = delta, size_acc = size_acc, method = method, better = better)["crit.val.mid"]
      
      # Calculate exact power
      df %>%
        mutate(reject = stat <= crit.val) %>%
        power(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA, better = better) ->
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
      teststat(n_C = n_C, n_E = n_E, delta = delta, method = method) ->
      df
    
    # Calculate raised nominal level
    crit.val <- critval(alpha = alpha, n_C = n_C, n_E = n_E, delta = delta, size_acc = size_acc, method = method)["crit.val.mid"]
    
    # Calculate exact power
    df %>%
      mutate(reject = stat <= crit.val) %>%
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

# function to compute p-values for a specific result
p_value <- function(x_E., x_C., n_E, n_C, method, delta, size_acc = 3){
  # Define grid for p_C
  p_C <- seq(10^-size_acc, 1-10^-size_acc, by = 10^-size_acc)
  
  # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
  # the null hypothesis
  p_E <- p_C.to.p_E(p_C, method, delta)
  
  p_C <- p_C[p_E >= 0 & p_E <= 1]
  p_E <- p_E[p_E >= 0 & p_E <= 1]
  
  df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
  
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat(n_E, n_C, delta, method) %>%
    mutate(
      reject = stat <= stat[x_E == x_E. & x_C == x_C.]
    ) %>%
    power(n_C, n_E, p_C, p_E) ->      # function power actually computes the rejection probability which in this case is the p-value
    p.values
  
  # return vector of p-values. The "one" p-value would be max(p.values).
  return(p.values)
}




# Funktionsaufruf:

n_E = 22
n_C = 18

df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
method = "RR"
delta = 0.7

method = "RD"
delta = 0.1
alpha = 0.05

cv <- critval(alpha = alpha, n_C = n_C, n_E = n_E, method = method, delta = delta, size_acc = 4)

p_CA = c(0.5)
p_EA = c(0.7)

reject <- teststat(df, n_E, n_C, delta, method)$stat < cv$crit.val.mid
df <- data.frame(df, reject = reject)

power(df, n_C, n_E, p_CA, p_EA)


samplesize_appr(p_EA, p_CA, delta = 0.1, alpha, beta=0.2, r=1, method = "RD")
samplesize_exact(p_EA, p_CA, delta = 0.1, alpha, beta = 0.2, r = 1, size_acc = 3, method = "RD")  


samplesize_appr(p_EA, p_CA, delta = 0.8, alpha, beta=0.2, r = 1, method = "RR")
samplesize_exact(p_EA, p_CA, delta = 0.8, alpha, beta = 0.2, r = 1, size_acc = 3, method = "RR")
  

max(p_value(15, 15, 50, 50, method = "RD", delta = 0.1, size_acc = 3))
max(p_value(20, 15, 50, 50, method = "RR", delta = 0.8, size_acc = 3))


# Example 8.6
p_CA = 0.4
p_EA = 0.4
delta = 0.105
r = 2
alpha = 0.025
beta = 0.15
samplesize_appr(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = r, method = "RD")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = r, size_acc = 3, method = "RD")


p_CA = 0.32
p_EA = 0.32
delta = 0.095
r = 2
alpha = 0.025
beta = 0.2
samplesize_appr(p_EA, p_CA, delta = delta, alpha = alpha, beta=beta, r=r, method = "RD")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = r, size_acc = 3, method = "RD")


# stimmt beides für r = 1, r = 2

# Example 8.7
p_CA = 0.0385
p_EA = 0.155
delta = 0.8
r = 1
alpha = 0.025
beta = 0.2
samplesize_appr(p_EA, p_CA, delta = delta, alpha = alpha, beta=beta, r=r, method = "RR")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = r, size_acc = 3, method = "RR", better = "high")

teststat(df, n_E, n_C, delta, method = "RR", better = "low")

