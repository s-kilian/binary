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

# function to compute test statistic for every response pair
teststat <- function(df, n_E, n_C, delta, method){
  # df is a data frame with all possible response pairs (x_E, x_C), i.e.
  # df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
  # n_E, n_C are the sample sizes in the groups
  # delta is the NI-margin.
  # method defines the specific test (statistic), i.e. "Risk ratio".
  # High values of stat favor the alternative.
  
  if (method == "RR") {
    df %>%
      mutate(
        stat = #... 
      ) %>%
      return()
  }
  
  if (method == "RD") {
    df %>%
      mutate(
        stat = #... 
      ) %>%
      return()
  }
}

# function to compute p_E from p_C and NI-margin delta s.t. effect(p_E, p_C) = delta
p_C.to.p_E <- function(p_C, method, delta){
  # method defines the specific effect measure (risk ratio, risk difference, ...)
  p_E <- rep(NA, length(p_C))
  
  if (method == "RR") {
    p_E <- p_C * delta
  }
  if (method == "RD") {
    p_E == p_C + delta
  }
}

# function to compute the critical value of a specific test statistic
critval <- function(alpha, n_C, n_E, method, delta, size_acc = 4){
  # method defines the used test with corresponding statistic and, size_acc defines
  # the accuracy of the grid used for the nuisance parameter p_C
  
  # Create data frame of all test statistics ordered by test statistic
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat(n_E, n_C, delta, method) %>%
    arrange(stat) ->
    df.stat
  
  # Extract stat, x_C and x_E as vector
  stat <- df.stat$stat
  x_C <- df.stat$x_C
  x_E <- df.stat$x_E
  
  # Find starting value for the search of critical value. E.g. take the
  # quantile of the approximate distribution of stat
  start_value # <- qX(alpha)
  
  # Find row number of df.stat corresponding to starting value
  i # <- row of df.stat where stat is maximal with stat <= start_value
  # Special case with very small sample sizes can lead to stat > start_value 
  # for all rows. Then set i <- 1
  
  # Define rough grid for p_C
  acc <- 1
  p_C <- seq(10^-acc, 1-10^-acc, by = 10^-acc)
  
  # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
  # the null hypothesis
  p_E <- p_C.to.p_E(p_C, method, delta)
  
  # Calculate exact size for every pair (p_C, p_E)
  sapply(
    1:length(p_C),
    function(j) dbinom(x_C[1:i], n_C, p_C[j])*dbinom(x_E[1:i], n_E, p_E[j])
  ) ->
    size.vec
  
  # Increase index if maximal size is too low
  while (max(size.vec) <= alpha) {
    skip.first <- TRUE
    i <- i+1
    
    # Compute new sizes
    size.vec <- size.vec + dbinom(x_C[i], n_C, p_C)*dbinom(x_E[i], n_E, p_E)
  }
  
  # Decrease index if maximal size is too high and iterate grid accuracy
  for (acc in 1:size_acc) {
    # Downstep with acc = 1 is not needed if Upstep has been executed with this accuracy
    if(skip.first){next()}
    
    # Define grid for p_C
    p_C <- seq(10^-acc, 1-10^-acc, by = 10^-acc)
    
    # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
    # the null hypothesis
    p_E <- p_C.to.p_E(p_C, method, delta)
    
    # Decrease index if maximal size is too high
    while (max(size.vec) > alpha & i >= 1) {
      # Compute new sizes
      size.vec <- size.vec - dbinom(x_C[i], n_C, p_C)*dbinom(x_E[i], n_E, p_E)
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
      max.size = max(size.vec)
    )
  )
}

power <- function(df, n_C, n_E, p_CA, p_EA){
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

# Function to compute approximate sample size
samplesize_normal_appr <- function(p_EA, p_CA, delta, alpha, beta, r, method){
  #...
  return(
    list(
      n_C = n_C,
      n_E = n_E
    )
  )
}

# Function to compute exact sample size
samplesize_exact <- function(p_EA, p_CA, delta, alpha, beta, r, size_acc = 3, method){
  # Calculate exact sample size for method "X and specified
  # level alpha, beta, allocation ratio r = n_E/n_C, true rates p_CA, p_EA and
  # OR-NI-margin delta
  # Accuracy of calculating the critical value can be specified by size_acc.
  # Output: Sample sizes per group (n_C, n_E), critical value and exact power.
  
  # Check wheter p_EA, p_CA lie in the alternative hypothesis
  if (
    p_C.to.p_E(p_CA, method, delta) >= p_EA
  ) {
    stop("p_EA, p_CA has to belong to alternative hypothesis.")
  }
  
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
    teststat(n_C = n_C, n_E = n_E, delta = delta, method = method) ->
    df
  
  # Calculate critical value
  crit.val <- critval(alpha = alpha, n_C = n_C, n_E = n_E, delta = delta, size_acc = size_acc, method = method)["crit.val.mid"]
  
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
      mutate(reject = cond_p <= crit.val) %>%
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
  
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat(n_E, n_C, delta, method) %>%
    mutate(
      reject = stat >= stat[x_E == x_E. & x_C == x_C.]
    ) %>%
    power(n_C, n_E, p_C, p_E) ->      # function power actually computes the rejection probability which in this case is the p-value
    p.values
  
  # return vector of p-values. The "one" p-value would be max(p.values).
  return(p.values)
}