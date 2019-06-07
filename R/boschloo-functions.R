##..............................................................................
##               Project: Fisher-Boschloo test
##               Purpose: Provide functions around Fisher-Boschloo test for the
##                        comparison of two groups with binary endpoint.
##                 Input: none
##                Output: Functions:
##                        Add.cond.p: Compute conditional p-values for every result
##                        Raise.level: Find critical value for Fisher-Boschloo test
##                        Compute.custom.reject.prob: Compute exact probability
##                          of rejecting H0 for an arbitrary rejection region
##                        Calculate.approximate.sample.size: for chi-square test
##                        Calculate.exact.sample.size: iteratively for Fisher-Boschloo test
##                        Calculate.exact.Fisher.sample.size: iteratively for
##                          Fisher's exact test
##      Date of creation: 2019-04-04
##   Date of last update: 2019-05-23
##                Author: Samuel Kilian
##..............................................................................


## Functions ###################################################################
## Superiority #################################################################
# Test problem:
# H0: p1 <= p0
# H1: p1 > p0

Add.cond.p <- function(df, n0, n1){
  # Take data frame df with variable x0 and x1 representing all possible
  # response pairs for group sizes n0 and n1 and add conditional fisher p-values
  # for H0: p1 <= p0.
  if (
    n0+1 != df %>% pull(x0) %>% unique() %>% length() |
    n1+1 != df %>% pull(x1) %>% unique() %>% length()
  ) {
    stop("Values of x0 and x1 have to fit n0 and n1.")
  }
  
  # Compute p-values of Fisher's exact test from hypergeometric distribution
  # for every s
  df %>%
    mutate(s = x0+x1) %>%
    group_by(s) %>%
    do(
      .,
      mutate(
        .,
        cond.p = phyper(x0, n0, n1, s[1])
      )
    ) %>%
    return()
}

Raise.level <- function(alpha, n0, n1, acc = 4){
  # Compute raised nominal level for Fisher-Boschloo test for true level alpha 
  # and sample sizes n0 and n1.
  # Accuracy of obtaining maximum size (dependent on p) can be defined by acc.
  # Output: Nominal level (critical value) and exact size.
  
  # Total sample size
  n <- n0+n1
  
  # Possible values for the total number of responders s
  s.area <- 0:n
  
  # Initiate elements for loop
  # Create list of p.values (test statistic) for every s
  p.value.list <- list()
  for (s in s.area) {
    p.value.list[[s+1]] <- phyper(max(s-n1, 0):min(s, n0), n0, n1, s)
  }
  
  # Ordered data frame of p-values mapped to every s
  data.frame(
    p.value = unlist(p.value.list),
    s = rep(s.area, c(1:min(n0, n1), rep(min(n0, n1)+1, max(n0, n1)-min(n0, n1)+1), n+1-(max(n0, n1)+1):n))
  ) %>%
    arrange(p.value, s) ->
    ordered.p.values
  
  # Vector of p-values and vector of s
  p.values <- ordered.p.values$p.value
  s.vec <- ordered.p.values$s
  
  # Start with critical value = alpha and define the corresponding index
  start.index <- sum(p.values <= alpha)
  
  # Calculate boundaries c(s) of rejection region for every s for first critical 
  # value = alpha
  ordered.p.values %>%
    group_by(s) %>%
    summarise(c = suppressWarnings(max(p.value[p.value <= alpha]))) %>%
    arrange(s) %>%
    pull(c) ->
    start.bounds
  
  # Determine rough approximation of critical value iteratively
  size <- 0
  i <- start.index
  bounds <- start.bounds
  while (size <= alpha) {
    # Iterate critical value to next step
    i <- i+1
    # Determine s where boundary changes in this step
    new.s <- s.vec[i]
    # Determine the new boundary for specific s
    new.c <- p.values[i]
    # Update boundaries of rejection region
    bounds[new.s+1] <- new.c
    # Calculate values to find approximate maximum of size
    order.help <- choose(n, s.area)*bounds
    # Determine p-values where size is approximately maximal
    max.ps <- s.area[order.help >= 0.9*max(order.help)]/n
    # Determine maximal size for the specific p-values
    size <- 0
    for (p in max.ps) {
      sum(bounds[bounds != -Inf]*dbinom(s.area[bounds != -Inf], n, p)) %>%
        max(c(size, .)) ->
        size
    }
  }
  # Exit function if size of smallest possible rejection region is too high
  if (i <= 1) {
    warning("The rejection region of the test is empty.")
    return(c(nom.alpha.mid = 0, size = 0))
  }
  
  # If two or more possible results have the same p-values, they have to fall
  # in the same region. The rejection region is shrinked until this condition
  # is fulfilled.
  while(p.values[i-1] == p.values[i] & i > 1){
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max())
    i <- i-1
  }
  # Exit function if size of smallest possible rejection region is too high
  if (i <= 1) {
    warning("The rejection region of the test is empty.")
    return(c(nom.alpha.mid = 0, size = 0))
  }
  
  bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max())
  i <- i-1
  # Create grid for p with 51 points to compute more accurate maximum of size.
  p <- seq(0, 1, by = 0.02)
  # Compute size for every p in grid and take maximum
  sapply(
    p,
    function(x) sum(bounds[bounds != -Inf]*dbinom(s.area[bounds != -Inf], n, x))
  ) %>%
    max() ->
    max.size
  # If maximum size is too high, shrink rejection region and compute new maximum
  # size
  while (max.size > alpha) {
    # Exit function if size of smallest possible rejection region is too high
    if (i <= 1) {
      warning("The rejection region of the test is empty.")
      return(c(nom.alpha.mid = 0, size = 0))
    }
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max())
    i <- i-1
    while(p.values[i-1] == p.values[i] & i > 1){
      bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max())
      i <- i-1
    }
    sapply(
      p,
      function(x) sum(bounds[bounds != -Inf]*dbinom(s.area[bounds != -Inf], n, x))
    ) %>%
      max() ->
      max.size
  }
  # Creaste grid for p with specified accuracy to compute maximum size with
  # desired accuracy
  p <- seq(0, 1, by = 10^-acc)
  # Compute maximum size
  sapply(
    p,
    function(x) sum(bounds[bounds != -Inf]*dbinom(s.area[bounds != -Inf], n, x))
  ) %>%
    max() ->
    max.size
  # If maximum size is too high, shrink rejection region
  while (max.size > alpha) {
    # Exit function if size of smallest possible rejection region is too high
    if (i <= 1) {
      warning("The rejection region of the test is empty.")
      return(c(nom.alpha.mid = 0, size = 0))
    }
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max())
    i <- i-1
    while(p.values[i-1] == p.values[i] & i > 1){
      bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max())
      i <- i-1
    }
    sapply(
      p,
      function(x) sum(bounds[bounds != -Inf]*dbinom(s.area[bounds != -Inf], n, x))
    ) %>%
      max() ->
      max.size
  }
  # Define nominal alpha as mean of highest p-value in rejection region and
  # lowest p-value in acceptance region
  nom.alpha.mid <- (p.values[i] + p.values[i+1])/2
  
  return(c(nom.alpha.mid = nom.alpha.mid, size = max.size))
}

Compute.custom.reject.prob <- function(df, n0, n1, p0, p1){
  # Take data frame df with variable x0 and x1 representing all possible
  # response pairs for group sizes n0 and n1, variable reject indicating
  # whether coordinates belong to rejection region.
  # Compute exact prob. of rejection region for all pairs (p0, p1).
  
  if (
    n0+1 != df %>% pull(x0) %>% unique() %>% length() |
    n1+1 != df %>% pull(x1) %>% unique() %>% length()
  ) {
    stop("Values of x0 and x1 have to fit n0 and n1.")
  }
  
  if (
    length(p0) != length(p1) |
    !all(p0 >= 0 & p0 <= 1 & p1 >= 0 & p1 <= 1)
  ) {
    stop("p0 and p1 must have same length and values in [0, 1].")
  }
  
  
  # compute uncond. size for every p
  sapply(
    1:length(p0),
    function(i) {
      df %>%
        filter(reject) %>%
        mutate(prob = dbinom(x0, n0, p0[i])*dbinom(x1, n1, p1[i])) %>%
        pull(prob) %>%
        sum()
    }
  ) ->
    result
  names(result) <- paste(p0, p1, sep = ", ")
  return(result)
}

Calculate.approximate.sample.size <- function(alpha, power, r = 1, p0, p1){
  # Calculate approximate sample size for normal approximation test for specified
  # level alpha, power, allocation ratio r = n.1/n.0 and true rates p0, p1.
  # Output: Sample sizes per group (n.0, n.1).
  p.0 <- (p0 + r*p1)/(1+r)
  Delta.A <- p1 - p0
  n.0 <- ceiling(1/r*(qnorm(1-alpha)*sqrt((1+r)*p.0*(1-p.0)) + qnorm(power)*sqrt(r*p0*(1-p0) + p1*(1-p1)))^2  / Delta.A^2)
  n.1 <- r*n.0 %>% ceiling()
  
  return(
    list(n.0, n.1)
  )
}

Calculate.exact.sample.size <- function(alpha, power, r = 1, p0, p1, size.acc = 4){
  # Calculate exact sample size for Fisher-Boschloo test and specified
  # level alpha, power, allocation ratio r = n1/n0 and true rates p0, p1.
  # Accuracy of calculating the critical value can be specified by size.acc.
  # Output: Sample sizes per group (n0, n1), nominal alpha and exact power.
  if (p0 >= p1) {
    stop("p1 has to be greater than p0.")
  }
  
  # Estimate sample size with approximate formula
  n.approx <- Calculate.approximate.sample.size(alpha, power, r, p0, p1)
  
  # Use estimates as starting values
  n0 <- n.approx[[1]]
  n1 <- n.approx[[2]]
  
  # Initiate data frame for starting sample size
  expand.grid(
    x0 = 0:n0,
    x1 = 0:n1
  ) %>%
    Add.cond.p(n0 = n0, n1 = n1) ->
    df
  
  # Calculate raised nominal level for starting values
  nom.alpha <- Raise.level(alpha, n0, n1, size.acc)["nom.alpha.mid"]
  
  # Calculate exact power for starting values
  df %>%
    mutate(reject = cond.p <= nom.alpha) %>%
    Compute.custom.reject.prob(n0, n1, p0, p1) ->
    exact.power
  
  # Decrease sample size if power is too high
  if(exact.power > power){
    while(exact.power > power){
      # Store power and nominal level of last iteration
      last.power <- exact.power
      last.alpha <- nom.alpha
      
      # Decrease sample size by minimal amount possible with allocation ratio r
      if (r >= 1) {
        n0 <- n0 - 1
        n1 <- ceiling(r*n0)
      } else {
        n1 <- n1 - 1
        n0 <- ceiling(1/r*n1)
      }
  
      # Initiate data frame
      expand.grid(
        x0 = 0:n0,
        x1 = 0:n1
      ) %>%
        Add.cond.p(n0 = n0, n1 = n1) ->
        df
      
      # Calculate raised nominal level
      nom.alpha <- Raise.level(alpha, n0, n1, size.acc)["nom.alpha.mid"]
      
      # Calculate exact power
      df %>%
        mutate(reject = cond.p <= nom.alpha) %>%
        Compute.custom.reject.prob(n0, n1, p0, p1) ->
        exact.power
    }
    # Go one step back
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    exact.power <- last.power
    nom.alpha <- last.alpha
  }
  
  # If power is too low: increase sample size until power is achieved
  while (exact.power < power) {
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    
    # Initiate data frame
    expand.grid(
      x0 = 0:n0,
      x1 = 0:n1
    ) %>%
      Add.cond.p(n0 = n0, n1 = n1) ->
      df
    
    # Calculate raised nominal level
    nom.alpha <- Raise.level(alpha, n0, n1, size.acc)["nom.alpha.mid"]
    
    # Calculate exact power
    df %>%
      mutate(reject = cond.p <= nom.alpha) %>%
      Compute.custom.reject.prob(n0, n1, p0, p1) ->
      exact.power
  }
  
  return(
    list(
      n0 = n0,
      n1 = n1,
      nom.alpha = nom.alpha,
      exact.power = exact.power
    )
  )
}

Calculate.exact.Fisher.sample.size <- function(alpha, power, r = 1, p0, p1){
  # Calculate exact sample size for Fisher's exact test and specified
  # level alpha, power, allocation ratio r = n1/n0 and true rates p0, p1.
  # Accuracy of calculating the critical value can be specified by size.acc.
  # Output: Sample sizes per group (n0, n1) and exact power.
  if (p0 >= p1) {
    stop("p1 has to be greater than p0.")
  }
  
  # Estimate sample size with approximate formula
  n.approx <- Calculate.approximate.sample.size(alpha, power, r, p0, p1)
  
  # Use estimates as starting values
  n0 <- n.approx[[1]]
  n1 <- n.approx[[2]]
  
  # Initiate data frame
  expand.grid(
    x0 = 0:n0,
    x1 = 0:n1
  ) %>%
    Add.cond.p(n0 = n0, n1 = n1) ->
    df
  
  # Calculate exact power
  df %>%
    mutate(reject = cond.p <= alpha) %>%
    Compute.custom.reject.prob(n0, n1, p0, p1) ->
    exact.power
  
  # Decrease sample size if power is too high
  if(exact.power > power){
    while(exact.power > power){
      # Store power and nominal level of last iteration
      last.power <- exact.power

      # Decrease sample size by minimal amount possible with allocation ratio r
      if (r >= 1) {
        n0 <- n0 - 1
        n1 <- ceiling(r*n0)
      } else {
        n1 <- n1 - 1
        n0 <- ceiling(1/r*n1)
      }
      
      # Initiate data frame
      expand.grid(
        x0 = 0:n0,
        x1 = 0:n1
      ) %>%
        Add.cond.p(n0 = n0, n1 = n1) ->
        df
      
      # Calculate exact power
      df %>%
        mutate(reject = cond.p <= alpha) %>%
        Compute.custom.reject.prob(n0, n1, p0, p1) ->
        exact.power
    }
    # Go one step back
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    exact.power <- last.power
  }
  
  # If power is too low: increase sample size until power is achieved
  while (exact.power < power) {
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    
    # Initiate data frame
    expand.grid(
      x0 = 0:n0,
      x1 = 0:n1
    ) %>%
      Add.cond.p(n0 = n0, n1 = n1) ->
      df
    
    # Calculate exact power
    df %>%
      mutate(reject = cond.p <= alpha) %>%
      Compute.custom.reject.prob(n0, n1, p0, p1) ->
      exact.power
  }
  
  return(
    list(
      n0 = n0,
      n1 = n1,
      exact.power = exact.power
    )
  )
}


# Non-Inferiority ##############################################################
# Test problem:
# H0: OR(p1, p0) <= delta
# H1: OR(p1, p0) > delta
# with 0 < delta < 1
Add.cond.p.NI <- function(df, n0, n1, delta){
  # Take data frame df with variable x0 and x1 representing all possible
  # response pairs for group sizes n0 and n1 and add conditional fisher p-values
  # for H0: OR(p1, p0) <= delta.
  if (
    n0+1 != df %>% pull(x0) %>% unique() %>% length() |
    n1+1 != df %>% pull(x1) %>% unique() %>% length()
  ) {
    stop("Values of x0 and x1 have to fit n0 and n1.")
  }
  
  # Compute p-values of Fisher's exact test from Fisher's noncentral 
  # hypergeometric distribution for every s
  df %>%
    mutate(s = x0+x1) %>%
    group_by(s) %>%
    do(
      .,
      mutate(
        .,
        cond.p = BiasedUrn::pFNCHypergeo(x0, n0, n1, s[1], 1/delta)
      )
    ) %>%
    return()
}

Raise.level.NI <- function(alpha, n0, n1, delta, acc = 3){
  # Compute raised nominal level for Fisher-Boschloo test for true level alpha 
  # and sample sizes n0 and n1.
  # Accuracy of obtaining maximum size (dependent on p) can be defined by acc.
  # Output: Nominal level (critical value) and exact size.
  
  # Total sample size
  n <- n0+n1
  
  # Possible values for the total number of responders s
  s.area <- 0:n
  
  # Initiate elements for loop
  # Create list of p.values (test statistic) for every s
  p.value.list <- list()
  for (s in s.area) {
    p.value.list[[s+1]] <- BiasedUrn::pFNCHypergeo(max(s-n1, 0):min(s, n0), n0, n1, s, 1/delta)
  }
  
  # Ordered data frame of p-values mapped to every s
  data.frame(
    p.value = unlist(p.value.list),
    s = rep(s.area, c(1:min(n0, n1), rep(min(n0, n1)+1, max(n0, n1)-min(n0, n1)+1), n+1-(max(n0, n1)+1):n))
  ) %>%
    arrange(p.value, s) ->
    ordered.p.values
  
  # Vector of p-value and vector of s
  p.values <- ordered.p.values$p.value
  s.vec <- ordered.p.values$s
  
  # Start with critical value = alpha and define the corresponding index
  start.index <- sum(p.values <= alpha)
  
  # Calculate boundaries c(s) of rejection region for every s for first critical 
  # value = alpha
  ordered.p.values %>%
    group_by(s) %>%
    summarise(c = max(p.value[p.value <= alpha])) %>%
    arrange(s) %>%
    pull(c) ->
    start.bounds
  
  # Determine maximal nominal alpha iteratively
  max.size <- 0
  i <- start.index
  bounds <- start.bounds
  
  # Help function to compute P(S=s) under constant odds ratio delta
  Compute.s.prob.vec <- function(p0){
    p1 <- 1/(1+(1-p0)/(delta*p0))
    sapply(
      s.area,
      function(x) {
        k <- max(x-n1, 0):min(x, n0)
        sum(choose(n0, k) * choose(n1, x - k) * 1/delta^k) * (1-p0)^n0 * p1^x * (1-p1)^(n1-x)
      }
    )
  }
  
  # Create grid with 9 points fo p0 (must not contain 0 or 1)
  p0 <- seq(0.1, 0.9, by = 10^-1)
  
  # Create list of probabilites P(S=s) for every p0 in grid
  lapply(
    p0,
    Compute.s.prob.vec
  ) ->
    s.prob.vec.list
  
  # Increase nominal alpha
  while (max.size <= alpha) {
    # Iterate critical value to next step
    i <- i+1
    # Determine s where boundary changes in this step
    new.s <- s.vec[i]
    # Determine the new boundary for specific s
    new.c <- p.values[i]
    # Update boundaries of rejection region
    bounds[new.s+1] <- new.c
    # Compute size for every p in grid and take maximum
    sapply(
      1:length(p0),
      function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
    ) %>%
      max() ->
      max.size
  }
  # Go one step back
  bounds[s.vec[i]+1] <- p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max()
  i <- i-1
  
  # If two or more possible results have the same p-values, they have to fall
  # in the same region. The rejection region is shrinked until this condition
  # is fulfilled.
  while(p.values[i-1] == p.values[i]){
    bounds[s.vec[i]+1] <- p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max()
    i <- i-1
  }
  
  # Compute maximal size with increasing accuracy
  for (grid.acc in 2:acc) {
    # Define grid
    p0 <- seq(10^-grid.acc, 1-10^-grid.acc, by = 10^-grid.acc)
    # Compute probabilities P(S=s)
    lapply(
      p0,
      Compute.s.prob.vec
    ) ->
      s.prob.vec.list
    # Compute maximum size
    sapply(
      1:length(p0),
      function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
    ) %>%
      max() ->
      max.size
    # Shrink rejection region if size is too high
    while (max.size > alpha) {
      # Reduce rejection region
      bounds[s.vec[i]+1] <- p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max()
      i <- i-1
      while(p.values[i-1] == p.values[i]){
        bounds[s.vec[i]+1] <- p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max()
        i <- i-1
      }
      # Compute maximum size
      sapply(
        1:length(p0),
        function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
      ) %>%
        max() ->
        max.size
    }
  }
  # If two or more possible results have the same p-values, they have to fall
  # in the same region. The rejection region is shrinked until this condition
  # is fulfilled.
  while(p.values[i-1] == p.values[i]){
    bounds[s.vec[i]+1] <- p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% tail(1) %>% max()
    i <- i-1
  }
  
  # Define nominal alpha as mean of highest p-value in rejection region and
  # lowest p-value in acceptance region
  nom.alpha.mid <- (p.values[i] + p.values[i+1])/2
  
  return(c(nom.alpha.mid = nom.alpha.mid, size = max.size))
}

Calculate.approximate.sample.size.NI <- function(alpha, power, r = 1, p0, p1, delta){
  # Calculate approximate sample size for approximate test for specified
  # level alpha, power, allocation ratio r = n.1/n.0, true rates p0, p1 and
  # OR-NI.margin delta.
  # Output: Sample sizes per group (n.0, n.1).
  theta.A <- p1*(1-p0)/(p0*(1-p1))
  n.0 <- ceiling(1/r*(qnorm(1-alpha) + qnorm(power))^2 * (1/(p1*(1-p1)) + r/(p0*(1-p0))) / (log(theta.A) - log(delta))^2)
  n.1 <- r*n.0 %>% ceiling()
  
  return(
    list(n.0, n.1)
  )
}

Calculate.exact.sample.size.NI <- function(alpha, delta, power, r = 1, p0, p1, size.acc = 3){
  # Calculate exact sample size for Fisher-Boschloo test and specified
  # level alpha, power, allocation ratio r = n1/n0, true rates p0, p1 and
  # OR-NI-margin delta.
  # Accuracy of calculating the critical value can be specified by size.acc.
  # Output: Sample sizes per group (n0, n1), nominal alpha and exact power.
  
  if (p0 >= p1) {
    stop("p1 has to be greater than p0.")
  }
  
  # Estimate sample size with approximate formula
  n.approx <- Calculate.approximate.sample.size.NI(alpha, power, r, p0, p1, delta)
  
  # Use estimates as starting values
  n0 <- n.approx[[1]]
  n1 <- n.approx[[2]]
  
  # Initiate data frame
  expand.grid(
    x0 = 0:n0,
    x1 = 0:n1
  ) %>%
    Add.cond.p.NI(n0 = n0, n1 = n1, delta = delta) ->
    df
  
  # Calculate raised nominal level
  nom.alpha <- Raise.level.NI(alpha, n0, n1, delta, size.acc)["nom.alpha.mid"]
  
  # Calculate exact power
  df %>%
    mutate(reject = cond.p <= nom.alpha) %>%
    Compute.custom.reject.prob(n0, n1, p0, p1) ->
    exact.power
  
  # Decrease sample size if power is too high
  if(exact.power > power){
    while(exact.power > power){
      # Store power and nominal level of last iteration
      last.power <- exact.power
      last.alpha <- nom.alpha
      
      # Decrease sample size by minimal amount possible with allocation ratio r
      if (r >= 1) {
        n0 <- n0 - 1
        n1 <- ceiling(r*n0)
      } else {
        n1 <- n1 - 1
        n0 <- ceiling(1/r*n1)
      }
      
      # Initiate data frame
      expand.grid(
        x0 = 0:n0,
        x1 = 0:n1
      ) %>%
        Add.cond.p.NI(n0 = n0, n1 = n1, delta = delta) ->
        df
      
      # Calculate raised nominal level
      nom.alpha <- Raise.level.NI(alpha, n0, n1, delta, size.acc)["nom.alpha.mid"]
      
      # Calculate exact power
      df %>%
        mutate(reject = cond.p <= nom.alpha) %>%
        Compute.custom.reject.prob(n0, n1, p0, p1) ->
        exact.power
    }
    # Go one step back
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    exact.power <- last.power
    nom.alpha <- last.alpha
  }
  
  # If power is too low: increase sample size until power is achieved
  while (exact.power < power) {
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    
    # Initiate data frame
    expand.grid(
      x0 = 0:n0,
      x1 = 0:n1
    ) %>%
      Add.cond.p.NI(n0 = n0, n1 = n1, delta = delta) ->
      df
    
    # Calculate raised nominal level
    nom.alpha <- Raise.level.NI(alpha, n0, n1, delta, size.acc)["nom.alpha.mid"]
    
    # Calculate exact power
    df %>%
      mutate(reject = cond.p <= nom.alpha) %>%
      Compute.custom.reject.prob(n0, n1, p0, p1) ->
      exact.power
  }
  
  return(
    list(
      n0 = n0,
      n1 = n1,
      nom.alpha = nom.alpha,
      exact.power = exact.power
    )
  )
}

Calculate.exact.Fisher.sample.size.NI <- function(alpha, delta, power, r = 1, p0, p1, size.acc = 3){
  # Calculate exact sample size for Fisher's exact test and specified
  # level alpha, power, allocation ratio r = n1/n0 and true rates p0, p1.
  # Accuracy of calculating the critical value can be specified by size.acc.
  # Output: Sample sizes per group (n0, n1) and exact power.
  if (p0 >= p1) {
    stop("p1 has to be greater than p0.")
  }
  
  # Estimate sample size with approximate formula
  n.approx <- Calculate.approximate.sample.size.NI(alpha, power, r, p0, p1, delta)
  
  # Use estimates as starting values
  n0 <- n.approx[[1]]
  n1 <- n.approx[[2]]
  
  # Initiate data frame
  expand.grid(
    x0 = 0:n0,
    x1 = 0:n1
  ) %>%
    Add.cond.p.NI(n0 = n0, n1 = n1, delta = delta) ->
    df
  
  # Calculate exact power
  df %>%
    mutate(reject = cond.p <= alpha) %>%
    Compute.custom.reject.prob(n0, n1, p0, p1) ->
    exact.power
  
  # Decrease sample size if power is too high
  if(exact.power > power){
    while(exact.power > power){
      # Store power and nominal level of last iteration
      last.power <- exact.power
      
      # Decrease sample size by minimal amount possible with allocation ratio r
      if (r >= 1) {
        n0 <- n0 - 1
        n1 <- ceiling(r*n0)
      } else {
        n1 <- n1 - 1
        n0 <- ceiling(1/r*n1)
      }
      
      # Initiate data frame
      expand.grid(
        x0 = 0:n0,
        x1 = 0:n1
      ) %>%
        Add.cond.p.NI(n0 = n0, n1 = n1, delta = delta) ->
        df
      
      # Calculate exact power
      df %>%
        mutate(reject = cond.p <= alpha) %>%
        Compute.custom.reject.prob(n0, n1, p0, p1) ->
        exact.power
    }
    # Go one step back
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    exact.power <- last.power
  }
  
  # If power is too low: increase sample size until power is achieved
  while (exact.power < power) {
    if (r >= 1) {
      n0 <- n0 + 1
      n1 <- ceiling(r*n0)
    } else {
      n1 <- n1 + 1
      n0 <- ceiling(1/r*n1)
    }
    
    # Initiate data frame
    expand.grid(
      x0 = 0:n0,
      x1 = 0:n1
    ) %>%
      Add.cond.p.NI(n0 = n0, n1 = n1, delta = delta) ->
      df
    
    # Calculate exact power
    df %>%
      mutate(reject = cond.p <= alpha) %>%
      Compute.custom.reject.prob(n0, n1, p0, p1) ->
      exact.power
  }
  
  return(
    list(
      n0 = n0,
      n1 = n1,
      exact.power = exact.power
    )
  )
}

