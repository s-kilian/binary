##..............................................................................
##               Project: Fisher-Boschloo test
##               Purpose: Provide functions around Fisher-Boschloo test for the
##                        comparison of two groups with binary endpoint.
##                 Input: none
##                Output: Functions:
##                        teststat_boschloo: Compute test statistic
##                        critval_boschloo: Find critical value
##                        power_boschloo: Compute exact probability
##                          of rejecting H_0 for an arbitrary rejection region
##                        samplesize_normal_appr: Compute sample size for
##                          normal approximation test
##                        samplesize_exact_boschloo: COmpute exact sample size
##                        samplesize_exact_Fisher: Compute sample size for
##                          Fisher's exact test
##                        Same functions for non-inferiority:
##                        teststat_boschloo_NI: Compute test statistic
##                        critval_boschloo_NI: Find critical value
##                        samplesize_Wang: Compute sample size for a
##                          normal approximation test by Wang
##                        samplesize_exact_boschloo_NI: COmpute exact sample size
##                        samplesize_exact_Fisher_NI: Compute sample size for
##                          Fisher's exact test
##                        p_value: Calculate p-value for NI (superiority is
##                          special case with gamma = 1)
##      Date of creation: 2019-04-04
##   Date of last update: 2020-02-18
##                Author: Samuel Kilian
##..............................................................................


## Superiority #################################################################
# Default test problem (alternative = "greater"):
# H_0: p_E <= p_C
# H_1: p_E > p_C

teststat_boschloo <- function(df, n_C, n_E){
  # Take data frame df with variable x_C and x_E representing all possible
  # response pairs for group sizes n_C and n_E and add test statistic (conditional
  # fisher p-values for H_0: p_E <= p_C).
  if (
    n_C+1 != df %>% dplyr::pull(x_C) %>% unique() %>% length() |
    n_E+1 != df %>% dplyr::pull(x_E) %>% unique() %>% length()
  ) {
    stop("Values of x_C and x_E have to fit n_C and n_E.")
  }
  
  # Compute p-values of Fisher's exact test from hypergeometric distribution
  # for every s
  df %>%
    dplyr::mutate(s = x_C+x_E) %>%
    dplyr::group_by(s) %>%
    dplyr::do(
      .,
      dplyr::mutate(
        .,
        cond_p = stats::phyper(x_C, n_C, n_E, s[1])
      )
    ) %>%
    return()
}

critval_boschloo <- function(alpha, n_C, n_E, size_acc = 4){
  # Compute raised nominal level for Fisher-Boschloo test for true level alpha 
  # and sample sizes n_C and n_E.
  # Accuracy of obtaining maximum size (dependent on p) can be defined by size_acc.
  # Output: Nominal level (critical value) and exact size.
  
  # Total sample size
  n <- n_C+n_E
  
  # Possible values for the total number of responders s
  s.area <- 0:n
  
  # Initiate elements for loop
  # Create list of p.values (test statistic) for every s
  p.value.list <- list()
  for (s in s.area) {
    p.value.list[[s+1]] <- stats::phyper(max(s-n_E, 0):min(s, n_C), n_C, n_E, s)
  }
  
  # Ordered data frame of p-values mapped to every s
  data.frame(
    p.value = unlist(p.value.list),
    s = rep(s.area, c(1:min(n_C, n_E), rep(min(n_C, n_E)+1, max(n_C, n_E)-min(n_C, n_E)+1), n+1-(max(n_C, n_E)+1):n))
  ) %>%
    dplyr::arrange(p.value, s) ->
    ordered.p.values
  
  # Vector of p-values and vector of s
  p.values <- ordered.p.values$p.value
  s.vec <- ordered.p.values$s
  
  # Start with critical value = alpha and define the corresponding index
  start.index <- sum(p.values <= alpha)
  
  # Calculate boundaries c(s) of rejection region for every s for first critical 
  # value = alpha
  ordered.p.values %>%
    dplyr::group_by(s) %>%
    dplyr::summarise(c = suppressWarnings(max(p.value[p.value <= alpha]))) %>%
    dplyr::arrange(s) %>%
    dplyr::pull(c) ->
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
      sum(bounds[bounds != -Inf]*stats::dbinom(s.area[bounds != -Inf], n, p)) %>%
        max(c(size, .)) ->
        size
    }
  }
  # Exit function if size of smallest possible rejection region is too high
  if (i <= 1) {
    warning("The rejection region of the test is empty.")
    return(c(nom_alpha_mid = 0, size = 0))
  }
  
  # If two or more possible results have the same p-values, they have to fall
  # in the same region. The rejection region is shrinked until this condition
  # is fulfilled.
  while(p.values[i-1] == p.values[i] & i > 1){
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
    i <- i-1
  }
  # Exit function if size of smallest possible rejection region is too high
  if (i <= 1) {
    warning("The rejection region of the test is empty.")
    return(c(nom_alpha_mid = 0, size = 0))
  }
  
  bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
  i <- i-1
  # Create grid for p with 51 points to compute more accurate maximum of size.
  p <- seq(0, 1, by = 0.02)
  # Compute size for every p in grid and take maximum
  sapply(
    p,
    function(x) sum(bounds[bounds != -Inf]*stats::dbinom(s.area[bounds != -Inf], n, x))
  ) %>%
    max() ->
    max.size
  # If maximum size is too high, shrink rejection region and compute new maximum
  # size
  while (max.size > alpha) {
    # Exit function if size of smallest possible rejection region is too high
    if (i <= 1) {
      warning("The rejection region of the test is empty.")
      return(c(nom_alpha_mid = 0, size = 0))
    }
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
    i <- i-1
    while(p.values[i-1] == p.values[i] & i > 1){
      bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
      i <- i-1
    }
    sapply(
      p,
      function(x) sum(bounds[bounds != -Inf]*stats::dbinom(s.area[bounds != -Inf], n, x))
    ) %>%
      max() ->
      max.size
  }
  # Create grid for p with specified accuracy to compute maximum size with
  # desired accuracy
  p <- seq(0, 1, by = 10^-size_acc)
  # Compute maximum size
  sapply(
    p,
    function(x) sum(bounds[bounds != -Inf]*stats::dbinom(s.area[bounds != -Inf], n, x))
  ) %>%
    max() ->
    max.size
  # If maximum size is too high, shrink rejection region
  while (max.size > alpha) {
    # Exit function if size of smallest possible rejection region is too high
    if (i <= 1) {
      warning("The rejection region of the test is empty.")
      return(c(nom_alpha_mid = 0, size = 0))
    }
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
    i <- i-1
    while(p.values[i-1] == p.values[i] & i > 1){
      bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
      i <- i-1
    }
    sapply(
      p,
      function(x) sum(bounds[bounds != -Inf]*stats::dbinom(s.area[bounds != -Inf], n, x))
    ) %>%
      max() ->
      max.size
  }
  # Define nominal alpha as mean of highest p-value in rejection region and
  # lowest p-value in acceptance region
  nom_alpha_mid <- (p.values[i] + p.values[i+1])/2
  
  return(c(nom_alpha_mid = nom_alpha_mid, size = max.size))
}

power_boschloo <- function(df, n_C, n_E, p_CA, p_EA){
  # Take data frame df with variable x_C and x_E representing all possible
  # response pairs for group sizes n_C and n_E, variable reject indicating
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
  sapply(
    1:length(p_CA),
    function(i) {
      df %>%
        dplyr::filter(reject) %>%
        dplyr::mutate(prob = stats::dbinom(x_C, n_C, p_CA[i])*stats::dbinom(x_E, n_E, p_EA[i])) %>%
        dplyr::pull(prob) %>%
        sum()
    }
  ) ->
    result
  names(result) <- paste(p_CA, p_EA, sep = ", ")
  return(result)
}


#' Calculate approximate sample size for normal approximation test
#' 
#' Calculate approximate sample size for normal approximation test for specified
#' level alpha, power, allocation ratio r = n_E/n_C and true rates p_CA, p_EA.
#' 
#'@return Sample sizes per group (n_C, n_E).
#' 
#' @template p_A
#' @template error_rate
#' @template r
#' 
#' @export
samplesize_normal_appr <- function(p_EA, p_CA, alpha, beta, r){
  p_0 <- (p_CA + r*p_EA)/(1+r)
  Delta_A <- p_EA - p_CA
  n_C <- ceiling(1/r*(stats::qnorm(1-alpha)*sqrt((1+r)*p_0*(1-p_0)) + stats::qnorm(1-beta)*sqrt(r*p_CA*(1-p_CA) + p_EA*(1-p_EA)))^2  / Delta_A^2)
  n_E <- r*n_C %>% ceiling()
  
  return(
    list(n_C = n_C, n_E = n_E)
  )
}


#' Calculate exact sample size for Fisher-Boschloo test 
#' 
#' Calculate exact sample size for Fisher-Boschloo test and specified
#' level alpha, power, allocation ratio r = n_E/n_C and true rates p_CA, p_EA.
#' Accuracy of calculating the critical value can be specified by size_acc.
#' 
#'@return Sample sizes per group (n_C, n_E), nominal alpha and exact power.
#' 
#' @template p_A
#' @template error_rate
#' @template r
#' @template size_acc
#' @template alternative 
#' 
#' @export
samplesize_exact_boschloo <- function(p_EA, p_CA, alpha, beta, r, size_acc = 4, alternative = "greater"){
  # Check if input is correctly specified
  check.0.1(
    c(p_EA, p_CA, alpha, beta),
    "p_EA, p_CA, alpha, beta have to lie in interval (0,1)."
  )
  check.pos.int(
    size_acc,
    "size_acc has to be positive integer."
  )
  
  if (
    any(
      sapply(
        list(p_EA, p_CA, alpha, beta, r, size_acc, alternative),
        length
      ) != 1
    )
  ) {
    stop("Input values have to be single values.")
  }
  
  if (r <= 0) {
    stop("r has to be greater than 0.")
  }
  
  if (!(alternative %in% c("less", "greater"))) {
    warning("alternative has to be \"less\" or \"greater\". Will be treated as \"greater\".")
    alternative <- "greater"
  }
  
  if (p_CA >= p_EA & alternative == "greater") {
    stop("p_EA has to be greater than p_CA.")
  }
  
  if (alternative == "less") {
    if (p_EA >= p_CA) {
      stop("p_CA has to be greater than p_EA.")
    }
    # Inverse p_EA and p_CA to test the switched hypothesis
    p_EA <- 1-p_EA
    p_CA <- 1-p_CA
  }
  
  # Estimate sample size with approximate formula
  n_appr <- samplesize_normal_appr(
    p_EA = p_EA,
    p_CA = p_CA,
    alpha = alpha,
    beta = beta,
    r = r
  )
  
  # Use estimates as starting values
  n_C <- n_appr[["n_C"]]
  n_E <- n_appr[["n_E"]]
  
  # Initiate data frame for starting sample size
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat_boschloo(n_C = n_C, n_E = n_E) ->
    df
  
  # Calculate raised nominal level for starting values
  nom_alpha <- critval_boschloo(alpha = alpha, n_C = n_C, n_E = n_E, size_acc = size_acc)["nom_alpha_mid"]
  
  # Calculate exact power for starting values
  df %>%
    dplyr::mutate(reject = cond_p <= nom_alpha) %>%
    power_boschloo(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
    exact_power
  
  # Decrease sample size if power is too high
  if(exact_power > 1-beta){
    while(exact_power > 1-beta){
      # Store power and nominal level of last iteration
      last_power <- exact_power
      last_alpha <- nom_alpha
      
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
        teststat_boschloo(n_C = n_C, n_E = n_E) ->
        df
      
      # Calculate raised nominal level
      nom_alpha <- critval_boschloo(alpha = alpha, n_C = n_C, n_E = n_E, size_acc = size_acc)["nom_alpha_mid"]
      
      # Calculate exact power
      df %>%
        dplyr::mutate(reject = cond_p <= nom_alpha) %>%
        power_boschloo(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
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
    nom_alpha <- last_alpha
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
      teststat_boschloo(n_C = n_C, n_E = n_E) ->
      df
    
    # Calculate raised nominal level
    nom_alpha <- critval_boschloo(alpha = alpha, n_C = n_C, n_E = n_E, size_acc = size_acc)["nom_alpha_mid"]
    
    # Calculate exact power
    df %>%
      dplyr::mutate(reject = cond_p <= nom_alpha) %>%
      power_boschloo(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
      exact_power
  }
  
  return(
    list(
      n_C = n_C,
      n_E = n_E,
      nom_alpha = nom_alpha,
      exact_power = exact_power
    )
  )
}


# Non-Inferiority ##############################################################
# Default test problem (alternative = "greater"):
# H_0: OR(p_E, p_A) <= gamma
# H_1: OR(p_E, p_A) > gamma
# with 0 < gamma < 1

teststat_boschloo_NI <- function(df, n_C, n_E, gamma){
  # Take data frame df with variable x_C and x_E representing all possible
  # response pairs for group sizes n_C and n_E and add conditional fisher p-values
  # for H_0: OR(p_E, p_A) <= gamma.
  if (
    n_C+1 != df %>% dplyr::pull(x_C) %>% unique() %>% length() |
    n_E+1 != df %>% dplyr::pull(x_E) %>% unique() %>% length()
  ) {
    stop("Values of x_C and x_E have to fit n_C and n_E.")
  }
  
  # Compute p-values of Fisher's exact test from Fisher's noncentral 
  # hypergeometric distribution for every s
  df %>%
    dplyr::mutate(s = x_C+x_E) %>%
    dplyr::group_by(s) %>%
    dplyr::do(
      .,
      dplyr::mutate(
        .,
        cond_p = BiasedUrn::pFNCHypergeo(x_C, n_C, n_E, s[1], 1/gamma)
      )
    ) %>%
    return()
}

critval_boschloo_NI <- function(alpha, n_C, n_E, gamma, size_acc = 3){
  # Compute raised nominal level for Fisher-Boschloo test for true level alpha 
  # and sample sizes n_C and n_E.
  # Accuracy of obtaining maximum size (dependent on p) can be defined by size_acc.
  # Output: Nominal level (critical value) and exact size.
  
  # Total sample size
  n <- n_C+n_E
  
  # Possible values for the total number of responders s
  s.area <- 0:n
  
  # Initiate elements for loop
  # Create list of p.values (test statistic) for every s
  p.value.list <- list()
  for (s in s.area) {
    p.value.list[[s+1]] <- BiasedUrn::pFNCHypergeo(max(s-n_E, 0):min(s, n_C), n_C, n_E, s, 1/gamma)
  }
  
  # Ordered data frame of p-values mapped to every s
  data.frame(
    p.value = unlist(p.value.list),
    s = rep(s.area, c(1:min(n_C, n_E), rep(min(n_C, n_E)+1, max(n_C, n_E)-min(n_C, n_E)+1), n+1-(max(n_C, n_E)+1):n))
  ) %>%
    dplyr::arrange(p.value, s) ->
    ordered.p.values
  
  # Vector of p-value and vector of s
  p.values <- ordered.p.values$p.value
  s.vec <- ordered.p.values$s
  
  # Start with critical value = alpha and define the corresponding index
  start.index <- sum(p.values <= alpha)
  
  # Calculate boundaries c(s) of rejection region for every s for first critical 
  # value = alpha
  ordered.p.values %>%
    dplyr::group_by(s) %>%
    dplyr::summarise(c = suppressWarnings(max(p.value[p.value <= alpha]))) %>%
    dplyr::arrange(s) %>%
    dplyr::pull(c) ->
    start.bounds
  
  # Determine maximal nominal alpha iteratively
  max.size <- 0
  i <- start.index
  bounds <- start.bounds
  
  # Help function to efficiently compute logarithmic binomial coefficient
  logchoose <- function(o, u){
    if(u > o){stop("u cannot be greater than o!")}
    sum(log(seq_len(o-max(o-u, u))+max(o-u, u))) - sum(log(seq_len(min(o-u, u))))
  }
  
  # Help function to compute P(S=s) under constant odds ratio gamma
  Compute.s.prob.vec <- function(p_CA){
    p_EA <- 1/(1+(1-p_CA)/(gamma*p_CA))
    k.range <- 0:n_C
    sapply(
      k.range,
      function(y) logchoose(n_C, y) - y*log(gamma)
    ) ->
      add.1
    s.minus.k.range <- 0:n_E
    sapply(
      s.minus.k.range,
      function(y) logchoose(n_E, y)
    ) ->
      add.2
    sapply(
      s.area,
      function(x){
        k <- max(x-n_E, 0):min(x, n_C)
        help.val <- add.1[k+1] + add.2[x-k+1]
        help.val <- help.val + n_C*log(1-p_CA) + (n_E-x)*log(1-p_EA) + x*log(p_EA)
        sum(exp(help.val))
      }
    )
  }
  
  # Create grid with 9 points fo p_CA (must not contain 0 or 1)
  p_CA <- seq(0.1, 0.9, by = 10^-1)
  
  # Create list of probabilites P(S=s) for every p_CA in grid
  lapply(
    p_CA,
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
      1:length(p_CA),
      function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
    ) %>%
      max() ->
      max.size
  }
  # Go one step back
  bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
  i <- i-1
  
  # If two or more possible results have the same p-values, they have to fall
  # in the same region. The rejection region is shrinked until this condition
  # is fulfilled.
  while(p.values[i+1] == p.values[i]){
    bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
    i <- i-1
  }
  
  # Compute maximal size with increasing accuracy
  for (grid.acc in 2:size_acc) {
    # Define grid
    p_CA <- seq(10^-grid.acc, 1-10^-grid.acc, by = 10^-grid.acc)
    # Compute probabilities P(S=s)
    lapply(
      p_CA,
      Compute.s.prob.vec
    ) ->
      s.prob.vec.list
    # Compute maximum size
    sapply(
      1:length(p_CA),
      function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
    ) %>%
      max() ->
      max.size
    # Shrink rejection region if size is too high
    while (max.size > alpha) {
      # Reduce rejection region
      bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
      i <- i-1
      while(p.values[i+1] == p.values[i]){
        bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
        i <- i-1
      }
      # Compute maximum size
      sapply(
        1:length(p_CA),
        function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
      ) %>%
        max() ->
        max.size
    }
  }
  # # If two or more possible results have the same p-values, they have to fall
  # # in the same region. The rejection region is shrinked until this condition
  # # is fulfilled.
  # while(p.values[i-1] == p.values[i]){
  #   bounds[s.vec[i]+1] <- suppressWarnings(p.values[1:(i-1)][s.vec[1:(i-1)] == s.vec[i]] %>% utils::tail(1) %>% max())
  #   i <- i-1
  # }
  # # Recalculate maximal size
  # sapply(
  #   1:length(p_CA),
  #   function(x) sum(bounds[bounds != -Inf]*s.prob.vec.list[[x]][bounds != -Inf])
  # ) %>%
  #   max() ->
  #   max.size
  
  # Define nominal alpha as mean of highest p-value in rejection region and
  # lowest p-value in acceptance region
  nom_alpha_mid <- (p.values[i] + p.values[i+1])/2
  
  return(c(nom_alpha_mid = nom_alpha_mid, size = max.size))
}


#' Calculate approximate sample size for approximate test
#' 
#' alculate approximate sample size for approximate test for specified
#' level alpha, power, allocation ratio r = n_E/n_C, true rates p_CA, p_EA and
#' OR-NI.margin gamma.
#' 
#'@return Sample sizes per group (n_C, n_E).
#' 
#' @template p_A
#' @template gamma
#' @template error_rate
#' @template r
#' 
#' @export
samplesize_Wang <- function(p_EA, p_CA, gamma, alpha, beta, r){
  theta_A <- p_EA*(1-p_CA)/(p_CA*(1-p_EA))
  n_C <- ceiling(1/r*(stats::qnorm(1-alpha) + stats::qnorm(1-beta))^2 * (1/(p_EA*(1-p_EA)) + r/(p_CA*(1-p_CA))) / (log(theta_A) - log(gamma))^2)
  n_E <- r*n_C %>% ceiling()
  
  return(
    list(n_C = n_C, n_E = n_E)
  )
}


#' Calculate exact sample size for Fisher-Boschloo test
#' 
#' Calculate exact sample size for Fisher-Boschloo test and specified
#' level alpha, power, allocation ratio r = n_E/n_C, true rates p_CA, p_EA and
#' OR-NI-margin gamma.
#' 
#' @details Accuracy of calculating the critical value can be specified by size_acc.
#' 
#' @template p_A
#' @template gamma
#' @template error_rate
#' @template r
#' @template size_acc
#' @template alternative
#' 
#'@return Sample sizes per group (n_C, n_E), nominal alpha and exact power.
#' 
#' @export 
samplesize_exact_boschloo_NI <- function(p_EA, p_CA, gamma, alpha, beta, r, size_acc = 3, alternative = "greater"){
  # Check if input is correctly specified
  check.0.1(
    c(p_EA, p_CA, alpha, beta),
    "p_EA, p_CA, alpha, beta have to lie in interval (0,1)."
  )
  check.pos.int(
    size_acc,
    "size_acc has to be positive integer."
  )
  
  if (
    any(
      sapply(
        list(p_EA, p_CA, gamma, alpha, beta, r, size_acc, alternative),
        length
      ) != 1
    )
  ) {
    stop("Input values have to be single values.")
  }
  
  if (any(c(r, gamma) <= 0)) {
    stop("r and gamma have to be positive.")
  }
  
  if (!(alternative %in% c("less", "greater"))) {
    warning("alternative has to be \"less\" or \"greater\". Will be treated as \"greater\".")
    alternative <- "greater"
  }

  if (p_EA*(1-p_CA)/(p_CA*(1-p_EA)) <= gamma & alternative == "greater") {
    stop("OR(p_EA, p_CA) has to be greater than gamma.")
  }
  
  if (alternative == "less") {
    if (p_EA*(1-p_CA)/(p_CA*(1-p_EA)) >= gamma) {
      stop("OR(p_EA, p_CA) has to be smaller than gamma.")
    }
    # Inverse p_EA, p_CA and gamma to test the switched hypothesis
    p_EA <- 1-p_EA
    p_CA <- 1-p_CA
    gamma <- 1/gamma
  }

  
  # Estimate sample size with approximate formula
  n_appr <- samplesize_Wang(
    p_EA = p_EA,
    p_CA = p_CA,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    r = r
  )
  
  # Use estimates as starting values
  n_C <- n_appr[["n_C"]]
  n_E <- n_appr[["n_E"]]
  
  # Initiate data frame
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat_boschloo_NI(n_C = n_C, n_E = n_E, gamma = gamma) ->
    df
  
  # Calculate raised nominal level
  nom_alpha <- critval_boschloo_NI(alpha = alpha, n_C = n_C, n_E = n_E, gamma = gamma, size_acc = size_acc)["nom_alpha_mid"]
  
  # Calculate exact power
  df %>%
    dplyr::mutate(reject = cond_p <= nom_alpha) %>%
    power_boschloo(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
    exact_power
  
  # Decrease sample size if power is too high
  if(exact_power > 1-beta){
    while(exact_power > 1-beta){
      # Store power and nominal level of last iteration
      last_power <- exact_power
      last_alpha <- nom_alpha
      
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
        teststat_boschloo_NI(n_C = n_C, n_E = n_E, gamma = gamma) ->
        df
      
      # Calculate raised nominal level
      nom_alpha <- critval_boschloo_NI(alpha = alpha, n_C = n_C, n_E = n_E, gamma = gamma, size_acc = size_acc)["nom_alpha_mid"]
      
      # Calculate exact power
      df %>%
        dplyr::mutate(reject = cond_p <= nom_alpha) %>%
        power_boschloo(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
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
    nom_alpha <- last_alpha
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
      teststat_boschloo_NI(n_C = n_C, n_E = n_E, gamma = gamma) ->
      df
    
    # Calculate raised nominal level
    nom_alpha <- critval_boschloo_NI(alpha = alpha, n_C = n_C, n_E = n_E, gamma = gamma, size_acc = size_acc)["nom_alpha_mid"]
    
    # Calculate exact power
    df %>%
      dplyr::mutate(reject = cond_p <= nom_alpha) %>%
      power_boschloo(n_C = n_C, n_E = n_E, p_CA = p_CA, p_EA = p_EA) ->
      exact_power
  }
  
  return(
    list(
      n_C = n_C,
      n_E = n_E,
      nom_alpha = nom_alpha,
      exact_power = exact_power
    )
  )
}

p_value_boschloo <- function(x_E., x_C., n_E, n_C, gamma, size_acc = 3, alternative = "greater"){
  # Calculate p-values for a specific result (x_E., x_C.)
  
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
        list(x_E., x_C., n_E, n_C, gamma, size_acc, alternative),
        length
      ) != 1
    )
  ) {
    stop("Input values have to be single values.")
  }
  
  if (gamma <= 0) {
    stop("gamma has to be positive.")
  }
  
  if (!(alternative %in% c("less", "greater"))) {
    warning("alternative has to be \"less\" or \"greater\". Will be treated as \"greater\".")
    alternative <- "greater"
  }
  
  if (alternative == "less") {
    # Inverse values for other-sided hypothesis
    x_E. <- n_E - x_E.
    x_C. <- n_C - x_C.
    gamma <- 1/gamma
  }
  
  # Define grid for p_C
  p_C <- seq(10^-size_acc, 1-10^-size_acc, by = 10^-size_acc)
  
  # Find corresponding values of p_E such that (p_C, p_E) lie on the border of
  # the null hypothesis
  p_E <- 1/(1+(1-p_C)/(gamma*p_C))
  
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat_boschloo_NI(n_E = n_E, n_C = n_C, gamma = gamma) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      reject = cond_p <= cond_p[x_E == x_E. & x_C == x_C.]
    ) %>%
    power_boschloo(n_C, n_E, p_C, p_E) ->      # function power actually computes the rejection probability which in this case is the p-value
    p.values
  
  return(p.values)
}
