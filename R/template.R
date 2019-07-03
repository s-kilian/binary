##..............................................................................
##               Project: Package binary (working title)
##               Purpose: Provide template to implement exact test and sample
##                        size calculation for arbitrary test statistic
##                 Input: None
##                Output: None
##      Date of creation: 2019-07-03
##   Date of last update: 2019-07-03
##                Author: Samuel Kilian
##..............................................................................

## Functions ###################################################################

# function to compute a specific test statistic for every response pair
teststat_X <- function(df, n_E, n_C){
  # df is a data frame with all possible response pairs (x_E, x_C), i.e.
  # df <- expand.grid(x_E = 0:n_E, x_C = 0:n_C)
  # n_E, n_C are the sample sizes in the groups
  df %>%
    mutate(
      stat = ... 
    ) %>%
    return()
}

# function to compute the critical value of a specific test statistic
critval <- function(alpha, n_C, n_E, teststat = teststat_X, size_acc = 4){
  # teststat is the function to compute the test statistic, size_acc defines
  # the accuracy of the grid used for the nuisance parameter p_C
  
  # Create data frame of all test statistics order by test statistic
  expand.grid(
    x_C = 0:n_C,
    x_E = 0:n_E
  ) %>%
    teststat_X(n_E, n_C) %>%
    arrange(stat)
  
  # Find starting value for the search of critical value. E.g. take the
  # quantile of the approximate distribution of stat
  start_value <- qX(alpha)
}