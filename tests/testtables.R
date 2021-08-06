

# damit R alle FUnktionen im Paket/Projekt kennt:
library(devtools)
load_all()



# Testtable for testing teststat and p_value

# Situationen noch hinzufügen, die Fehlermeldung erzeugen sollen

data.frame(
  x_E = rep(c(0, 41, 13, 0, 41),3),
  x_C = rep(c(0, 50,  2, 0, 37),3),
  n_E = rep(c(5, 50, 25, 5, 50),3),
  n_C = rep(c(1, 50, 30, 1, 50),3),
  delta = c(-0.1, 0, 0.9, -0.3, 0.05, 
            0.1, 1, 1.3, 0.9, 1.1, 
            0.1, 1.3, 1, 0.9, 1.1), 
  method = c(rep("RD", 5), rep("RR", 5), rep("OR", 5)), 
  better = rep(c("high", "high", "high", "low", "low"), 3)
) ->
  df

df = df[c(5, 7),]

ts <- as.numeric()
p_v <- as.numeric()

for(i in 1:dim(df)[1]){
  ts[i] <- teststat(data.frame(x_E = df$x_E[i], x_C = df$x_C[i]), df$n_E[i], df$n_C[i], df$delta[i], df$method[i], df$better[i])$stat
  
  p_v[i] <- p_value(
    x_E. = df$x_E[i],
    x_C. = df$x_C[i],
    n_E = df$n_E[i],
    n_C = df$n_C[i],
    method = df$method[i] ,
    delta = df$delta[i],
    size_acc = 3,
    better = df$better[i],
    calc_method = "uniroot")$p_max
  
  #p_value(
  #  x_E. = df$x_E[i],
  #  x_C. = df$x_C[i],
  #  n_E = df$n_E[i],
  #  n_C = df$n_C[i],
  #  method = df$method[i] ,
  #  delta = df$delta[i],
  #  size_acc = 3,
  #  better = df$better[i],
  #  calc_method = "grid search")
  
}

cbind(df,ts, p_v)


# Fehler in FUnktion: find_max_prob_uniroot

power(
  df = df,
  n_C = n_C,
  n_E = n_C,  # hier ist der Fehler!
  p_CA = p_CA,
  p_EA = p_EA
  
  
  
  
  
# Testtable for testing critval and power
# Situationen noch hinzufügen, die Fehlermeldung erzeugen sollen
  

data.frame(
    alpha = rep(c(0.01, 0.05, 0.9, 0.5, 0.005),3),
    p_EA = rep(c(0.7, 0.7, 0.7, 0.3, 0.2), 3),
    p_CA = rep(c(0.6, 0.5, 0.4, 0.4, 0.9), 3),
    n_E = rep(c(5, 50, 25, 5, 50),3),
    n_C = rep(c(1, 50, 30, 1, 50),3),
    delta = c(-0.1, 0, 0.9, -0.3, 0.05, 
              0.1, 1, 1.3, 0.9, 1.1, 
              0.1, 1.3, 1, 0.9, 1.1), 
    method = c(rep("RD", 5), rep("RR", 5), rep("OR", 5)), 
    better = rep(c("high", "high", "high", "low", "low"), 3)
  ) ->
    df
  
df = df[-2,]


pow <- as.numeric()
for(i in 1:dim(df)[1]){
  crit.val <- critval(alpha = df$alpha[i], n_C = df$n_C[i], n_E = df$n_E[i], 
                    method = df$method[i], delta = df$delta[i], size_acc = 2, better = df$better[i])["crit.val.mid"]

  expand.grid(
    x_C = 0:df$n_C[i],
    x_E = 0:df$n_E[i]
  ) %>%
    teststat(n_C = df$n_C[i], n_E = df$n_E[i], delta = df$delta[i], method = df$method[i], better = df$better[i]) ->
    df2
  
  # Calculate exact power
  df2 %>%
    dplyr::mutate(reject = stat >= crit.val) %>%
    power(n_C = df$n_C[i], n_E = df$n_E[i], p_CA = df$p_CA[i], p_EA = df$p_EA[i]) ->
    exact_power
  pow[i] = exact_power
  
  
}

pow




# für i = 2 gibt es Fehlermeldung, warum???
# Warnmeldungen:
#  1: Problem with `mutate()` input `stat`.
# i NaNs wurden erzeugt
# i Input `stat` is `test_RD(x_E, x_C, n_E, n_C, delta, better)`. 
# 2: Problem with `mutate()` input `stat`.
# i NaNs wurden erzeugt
# i Input `stat` is `test_RD(x_E, x_C, n_E, n_C, delta, better)`. 

# Problem mit x_C = 0, x_E = 0

alpha = df$alpha[i]
n_C = df$n_C[i]
n_E = df$n_E[i]
method = df$method[i]
delta = df$delta[i]
size_acc = 2
better = df$better[i]

df2 = expand.grid(
  x_C = 0:n_C,
  x_E = 0:n_E
)


t2 = teststat(df2, n_E, n_C, delta, method, better)

t = as.numeric()
for(i in 1:dim(df2)[1]){
t[i] = test_RD(df2$x_E[i], df2$x_C[i], n_E, n_C, delta, better)
}


return = df2 %>%
  dplyr::mutate(stat = test_RD(x_E, x_C, n_E, n_C, delta, better))


df2 %>%
  dplyr::mutate(stat = t)



# Testtable for samplesize_approx und samplesize_exact
# Situationen noch hinzufügen, die Fehlermeldung erzeugen sollen

data.frame(
  alpha = rep(c(0.01, 0.05, 0.98, 0.5, 0.005),3),
  p_EA = rep(c(0.7, 0.7, 0.7, 0.4, 0.2), 3),
  p_CA = rep(c(0.6, 0.5, 0.4, 0.4, 0.9), 3),
  beta = rep(c(0.2, 0.01, 0.4, 0.3, 0.2),3),
  r = rep(c(1, 5, 0.1, 2, 0.5),3),
  delta = c(-0.1, 0, -0.9, 0.3, 0.05, 
            0.1, 1, 1.3, 5, 1.1, 
            0.1, 1.3, 1, 5, 1.1), 
  method = c(rep("RD", 5), rep("RR", 5), rep("OR", 5)), 
  better = rep(c("high", "high", "high", "low", "low"), 3)
) ->
  df


# Fehler bei 3, 8, 13: Hier ist alpha > 0.5 (bei approx kommt was anderes raus als bei PASS und exact gibt Fehlermeldung)
# Fehler bei 2: Problem with `mutate()` input `stat`.
# Fehler bei 14, 15: vermutlich liegt es an low (crit.val ist immer Inf, egal wie hoch das n geht)


df = df[-2,3, 8, 13; 14 und 15 dauert super lange)]
sa = as.numeric()
se = as.numeric()
for(i in 1:dim(df)[1]){
  sa = rbind(sa, samplesize_appr(df$p_EA[i], df$p_CA[i], df$delta[i], df$alpha[i], df$beta[i], df$r[i], df$method[i], df$better[i]))
  # se = rbind(se, samplesize_exact(df$p_EA[i], df$p_CA[i], df$delta[i], df$alpha[i], df$beta[i], df$r[i], size_acc = 1, df$method[i], df$better[i]))
  
}
sa

cbind(df,sa)

samplesize_exact(df$p_EA[i], df$p_CA[i], df$delta[i], df$alpha[i], df$beta[i], df$r[i], size_acc = 1, df$method[i], df$better[i])


samplesize_appr(df$p_EA[i], df$p_CA[i], df$delta[i], df$alpha[i], df$beta[i], df$r[i], df$method[i], df$better[i])


i = 14
p_EA = df$p_EA[i]
p_CA = df$p_CA[i]
delta = df$delta[i]
alpha = df$alpha[i]
beta = df$beta[i]
r = df$r[i]
method = df$method[i]
better = df$better[i]

size_acc = 2


# i = 3: folgendes sollte nur ausgeführt werden, wenn n_C>1:         n_C <- n_C - 1
# i = 8: 

samplesize_exact(df$p_EA[i], df$p_CA[i], df$delta[i], df$alpha[i], df$beta[i], df$r[i], size_acc = 1, df$method[i], df$better[i])


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




or = p_E[k]/(1-p_E[k])/(p_C[k]/(1-p_C[k]))



