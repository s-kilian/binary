# Parameter definition for function development
p_EA = 0.7
p_CA = 0.4
delta = 0
alpha = 0.025
beta = 0.2
r = 1
size_acc = 3
eff_meas = "RD"
test_stat = NULL
better = "high"
start_value = NULL


ts_fun <- test_stat_FM_RD
ts_fun_args <- list(n_E = 10)
quant_fun <- qnorm
quant_fun_args = list()

arguments <- list(
  n_E = n_E,
  n_C = n_C,
  delta = delta
)
