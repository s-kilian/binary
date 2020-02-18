
# Example 8.6
p_CA = 0.4
p_EA = 0.4
delta = -0.105
alpha = 0.025
beta = 0.15
#samplesize_appr(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = r, method = "RD")
#samplesize_FM_diff_rates(p_EA, p_CA, delta = -delta, better = "low", alpha, beta, r)
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = 2, size_acc = 3, method = "RD", better = "low")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = 1, size_acc = 3, method = "RD", better = "low")


p_CA = 0.32
p_EA = 0.32
delta = -0.095
alpha = 0.025
beta = 0.2
#samplesize_appr(p_EA, p_CA, delta = delta, alpha = alpha, beta=beta, r=r, method = "RD")
#samplesize_FM_diff_rates(p_EA, p_CA, delta = -delta, better = "low", alpha, beta, r)
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = 2, size_acc = 3, method = "RD", better = "low")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = 1, size_acc = 3, method = "RD", better = "low")


# Example 8.7
p_EA = 0.0385
p_CA = 0.055
delta = 1/1.25
alpha = 0.025
beta = 0.08
#samplesize_appr(p_EA, p_CA, delta = delta, alpha = alpha, beta=beta, r=r, method = "RR")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = 1, size_acc = 3, method = "RR", better = "low")
samplesize_exact(p_EA, p_CA, delta = delta, alpha = alpha, beta = beta, r = 2, size_acc = 3, method = "RR", better = "low")

