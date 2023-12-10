# rm(list = ls())
# set.seed(123L)
#
# N <- 1000L
# n <- N / 10L
# simulation <- lin_sim(pX = 2L, pZ = 1L, Mbar = 2L)
#
# X <- extract_X_from_sim(simulation)
# Y_tab <- extract_Y_tab_from_sim(simulation)
# phi_tab <- simulation$phi_tab
# probaModes <- simulation$probsM
#
# sampler <- MMSRS(n = n, N = N, X = X,
#                  phi_tab = phi_tab, probaModes = probaModes)
#
# res <- sampler$sampling()
# res$Ytilde <- Y_tab[, 1L]
