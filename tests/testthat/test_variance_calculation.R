# rm(list = ls())
#
# N <- 1000L
# n <- N / 10.0
# set.seed(123L)
# simulation <- lin_sim(N = N, pX = 3L, pZ = 2L)
# list2env(simulation, envir = globalenv())
# rm(simulation)
# X <- extract_X_from_sim(data)
# Y_tab <- extract_Y_tab_from_sim(data)
# M <- colnames(Y_tab)
# phi_tab <- simulation$phi_tab
# sampler <- MMSRS$new(n, N, X = X, probaModes = probsM, Y_tab = Y_tab,
#                      phi_tab = phi_tab, M = M)
