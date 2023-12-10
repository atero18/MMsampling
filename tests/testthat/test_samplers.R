# rm(list = ls())
# set.seed(123L)
#
# library(ggplot2)
# library(tibble)
#
# N <- 1000L
# n <- 100L
# pi <- n / N
# simulation <- lin_sim(N, pX = 2L, pZ = 0L)
# X <- extract_X_from_sim(simulation)
# Y_tab <- extract_Y_tab_from_sim(simulation)
# phi_tab <- simulation$phi_tab
# probaModes <- simulation$probsM
# sampler <- MMSRS$new(n = n, N = N, X = X, Y_tab = Y_tab, phi_tab = phi_tab, probaModes = probaModes)
#
# B <- 1500L
# samples <- sampler$samplings(B = B)
#
# # Each sample is of size n
# expect_true(all(vapply(samples, function(sample) sum(sample$I) == n, integer(1L))))
#
# firstPlan <-  vapply(samples, function(sample) sample$I, integer(N))
# propInd <- rowSums(firstPlan) / B
#
# abs(prodInd - pi) <= qnorm(1.0 - 0.05 / 2.0) * sqrt(pi * (1.0 - pi) / B)
#
# summary(propInd)
#
# propInd <- tibble(i = seq_len(N), proba = propInd)
# propInd %>% ggplot(aes(x = proba)) + geom_histogram()
#
# pi <- n / N
# for (m in colnames(probaModes))
# {
#   maskMode <-
#     vapply(samples, function(sample) sample$respondents_mode(m), logical(N))
#   propIndm <- rowSums(maskMode) / B
#   summary(propIndm - pi * probaModes[, m])
# }
#
# sample <- samples[[1L]]
# sample$show()
