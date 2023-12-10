sim <- lin_sim(N = 20L, pX = 3L, pZ = 1L, seed = 123L)

res <- grid_sim(seed = 123L)

bias_MC(res$HT, res$totY[1L])

CV_MC(res$HT, res$totY[1L])
