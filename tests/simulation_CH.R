library(dplyr)

N <- 1000L

listeZ <- 0L:1L
listeDeltas <- c(-0.25, 0.0, 0.25)
listeBetaX <- c(-5.0, 5.0)

#' @importFrom tidyr expand_grid
gridParams <- expand_grid(pZ = listeZ, delta = listeDeltas, betaX = listeBetaX)

simulations <- list()

for (i in seq_len(nrow(gridParams)))
{
  simulations[[i]] <- sim_CH(N = N, pZ = gridParams$pZ[i], delta = gridParams$delta[i], betaXtel = gridParams$betaX[i])
}

set.seed(123L)
sim1 <- sim_CH(N = N, pZ = 0L)
sim2 <- sim_CH(N = N, pZ = 1L)
gridIntTel <- make_grid_int_tel(N, seed = 123L)

gridIntTel <- gridIntTel %>%
  filter(pmEstim == "true_values") %>%
  #filter(deltaEstim == "CF") %>%
  #filter(imputation == "MCO") %>%
  distinct()

B <- 100L
res <- grid_sim(B = B, grid = gridIntTel, simulations = list(sim1, sim2))

res$models <- res$models %>% select(-checkEquality, -checkNullityBias)


res$MC %>% group_by(idModel) %>% summarize(meanHT = mean(HT)) %>% print(n = 40)

res$models %>% group_by(idProblem) %>%  arrange(CV_tot, .by_group = TRUE) %>% select(imputation, deltaEstim, CV_tot)
