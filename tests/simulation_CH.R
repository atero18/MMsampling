rm(list = ls())

library(dplyr)

N <- 1000L

listeZ <- 0L:1L
listeDeltas <- c(-0.25, 0.0, 0.25)
listeBetaX <- c(-5.0, 5.0)

#' @importFrom tidyr expand_grid
gridParams <- expand_grid(pZ = listeZ, delta = listeDeltas, betaX = listeBetaX)

simulations <- list()

set.seed(123L)

for (i in seq_len(nrow(gridParams)))
{
  simulations[[i]] <- sim_CH(N = N, pZ = gridParams$pZ[i],
                             delta = gridParams$delta[i],
                             betaXtel = gridParams$betaX[i])
}


gridIntTel <- make_grid_int_tel(N, seed = 123L)

gridIntTel <- gridIntTel %>%
  filter(pmEstim == "true_values", pZ == 0L) %>%
  #filter(deltaEstim == "CF") %>%
  #filter(imputation == "MCO") %>%
  distinct()

B <- 10L
res <- grid_sim(B = B, grid = gridIntTel, simulations = simulations)

res$models <- res$models %>% select(-checkEquality, -checkNullityBias)


res$MC %>% group_by(idModel) %>% summarize(meanHT = mean(HT)) %>% print(n = 40)

res$models %>% group_by(idProblem) %>%
  arrange(CV_tot, .by_group = TRUE) %>%
  select(imputation, deltaEstim, CV_tot)

res$models %>%
  filter(imputation != "true_values") %>%
  group_by(idSample) %>%
  arrange(CV_tot, .by_group = TRUE) %>%
  slice_min(CV_tot, n = 1L)
