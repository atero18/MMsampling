N <- 1000L
gridIntTel <- make_grid_int_tel(N, seed = 123L)

gridIntTel <- gridIntTel %>%
  filter(pmEstim == "true_values") %>%
  #filter(deltaEstim == "CF") %>%
  #filter(imputation == "MCO") %>%
  distinct()

B <- 10L
res <- grid_sim(B = B, grid = gridIntTel)


res$MC %>% group_by(idModel) %>% summarize(meanHT = mean(HT)) %>% print(n = 40)
