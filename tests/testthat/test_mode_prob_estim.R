set.seed(245L)

N <- 1000L

p <- 2L

varMat <- genVarMatrix(p)

library(MASS)
X <- mvrnorm(N, rep(0.0, p), varMat)

I <- rep(TRUE, N)
beta1 <- rnorm(n = p, mean = 0, sd = .5)
p1 <- 1.0 / (1 + exp(-X %*% beta1)) %>% as.vector()
R1 <- runif(N) <= p1

# With MAR

beta2 <- rnorm(n = p, mean = 0, sd = .5)
p2 <- 1.0 / (1 + exp(-X %*% beta2)) %>% as.vector()

R2 <- runif(N) <= p2

modes <- rep(NA_character_, N)
modes[R1] <- "1"
modes[!R1 & R2] <- "2"

estimProbs <-
  estim_response_prob_sequential(I, X, modes, c("1", "2"), chosenOnly = FALSE)

glue("MSE between p1 and estimated conditional :")
mean((estimProbs$conditional[, "1"] - p1)^2L)
glue("MSE between p1 and estimated unconditional :")
mean((estimProbs$unconditional[, "1"] - p1)^2L)

glue("MSE between p2 and estimated conditional :")
mean((estimProbs$conditional[, "2"] - p2)^2L)
glue("MSE between p2 and estimated unconditional :")
mean((estimProbs$unconditional[, "2"] - p2)^2L)

# With NMAR
betaNMAR <- sample(c(-10L, 10L), size = p)
p2 <- p2 <- 1.0 / (1 + exp(-X %*% beta2 + ifelse(R1, -X %*% betaNMAR, 0))) %>% as.vector()

R2 <- runif(N) <= p2

modes <- rep(NA_character_, N)
modes[R1] <- "1"
modes[!R1 & R2] <- "2"

estimProbs <- estim_response_prob_sequential(I, X, modes, c("1", "2"), chosenOnly = FALSE)

glue("MSE between p1 and estimated conditional :")
mean((estimProbs$conditional[, "1"] - p1)^2L)
glue("MSE between p1 and estimated unconditional :")
mean((estimProbs$unconditional[, "1"] - p1)^2L)

glue("MSE between p2 and estimated conditional :")
mean((estimProbs$conditional[, "2"] - p2)^2L)
glue("MSE between p2 and estimated unconditional :")
mean((estimProbs$unconditional[, "2"] - p2)^2L)
