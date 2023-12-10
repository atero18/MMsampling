# N <- 1000L
#
# set.seed(123L)
# X <- cbind(seq_len(N), seq_len(N)^2L)
# eps1 <- rnorm(N)
# beta0 <- 10.0
# beta <- c(3.0, 5.0)
# Y1 <- beta0 + X %*% beta + eps1
# eps2 <- rnorm(N, sd = 3.0)
# Y2 <- Y1 + eps2
# eps3 <- rnorm(N, sd = 5.0)
# Y3 <- Y1 + eps3
#
# Y <- cbind(Y1, Y2, Y3)
# colnames(Y) <- c("Y1", "Y2", "Y3")
#
# mode <- sample(c("Y1", "Y2", "Y3"), size = N, replace = TRUE)
# n <- ceiling(N / 2.0)
# I <- srswor(n, N) == 1L
#
# plan <-
#   data.frame(I = I, R = I, mode = mode, Ytilde = get_value_by_mode(Y, mode))
