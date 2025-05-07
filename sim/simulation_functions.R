library(Matrix)
library(mc2d)

# Generation of a vector of probabilities under a logistic model
gen_selection_probabilities <- function(Z, alpha)
{
  expit <- function(x) 1.0 / (1.0 + exp(-x))

  expit(Z %*% alpha) %>%
    as.vector()
}


# Randomly affect a mode (or non-response) to each unit
# independency within and between the modes
gen_choice_bimode <- function(I, p1, p2, Z, checkRank = FALSE,
                              modesName = c("m1", "m2"),
                              seed = NULL)
{
  if (!is.null(seed))
    set.seed(seed)

  N <- length(p1)

  R1 <- runif(N) <= p1
  R2 <- runif(N) <= p2

  modes <- rep("nr", N)

  modes[I & R1] <- modesName[1L]
  modes[I & !R1 & R2] <- modesName[2L]

  q <- ncol(Z)

  # While the rank for the submatrices are not enough
  # we rerun the simulation
  if (checkRank &&
      (Matrix::rankMatrix(Z[R1, , drop = FALSE]) < q ||
       Matrix::rankMatrix(Z[!R1 & R2, , drop = FALSE]) < q))
  {
    warning("Rank for submatrices are not enough")
    return(gen_choice_bimode(I, p1, p2, checkRank = TRUE, Z))
  }

  # probSelect : probability of selecting the chosen mode
  # if selected in sample S
  probSelect <- rep(NA_real_, N)

  probSelect[R1] <- p1[R1]
  probSelect[!R1 & R2] <- (1.0 - p1[!R1 & R2]) * p2[!R1 & R2]

  tibble(p1 = p1, r1 = R1, p2 = p2, r2 = R2,
         mode = modes, probSelect = probSelect)
}

# Generate the conditional expectations of the counterfactuals
gen_expY <- function(X,
                     beta1, Y1Law = "gaussian", sd1,
                     beta2, Y2Law = "gaussian", sd2 = sd1,
                     seed = NULL)
{
  if (!is.null(seed))
    set.seed(seed)

  if (Y1Law == "gaussian")
  {

    Y1exp <- as.vector(X %*% beta1)

    sd1 <- as.numeric(sd1)
    sdY1 <- rep(sd1, N)
  }

  else if (Y1Law == "exponential")
  {
    sd1 <- "null"
    Y1exp <- as.vector(abs(X %*% beta1))

    sdY1 <- Y1exp
  }

  ## UTILE ?
  else if (Y1Law == "split3")
  {
    if (is.null(Y1exp))
      Y1exp <- as.vector(abs(X %*% beta1))

    quantiles1 <- quantile(Y1exp, probs = c(1.0 / 3.0, 2.0 / 3.0))

    sd1 <- as.numeric(sd1)
    sdY1 <- rep(as.numeric(sd1), N)

    # Between first and second quantile of 3rd degree we add one
    Y1exp[Y1exp > quantiles1[1L]] <-
      Y1exp[Y1exp > quantiles1[1L]] + 1.0

    # After the second quantile of 3rd degree we add one more
    Y1exp[Y1exp > quantiles1[2L]] <-
      Y1exp[Y1exp > quantiles1[2L]] + 1.0
  }

  if (Y2Law == "gaussian")
  {

    Y2exp <- as.vector(X %*% beta2)

    sd2 <- as.numeric(sd2)
    sdY2 <- rep(sd2, N)
  }
  else if (Y2Law == "exponential")
  {
    sd2 <- "null"

    Y2exp <- as.vector(abs(X %*% beta2))

    sdY2 <- Y2exp
  }

  list(Y1exp = Y1exp, sdY1 = sdY1,
       Y2exp = Y2exp, sdY2 = sdY2)
}


# Generate the counterfactuals with the conditional distributions
gen_Y <- function(X,
                  beta1, sd1, Y1Law = "gaussian",
                  beta2, sd2 = sd1, Y2Law = "gaussian",
                  Y1exp = NULL, Y2exp = NULL,
                  rho = 0.0)
{

  N <- nrow(X)

  # If the expectations have not been generated yet
  recExpY <- is.null(Y1exp) || is.null(Y2exp)
  if (recExpY)
  {
    dataY <- gen_expY(X,
                      beta1, Y1Law, sd1,
                      beta2, Y2Law, sd2)

    Y1exp <- dataY$Y1exp
    sdY1 <- dataY$sdY1

    Y2exp <- dataY$Y2exp
    sdY2 <- dataY$sdY2
  }
  else
  {
    if (Y1Law == "gaussian")
      sdY1 <- rep(sd1, N)
    else if (Y1Law == "exponential")
      sdY1 <- Y1exp

    if (Y2Law == "gaussian")
      sdY2 <- rep(sd2, N)
    else if (Y2Law == "exponential")
      sdY2 <- Y2exp
  }

  # Affecting values to Y1 and Y2
  if (Y1Law %in% c("gaussian", "split3"))
    Y1 <- Y1exp + rnorm(n = N, sd = sdY1)

  else if (Y1Law == "exponential")
    Y1 <- rexp(n = N, rate = Y1exp^-1L)

  if (Y2Law == "gaussian")
    Y2 <- Y2exp + rnorm(n = N, sd = sdY2)

  else if (Y2Law == "exponential")
    Y2 <- rexp(n = N, rate = Y2exp^-1L)

  # Consider the case of inter-mode correlation
  # (use of the Iman and Conover method)
  # Note : the correlation is approximated
  if (rho != 0.0)
  {
    results <- cbind(Y1 = Y1, Y2 = Y2) %>%
      mc2d::cornode(target = rho, result = FALSE) %>%
      as_tibble()

    rhoExper <- cor(results$Y1, results$Y2)

    results <- results %>%
      mutate(rho = rhoExper)
  }
  else
    results <- tibble(Y1 = Y1, Y2 = Y2, rho = 0.0)

  if (recExpY)
  {
    results <- results %>%
      mutate(Y1exp = Y1exp, .before = "Y1") %>%
      mutate(Y2exp = Y2exp, .before = "Y2")
  }

  results
}

set_Yobs <- function(Y1, Y2, modes, modesName = c("m1", "m2"))
{
  N <- length(Y1)

  Yobs <- rep(NA_real_, N)
  m1 <- modesName[1L]
  Yobs[modes == m1] <- Y1[modes == m1]
  m2 <- modesName[2L]
  Yobs[modes == m2] <- Y2[modes == m2]

  Yobs
}

# Generates the measure-weights vector phi
set_phi <- function(phi = "eq", N)
{
  if (is.numeric(phi) && length(phi) == 1L)
    phi <- rep(phi, N)
  else if (is.numeric(phi))
    return(phi)
  else if (phi == "eq")
    phi <- rep(0.5, N)
  else
    stop("Unknown value for phi")

  ## SUITE UTILE ?
  #
  # else if (phi == "1/3")
  #   phi <- rep(1.0 / 3.0, N)
  #
  # else if (phi == "2/3")
  #   phi <- rep(2.0 / 3.0, N)
  #
  # else if (phi == "var")
  #   phi <- seq_len(N) / N

  # When there is no bias, phi = 0.5
  # When the bias increase (in absolute), phi tends to zero
  # (it gives more weights to the unbiased value)
  # else if (phi == "linear")
  # {
  #   expDeltas <- expY1s - expY2s
  #   maxBias <- max(abs(expDeltas))
  #   phi <- 0.5 - 0.5 * abs(expDeltas) / maxBias
  # }
  # else if (phi == "split3")
  # {
  #   expDeltas <- expY1s - expY2s
  #   absExpDeltas <- abs(expDeltas)
  #   quantiles <- quantile(absExpDeltas, c(1.0 / 3.0, 2.0 / 3.0))
  #   phi <- rep(0.5, N)
  #   phi[absExpDeltas >= quantiles[1L]] <- 1.0 / 3.0
  #   phi[absExpDeltas >= quantiles[2L]] <- 1.0 / 6.0
  # }

  phi
}

# Generates the second order inclusion probabilities matrix of a
# Simple Random Sampling (SRS) design
piMat_SRS <- function(n, N)
{
  piMat <- matrix(n * (n - 1L) / (N * (N - 1L)), nrow = N, ncol = N)
  diag(piMat) <- n / N

  piMat
}

# Generates the covariance matrix of a Simple Random Sampling (SRS) design
covar_SRS <- function(n, N)
{
  ## pi2_to_covarInc DÉFINIE OÙ ?
  piMat_SRS(n, N) %>% pi2_to_covarInc()
}
