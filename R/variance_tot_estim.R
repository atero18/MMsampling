covar_difference_HT <- function(Y1exp, Y2exp, I,
                                pi_mat,
                                p1_mat, p2_mat,
                                Z,
                                biasedMode,
                                refMode,
                                modes,
                                phi = rep(1.0, length(Y1exp)),
                                correcEstimWeights = FALSE,
                                constY = FALSE,
                                covarY1 = NULL,
                                covarY2 = NULL)
{
  if (all(phi == 0.0))
    return(0.0)

  weightedY1 <- phi * Y1exp
  weightedY2 <- phi * Y2exp

  # Sampling design variability (p, S)
  covarPi <- pi2_to_covarInc(pi_mat)
  if (constY)
  {
    ## TO DO
  }
  else
  {
    varp <- weightedY1 / pi %*% t(weightedY2 / pi) * covarPi
  }

  # m1 selection variability (q1, R1)

  covarq1 <- pi2_to_covarInc(p1_mat)
  if (constY)
  {
    ## TO DO
  }
  else
  {

    correctedY1 <- weightedY1 / p1

    if (correcEstimWeights)
    {
      estbPhi11 <- .estim_bPhi11(Y1exp, I, maskSr, p1, Z, phi, constY)
      correctedY1 <- correctedY1Sr - crossprod(Z, estbPhi11)
    }

    correctedY2 <- weightedY2 / p1bar

    if (correcEstimWeights)
    {
      estbPhi21 <- .estim_bPhi21(Y2exp, I, p1, R2, p2, Z, phi, constY)
      correctedY2 <- correctedY2 + crossprod(Z, estbPhi21)
    }

    varq1 <- -(correctedY1 / pi) %*% t(correctedY2 / pi) * pi_mat * covarq1
  }


  sum(varp + varq1)
}

var_difference_HT <- function(Y1exp, Y2exp, I,
                              pi_mat,
                              p1_mat, p2_mat,
                              Z,
                              biasedMode,
                              refMode,
                              modes,
                              phi = rep(1.0, length(Y2exp)),
                              correcEstimWeights = FALSE,
                              constY = FALSE,
                              covarY1 = NULL,
                              covarY2 = NULL)
{
  var_expansion_seq_phi1(Y1exp, I, pi_mat, p1_mat, Z, biasedMode,
                  modes, phi, correcEstimWeights, constY1, covarY1) -
    2.0 * covar_difference_HT(Y1exp, Y2exp, I,
                              pi_mat,
                              p1_mat, p2_mat,
                              Z,
                              biasedMode,
                              refMode,
                              modes,
                              phi,
                              correcEstimWeights,
                              constY,
                              covarY1,
                              covarY2) +
    var_expansion_seq_phi2(Y2exp, pi_mat, p1_mat, p2_mat, biasedMode,
                    refMode, modes, phi, correcEstimWeights,
                    constY2, covarY2)
}

#' @export
covar_HT_seq_phi1_phi2 <- function(Y1exp, Y2exp,
                                   pi_mat,
                                   p1_mat, p2_mat,
                                   phi = numeric(length(Y2exp)))
{

  if (all(phi == 0.0) || all(phi == 1.0))
    return(0.0)

  pi <- diag(pi_mat)
  p1 <- diag(p1_mat)

  prodexpY12Mat <- Y1exp %*% t(Y2exp)

  # mode 1 selection variability
  covarq1 <- pi2_to_covarInc(p1_mat)
  varq1 <- -prodexpY12Mat * pi_mat * covarq1 / (p1 %*% t(1.0 - p1))

  # sampling variability
  covarPi <- pi2_to_covarInc(pi_mat)
  varp <- prodexpY12Mat * covarPi

  vMat <- varq1 + varp

  sum((phi / pi) %*% t((1.0 - phi) / pi) * vMat)
}


## NEEDS TO BE COMPLETELY UPDATED
#' @export
var_estim_tot_BM <- function(modeTotBiased = "HT", modeTotRef = "HT",
                             calculTotal = "population",
                             Y1exp, Y2exp,
                             covarY1, covarY2,
                             pi_mat,
                             p1_mat, p2_mat,
                             phi = numeric(length(Y2exp)),
                             subResults = FALSE)
{

  return(NA_real_)
  N <- length(Y2exp)

  pi <- diag(pi_mat)
  covarPi <- pi2_to_covarInc(pi_mat)

  p1 <- diag(p1_mat)
  p1bar <- 1.0 - p1
  pq1BarMat <- 1.0 -
    matrix(p1, nrow = N, ncol = N, byrow = TRUE) -
    matrix(p1, nrow = N, ncol = N, byrow = FALSE) +
    p1_mat
  covarq1 <- pi2_to_covarInc(p1_mat)
  p1Mat <- p1bar %*% t(p1bar)
  invp1Mat <- p1Mat^-1L
  invp1barMat <- (p1bar %*% t(p1bar))^-1L

  p2 <- diag(p2_mat)
  covarq2 <- pi2_to_covarInc(p2_mat)


  invProbsMatSelecMat <- invp1barMat * (p2 %*% t(p2))^-1L

  prodexpY1Mat <- Y1exp %*% t(Y1exp)

  prodexpY2Mat <- Y2exp %*% t(Y2exp)

  expDeltas <- Y1exp - Y2exp

  phiBar <- 1.0 - phi

  ## Functions to update
  varPhi1 <- var_expansion_seq_phi1(Y1exp, covarY1, pi_mat, p1_mat, SD1, phi)
  varPhi2 <- var_expansion_seq_phi2(Y2exp, covarY2, pi_mat, p1_mat, p2_mat, SD2, phi)
  covarPhi12 <- covar_HT_seq_phi1_phi2(Y1exp, Y2exp, pi_mat, p1_mat, p2_mat, COV12, phi)

  # Variance of the total MB estimator
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY2 <- pi_mat * pq1BarMat * p2_mat * covarY2 * invProbsMatSelecMat

    varY1 <- pi_mat * p1_mat * covarY1 * invp1Mat

    varq2 <- prodexpY2Mat * pi_mat * pq1BarMat * covarq2 * invProbsMatSelecMat

    sumY1Y2 <- Y1exp / p1 + Y2exp / p1bar
    varq1 <- sumY1Y2 %*% t(sumY1Y2) * pi_mat * covarq1

    varp <- expDeltas %*% t(expDeltas) * covarPi

    vDelta <- varY2 + varY1 + varq2 + varq1 + varp

    b <- t(t(phi) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) / pi)

    varPhiDelta <- sum(b %*% t(b) * vDelta)
  }

  # Covariance between ^t_phi1 and ^t_phiDelta
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY1 <- pi_mat * p1_mat * covarY1 * invp1Mat

    weightedY1 <- Y1exp / p1
    varq1 <- weightedY1 %*% t(weightedY1 + Y2exp / p1bar) * pi_mat * covarq1

    varp <- Y1exp %*% t(expDeltas) * covarPi

    v1Delta <- varY1 + varq1 + varp

    #b <- t(phi) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) / pi

    covarPhi1Delta <- sum((phi / pi) %*% t(b) * v1Delta)
  }

  # Covariance between ^t_phi2 and ^t_phiDelta
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY2 <- pi_mat * pq1BarMat * p2_mat * covarY2 * invProbsMatSelecMat

    varq2 <- prodexpY2Mat * pi_mat * pq1BarMat * covarq2 * invProbsMatSelecMat

    weightedY2 <- Y2exp / p1bar
    varq1 <- weightedY2 %*% t(Y1exp / p1 + weightedY2) * pi_mat * covarq1

    varp <- -Y2exp %*% t(expDeltas) * covarPi

    v2Delta <- varY2 + varq2 + varq1 + varp

    #b <- t(phi) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) / pi

    covarPhi2Delta <- -sum((phiBar / pi) %*% t(b) * v2Delta)
  }

  varEstim <- varPhi1 +
    varPhi2 +
    varPhiDelta +
    2.0 * covarPhi12 -
    2.0 * covarPhi1Delta -
    2.0 * covarPhi2Delta

  if (!subResults)
    return(varEstim)
  else
  {
    c(expVar2 = varEstim,
      expVarPhi1 = varPhi1,
      expVarPhi2 = varPhi2,
      expCovarPhi12 = covarPhi12,
      expVarPhiDelta = varPhiDelta,
      expCovarPhi1Delta = covarPhi1Delta,
      expCovarPhi2Delta = covarPhi2Delta)
  }

}
