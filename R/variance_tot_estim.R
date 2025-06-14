.estim_bPhi11 <- function(Yobs, I, pi, p1, maskSr, Z,
                    phi = rep(1.0, length(Yobs)))
{
  p1Sr <- p1[maskSr]

  # Estimation of the expected value of the derivative of the total estimator
  estVPhi11 <- -crossprod(Z[maskSr, , drop = FALSE],
                          (pi[maskSr] * p1Sr)^-1L *
                            phi[maskSr] * Yobs[maskSr] * (1.0 - p1Sr))


  solve(Fisher_Information_Matrix(p1, Z, I)) %*% estVPhi11
}


.estim_bPhi21 <- function(Yobs, I, pi, p1, p2, maskSmr, Z,
                          phi = rep(1.0, length(Y1exp)))
{

  p1Smr <- p1[maskSmr]

  estVPhi21 <- crossprod(Z[maskSmr, , drop = FALSE],
                         (pi[maskSmr] * (1.0 - p1Smr) * p2[maskSmr])^-1L *
                           phi[maskSmr] * Yobs[maskSmr] * p1Smr)



  solve(Fisher_Information_Matrix(p1, Z, I)) %*% estVPhi21
}


.estim_bPhi22 <- function(Yobs, pi, p1, p2, maskSm, maskSmr, Z,
                          phi = rep(1.0, length(Y1exp)))

{

  p2Smr <- p2[maskSmr]

  estVPhi22 <- -crossprod(Z[maskSmr, , drop = FALSE],
                          (pi[maskSmr] * (1.0 - p1[maskSmr]) * p2Smr)^-1L *
                            (1.0 - p2Smr) * phi[maskSmr] * Yobs[maskSmr])


  solve(Fisher_Information_Matrix(p2, Z, maskSm)) %*% estVPhi22
}

#' Estimates the approximate variance of the HT estimator of t_phi1 with known
#' or estimated mode selection probabilities. Homoscedasticity and independence
#' in the m1 mode selection mechanism assumed.
#' @export
estim_appr_var_seq_phi1 <- function(Yobs,
                                    modes,
                                    I,
                                    pi_mat,
                                    p1,
                                    sd1 = NULL,
                                    phi = rep(1.0, length(Yobs)),
                                    correcEstimWeights = FALSE,
                                    Z = matrix(1.0, nrow = len(Yobs),
                                               ncol = 1L))
{

  if (all(phi == 0.0))
    return(0.0)

  maskSr <- modes == "m1"

  if (is.null(sd1))
    sd1 <- sd(Yobs[maskSr], na.rm = FALSE)

  pi <- diag(pi_mat)
  piSr <- pi[maskSr]
  p1Sr <- p1[maskSr]

  # The y_1k are weighted with the phi_k

  weightedY1Sr <- phi[maskSr] * Yobs[maskSr]

  # Sampling design variability (p, S)
  # There is no correction needed for probabilities estimation
  piSr_mat <- pi_mat[maskSr, maskSr, drop = FALSE]
  covarpSr <- pi2_to_covarInc(piSr_mat)

  correctedY1Srp <- (piSr * p1Sr)^-1L * weightedY1Sr
  varpEst <-
    t(correctedY1Srp) %*%
    (covarpSr / piSr_mat) %*%
    correctedY1Srp %>%
    as.numeric()
  varpEst <- varpEst +
    sum((1.0 - piSr) * piSr^-2L * p1Sr^-1L * (1.0 - p1Sr^-1L) *
          weightedY1Sr^2L)

  # m1 selection variability (q1, R1)
  # If we use estimated p_1k weights we have to make a correction
  correctedY1Srq1 <- correctedY1Srp

  #   If the probabilities p_1k are estimations we add a linearisation term
  if (correcEstimWeights)
  {
    estbPhi11 <- .estim_bPhi11(Yobs, I, pi, p1, maskSr, Z, phi)
    correctedY1Srq1 <- correctedY1Srq1 + Z[maskSr, , drop = FALSE] %*% estbPhi11
  }

  varq1Est <- sum((1.0 - p1Sr) * correctedY1Srq1^2L)


  # Y1 variability
  varY1Est <- sum(phi^2L) * sd1^2L

  varpEst + varq1Est + varY1Est
}

#' Calculate the true variance of the HT estimator of t_phi1 in the case
#' of known selection probabilities. Using estimated values give an estimation
#' The expected values of the m1 counterfactuals have to be given.
#' The calculation will be correct under the additional strong assumption
#' of conditional independence between the sampling design p
#' and the mode selection mechanism q1.
#' @export
var_expansion_seq_phi1 <- function(Y1exp,
                                   pi_mat,
                                   p1,
                                   phi = rep(1.0, length(Y1exp)),
                                   sd1 = 0.0)
{

  if (all(phi == 0.0))
    return(0.0)

  pi <- diag(pi_mat)

  # Sampling design variability (p, S)
  covarPi <- pi2_to_covarInc(pi_mat)
  correctedY1p <- pi^-1L * phi * Y1exp
  varp <- t(correctedY1p) %*% covarPi %*% correctedY1p +
    sum(pi^-1L * (1.0 - pi) * phi^2L) * sd1^2L

  # m1 selection variability (q1, R1)
  #  With the independence between p and q1
  varq1 <- sum((pi * p1)^-1L * (1.0 - p1) * phi^2L * (sd1^2L + Y1exp^2L))


  # Y1 variability
  varY1 <- sum(phi^2L) * sd1^2L

  varp + varq1 + varY1
}

#' Estimates the approximate variance of the HT estimator of t_phi2 with known
#' or estimated mode selection probabilities. Homoscedasticity and independence
#' in the m1 and m2 mode selection mechanisms assumed.
#' @export
estim_appr_var_seq_phi2 <- function(Yobs,
                                    modes,
                                    I,
                                    pi_mat,
                                    p1,
                                    p2,
                                    sd2 = NULL,
                                    phi = rep(1.0, length(Yobs)),
                                    correcEstimWeights = FALSE,
                                    Z = matrix(1.0, nrow = len(Yobs),
                                               ncol = 1L))
{
  if (all(phi == 0.0))
    return(0.0)

  maskSmr <- modes == "m2"

  if (is.null(sd2))
    sd2 <- sd(Yobs[maskSmr])

  pi <- diag(pi_mat)
  piSmr <- pi[maskSmr]

  p1Smr <- p1[maskSmr]
  p1barSmr <- 1.0 - p1Smr

  p2Smr <- p2[maskSmr]

  # The y_2k are weighted with the phi_k
  weightedY2Smr <- phi[maskSmr] * Yobs[maskSmr]

  # Sampling design variability (p, S)
  piSmr_mat <- pi_mat[maskSmr, maskSmr, drop = FALSE]
  covarpSmr <- pi2_to_covarInc(piSmr_mat)

  correctedY2Smrp <- (piSmr * p1barSmr * p2Smr)^-1L * weightedY2Smr
  varpEst <-
    t(correctedY2Smrp) %*%
    (covarpSmr / piSmr_mat) %*%
    correctedY2Smrp %>%
    as.numeric()
  varpEst <- varpEst +
    sum((1.0 - piSmr) * piSmr^-2L * (p1barSmr * p2Smr)^-1L *
          (1.0 - (p1barSmr * p2Smr)^-1L) * weightedY2Smr^2L)


  # m1 selection variability (q1, R1)
  correctedY2Smrq1 <- (piSmr * p1barSmr)^-1L * weightedY2Smr

  #   If the probabilities p_1k are estimations we add a term
  #   that uses the covariates used by the regression estimators of
  #   alpha1 (the parameter of the logistic model for the m1 response)
  if (correcEstimWeights)
  {
    estbPhi21 <- .estim_bPhi21(Yobs, I, pi, p1, p2, maskSmr, Z, phi)
    correctedY2Smrq1 <- correctedY2Smrq1 - Z[maskSmr, ] %*% estbPhi21
  }

  varq1Est <- sum(p1Smr * p2Smr^-1L * correctedY2Smrq1^2L)


  # m2 selection variability (q2, R2)
  correctedY2Smrq2 <- correctedY2Smrp

  #   If the probabilities p_1k are estimations we add a term
  #   that uses the covariates used by the regression estimators of
  #   alpha2 (the parameter of the logistic model for the mode-2 response)
  if (correcEstimWeights)
  {
    estbPhi22 <- .estim_bPhi22(Yobs, pi, p1, p2,
                            I & modes != "m1", maskSmr, Z, phi)
    correctedY2Smrq2 <- correctedY2Smrq2 + Z[maskSmr, ] %*% estbPhi22
  }

  varq2Est <- sum((1.0 - p2Smr) * correctedY2Smrq2^2L)


  # Y2 variability
  varY2Est <- sum(phi^2L) * sd2^2L

  varpEst + varq1Est + varq2Est + varY2Est
}

#' Calculate the true variance of the HT estimator of t_phi2 in the case
#' of known selection probabilities. Using estimated values give an estimation
#' The expected values of the m2 counterfactuals have to be given.
#' The calculation will be correct under the additional strong assumptions
#' of conditional independence between the sampling design p
#' and the mode selection mechanism q1 and q2, plus the conditional independence
#' between the mode selection mechanisms.
#' @export
var_expansion_seq_phi2 <- function(Y2exp,
                                   pi_mat,
                                   p1, p2,
                                   sd2 = 0.0,
                                   phi = rep(1.0, length(Y2exp)))
{

  if (all(phi == 0.0))
    return(0.0)

  pi <- diag(pi_mat)
  p1bar <- 1.0 - p1

  # Sampling design variability (p, S)
  covarPi <- pi2_to_covarInc(pi_mat)
  correctedY2p <- pi^-1L * phi * Y2exp
  varp <- t(correctedY2p) %*% covarPi %*% correctedY2p +
    sum(pi^-1L * (1.0 - pi) * phi^2L) * sd2^2L

  # m1 selection variability (q1, R1)
  varPhiY2 <- phi^2L * (sd2^2L + Y2exp^2L) # Variance of each phi_k y_2k
  varq1 <- sum((pi * p1bar)^-1L * p1 * varPhiY2)

  # m2 selection variability (q2, R2)
  varq2 <- sum((pi * p1bar * p2)^-1L * (1.0 - p2) * varPhiY2)

  # Y2 variability
  varY2 <- sum(phi^2L) * sd2^2L

  varp + varq1 + varq2 + varY2
}

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
