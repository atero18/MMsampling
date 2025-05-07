#' Return the Fisher Information Matrix ("FIM") of a logistic model
#' @param prob probability for each unit to be sampled, considered known
#'  (numeric vector).
#' @param Z covariates matrix (numeric matrix).
#' @param maskSubset logical vector indicating which unit will be use in the
#' model (logical vector).
Fisher_Information_Matrix <- function(prob, Z, maskSubset = !logical(nrow(Z)))
{
  probSubset <- prob[maskSubset]
  ZSubset <- Z[maskSubset, , drop = FALSE]

  crossprod(ZSubset, probSubset * (1.0 - probSubset) * ZSubset)
}


.estim_b1 <- function(estV, I, p1, Z)
{
  solve(Fisher_Information_Matrix(p1, Z, I)) %*% estV
}

.estim_bPhi11 <- function(Yobs, I, pi, p1, maskSr, Z,
                    phi = rep(1.0, length(Yobs)))
{
  p1Sr <- p1[maskSr]

  # Estimation of the expected value of the derivative of the total estimator
  estVPhi11 <- -crossprod(Z[maskSr, , drop = FALSE],
                          (pi[maskSr] * p1Sr)^-1L *
                            phi[maskSr] * Yobs[maskSr] * (1.0 - p1Sr))


  .estim_b1(estVPhi11, I, p1, Z)
}


.estim_bPhi21 <- function(Yobs, I, pi, p1, p2, maskSmr, Z,
                          phi = rep(1.0, length(Y1exp)))
{

  p1Smr <- p1[maskSmr]

  estVPhi21 <- crossprod(Z[maskSmr, , drop = FALSE],
                         (pi[maskSmr] * (1.0 - p1Smr) * p2[maskSmr])^-1L *
                           phi[maskSmr] * Yobs[maskSmr] * p1Smr)



  .estim_b1(estVPhi21, I, p1, Z)
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
                                    I, piMat,
                                    p1,
                                    Z,
                                    phi = rep(1.0, length(Yobs)),
                                    estSD1 = 0.0,
                                    correcEstimWeights = FALSE,
                                    ...)
{

  if (all(phi == 0.0))
    return(0.0)

  args <- list(...)

  maskSr <- modes == "m1"

  pi <- diag(piMat)
  piSr <- pi[maskSr]
  p1Sr <- p1[maskSr]

  # The y_1k are weighted with the phi_k

  weightedY1Sr <- phi[maskSr] * Yobs[maskSr]

  # Sampling variability (S)
  # There is no correction needed for probabilities estimation

  if ("piMatSr" %in% names(args))
    piMatSr <- args[["piMatSr"]]
  else
    piMatSr <- piMat[maskSr, maskSr]

  if ("covarpSr" %in% names(args))
    covarpSr <- args[["covarpSr"]]
  else
    covarpSr <- pi2_to_covarInc(piMatSr)

  correctedY1Srp <- (piSr * p1Sr)^-1L * weightedY1Sr
  varSEst <-
    t(correctedY1Srp) %*%
    (covarpSr / piMatSr) %*%
    correctedY1Srp %>%
    as.numeric()
  varSEst <- varSEst +
    sum((1.0 - piSr) * piSr^-2L * p1Sr^-1L * (1.0 - p1Sr^-1L) *
          weightedY1Sr^2L)

  # q1 variability (R1)
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
  # We suppose the potential outcomes homoscedastic and
  # we have an estimator of the variance of the Y1
  varY1Est <- sum(phi^2L) * estSD1^2L

  varSEst + varq1Est + varY1Est
}

#' Calculate the true variance of the HT estimator of t_phi1 in the case
#' of known selection probabilities. Using estimated values give an estimation
#' The expected values of the m1 counterfactuals have to be given.
#' The calculation will be correct under the additional strong assumption
#' of conditional independence between the sampling design p
#' and the mode selection mechanism q1.
#' @export
var_expansion_seq_phi1 <- function(Y1exp,
                                   piMat,
                                   p1,
                                   phi = rep(1.0, length(Y1exp)),
                                   sd1 = 0.0)
{

  if (all(phi == 0.0))
    return(0.0)

  pi <- diag(piMat)

  # Sampling variability (S)
  covarPi <- pi2_to_covarInc(piMat)
  correctedY1p <- pi^-1L * phi * Y1exp
  varS <- t(correctedY1p) %*% covarPi %*% correctedY1p +
    sum(pi^-1L * (1.0 - pi) * phi^2L) * sd1^2L

  # q1 variability (R1)
  #  With the independence between p and q1
  varq1 <- sum((pi * p1)^-1L * (1.0 - p1) * phi^2L * (sd1^2L + Y1exp^2L))


  # Y1 variability
  varY1 <- sum(phi^2L) * sd1^2L

  varS + varq1 + varY1
}

#' Estimates the approximate variance of the HT estimator of t_phi2 with known
#' or estimated mode selection probabilities. Homoscedasticity and independence
#' in the m1 and m2 mode selection mechanisms assumed.
#' @export
estim_appr_var_seq_phi2 <- function(Yobs,
                                    modes,
                                    I, piMat,
                                    p1,
                                    p2,
                                    Z,
                                    phi = rep(1.0, length(Yobs)),
                                    estSD2 = 0.0,
                                    correcEstimWeights = FALSE,
                                    ...)
{
  if (all(phi == 0.0))
    return(0.0)

  args <- list(...)

  maskSmr <- modes == "m2"

  pi <- diag(piMat)
  piSmr <- pi[maskSmr]

  piMatSmr <- piMat[maskSmr, maskSmr]

  p1Smr <- p1[maskSmr]
  p1Smrbar <- 1.0 - p1Smr

  p2Smr <- p2[maskSmr]

  # The y_2k are weighted with the phi_k
  weightedY2Smr <- phi[maskSmr] * Yobs[maskSmr]

  # Sampling variability (S)
  if ("covarpSmr" %in% names(args))
    covarpSmr <- args[["covarpSmr"]]
  else
    covarpSmr <- pi2_to_covarInc(piMatSmr)

  correctedY2Smrp <- (piSmr * p1Smrbar * p2Smr)^-1L * weightedY2Smr
  varSEst <-
    t(correctedY2Smrp) %*%
    (covarpSmr / piMatSmr) %*%
    correctedY2Smrp %>%
    as.numeric()
  varSEst <- varSEst +
    sum((1.0 - piSmr) * piSmr^-2L * (p1Smrbar * p2Smr)^-1L *
          (1.0 - (p1Smrbar * p2Smr)^-1L) * weightedY2Smr^2L)


  # q1 variability (R1)
  correctedY2Smrq1 <- (piSmr * p1Smrbar)^-1L * weightedY2Smr

  #   If the probabilities p_1k are estimations we add a term
  #   that uses the covariates used by the regression estimators of
  #   alpha1 (the parameter of the logistic model for the m1 response)
  if (correcEstimWeights)
  {
    estbPhi21 <- .estim_bPhi21(Yobs, I, pi, p1, p2, maskSmr, Z, phi)
    correctedY2Smrq1 <- correctedY2Smrq1 - Z[maskSmr, ] %*% estbPhi21
  }

  varq1Est <- sum(p1Smr * p2Smr^-1L * correctedY2Smrq1^2L)


  # q2 variability (R2)
  correctedY2Smrq2 <- correctedY2Smrp

  #   If the probabilities p_1k are estimations we add a term
  #   that uses the covariates used by the regression estimators of
  #   alpha2 (the parameter of the logistic model for the mode-2 response)
  if (correcEstimWeights)
  {
    bPhi22 <- .estim_bPhi22(Yobs, pi, p1, p2,
                            I & modes != "m1", maskSmr, Z, phi)
    correctedY2Smrq2 <- correctedY2Smrq2 + Z[maskSmr, ] %*% bPhi22
  }

  varq2Est <- sum((1.0 - p2Smr) * correctedY2Smrq2^2L)


  # Y2 variability
  # We suppose we have an estimator of the variance
  # of the Y2
  varY2Est <- sum(phi^2L) * estSD2^2L

  varSEst + varq1Est + varq2Est + varY2Est
}

#' Calculate the true variance of the HT estimator of t_phi2 in the case
#' of known selection probabilities. Using estimated values give an estimation
#' The expected values of the m2 counterfactuals have to be given.
#' The calculation will be correct under the additional strong assumptions
#' of conditional independence between the sampling design p
#' and the mode selection mechanism q1 and q2, plus the conditional independence
#' between the mode selection mechanisms.
#' @export
var_expansion_seq_phi2 <- function(Y2exp, I,
                                   piMat,
                                   p1, p2,
                                   Z,
                                   phi = rep(1.0, length(Y2exp)),
                                   sd2 = 0.0)
{

  if (all(phi == 0.0))
    return(0.0)

  pi <- diag(piMat)
  p1Bar <- 1.0 - p1

  # Sampling variability (S)
  covarPi <- pi2_to_covarInc(piMat)
  correctedY2p <- pi^-1L * phi * Y2exp
  varS <- t(correctedY2p) %*% covarPi %*% correctedY2p +
    sum(pi^-1L * (1.0 - pi) * phi^2L) * sd2^2L

  # q1 variability (R1)
  varPhiY2 <- phi^2L * (sd2^2L + Y2exp^2L) # Variance of each phi_k y_2k
  varq1 <- sum((pi * p1Bar)^-1L * p1 * varPhiY2)

  # q2 variability (R2)
  varq2 <- sum((pi * p1Bar * p2)^-1L * (1.0 - p2) * varPhiY2)

  # Y2 variability
  varY2 <- sum(phi^2L) * sd2^2L

  varS + varq1 + varq2 + varY2
}

covar_difference_HT <- function(Y1exp, Y2exp, I,
                                piMat,
                                pq1Mat, pq2Mat,
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

  # Sampling variability (S)
  covarPi <- pi2_to_covarInc(piMat)
  if (constY)
  {
    ## TO DO
  }
  else
  {
    varS <- weightedY1 / pi %*% t(weightedY2 / pi) * covarPi
  }

  # q1 variability (R1)

  covarq1 <- pi2_to_covarInc(pq1Mat)
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

    correctedY2 <- weightedY2 / p1Bar

    if (correcEstimWeights)
    {
      estbPhi21 <- .estim_bPhi21(Y2exp, I, p1, R2, p2, Z, phi, constY)
      correctedY2 <- correctedY2 + crossprod(Z, estbPhi21)
    }

    varq1 <- -(correctedY1 / pi) %*% t(correctedY2 / pi) * piMat * covarq1
  }


  sum(varS + varq1)
}

var_difference_HT <- function(Y1exp, Y2exp, I,
                              piMat,
                              pq1Mat, pq2Mat,
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
  var_expansion_seq_phi1(Y1exp, I, piMat, pq1Mat, Z, biasedMode,
                  modes, phi, correcEstimWeights, constY1, covarY1) -
    2.0 * covar_difference_HT(Y1exp, Y2exp, I,
                              piMat,
                              pq1Mat, pq2Mat,
                              Z,
                              biasedMode,
                              refMode,
                              modes,
                              phi,
                              correcEstimWeights,
                              constY,
                              covarY1,
                              covarY2) +
    var_expansion_seq_phi2(Y2exp, I, piMat, pq1Mat, pq2Mat, Z, biasedMode,
                    refMode, modes, phi, correcEstimWeights,
                    constY2, covarY2)
}

#' @export
covar_HT_seq_phi1_phi2 <- function(Y1exp, Y2exp,
                                   piMat,
                                   pq1Mat, pq2Mat,
                                   phi = numeric(length(Y2exp)))
{

  if (all(phi == 0.0) || all(phi == 1.0))
    return(0.0)

  pi <- diag(piMat)
  p1 <- diag(pq1Mat)

  prodexpY12Mat <- Y1exp %*% t(Y2exp)

  # mode 1 selection variability
  covarq1 <- pi2_to_covarInc(pq1Mat)
  varq1 <- -prodexpY12Mat * piMat * covarq1 / (p1 %*% t(1.0 - p1))

  # sampling variability
  covarPi <- pi2_to_covarInc(piMat)
  varS <- prodexpY12Mat * covarPi

  vMat <- varq1 + varS

  sum((phi / pi) %*% t((1.0 - phi) / pi) * vMat)
}


#' @export
var_estim_tot_BM <- function(modeTotBiased = "HT", modeTotRef = "HT",
                             calculTotal = "population",
                             Y1exp, Y2exp,
                             covarY1, covarY2,
                             piMat,
                             pq1Mat, pq2Mat,
                             phi = numeric(length(Y2exp)),
                             subResults = FALSE)
{

  return(NA_real_)
  N <- length(Y2exp)

  pi <- diag(piMat)
  covarPi <- pi2_to_covarInc(piMat)

  p1 <- diag(pq1Mat)
  p1Bar <- 1.0 - p1
  pq1BarMat <- 1.0 -
    matrix(p1, nrow = N, ncol = N, byrow = TRUE) -
    matrix(p1, nrow = N, ncol = N, byrow = FALSE) +
    pq1Mat
  covarq1 <- pi2_to_covarInc(pq1Mat)
  p1Mat <- p1Bar %*% t(p1Bar)
  invp1Mat <- p1Mat^-1L
  invp1BarMat <- (p1Bar %*% t(p1Bar))^-1L

  p2 <- diag(pq2Mat)
  covarq2 <- pi2_to_covarInc(pq2Mat)


  invProbsMatSelecMat <- invp1BarMat * (p2 %*% t(p2))^-1L

  prodexpY1Mat <- Y1exp %*% t(Y1exp)

  prodexpY2Mat <- Y2exp %*% t(Y2exp)

  expDeltas <- Y1exp - Y2exp

  phiBar <- 1.0 - phi

  varPhi1 <- var_expansion_seq_phi1(Y1exp, covarY1, piMat, pq1Mat, phi)
  varPhi2 <- var_expansion_seq_phi2(Y2exp, covarY2, piMat, pq1Mat, pq2Mat, phi)
  covarPhi12 <- covar_HT_seq_phi1_phi2(Y1exp, Y2exp, piMat, pq1Mat, pq2Mat, phi)

  # Variance of the total MB estimator
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY2 <- piMat * pq1BarMat * pq2Mat * covarY2 * invProbsMatSelecMat

    varY1 <- piMat * pq1Mat * covarY1 * invp1Mat

    varq2 <- prodexpY2Mat * piMat * pq1BarMat * covarq2 * invProbsMatSelecMat

    sumY1Y2 <- Y1exp / p1 + Y2exp / p1Bar
    varq1 <- sumY1Y2 %*% t(sumY1Y2) * piMat * covarq1

    varS <- expDeltas %*% t(expDeltas) * covarPi

    vDelta <- varY2 + varY1 + varq2 + varq1 + varS

    b <- t(t(phi) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) / pi)

    varPhiDelta <- sum(b %*% t(b) * vDelta)
  }

  # Covariance between ^t_phi1 and ^t_phiDelta
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY1 <- piMat * pq1Mat * covarY1 * invp1Mat

    weightedY1 <- Y1exp / p1
    varq1 <- weightedY1 %*% t(weightedY1 + Y2exp / p1Bar) * piMat * covarq1

    varS <- Y1exp %*% t(expDeltas) * covarPi

    v1Delta <- varY1 + varq1 + varS

    #b <- t(phi) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) / pi

    covarPhi1Delta <- sum((phi / pi) %*% t(b) * v1Delta)
  }

  # Covariance between ^t_phi2 and ^t_phiDelta
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY2 <- piMat * pq1BarMat * pq2Mat * covarY2 * invProbsMatSelecMat

    varq2 <- prodexpY2Mat * piMat * pq1BarMat * covarq2 * invProbsMatSelecMat

    weightedY2 <- Y2exp / p1Bar
    varq1 <- weightedY2 %*% t(Y1exp / p1 + weightedY2) * piMat * covarq1

    varS <- -Y2exp %*% t(expDeltas) * covarPi

    v2Delta <- varY2 + varq2 + varq1 + varS

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
