.bPhi11 <- function(Yobs, I, pi, maskSr, p1, Z,
                    phi = rep(1.0, length(Yobs)))
{

  p1Sr <- p1[maskSr]

  # Estimation of the expected value of the derivative of the total estimator
  muPhi11 <- 2.0 * crossprod(Z[maskSr, , drop = FALSE],
                       phi[maskSr] * Yobs[maskSr] *
                         (1.0 - p1Sr) / (pi[maskSr] * p1Sr^2L))



  # Estimation of the partial derivative expectation of the logistic score
  # function
  n <- sum(I)
  p1S <- p1[I]

  lambda1S <- diag(p1S * (1.0 - p1S) / pi[I], nrow = n, ncol = n)
  minusPartialW1 <-
    t(Z[I, , drop = FALSE]) %*%
    lambda1S %*%
    Z[I, , drop = FALSE]

  bPhi11 <- -solve(minusPartialW1) %*% muPhi11

  return(bPhi11)
}


.bPhi21 <- function(expY2, I, p1, R2, p2, Z,
                    phi = rep(1.0, length(expY1)),
                    constY = FALSE)
{

  weightedY2 <- phi * expY2

  p1Bar <- 1.0 - p1

  # muPhi21 (expected value of the derivative of the total
  # with alpha_1 = alpha_1^*) depends on Y2, so if the y_2k are
  # considered as constant then we can only use Smr.
  if (constY2)
  {
    p1BarSmr <- p1Bar[R2]
    p2Smr <- p2[R2]
    weightedY2Smr <- weightedY2[R2]

    muPhi21 <- crossprod(2.0 * weightedY2Smr *
                           p1Smr / (piSmr * p1BarSmr^2L * p2Smr),
                         Z[R2, , drop = FALSE])
  }
  else
  {
    muPhi21 <- crossprod(2.0 * weightedY2 * p1 / p1Bar, Z)
  }


  # diagonal matrix of ^p_1k * (1-^p1k)
  lambda1S <- diag(p1[I] * (1.0 - p1[I]) / pi[I])
  partialW1 <- -t(Z[I, , drop = FALSE]) %*%
    lambda1S %*%
    Z[I, , drop = FALSE]

  bPhi21 <- t(muPhi21) %*% solve(partialW1)

  return(bPhi21)
}


.bPhi22 <- function(Yobs, I, pi, maskSr, p1, R2, p2, Z,
                    phi = rep(1.0, length(Yobs)),
                    constY = FALSE)

{

  weightedY2 <- phi * expY2

  p1Bar <- 1.0 - p1

  # muPhi22 (expected value of the derivative of the total
  # with alpha_2 = alpha_2^*) depends on Y2, so if the y_2k are
  # considered as constant then we can only use Smr.
  if (constY2)
  {
    p1BarSmr <- p1Bar[R2]
    p2Smr <- p2[R2]
    weightedY2Smr <- weightedY2[R2]

    muPhi22 <- crossprod(2.0 * weightedY2Smr *
                           (1.0 - p2Smr) / (piSmr * p1BarSmr * p2Smr^2L),
                         Z[R2, , drop = FALSE])
  }
  else
  {
    muPhi22 <- crossprod(2.0 * weigtedY2 * (1.0 - p2) / p2, Z)
  }

  indicatorSm <- I & !maskSr
  # diagonal matrix of ^p_1k * (1-^p1k)
  lambda2Sm <- diag(p2[indicatorSm] * (1.0 - p2[indicatorSm]) /
                      (pi[indicatorSm] * p1Bar[indicatorSm]))
  partialW2 <- -t(Z[indicatorSm, , drop = FALSE]) %*%
    lambda2Sm %*%
    Z[indicatorSm, , drop = FALSE]

  bPhi22 <- t(muPhi22) %*% solve(partialW2)

}

#' Estimates the approximate variance of the HT estimator of t_phi1 with known
#' or estimated mode selection probabilities. Homoscedasticity is assumed.
#' @export
estim_appr_var_seq_phi1 <- function(Yobs,
                                    modes,
                                    I, piMat,
                                    pq1Mat,
                                    Z,
                                    phi = rep(1.0, length(Yobs)),
                                    sd1 = 0.0,
                                    correcEstimWeights = FALSE)
{
  if (all(phi == 0.0))
    return(0.0)

  maskSr <- modes == "int"

  pi <- diag(piMat)
  piSr <- pi[maskSr]
  piMatSr <- piMat[maskSr, maskSr]

  p1 <- diag(pq1Mat)
  pq1MatSr <- pq1Mat[maskSr, maskSr]

  # The y_1k are weighted with the phi_k
  weightedY1Sr <- phi[maskSr] * Yobs[maskSr]

  # Sampling variability (S)
  # There is no correction needed for probabilities estimation

  covarpSr <- pi2_to_covarInc(piMatSr)
  varSEst <-
    t(weightedY1Sr / piSr) %*%
    (covarpSr / (piMatSr * pq1MatSr)) %*%
    (weightedY1Sr / piSr) %>%
    as.numeric()

  # q1 variability (R1)
  # If we use estimated p_1k weights we have to make a correction
  correctedY1Sr <- weightedY1Sr / p1[maskSr]
  covarq1Sr <- pi2_to_covarInc(pq1MatSr)

  #   If the probabilities p_1k are estimations we add a term
  #   that uses the covariates used by the regression estimators
  if (correcEstimWeights)
  {
    bPhi11 <- .bPhi11(Yobs, I, pi, maskSr, p1, Z, phi)
    correctedY1Sr <- correctedY1Sr - Z[maskSr, ] %*% bPhi11
  }

  varq1Est <- t(correctedY1Sr / piSr) %*%
    (covarq1Sr / pq1MatSr) %*%
    (correctedY1Sr / piSr) %>%
    as.numeric()

  # Y1 variability
  # We suppose we have an estimator of the variance
  # of the Y1
  varY1Est <- sum(phi^2L) * sd1^2L

  varSEst + varq1Est + varY1Est
}

#' Calculate the true variance of the HT estimator of t_phi1 in the case
#' of known selection probabilities. Using estimated values give an approximate
#' variance.
#' @export
var_HT_seq_phi1 <- function(expY1,
                            I,
                            piMat,
                            pq1Mat,
                            Z,
                            biasedMode,
                            modes,
                            phi = rep(1.0, length(expY1)),
                            correcEstimWeights = FALSE,
                            constY1 = FALSE,
                            covarY1 = NULL)
{

  if (all(phi == 0.0))
    return(0.0)

  pi <- diag(piMat)
  piMatSr <- piMat[maskSr, maskSr]

  p1 <- diag(pq1Mat)
  invP1Mat <- (p1 %*% t(p1))^-1L

  maskSr <- as.numeric(modes == biasedMode)

  covarPi <- pi2_to_covarInc(piMat)

  pq1MatSr <- pq1Mat[maskSr, maskSr]

  weightedY1 <- phi * expY1


  ## constY1 : important sachant que l'on a une fonction qui estime la variance approchée?
  # Sampling variability (S)

  # If Y1 is considered as constant we can only estimate the sub-variance
  # in Sr
  if (constY1)
  {
    piSr <- pi[maskSr]
    weightedY1Sr <- weightedY1[maskSr]
    varS <- weightedY1Sr / piSr %*% t(weightedY1Sr / piSr) *
      covarPi[maskSr, maskSr] /
      (piMatSr * pq1MatSr)
  }
  # Otherwise we can do it on U
  else
    varS <- weightedY1 / pi %*% t(weightedY1 / pi) * covarPi

  # q1 variability (R1)
  # If we use estimated p_1k weights we have to do a correction.

  covarq1 <- pi2_to_covarInc(pq1Mat)

  if (constY1)
  {
    correctedY1Sr <- weightedY1Sr / p1Sr

    if (correcEstimWeights)
    {
      bPhi11 <- .bPhi11(expY1, I, maskSr, p1, Z, phi, constY)
      correctedY1Sr <- correctedY1Sr -
        crossprod(Z[maskSr, , drop = FALSE], bPhi11)
    }

    varq1 <- (correctedY1Sr / piSr) %*% t(correctedY1Sr / piSr) *
      covarq1[maskSr, maskSr] / pq1MatSr
  }
  else
  {
    correctedY1 <- weightedY1 / p1

    if (correcEstimWeights)
    {
      bPhi11 <- .bPhi11(expY1, I, maskSr, p1, Z, phi, constY)
      correctedY1 <- correctedY1Sr - crossprod(Z, bPhi11)
    }

    varq1 <- (correctedY1 / pi) %*% t(correctedY1 / pi) *
      piMat * covarq1
  }

  var <- sum(varS + varq1)

  # If the y_1k are constants the complete
  # variance has been calculated
  # Otherwise is missing the part of variability from y_1k
  if (constY1)
    return(var)

  # Y1 variability

  # Can be calculated on the entire population U
  variables <- phi / (pi * p1)

  varY1 <- variables %*% t(variables) * piMat * pq1Mat * covarY1


  var + sum(varY1)
}


#' Calculate the true variance of the HT estimator of t_phi2
#' (or t_2 as a special case) in the case of known selection probabilities.
#' Using estimated values give an approximate variance in the case of
#'  known selection probabilities.
#' @export
var_HT_seq_phi2 <- function(expY2, I,
                            piMat,
                            pq1Mat, pq2Mat,
                            Z,
                            biasedMode,
                            refMode,
                            modes,
                            phi = rep(1.0, length(expY2)),
                            correcEstimWeights = FALSE,
                            constY2 = FALSE,
                            covarY2 = NULL)
{

  if (all(phi == 0.0))
    return(0.0)


  pi <- diag(piMat)

  p1 <- diag(pq1Mat)

  pq1BarMat <- 1.0 -
    matrix(p1, nrow = N, ncol = N, byrow = TRUE) -
    matrix(p1, nrow = N, ncol = N, byrow = FALSE) +
    pq1Mat

  p1Bar <- 1.0 - p1

  p2 <- diag(pq2Mat)

  maskSr <- as.numeric(modes == biasedMode)
  R2 <- as.numeric(modes == refMode)

  weightedY2 <- phi * expY2


  # Sampling variability (S)

  # If Y2 is considered as constant we can only estimate the sub-variance
  # in Sr
  if (constY2)
  {
    piSmr <- pi[R2]
    weightedY2Smr <- weightedY2[R2]
    piMatSmr <- piMat[R2, R2]
    pq1BarMatSmr <- pq1BarMat[R2, R2]
    pq2MatSmr <- pq2Mat[R2, R2]

    varS <- weightedY2Smr / piSmr %*% t(weightedY2Smr / piSmr) *
      covarPi[R2, R2] /
      (piMatSmr * pq1BarMatSmr * pq2MatSmr)
  }
  # Otherwise we can do it on U
  else
  {
    varS <- weightedY1 / pi %*% t(weightedY1 / pi) * covarPi
  }


  # q1 variability (R1)
  # If we use estimated p_1k weights we have to do a correction.

  covarq1 <- pi2_to_covarInc(pq1Mat)

  if (constY2)
  {
    correctedY2Smr <- weightedY2Smr / p1BarSmr


    if (correcEstimWeights)
    {
      bPhi21 <- .bPhi21(expY2, I, p1, R2, p2, Z, phi, constY)
      correctedY2Smr <- correctedY2Smr +
        crossprod(Z[R2, , drop = FALSE], bPhi21)
    }

    varq1 <- (correctedY2Smr / piSmr) %*% t(correctedY2Smr / piSr) *
      covarq1[R2, R2] /
      pq1Mat[R2, R2]
  }
  else
  {
    correctedY2 <- weightedY2 / p1Bar

    if (correcEstimWeights)
    {
      bPhi21 <- .bPhi21(expY2, I, p1, R2, p2, Z, phi, constY)
      correctedY2 <- correctedY2 + crossprod(Z, bPhi21)
    }

    varq1 <- (correctedY2 / pi) %*% t(correctedY2 / pi) *
      piMat * covarq1
  }

  # q2 variability (R2)
  # If we use estimated p_2k weights we have to do a correction.


  covarq2 <- pi2_to_covarInc(pq2Mat)

  if (constY2)
  {
    correctedY2Smr <- weightedY2Smr / p2Smr

    if (correcEstimWeights)
    {
      bPhi22 <- .bPhi22(expY2, I, maskSr, p1, R2, p2, Z, phi, constY)
      correctedY2Smr <- correctedY2Smr -
        crossprod(Z[R2, , drop = FALSE], bPhi22)
    }

    varq2 <- (correctedY2Smr / (piSmr * p1BarSmr)) %*%
      t(correctedY2Smr / (piSmr * p1BarSmr)) *
      pq1BarMat[R2, R2] *
      covarq2[R2, R2] /
      pq1Mat[R2, R2]
  }
  else
  {
    correctedY2 <- weightedY2 / p2

    if (correcEstimWeights)
    {
      bPhi22 <- .bPhi22(expY2, I, maskSr, p1, R2, p2, Z, phi, constY)
      correctedY2 <- correctedY2 - crossprod(Z, bPhi22)
    }

    varq2 <- (correctedY2 / (pi * p1Bar)) %*% t(correctedY2 / (pi * p1Bar)) *
      piMat * pq1BarMat * covarq2
  }

  var <- sum(varS + varq1 + varq2)

  # If the y_2k are constants the complete
  # variance has been calculated
  # Otherwise is missing the part of variability from y_2k
  if (constY2)
    return(var)

  # Y2 variability

  variables <- phi / (pi * p1Bar * p2)

  varY2 <- variables %*% t(variables) * piMat * pq1Mat * pq2Mat * covarY1


  var + sum(varY2)

}

covar_difference_HT <- function(expY1, expY2, I,
                                piMat,
                                pq1Mat, pq2Mat,
                                Z,
                                biasedMode,
                                refMode,
                                modes,
                                phi = rep(1.0, length(expY1)),
                                correcEstimWeights = FALSE,
                                constY = FALSE,
                                covarY1 = NULL,
                                covarY2 = NULL)
{
  if (all(phi == 0.0))
    return(0.0)

  weightedY1 <- phi * expY1
  weightedY2 <- phi * expY2

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
      bPhi11 <- .bPhi11(expY1, I, maskSr, p1, Z, phi, constY)
      correctedY1 <- correctedY1Sr - crossprod(Z, bPhi11)
    }

    correctedY2 <- weightedY2 / p1Bar

    if (correcEstimWeights)
    {
      bPhi21 <- .bPhi21(expY2, I, p1, R2, p2, Z, phi, constY)
      correctedY2 <- correctedY2 + crossprod(Z, bPhi21)
    }

    varq1 <- -(correctedY1 / pi) %*% t(correctedY2 / pi) * piMat * covarq1
  }


  sum(varS + varq1)
}

var_difference_HT <- function(expY1, expY2, I,
                              piMat,
                              pq1Mat, pq2Mat,
                              Z,
                              biasedMode,
                              refMode,
                              modes,
                              phi = rep(1.0, length(expY2)),
                              correcEstimWeights = FALSE,
                              constY = FALSE,
                              covarY1 = NULL,
                              covarY2 = NULL)
{
  var_HT_seq_phi1(expY1, I, piMat, pq1Mat, Z, biasedMode,
                  modes, phi, correcEstimWeights, constY1, covarY1) -
    2.0 * covar_difference_HT(expY1, expY2, I,
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
    var_HT_seq_phi2(expY2, I, piMat, pq1Mat, pq2Mat, Z, biasedMode,
                    refMode, modes, phi, correcEstimWeights,
                    constY2, covarY2)
}

#' @export
covar_HT_seq_phi1_phi2 <- function(expY1, expY2,
                                   piMat,
                                   pq1Mat, pq2Mat,
                                   phi = numeric(length(expY2)))
{

  if (all(phi == 0.0) || all(phi == 1.0))
    return(0.0)

  pi <- diag(piMat)
  p1 <- diag(pq1Mat)

  prodexpY12Mat <- expY1 %*% t(expY2)

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
                             expY1, expY2,
                             covarY1, covarY2,
                             piMat,
                             pq1Mat, pq2Mat,
                             phi = numeric(length(expY2)),
                             subResults = FALSE)
{

  return(NA_real_)
  N <- length(expY2)

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

  prodexpY1Mat <- expY1 %*% t(expY1)

  prodexpY2Mat <- expY2 %*% t(expY2)

  expDeltas <- expY1 - expY2

  phiBar <- 1.0 - phi

  varPhi1 <- var_HT_seq_phi1(expY1, covarY1, piMat, pq1Mat, phi)
  varPhi2 <- var_HT_seq_phi2(expY2, covarY2, piMat, pq1Mat, pq2Mat, phi)
  covarPhi12 <- covar_HT_seq_phi1_phi2(expY1, expY2, piMat, pq1Mat, pq2Mat, phi)

  # Variance of the total MB estimator
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY2 <- piMat * pq1BarMat * pq2Mat * covarY2 * invProbsMatSelecMat

    varY1 <- piMat * pq1Mat * covarY1 * invp1Mat

    varq2 <- prodexpY2Mat * piMat * pq1BarMat * covarq2 * invProbsMatSelecMat

    sumY1Y2 <- expY1 / p1 + expY2 / p1Bar
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

    weightedY1 <- expY1 / p1
    varq1 <- weightedY1 %*% t(weightedY1 + expY2 / p1Bar) * piMat * covarq1

    varS <- expY1 %*% t(expDeltas) * covarPi

    v1Delta <- varY1 + varq1 + varS

    #b <- t(phi) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) / pi

    covarPhi1Delta <- sum((phi / pi) %*% t(b) * v1Delta)
  }

  # Covariance between ^t_phi2 and ^t_phiDelta
  if (modeTotBiased == "HT" && modeTotRef == "HT" && calculTotal == "population")
  {
    varY2 <- piMat * pq1BarMat * pq2Mat * covarY2 * invProbsMatSelecMat

    varq2 <- prodexpY2Mat * piMat * pq1BarMat * covarq2 * invProbsMatSelecMat

    weightedY2 <- expY2 / p1Bar
    varq1 <- weightedY2 %*% t(expY1 / p1 + weightedY2) * piMat * covarq1

    varS <- -expY2 %*% t(expDeltas) * covarPi

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
  {
    return(varEstim)
  }
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
