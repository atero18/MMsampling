# @importFrom checkmate assertFlag
# HT_mm_weights_with_control <-
#   function(pi, probaMode, control = logical(n),
#            nu = NULL, nuBar = NULL, q = NULL, inv = TRUE)
# {
#
#   ## ajouter validations
#   assertFlag(inv)
#
#   n <- length(pi)
#
#   if (any(pi == 0.0))
#     stop("At least one pi probability is equal to zero")
#
#   inclusionWeights <- pi
#
#   existingControlSet <- any(control)
#   missingDataControl <- is.null(nu) || is.null(nuBar) || is.null(q)
#   if (existingControlSet && missingDataControl)
#     stop("controls are present but probabilities are missing")
#
#   else if (!existingControlSet)
#   {
#     nu <- q <- numeric(n)
#     nuBar <- rep(1.0, n)
#   }
#
#   modeAndNRWeights <- numeric(n)
#
#   if (existingControlSet)
#     modeAndNRWeights[control] <- nu * q
#   if (!existingControlSet || !all(control))
#     modeAndNRWeights[!control] <- nuBar * probaMode
#
#   weights <- inclusionWeights * modeAndNRWeights
#
#   if (!inv)
#     weigths <- 1.0 / weigths
#
#   as.double(weights)
# }

# @importFrom checkmate testSubset assertInt testInt
# HT_mm_with_control <- function(sample, delta,
#                                modesRef = NULL,
#                                p, K = nrow(p), phi = NULL, nu = NULL, q = NULL)
# {
#   masqueRepondants <- sample$R
#   masqueRepondantsControle <- masque_repondants_controle(plan)
#   masqueRepondantsNonControle <- masque_repondants_non_controle(plan)
#
#   nr <- sum(masqueRepondants)
#
#   if (is.null(nu) && any(plan$C))
#     abort("Probas for control are not gave")
#
#
#   plan <- plan[masqueRepondants, , drop = FALSE]
#
#   assertPlan(plan)
#   plan <- preprocess_plan(plan)
#
#   Ytilde <- Ytilde[masqueRepondants]
#   assertY(Ytilde, nr)
#
#   usedModes <- modes_utilises(plan)
#   usedBiasedModes <- modes_biaises_utilises(plan, modesRef)
#
#   if (length(usedBiasedModes) > 0L)
#   {
#     if (is.vector(delta))
#       deltaRepondants <- delta[masqueRepondants]
#
#     else
#     {
#       if (!testSubset(usedBiasedModes, colnames(delta)))
#         stop("At least one used biased mode doesn't have delta estimations")
#
#       delta <- delta[masqueRepondants, usedBiasedModes, drop = FALSE]
#
#       deltaRepondants <- numeric(nr)
#       inBiasModes <- which(plan$mode %in% usedBiasedModes)
#       deltaRepondants[inBiasModes] <-
#         get_value_by_mode(delta[inBiasModes, ], plan$mode)
#     }
#
#     Y <- Ytilde - deltaRepondants
#
#   }
#   else
#     Y <- Ytilde
#
#   if (any(masqueRepondantsNonControle))
#   {
#     if (is.vector(p))
#       p <- p[masqueRepondantsNonControle]
#     else
#     {
#       p <- p[masqueRepondantsNonControle, , drop = FALSE]
#       assertProbTable(p, nr)
#
#       if (!testSubset(usedModes, colnames(p)))
#         stop("Some modes are not indicated in probabilities")
#
#       p <- get_value_by_mode(p, plan$mode)
#     }
#
#     if (!all(p > 0.0))
#       stop("Some responding probabilities are equal to 0")
#   }
#
#   if (is.null(nu) && any(masqueRepondantsControle))
#     stop("Control sample is not empty and control probabilities aren't given")
#
#   if (is.null(nu))
#     nuBar <- NULL
#
#   if (!is.null(nu))
#   {
#
#     if (is.vector(nu))
#     {
#       nu <- nu[masqueRepondants]
#     }
#
#     else
#     {
#       nu <- nu[masqueRepondants, , drop = FALSE]
#       assertProbTable(nu)
#       nu <- add_nr_prob(nu, "nC")
#       nuBar <- nu$nC
#
#       if (!testSubset(usedModesControl, colnames(nu)))
#         stop("Some modes are not indicated in probabilities")
#
#       nu <- get_value_by_mode(nu, plan$mode)
#     }
#
#     if (!all(nu > 0.0))
#       stop("Some control probabilities are equal to 0")
#
#
#     if (any(masqueRepondantsControle))
#     {
#
#       if (is.vector(q))
#         q <- q[masqueRepondantsControle]
#
#       else
#       {
#         usedModesControl <- modes_utilises_controle(plan)
#         q <- q[masqueRepondantsControle, usedModesControl, drop = FALSE]
#         assertProbTable(q, probVec = FALSE)
#
#         if (!testSubset(usedModes(usedModesControl), colnames(q)))
#           stop("Some modes are not indicated in probabilities")
#
#         q <- get_value_by_mode(q, plan$mode)
#       }
#
#       if (!all(q > 0.0))
#         stop("Some control responding probabilities are equal to 0")
#     }
#   }
#
#   invWeights <- HT_mm_weights(pi, p, plan$C, nu, nuBar, q)
#
#
#   if (is.null(phi) || (testInt(phi) && phi == 0.0))
#   {
#     assertInt(K, lower = length(usedModes))
#
#     # Same weight for each answer
#     if (is.null(phi))
#       phi <- rep(1.0 / K, nr)
#
#     # More important weight for reference modes
#     else
#     {
#       phi <- rep(1.0 / (K - length(modesRef)))
#       masqueBiaises <- masque_repondants_biaises(plan, modesRef)
#       phi[masqueBiaises] <- 1.0 / length(modesRef)
#     }
#
#   }
#   else
#   {
#     if (is.vector)
#       phi <- phi[masqueRepondants]
#     else
#     {
#       phi <- phi[masqueRepondants, usedModes, drop = FALSE]
#       phi <- get_value_by_mode(phi, plan$modes)
#     }
#
#     assertProbabilityVec(phi)
#   }
#
#   Y <- Y * phi
#
#   crossprod(Y, 1.0 / invWeights)
# }


.bPhi11 <- function(expY1, I, R1, p1, Z,
                    phi = rep(1.0, length(expY1)),
                    constY = FALSE)
{

  weightedY1 <- phi * expY1


  # muPhi11 (expected value of the derivative of the total
  # with alpha_1 = alpha_1^*) depends on Y1, so if the y_1k are
  # considered as constant then we can only use Sr.
  if (constY1)
  {
    p1Sr <- p1[R1]
    weightedY1Sr <- weightedY1[R1]
    muPhi11 <- crossprod(2.0 * weightedY1Sr *
                           (1.0 - p1Sr) / (piSr * p1Sr^2L),
                         Z[R1, , drop = FALSE])
  }
  else
  {
    muPhi11 <- crossprod(2.0 * weightedY1 *
                           (1.0 - p1) / p1, Z)
  }


  # diagonal matrix of ^p_1k * (1-^p1k)
  lambda1S <- diag(p1[I] * (1.0 - p1[I]) / pi[I])
  partialW1 <- t(Z[I, , drop = FALSE]) %*%
    lambda1S %*%
    Z[I, , drop = FALSE]

  bPhi11 <- t(muPhi11) %*% solve(partialW1)

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
  partialW1 <- t(Z[I, , drop = FALSE]) %*%
    lambda1S %*%
    Z[I, , drop = FALSE]

  bPhi21 <- t(muPhi21) %*% solve(partialW1)

  return(bPhi21)
}


.bPhi22 <- function(expY2, I, R1, p1, R2, p2, Z,
                    phi = rep(1.0, length(expY1)),
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

  indicatorSm <- I & !R1
  # diagonal matrix of ^p_1k * (1-^p1k)
  lambda2Sm <- diag(p2[indicatorSm] * (1.0 - p2[indicatorSm]) /
                      (pi[indicatorSm] * p1Bar[indicatorSm]))
  partialW2 <- t(Z[indicatorSm, , drop = FALSE]) %*%
    lambda2Sm %*%
    Z[indicatorSm, , drop = FALSE]

  bPhi22 <- t(muPhi22) %*% solve(partialW2)

}

#' Calculate the true variance of the HT estimator of t_phi1 in the case
#' of known selection probabilities. Using estimated values give an approximate
#' variance in the case of known selection probabilities
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

  p1 <- diag(pq1Mat)
  invP1Mat <- (p1 %*% t(p1))^-1L

  R1 <- as.numeric(modes == biasedMode)

  covarPi <- pi2_to_covarInc(piMat)

  weightedY1 <- phi * expY1


  # Sampling variability (S)

  # If Y1 is considered as constant we can only estimate the sub-variance
  # in Sr
  if (constY1)
  {
    piSr <- pi[R1]
    weightedY1Sr <- weightedY1[R1]
    varS <- weightedY1Sr / piSr %*% t(weightedY1Sr / piSr) *
      covarPi[R1, R1] /
      (piMat[R1, R1] * pq1Mat[R1, R1])
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
      bPhi11 <- .bPhi11(expY1, I, R1, p1, Z, phi, constY)
      correctedY1Sr <- correctedY1Sr -
        crossprod(Z[R1, , drop = FALSE], bPhi11)
    }

    varq1 <- (correctedY1Sr / piSr) %*% t(correctedY1Sr / piSr) *
      covarq1[R1, R1] / pq1Mat[R1, R1]
  }
  else
  {
    correctedY1 <- weightedY1 / p1

    if (correcEstimWeights)
    {
      bPhi11 <- .bPhi11(expY1, I, R1, p1, Z, phi, constY)
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

  R1 <- as.numeric(modes == biasedMode)
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
      bPhi22 <- .bPhi22(expY2, I, R1, p1, R2, p2, Z, phi, constY)
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
      bPhi22 <- .bPhi22(expY2, I, R1, p1, R2, p2, Z, phi, constY)
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
      bPhi11 <- .bPhi11(expY1, I, R1, p1, Z, phi, constY)
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

test_absence_measure_bias <- function(expY1, expY2, I,
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
                                      covarY2 = NULL,
                                      alpha = 0.05)
{

  pi <- diag(piMat)


  R1 <- as.numeric(modes == biasedMode)
  p1 <- diag(pq1Mat)

  tPhi1 <- sum(phi[R1] * expY1[R1] / (pi[R1] * p1[R1]))

  R2 <- as.numeric(modes == refMode)
  p2 <- diag(pq2Mat)

  tPhi2 <- sum(phi[R2] * expY2[R2] / (pi[R2] * (1.0 - p1[R2]) * p2[R2]))

  delta <- tPhi1 - tPhi2

  varDelta <- var_difference_HT(expY1, expY2, I,
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
                                covarY2)

  Zstat <- delta / sqrt(varDelta)

  pValue <- 2.0 * (1.0 - pnorm(Zstat))

  uAlpha <- qnorm(1.0 - alpha / 2.0)

  if (pValue <= alpha)
    finalHyp <- "H0 (no measure bias)"
  else
    finalHyp <- "H1 (at least one measure bias)"

  results <- c(Z = Zstat, alpha = alpha,
               pValue = pValue, uAlpha = uAlpha,
               hypothesis = finalHyp)

  results
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



#' @importFrom stats fitted
#' @importFrom checkmate assertChoice
# HT_mm_with_fitting <- function(sample,
#                                imputationY = NULL,
#                                estimMesBias = NULL, estimPm = NULL,
#                                checkEquality = "no", checkNullityBias = "no")
# {
#
#
#   # Size of the population
#   N <- sample$N
#
#   # Mask of respondents
#   maskResp <- sample$R
#
#   maskRespRef <- sample$respondents_ref()
#
#   # If some values in the sample might be changed we
#   # make a copy of it
#   if (checkEquality != "no" || checkNullityBias != "no")
#     sample <- sample$copy(deep = TRUE)
#
#   # Get answers of each unit
#   # (NA when the unit haven't been selected and "nr" when the unit has been
#   # selected but didn't answer)
#   Ytilde <- sample$Ytilde()
#
#   if (checkEquality %in% c("MCO", "MCO_agreg"))
#   {
#     mergeEquivalents <- checkEquality == "MCO_agreg"
#     sample <-
#       check_modes_equality_MCO(sample,
#                                mergeEquivalents = mergeEquivalents)$plan
#   }
#
#
#   # Some measure biases estimators don't need
#   # calculation of countefactuals
#   if (is.null(estimMesBias) || estimMesBias %in% c("G-COMP", "MCO_tots"))
#     imputationY <- NULL
#   else
#   {
#     assertString(imputationY)
#
#     if (imputationY == "true_values")
#       Ycf <- sample$Yref()
#     else if (imputationY == "MCO")
#       Ycf <- MCO_Y(sample)
#     else # Use of MatchIt package
#       Ycf <- Y_matching(sample, imputationY)
#   }
#
#   # Estimation of measure biases
#   usedBiasedModes <- sample$used_biased_modes()
#
#   # Vector of measure biases
#   delta <- numeric(N)
#
#   if (is.null(estimMesBias) || estimMesBias == "true_values")
#     delta <- sample$Ytilde() - sample$Yref()
#   else if (estimMesBias == "CF")
#     delta <- sample$Ytilde() - Ycf
#   else if (estimMesBias == "G-COMP")
#     delta <- estim_delta_G_comp_MCO(sample)
#
#   # Estimation with totals
#   else if (estimMesBias %in% c("estimDeltaCF", "doubleHT", "totYcf"))
#   {
#     for (biasedMode in usedBiasedModes)
#     {
#       masqueMode <- sample$respondents_mode(biasedMode)
#       delta[masqueMode] <-
#         estim_MB_MCO_tot(sample, biasedMode,
#                          typeTot = estimMesBias, Ycf = Ycf)[masqueMode]
#     }
#   }
#
#   # Estimation with imputation
#   else
#   {
#     # Case when Y^m is replaced by its counterfactuel
#     if (estimMesBias == "CF")
#     {
#       delta <- Ytilde - Ycf
#       usedBiasedModes <- NULL
#     }
#     else if (estimMesBias == "MCO")
#     {
#       for (biasedMode in usedBiasedModes)
#       {
#         masqueMode <- sample$respondents_mode(biasedMode)
#         delta[masqueMode] <-
#           estim_MB_MCO_cf(sample, biasedMode, Ycf)[masqueMode]
#       }
#     }
#
#     # Case when we want to check if there is negligeable biases
#     else if (checkNullityBias != "no")
#     {
#       # MCO_Ym : Ym is maintaned if the bias is negligeable
#       # MCO_CF : Ym is replaced by its counterfactual if the bias is negligeable
#       if (checkNullityBias %in% c("MCO_Ym", "MCO_CF"))
#       {
#         replaceByCF <- checkNullityBias == "MCO_CF"
#         resCheck <-
#           check_nullity_bias_MCO(sample, Ycf, replaceByCF = replaceByCF)
#
#         delta <- resCheck$delta
#         usedBiasedModes <- resCheck$remainingBiasedModes
#
#       }
#     }
#
#   }
#
#   masqueRepRef <- sample$respondents_ref()
#   delta[masqueRepRef] <- 0.0
#
#
#   if (is.null(estimPm) || estimPm == "true_values")
#     pHat <- sample$probaModes
#
#   else if (estimPm == "multinomial")
#   {
#     pHat <- estim_response_prob_global(plan, problem, regroupModesRef = TRUE)
#   }
#
#   HT_mm(sample$pi[maskResp], Ytilde[maskResp],
#         sample$phi_tab[maskResp, , drop = FALSE],
#         delta[maskResp], pHat[maskResp, , drop = FALSE],
#         chosenModes = sample$mode[maskResp])
# }




