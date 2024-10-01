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

#' Calculate the true variance of the HT estimator of t_phi1 in the case
#' of known selection probabilities. Using estimated values give an approximate
#' variance in the case of known selection probabilities
#' @export
var_HT_seq_phi1 <- function(expY1, covarY1,
                            piMat,
                            pq1Mat,
                            phi = rep(1.0, length(expY2)))
{

  if (all(phi == 0.0))
    return(0.0)

  pi <- diag(piMat)

  p1 <- diag(pq1Mat)
  invP1Mat <- (p1 %*% t(p1))^-1L

  # Y_1 variability
  varY1 <- piMat * pq1Mat * covarY1 * invP1Mat

  # mode 1 selection variability
  prodexpY1Mat <- expY1 %*% t(expY1)
  covarq1 <- pi2_to_covarInc(pq1Mat)
  varq1 <- prodexpY1Mat * piMat * covarq1 * invP1Mat

  # sampling variability
  covarPi <- pi2_to_covarInc(piMat)
  varS <- prodexpY1Mat * covarPi

  vMat <- varY1 + varq1 + varS

  weightedPhis <- phi / pi

  sum(weightedPhis %*% t(weightedPhis) * vMat)
}


#' Calculate the true variance of the HT estimator of t_phi2
#' (or t_2 as a special case) in the case of known selection probabilities.
#' Using estimated values give an approximate variance in the case of
#'  known selection probabilities.
#' @export
var_HT_seq_phi2 <- function(expY2, covarY2,
                            piMat,
                            pq1Mat, pq2Mat,
                            phi = numeric(length(expY2)))
{

  if (all(phi == 1.0))
    return(0.0)

  N <- length(expY2)

  pi <- diag(piMat)

  p1 <- diag(pq1Mat)

  pq1BarMat <- 1.0 -
    matrix(p1, nrow = N, ncol = N, byrow = TRUE) -
    matrix(p1, nrow = N, ncol = N, byrow = FALSE) +
    pq1Mat

  p1Bar <- 1.0 - p1

  p2 <- diag(pq2Mat)

  invp1BarMat <- (p1Bar %*% t(p1Bar))^-1L
  invProbsMatSelecMat <- invp1BarMat * (p2 %*% t(p2))^-1L

  # Y_2 variability
  varY2 <- piMat *
    pq1BarMat *
    pq2Mat *
    covarY2 *
    invProbsMatSelecMat


  # mode 2 selection variability
  prodexpY2Mat <- expY2 %*% t(expY2)
  covarq2 <- pi2_to_covarInc(pq2Mat)
  varq2 <- prodexpY2Mat *
    piMat *
    pq1BarMat *
    covarq2 *
    invProbsMatSelecMat


  # mode 1 selection variability
  covarq1 <- pi2_to_covarInc(pq1Mat)
  varq1 <- prodexpY2Mat *
    piMat *
    covarq1 *
    invp1BarMat

  # sampling variability
  covarPi <- pi2_to_covarInc(piMat)
  varS <- prodexpY2Mat * covarPi

  vMat <- varY2 + varq2 + varq1 + varS

  weightedPhis <- (1.0 - phi) / pi

  sum(weightedPhis %*% t(weightedPhis) * vMat)
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




