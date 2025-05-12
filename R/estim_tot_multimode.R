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




test_absence_measure_bias <- function(Y1exp, Y2exp, I,
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
                                      covarY2 = NULL,
                                      alpha = 0.05)
{

  pi <- diag(pi_mat)


  maskInt <- as.numeric(modes == biasedMode)
  p1 <- diag(p1_mat)

  tPhi1 <- sum(phi[maskInt] * Y1exp[maskInt] / (pi[maskInt] * p1[maskInt]))

  R2 <- as.numeric(modes == refMode)
  p2 <- diag(p2_mat)

  tPhi2 <- sum(phi[R2] * Y2exp[R2] / (pi[R2] * (1.0 - p1[R2]) * p2[R2]))

  delta <- tPhi1 - tPhi2

  varDelta <- var_difference_HT(Y1exp, Y2exp, I,
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




