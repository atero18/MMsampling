#' @importFrom checkmate assertFlag assertSubset
#' @export
HT_mm_weights <- function(pi, probaModes, chosenModes = NULL, inv = FALSE)
{
  assertFlag(inv)

  assertProbabilityVec(pi, striclyPos = TRUE)

  n <- length(pi)

  inclusionWeights <- pi

  modeAndNRWeights <-
    data_proba_to_vec(probaModes, N = n, modes = chosenModes,
                      probVec = TRUE)

  weights <- inclusionWeights * modeAndNRWeights


  if (!inv)
    weights <- 1.0 / weights

  as.double(weights)
}


#' @importFrom checkmate assertFlag
HT_mm_weights_with_control <-
  function(pi, probaMode, control = logical(n),
           nu = NULL, nuBar = NULL, q = NULL, inv = TRUE)
{

  ## ajouter validations
  assertFlag(inv)

  n <- length(pi)

  if (any(pi == 0.0))
    stop("At least one pi probability is equal to zero")

  inclusionWeights <- pi

  existingControlSet <- any(control)
  missingDataControl <- is.null(nu) || is.null(nuBar) || is.null(q)
  if (existingControlSet && missingDataControl)
    stop("controls are present but probabilities are missing")

  else if (!existingControlSet)
  {
    nu <- q <- numeric(n)
    nuBar <- rep(1.0, n)
  }

  modeAndNRWeights <- numeric(n)

  if (existingControlSet)
    modeAndNRWeights[control] <- nu * q
  if (!existingControlSet || !all(control))
    modeAndNRWeights[!control] <- nuBar * probaMode

  weights <- inclusionWeights * modeAndNRWeights

  if (!inv)
    weigths <- 1.0 / weigths

  as.double(weights)
}

#' @export
HT_mm <- function(pi, Y, phi, delta = NULL,
                  probaModes, chosenModes = NULL,
                  initialWeights = NULL)
{
  n <- length(pi)
  weights <- HT_mm_weights(pi, probaModes, chosenModes, inv = FALSE)

  assertY(Y, N = n)

  if (!is.null(initialWeights))
  {
    assertNumericVector(initialWeights, finite = TRUE,
                        any.missing = FALSE, all.missing = FALSE,
                        len = n)

    if (any(initialWeights == 0.0))
      return("some initial Weight is equal to zero")

    weights <- initialWeights * weights
  }


  phi <- data_proba_to_vec(phi, N = n, modes = chosenModes,
                           probVec = TRUE)


  if (is.null(delta))
    delta <- numeric(n)
  else
  {
    delta <- get_value_by_mode(delta, chosenModes)
    assertBiasesMes(delta, N = n)

  }




  crossprod(phi * (Y - delta), weights)
}

#' @importFrom checkmate assertChoice
#' @export
HT_Ym <- function(pi, Y, mode, phi,
                  probaModes, chosenModes = NULL)

{
  assertChoice(mode, chosenModes)

  Ym <- Y
  Ym[chosenModes != mode] <- O.0

  HT_mm(pi, Ym, phi, delta = NULL, probaModes, chosenModes)
}

HT_Yref <- function(pi, Y, modesRef, phi, probaModes, chosenModes) NULL ## à faire

#' @importFrom checkmate testSubset assertInt testInt
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

.estim_var_Inc_mm <- function(pi2, Ytilde, delta)
{
  pi <- diag(pi2)
  Ypi <- (Ytilde - delta) / pi
  covarInc <- pi2_to_covarInc(pi2)

  sum(Ypi %*% t(Ypi) * covarInc / pi2)
}

.estim_var_mode_mm <- function(mode, pi, p_tab, Ytilde, delta, phi)
{
  probaUsed <- get_value_by_mode(p, mode)
  Ysqpip <- (Ytilde - delta)^2L / (pi^2L * probaUsed)
  p <- p[, colnames(p) != "nr", drop = FALSE]
  varPhi <- rowSums(phi^2L / p) - 1.0

  sum(Ysqpip * varPhi)

}

#' @export
estim_var_HT_mm <- function(mode, pi2, p_tab, phi, Ytilde, delta)
{
  masque <- mode[!is.na(mode) & mode != "nr"]
  p_tab <- p_tab[masque, , drop = FALSE]
  pi2 <- pi2[masque, masque, drop = FALSE]
  pi <- diag(pi2)
  Ytilde <- Ytilde[masque]
  delta <- delta[masque]
  phi <- phi[masque, , drop = FALSE]

  .estim_var_Inc_mm(pi2, Ytilde, delta) +
    .estim_var_mode_mm(mode, pi, p_tab, Ytilde, delta, phi)
}

#' @importFrom stats fitted
HT_mm_with_fitting <- function(sample,
                               imputationY = NULL,
                               estimMesBias = NULL, estimPm = NULL,
                               checkEquality = "no", checkNullityBias = "no")
{


  N <- sample$N

  maskResp <- sample$R

  if (checkEquality != "no" || checkNullityBias != "no")
    sample <- sample$copy(deep = TRUE)

  Ytilde <- sample$Ytilde()

  if (checkEquality %in% c("MCO", "MCO_agreg"))
  {
    mergeEquivalents <- checkEquality == "MCO_agreg"
    sample <-
      check_modes_equality_MCO(sample,
                               mergeEquivalents = mergeEquivalents)$plan
  }

  delta <- numeric(N)

  # Some measure biases estimators don't need
  # calculation of contrefactuals
  if (is.null(estimMesBias) || estimMesBias %in% c("G-COMP", "MCO_tots"))
    imputationY <- NULL
  else
  {
    if (imputationY == "true_values")
      Ycf <- sample$Yref()
    else if (imputationY == "MCO")
      Ycf <- MCO_Y(sample)
    else # Use of MatchIt package
      Ycf <- Y_matching(sample, imputationY)
    assertString(imputationY, null.ok = TRUE)

  }

  # Estimation of measure biases
  usedBiasedModes <- sample$used_biased_modes()

  if (is.null(estimMesBias) || estimMesBias == "true_values")
    delta <- sample$Ytilde() - sample$Yref()
  else if (estimMesBias == "CF")
    delta <- sample$Ytilde() - Ycf
  else if (estimMesBias == "G-COMP")
    delta <- estim_delta_G_comp_MCO(sample)

  # Estimation with totals
  else if (estimMesBias %in% c("MCO_tots", "MCO_cf", "MCO_tot_cf"))
  {
    typeTot <- switch(estimMesBias,
                      MCO_tots = "doubleHT",
                      MCO_cf = "estimDeltaCF",
                      MCO_tot_cf = "totYimp")

    for (biasedMode in usedBiasedModes)
    {
      masqueMode <- sample$respondents_mode(biasedMode)
      delta[masqueMode] <-
        estim_MB_MCO_tot(sample, biasedMode, typeTot = typeTot, Ycf = Ycf)[masqueMode]
    }
  }

  # Estimation with imputation
  else
  {
    if (imputationY == "true_value")
      Ycf <- sample$Yref()

    else if (imputationY == "MCO")
      Ycf <- MCO_Y(sample)

    else
      Ycf <- Y_matching(sample, imputationY)

    # Case when Y^m is replaced by its counterfactuel
    if (estimMesBias == "CF")
    {
      delta <- Ytilde - Ycf
      usedBiasedModes <- NULL
    }

    # Case when we want to check if there is negligeable biases
    else if (checkNullityBias != "no")
    {
      # MCO_Ym : Ym is maintaned if the bias is negligeable
      # MCO_CF : Ym is replaced by its counterfactual if the bias is negligeable
      if (checkNullityBias %in% c("MCO_Ym", "MCO_CF"))
      {
        replaceByCF <- checkNullityBias == "MCO_CF"
        resCheck <-
          check_nullity_bias_MCO(sample, Ycf, replaceByCF = replaceByCF)

        delta <- resCheck$delta
        usedBiasedModes <- resCheck$remainingBiasedModes

      }
    }

  }

  masqueRepRef <- sample$respondents_ref()
  delta[masqueRepRef] <- 0.0


  if (is.null(estimPm) || estimPm == "true_values")
    pHat <- sample$probaModes

  else if (estimPm == "multinomial")
  {
    pHat <- estim_response_prob_global(plan, problem, regroupModesRef = TRUE)
  }

  HT_mm(sample$pi[maskResp], Ytilde[maskResp],
        sample$phi_tab[maskResp, , drop = FALSE], delta[maskResp], pHat[maskResp, , drop = FALSE], chosenModes = sample$mode[maskResp])
}
