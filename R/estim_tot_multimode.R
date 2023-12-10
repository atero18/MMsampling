#' @importFrom checkmate assertFlag assertSubset
#' @export
HT_mm_weights <- function(pi, probaModes, chosenModes = NULL, inv = TRUE)
{
  assertFlag(inv)

  assertProbabilityVec(pi, striclyPos = TRUE)

  n <- length(pi)

  inclusionWeights <- pi

  modeAndNRWeights <-
    data_proba_to_vec(probaModes, N = n, modes = chosenModes,
                      probVec = TRUE, addNr = TRUE,
                      strictlyPos = TRUE, striclyUnsure = FALSE)

  weights <- inclusionWeights * modeAndNRWeights


  if (inv)
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
HT_mm <- function(pi, Y, phi, delta, probaModes, chosenModes = NULL,
                  initialWeights = NULL)
{
  n <- length(pi)
  weights <- HT_mm_weights(pi, probaModes, chosenModes, inv = TRUE)

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
                           probVec = TRUE, addNr = FALSE,
                           strictlyPos = FALSE, striclyUnsure = FALSE)

  assertBiasesMes(delta, N = N)
  delta <- get_value_by_mode(delta, chosenModes)



  crossprod(phi * (Y - delta), weights)
}

#' @importFrom checkmate testSubset assertInt testInt
#' @export
HT_mm_with_control <- function(sample, delta,
                               modesRef = NULL,
                               p, K = nrow(p), phi = NULL, nu = NULL, q = NULL)
{
  masqueRepondants <- sample$R
  masqueRepondantsControle <- masque_repondants_controle(plan)
  masqueRepondantsNonControle <- masque_repondants_non_controle(plan)

  nr <- sum(masqueRepondants)

  if (is.null(nu) && any(plan$C))
    stop("Probas for control are not gave")


  plan <- plan[masqueRepondants, , drop = FALSE]

  assertPlan(plan)
  plan <- preprocess_plan(plan)

  Ytilde <- Ytilde[masqueRepondants]
  assertY(Ytilde, nr)

  usedModes <- modes_utilises(plan)
  usedBiasedModes <- modes_biaises_utilises(plan, modesRef)

  if (length(usedBiasedModes) > 0L)
  {
    if (is.vector(delta))
      deltaRepondants <- delta[masqueRepondants]

    else
    {
      if (!testSubset(usedBiasedModes, colnames(delta)))
        stop("At least one used biased mode doesn't have delta estimations")

      delta <- delta[masqueRepondants, usedBiasedModes, drop = FALSE]

      deltaRepondants <- numeric(nr)
      inBiasModes <- which(plan$mode %in% usedBiasedModes)
      deltaRepondants[inBiasModes] <-
        get_value_by_mode(delta[inBiasModes, ], plan$mode)
    }

    Y <- Ytilde - deltaRepondants

  }
  else
    Y <- Ytilde

  if (any(masqueRepondantsNonControle))
  {
    if (is.vector(p))
      p <- p[masqueRepondantsNonControle]
    else
    {
      p <- p[masqueRepondantsNonControle, , drop = FALSE]
      assertProbTable(p, nr)

      if (!testSubset(usedModes, colnames(p)))
        stop("Some modes are not indicated in probabilities")

      p <- get_value_by_mode(p, plan$mode)
    }

    if (!all(p > 0.0))
      stop("Some responding probabilities are equal to 0")
  }

  if (is.null(nu) && any(masqueRepondantsControle))
    stop("Control sample is not empty and control probabilities aren't given")

  if (is.null(nu))
    nuBar <- NULL

  if (!is.null(nu))
  {

    if (is.vector(nu))
    {
      nu <- nu[masqueRepondants]
    }

    else
    {
      nu <- nu[masqueRepondants, , drop = FALSE]
      assertProbTable(nu)
      nu <- add_nr_prob(nu, "nC")
      nuBar <- nu$nC

      if (!testSubset(usedModesControl, colnames(nu)))
        stop("Some modes are not indicated in probabilities")

      nu <- get_value_by_mode(nu, plan$mode)
    }

    if (!all(nu > 0.0))
      stop("Some control probabilities are equal to 0")


    if (any(masqueRepondantsControle))
    {

      if (is.vector(q))
        q <- q[masqueRepondantsControle]

      else
      {
        usedModesControl <- modes_utilises_controle(plan)
        q <- q[masqueRepondantsControle, usedModesControl, drop = FALSE]
        assertProbTable(q, probVec = FALSE)

        if (!testSubset(usedModes(usedModesControl), colnames(q)))
          stop("Some modes are not indicated in probabilities")

        q <- get_value_by_mode(q, plan$mode)
      }

      if (!all(q > 0.0))
        stop("Some control responding probabilities are equal to 0")
    }
  }

  invWeights <- HT_mm_weights(pi, p, plan$C, nu, nuBar, q)


  if (is.null(phi) || (testInt(phi) && phi == 0.0))
  {
    assertInt(K, lower = length(usedModes))

    # Same weight for each answer
    if (is.null(phi))
      phi <- rep(1.0 / K, nr)

    # More important weight for reference modes
    else
    {
      phi <- rep(1.0 / (K - length(modesRef)))
      masqueBiaises <- masque_repondants_biaises(plan, modesRef)
      phi[masqueBiaises] <- 1.0 / length(modesRef)
    }

  }
  else
  {
    if (is.vector)
      phi <- phi[masqueRepondants]
    else
    {
      phi <- phi[masqueRepondants, usedModes, drop = FALSE]
      phi <- get_value_by_mode(phi, plan$modes)
    }

    assertProbabilityVec(phi)
  }

  Y <- Y * phi

  crossprod(Y, 1.0 / invWeights)
}

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
HT_mm_with_fitting <- function(plan, problem,
                               imputationY = NULL,
                               estimMesBias = NULL, estimPm = NULL,
                               checkEquality = "no", checkNullityBias = "no")
{

  assertPlan(plan)

  N <- problem$N

  if (checkEquality != "no" || checkNullityBias != "no")
    problem <- problem$copy(deep = TRUE)

  Ytilde <- problem$Ytilde(plan$mode)

  if (checkEquality %in% c("MCO", "MCO_agreg"))
  {
    mergeEquivalents <- checkEquality == "MCO_agreg"
    plan <-
      check_modes_equality_MCO(problem, plan,
                               mergeEquivalents = mergeEquivalents)$plan
  }

  delta <- numeric(N)

  usedBiasedModes <- modes_biaises_utilises(plan, problem$modesRef)

  # Estimation of measure bias without imputation
  if (estimMesBias == "G-COMP")
  {
    for (biasedMode in usedBiasedModes)
    {
      masqueMode <- masque_repondants_mode(plan, biasedMode)
      delta[masqueMode] <-
        estim_delta_G_comp_MCO(plan, problem, biasedMode)[masqueMode]
    }
  }

  # Estimation with imputation
  else
  {
    if (imputationY == "true_value")
      Ycf <- problem$Y

    else if (imputationY == "MCO")
      Ycf <- MCO_Y(plan, problem) %>% fitted()

    else
      Ycf <- Y_matching(plan, problem, imputationY)

    # Case when Y^m is replaced by its counterfactuel
    if (estimMesBias == "CF")
    {
      delta <- Ytilde - Ycf
      usedBiasedModes <- NULL
    }

    # Case when we want to check if there is negligeable biases
    else if (checkNullityBias == "no")
    {
      # MCO_Ym : Ym is maintaned if the bias is negligeable
      # MCO_CF : Ym is replaced by its counterfactual if the bias is negligeable
      if (checkNullityBias %in% c("MCO_Ym", "MCO_CF"))
      {
        replaceByCF <- checkNullityBias == "MCO_CF"
        resCheck <- check_nullity_bias_MCO(problem, plan,
                                           Ycf, replaceByCF = replaceByCF)
      }

      delta <- resCheck$delta

      usedBiasedModes <- resCheck$remainingBiasedModes
    }

    # For the remaining biased modes with biaises to determine
    # we do some estimations
    for (biasedMode in usedBiasedModes)
    {
      masqueMode <- masque_repondants_mode(plan, biasedMode)
      delta[masqueMode] <- estimMesBias(Xm = problem$Xm(masqueMode),
                                        Ym = Ytilde[masqueMode],
                                        Ycf = Ycf[masqueMode])
    }
  }

  masqueRepRef <-
    masque_repondants_reference(plan, problem$modesRef)
  delta[masqueRepRef] <- 0.0


  if (estimPm == "true_values")
    pHat <- problem$probaModes

  else if (estimPm == "multinomial")
  {
    pHat <- estim_response_prob_global(plan, problem, regroupModesRef = TRUE)
  }

  HT_mm(plan, Ytilde, delta,
        problem$modesRef, pHat, K = problem$K)
}
