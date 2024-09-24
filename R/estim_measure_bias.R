#' Parametric regression for modes
#'
#' Estimate delta_m MCO parameters for Ym
#' based on the subpopulation that used this mode
#' to answer
#' @importFrom stats lm predict.lm
#' @keywords internal
#' @name regression_Ym
estim_coefs_Ym_MCO <- function(Z, Yobs, m, modes)
{
  maskMode <- modes == m

  if (!any(maskMode))
    glue("Nobody answered with the mode {m}") %>% abort()

  Ym <- Yobs[maskMode] # nolint: object_usage_linter
  Zm <- Z[maskMode, , drop = FALSE] # nolint: object_usage_linter

  data <- cbind(Ym, Zm) %>% as.data.frame()

  lm(Ym ~ . -1L, data = data) %>% coef()
}



#' @export
estim_delta_MCO_by_totals <- function(Z, totalBiased, totalRef)
{
  solve(t(Z) %*% Z) %*% (totalBiased - totalRef)
}

#' @export
estim_MB_by_delta_MCO_HT <- function(delta, Z, pi,
                                     m, modes,
                                     probsSelect,
                                     phi = rep(1.0 / 2.0, nrow(Z)))
{
  mask <- modes == m
  Z <- Z[mask, , drop = FALSE]
  pi <- pi[mask]
  probsSelect <- probsSelect[mask]
  phi <- phi[mask]

  weights <- (pi * probSelect)^-1L

  estim_MB_by_MCO(delta, Z, weights, phi)
}

#' @export
estim_MB_by_MCO <- function(delta, Z,
                            weights = rep(1.0, nrow(Z)),
                            phi = rep(1/2, nrow(Z)),
                            mask = rep(TRUE, nrow(Z)))
{
  diag(weights[mask] * phi[mask]) %*% Z[mask, , drop = FALSE] %*% delta
}


#' @inheritParams regression_Ym
#' @importFrom stats predict.lm
#' @export
estim_delta_G_comp_MCO <- function(Z, Yobs, modes, biasedMode, refMode)
{
  coefsBiased <- estim_coefs_Ym_MCO(Z, Yobs, biasedMode, modes)
  coefsRef <- estim_coefs_Ym_MCO(Z, Yobs, refMode, modes)

  coefsBiased - coefsRef
}

#' @export
estim_delta_MCO <- function(Z, Yobs,
                            modes, biasedMode, refMode,
                            modeTotBiased = "HT", modeTotRef = "HT",
                            pi, probsSelect, sampleMatrix = FALSE)
{

  Z <- as.matrix(Z)

  N <- length(Yobs)
  if (modeTotBiased == "HT" || modeTotRef == "HT")
  {
    weights <- numeric(N)
    weights[!is.na(probsSelect)] <-
      (pi[!is.na(probsSelect)] * probsSelect[!is.na(probsSelect)])^-1L

    # Case when all weights are equal:
    if (sampleMatrix && all(weights == weights[1L]))
      modeTotBiased <- modeTotRef <- "MCO"
  }

  if ((modeTotBiased == "MCO" && modeTotRef == "MCO") ||
      modeTotBiased == "G-COMP" ||
      modeTotRef == "G-COMP")
  {
    return(estim_delta_G_comp_MCO(Z, Yobs, modes, biasedMode, refMode))
  }


  ZtZInv <- NULL
  if (modeTotBiased == "HT")
  {
    maskBiased <- modes == biasedMode
    Zbiased <- Z[maskBiased, , drop = FALSE]
    totBiased <- t(Zbiased) %*%
      diag(weights[maskBiased]) %*%
      Yobs[maskBiased]

    if (sampleMatrix)
    {
      coefsBiased <- solve(t(Zbiased) %*%
                             diag(weights[maskBiased]) %*%
                             Zbiased) %*%
        totBiased
    }
    else
    {
      ZtZInv <- solve(t(Z) %*% Z)
      coefsBiased <- ZtZInv %*% totBiased
    }

  }
  else if (modeTotBiased == "MCO")
    coefsBiased <- estim_coefs_Ym_MCO(Z, Yobs, biasedMode, modes)
  else
    abort("No method recognized for the total of the bias mode.")

  if (modeTotRef == "HT")
  {
    maskRef <- modes == refMode
    Zref <- Z[maskRef, , drop = FALSE]
    totRef <- t(Zref) %*%
      diag(weights[maskRef]) %*%
      Yobs[maskRef]

    if (sampleMatrix)
    {
      coefsRef <- solve(t(Zref) %*%
                          diag(weights[maskRef])
                        %*% Zref) %*%
      totRef
    }
    else
    {
      if (is.null(ZtZInv))
        ZtZInv <- solve(t(Z) %*% Z)

      coefsRef <- ZtZInv %*% totRef
    }

  }
  else if (modeTotRef == "MCO")
  {
    coefsRef <- estim_coefs_Ym_MCO(Z, Yobs, refMode, modes)
  }
  else
    abort("No method recognized for the total of the reference mode.")

  delta <- coefsBiased - coefsRef

  delta
}

# Équivalent G-COMP
# estim_delta_MCO_unique_model <- function(Z, Yobs,
#                                          modes, biasedMode, refMode)
# {
#   maskBiased <- modes == biasedMode
#
#   maskRef <- modes == refMode
#
#   Zbiased <- Z[maskBiased, , drop = FALSE]
#   Zref <- Z[maskRef, , drop = FALSE]
#
#   ZtZBiased <- t(Zbiased) %*% Zbiased
#   ZtZRef <- t(Zref) %*% Zref
#
#
#   Ybiased <- Yobs[maskBiased]
#   Yref <- Yobs[maskRef]
#
#   beta <- solve(ZtZRef) %*% t(Zref) %*% Yref
#
#   delta <- solve(ZtZBiased) %*% t(Zbiased) %*% Ybiased - beta
#
#   as.numeric(delta)
#
# }

#' @export
estim_delta_MCO_unique_model_const <- function(Z, Yobs,
                                               modes, biasedMode, refMode)
{
  maskBiased <- modes == biasedMode
  nBiased <- sum(maskBiased)

  maskRef <- modes == refMode

  Zbiased <- Z[maskBiased, , drop = FALSE]
  Zref <- Z[maskRef, , drop = FALSE]
  transZref <- t(Zref)

  ZtZBiased <- t(Zbiased) %*% Zbiased
  Ybiased <- Yobs[maskBiased]
  sumYBiased <- sum(Ybiased)
  Yref <- Yobs[maskRef]

  sumZBiased <- apply(Zbiased, 2L, sum)
  transSumZBiased <- t(sumZBiased)

  A <- ZtZBiased -
    sumZBiased %*% transSumZBiased / nBiased +
    transZref %*% Zref

  b <- t(Zbiased) %*% Ybiased -
    sumZBiased * sumYBiased / nBiased +
    transZref %*% Yref

  beta <- solve(A, b)


  delta <- sumYBiased - transSumZBiased %*% beta


  delta <- delta / nBiased

  as.numeric(delta)
}


## À vérifier
#' @importFrom stats qf
#' @importFrom checkmate assertFlag
#' @export
check_nullity_bias_MCO <- function(sample, Ycf, alpha = 0.01,
                                   replaceByCF = FALSE)
{
  assertFlag(replaceByCF)
  assertProbabilityScalar(alpha, striclyPos = TRUE, striclyUnsure = TRUE)
  assertModesRef(sample$modesRef)

  N <- sample$N
  usedBiasedModes <- sample$used_biased_modes()

  delta <- rep(NA_real_, N)
  masqueRepondantRef <- sample$respondents_ref()
  delta[masqueRepondantRef] <- 0.0

  modifs <- character(0L)
  Ytilde <- sample$Ytilde()

  unbiasedModes <- NULL
  for (mode in usedBiasedModes)
  {
    masqueMode <- sample$respondents_mode(mode)

    Xm <- sample$Xm(masqueMode)
    Ym <- Ytilde[masqueMode]
    Ycfm <- Ycf[masqueMode]

    model <- modelisation_delta_MCO(Xm, Ym, Ycfm)

    FstatData <- summary(model)$fstatistic
    Fstat <- FstatData[1L]
    df1 <- FstatData[2L]
    df2 <- FstatData[3L]
    qFisher <- qf(1.0 - alpha, df1, df2)

    if (Fstat <= qFisher)
    {
      if (replaceByCF)
        delta[masqueMode] <- Ym - Ycfm
      else
        delta[masqueMode] <- 0.0

      modifs <- c(modifs, paste0(mode, "->", sample$modesRef[1L]))
      unbiasedModes <- c(unbiasedModes, mode)
    }
  }

  list(changes = modifs,
       remainingBiasedModes = setdiff(usedBiasedModes, unbiasedModes),
       delta = delta)
}
