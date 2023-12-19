#' Parametric regression for modes
#'
#' Estimate Ym on the whole population based
#' on the subpopulation that used this mode
#' to answer
#' @importFrom stats lm predict.lm
#' @keywords internal
#' @name regression_Ym
estim_Ym_MCO <- function(sample, m)
{
  maskMode <- sample$respondents_mode(m)

  if (!any(maskMode))
    glue("Nobody answered with the mode {m}") %>% abort()

  Ym <- sample$Ym_resp(m) # nolint: object_usage_linter
  Xm <- sample$Xm(m) %>% as.matrix()# nolint: object_usage_linter


  lm(Ym ~ Xm) %>% predict.lm(newdata = sample$X)
}

#' @inheritParams regression_Ym
#' @importFrom stats predict.lm
estim_delta_G_comp_MCO <- function(sample)
{
  Ycf <- MCO_Y(sample)


  usedBiasedModes <- sample$used_biased_modes()

  if (length(usedBiasedModes) == 0L)
    return(numeric(sample$N))

  delta <- rep(NA_real_, sample$N)
  delta[sample$respondents_ref()] <- 0.0


  for (m in usedBiasedModes)
  {
    maskRespMode <- sample$respondents_mode(m)
    delta[maskRespMode] <- (estim_Ym_MCO(sample, m) - Ycf)[maskRespMode]
  }

  return(delta)
}

#' @importFrom checkmate assertFlag
#' @importFrom stats lm coef
estim_MB_MCO_cf <- function(sample, m, Ycf)
{
  maskMode <- sample$respondents_mode(m)
  n <- sum(maskMode)
  assertY(Ycf, N = n)

  lm(sample$Ytilde() - Ycf ~ sample$X, subset = maskMode) %>%
    coef(newdata = sample$X)
}

#' @importFrom checkmate assertChoice
estim_MB_MCO_tot <- function(sample, m, typeTot = "doubleHT", Ycf = NULL)
{
  assertChoice(typeTot, c("doubleHT", "totYimp", "estimDeltaCF"))
  maskMode <- sample$respondents_mode(m)

  if (!any(maskMode))
    return(numeric(sample$N))

  X <- sample$X %>% as.matrix()
  X <- cbind(rep(1.0, nrow(X)), X)
  phi <- sample$phi_tab[, m]
  inverseMat <- (t(X) %*% diag(phi) %*% X) %>% solve()
  Xm <- sample$Xm(m) %>% as.matrix()
  Xm <- cbind(rep(1.0, nrow(Xm)), Xm)
  Ym <- sample$Ym_resp(mode = m, answersOnly = TRUE)
  phim <- phi[maskMode]
  weightsm <- sample$weights(rep(m, sample$N), inv = FALSE)[maskMode]

  if (typeTot == "estimDeltaCF")
  {
    assertNumericVector(Ycf, len = sample$N, null.ok = FALSE)

    delta <- Ym - Ycf[maskMode]
    HTdelta <- apply(diag(delta * phim * weightsm) %*% Xm, 2L, sum)
    coefs <- inverseMat %*% HTdelta
  }
  else if (typeTot == "totYimp")
  {
    weightsref <- sample$weights_ref(inv = FALSE)
    totYcf <- apply(diag(Ycf * phi) %*% X, 2L, sum)
    HTm <- apply(diag(Ym * phim * weightsm) %*% Xm, 2L, sum)
    coefs <- inverseMat %*% (HTm - totYcf)
  }

  else if (typeTot == "doubleHT")
  {
    Xref <- sample$Xref() %>% as.matrix()
    Yref <- sample$Yref_resp(answersOnly = TRUE)
    weightsref <- sample$weights_ref(inv = FALSE)[sample$respondents_ref()]
    phiref <- phi[sample$respondents_ref()]

    HTm <- apply(diag(Ym * phim * weightsm) %*% Xm, 2L, sum)
    HTref <- apply(diag(Yref * phiref * weightsref) %*% Xref, 2L, sum)

    coefs <- inverseMat %*% (HTm - HTref)
  }

  coefs <- as.vector(coefs)

  X %*% coefs
}

#' @importFrom stats qf
#' @importFrom checkmate assertFlag
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
