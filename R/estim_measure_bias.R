#' Parametric regression for modes
#'
#' Estimate Ym on the whole population based
#' on the subpopulation that used this mode
#' to answer
#' @importFrom stats lm
#' @keywords internal
#' @name regression_Ym
MCO_Ym <- function(sample, m)
{
  maskMode <- sample$respondents_mode(m)

  if (!any(maskMode))
    glue("Nobody answered with the mode {m}") %>% abort()

  Ym <- sample$Ym_resp(m) # nolint: object_usage_linter
  Xm <- sample$Xm(m) # nolint: object_usage_linter


  lm(Ym ~ Xm)
}

#' @inheritParams regression_Ym
#' @importFrom stats predict.lm
#' @export
estim_delta_G_comp_MCO <- function(X, Ytilde, choixMode, modesRef, modes = NULL)
{
  N <- length(Ytilde)

  if (is.null(modes))
    modes <- setdiff(unique(choixMode), modesRef)
  else
    modes <- unique(modes)

  if (length(modes) == 0L)
    abort("No mode asked fort bias estimation")

  modeleY <- MCO_Y(X, Ytilde, modesRef, choixMode)
  predY <- predict.lm(modeleY, newdata = X)

  estimYm <- function(m)
  {
    if (m %in% modesRef)
      return(numeric(N))

    else
    {
      modeleYm <- MCO_Ym(Xm = X, Ym = Ytilde, m = m, choixMode = choixMode)
      predict.lm(modeleYm, newdata = X) - predY
    }
  }

  vapply(modes, estimYm, numeric(N))
}

#' @importFrom checkmate assertFlag
#' @importFrom stats lm coef
coefs_MB_MCO_cf <- function(sample, m, Ycf)
{
  maskMode <- sample$respondents_mode(m)
  n <- sum(maskMode)
  assertY(Ycf, N = n)

  lm(sample$Ytilde() - Ycf ~ sample$X, subset = maskMode) %>% coef()
}

#' @importFrom checkmate assertChoice
coefs_MB_MCO_tot <- function(sample, m, typeTot = "doubleHT", Ycf = NULL)
{
  assertChoice(typeTot, c("doubleHT", "totYcf", "doubleHT"))
  maskMode <- sample$respondents_mode(m)

  phi <- sample$phi_tab[, m]
  inverseMat <- (t(sample$X) %*% diag(phi) %*% sample$X) %>% solve()
  Xm <- sample$Xm(m)
  Ym <- sample$Ym_resp(answersOnly = TRUE)
  phim <- phi[maskMode]
  weightsm <- sample$weights(m, inv = FALSE)[maskMode]

  if (typeTot == "estimDeltaCF")
  {
    delta <- Ym - Ycf[maskMode]
    HTdelta <- Xm %*% delta * phim * weightsm
    inverseMat %*% HTdelta
  }
  else if (typeTot == "totYcf")
  {
    totYcf <- sum(Ycf)
    HTm <- Xm %*% Ym * phim * weightsm
    inverseMat %*% (HTm - totYcf)
  }

  else if (typeTot == "doubleHT")
  {
    Xref <- sample$Xref()
    Yref <- sample$Yref_resp(answersOnly = TRUE)
    weightsref <- sample$weights_ref(inv = FALSE)
    phiref <- phi[sample$respondents_ref()]

    HTm <- Xm %*% Ym * phim * weightsm
    HTref <- Xref %*% Yref * phiref * weightsref

    inverseMat %*% (HTm - HTref)
  }
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
