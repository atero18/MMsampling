#' Parametric regression for modes
#'
#' Estimate Ym on the whole population based
#' on the subpopulation that used this mode
#' to answer
#' @importFrom stats lm
#' @keywords internal
#' @name regression_Ym
MCO_Ym <- function(Xm, Ym, m = NULL, choixMode = NULL)
{
  if (is.null(m))
    mask <- rep(TRUE, length(Ym))
  else
  {
    mask <- choixMode == m
    if (!any(mask))
      glue("Nobody answered with the mode {m}") %>% abort()
  }

  lm(Ym ~ Xm, subset = mask)
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
#' @export
modelisation_delta_MCO <- function(Xm, Ym, Ycf, addYm = TRUE)
{
  assertY(Ym)
  assertFlag(addYm)

  n <- length(Ym)

  assertY(Ycf, N = n)

  ## Assert X

  if (addYm)
    Xm <- cbind(Xm, Ym)

  lm(Ym - Ycf ~ Xm)
}

#' @importFrom stats predict.lm
estim_delta_MCO <- function(Xm, Ym, Ycf)
{
  model <- modelisation_delta_MCO(Xm, Ym, Ycf)
  predict.lm(model, newdata = Xm)
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

      modifs <- c(modifs, paste0(mode, "->", problem$modesRef[1L]))
      unbiasedModes <- c(unbiasedModes, mode)
    }
  }

  list(changes = modifs,
       remainingBiasedModes = setdiff(usedBiasedModes, unbiasedModes),
       delta = delta)
}
