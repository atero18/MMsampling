#' Estimation of the counterfactuals of each respondent
#'
#' In this version considers only the case with two modes. Could be
#' generalised.
#'
#' @param Yobs vector of the observed outcomes. For the non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The first mode of the protocol is defined as "m1" and the second as "m2"
#' (character vector or factor of size N).
#' @param X design matrix. The rows corresponding to the non-respondents
#' can contain NA (numeric matrix with N rows and p columns).
#' @param ... arguments for the function `MatchIt::matchit`.
#' @return for each mode a vector of size N containing for each respondent
#' its outcome or its estimated counterfactual depending on if it answered
#' by the mode or the other. Equal to NA for the non-respondents
#' (numeric matrix of dimension (N,2)).
#' @export
estim_counterfactuals <- function(Yobs,
                                  modes,
                                  X,
                                  ...)
{
  ## Generalise to any number of modes

  requireNamespace(package = "MatchIt", versionCheck = NULL, quietly = FALSE)

  maskSr <- modes == "m1"
  maskSmr <- modes == "m2"
  maskSa <- maskSr | maskSmr

  # Vectors that contains "estimation" of the counterfactuals of the outcomes.
  Y1Est <- Y2Est <- Yobs
  YobsSa <- Yobs[maskSa]

  dataMatchingSa <- as.data.frame(X[maskSa, , drop = FALSE])
  dataMatchingSa$mode <- as.integer(modes[maskSa] == "m1")
  row.names(dataMatchingSa) <- seq_len(sum(maskSa))

  # Estimation of the y_1k on Smr
  matchingSmr <-
    MatchIt::matchit(mode ~ ., data = dataMatchingSa,
                     estimand = "ATC", replace = TRUE, ...)

  Y1Est[maskSmr] <- YobsSa[as.integer(matchingSmr$match.matrix[, 1L])]

  # Estimation of the y_2k on Sr
  matchingSr <-
    MatchIt::matchit(mode ~ ., data = dataMatchingSa,
                     estimand = "ATT", replace = TRUE, ...)

  Y2Est[maskSr] <- YobsSa[as.integer(matchingSr$match.matrix[, 1L])]

  cbind(Y1Est = Y1Est, Y2Est = Y2Est)
}

#' @importFrom stats predict
estim_tot_imput <- function(X, model)
{
  Ypred <- predict(model, newdata = X)
  sum(Ypred)
}

#' @importFrom stats lm predict
MCO_Y <- function(sample)
{
  maskRespRef <- sample$respondents_ref()

  Ytilde <- sample$Yref() # nolint: object_usage_linter
  X <- sample$X %>% as.matrix() # nolint: object_usage_linter

  lm(Ytilde ~ X, subset = maskRespRef) %>%
    predict(newdata = sample$X)
}

#' @importFrom checkmate assertChoice
#' @importFrom stats coef
#' @export
estim_tot_MCO_Y <- function(plan, modesRef, X, Y, type = "pred")
{
  type <- tolower(type)
  assertChoice(type, c("pred", "beta"))
  modele <- MCO_Y(plan, modesRef, X, Y)

  if (type == "pred")
    estim_tot_imput(X, modele)

  else if (type == "beta")
  {
    betaSR <- coef(modele)
    (t(X) %*% X %*% betaSR)[1L]
  }
}

#' @importFrom checkmate assertVector
estim_tot_GH <- function(sample, groups)
{
  assertVector(groups,
               any.missing = FALSE, all.missing = FALSE,
               len = length(Y))

  if (length(unique(groups)) == 1L)
    return(estim_tot_1GH(plan))

  masqueRepRef <- sample$respondents_ref()

  groupSizes <- table(groups)
  groupNames <- names(groupSizes)
  groups <- groups[masque_repondants_reference]

  if (!all(groupNames %in% groups))
      stop("At least one group is not represented")

  crossprod(tapply(sample$Yref_resp(), groups, mean), groupSizes)
}

