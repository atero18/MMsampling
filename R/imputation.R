#' Estimation of the counterfactuals of each respondent
#'
#' In this version considers only the case with two modes. Could be
#' generalised. Requires the package [MatchIt](https://cran.r-project.org/web/packages/MatchIt/).
#'
#' @inheritParams common_arguments
#' @param X design matrix. The rows corresponding to the non-respondents
#' can contain NA (numeric matrix with N rows and p columns).
#' @param ... arguments for the function `MatchIt::matchit`.
#' @return for each mode a vector of size N containing for each respondent
#' its outcome or its estimated counterfactual depending on if it answered
#' by the mode or the other. Equal to NA for the non-respondents
#' (numeric matrix of dimension (N,2)).
#' @seealso [MatchIt::matchit()]
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
