#' Estimation of the covariance between the potential outcomes of two modes
#' of response under homoscedasticity assumptions and linearity
#'
#' We suppose the homoscedasticity within the first and second
#' and between the modes of response potential outcomes.
#' We will have an unbiased estimator if the non-informativeness,
#' unconfoundedness and conditional mutual independence assumptions
#' are verified ; the measure effect y_1k - y_2k  follow a linear model
#' with x_k and for some units both potential outcomes are known
#' (which can't be under a sequential protocol). We must have at least p + 1
#' (p the number of covariates) units such that their measure effect are
#' known or well estimated. The standard deviations `sd1` and `sd2` must be
#' unbiased as well.
#' @param Y1 vector of the first mode outcomes
#' (numeric vector of size N the size of the population).
#' @param Y2 vector of the second mode outcomes
#' (numeric vector of size N).
#' @param sd1 true or estimated value of the standard deviation of the first
#' mode of response potential outcomes (positive scalar).
#' @param sd2 true or estimated value of the standard deviation of the
#' second mode of response potential outcomes (positive scalar).
#' @param X design matrix (numeric matrix with N rows and p columns).
#' @param Yobs vector of the observed outcomes. Optional. Useful if `Y1`
#' and `Y2` are not given. In that case the counterfactuals of the respondents
#' are estimated with the function `estim_counterfactuals`.
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit. Optional. Used if
#' the counterfactuals must be estimated
#' (character vector or factor of size N).
#' @param clamp TRUE if the estimation of the covariance must be clamped if
#' its absolute value is superior to `sd1` * `sd2` (boolean).
#' @param warnClamp TRUE if a warning must be sent when a clamp is made (boolean).
#' @param ... arguments for the function `MatchIt::matchit`.
#' @return an estimator of the covariance between the m1 and m2
#' potential outcomes, unbiased under the assumptions, if `sd1`^2 and `sd2`^2
#' are unbiased and the counterfactuals are the true values (scalar).
#' @export
estim_cov_12 <- function(Y1, Y2, sd1, sd2, X, Yobs = NULL, modes = NULL,
                         clamp = FALSE, warnClamp = TRUE, ...)
{
  if (sd1 == 0.0 || sd2 == 0.0)
    return(0.0)

  # If the observation vector Yobs is given we do the estimation
  # of the counterfactuals
  if (!is.null(Yobs))
  {
    Yest <- estim_counterfactuals(Yobs, modes, X, ...)
    Y1 <- Yest[, "Y1Est"]
    Y2 <- Yest[, "Y2Est"]
    maskResp <- modes %in% c("m1", "m2")
  }
  else
  {
    # The estimation is made with units such that
    # both potential outcomes are known (or generally estimated)
    maskResp <- !is.na(Y1) & !is.na(Y2)
  }

  # Assuming a linear model for the differences y_1k - y_2k
  sdDelta <- summary.lm(lm(Y1 - Y2 ~ X, subset = maskResp))$sigma

  c12 <- (sd1^2L - sdDelta^2L + sd2^2L) / 2.0

  # If |c12| > sd1 * sd2 because the model is not correct,
  # we clamp it if clamp == TRUE
  if (isTRUE(clamp) && abs(c12) > sd1 * sd2)
  {
    c12 <- sign(c12) * sd1 * sd2
    if (isTRUE(warnClamp))
    {
      warningClamp <- paste("Covariance between the two modes counterfactuals",
                            "not well estimated. Value clampd")
      warning(warningClamp)
    }
  }


  c12
}
