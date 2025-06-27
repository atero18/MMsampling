#' Estimation of t1
#'
#' Estimation of the total of the first mode of response potential
#' outcomes \eqn{t_1}.
#'
#' @inheritParams common_arguments
#' @return the estimation of the total \eqn{t_1}, possibly weighted (scalar).
#' @family total estimators
#' @export
estim_tot_m1 <- function(Yobs, modes, pi, p1, phi = rep(1.0, length(Yobs)))
{
  maskSr <- modes == "m1"

  sum((pi[maskSr] * p1[maskSr])^-1L * phi[maskSr] * Yobs[maskSr])
}

#' Estimation of t2
#'
#' Estimation of the total of the second mode of response
#' potential outcomes \eqn{t_2}.
#'
#' @inheritParams common_arguments

#' @return the estimation of the total \eqn{t_2}, possibly weighted (scalar).
#' @family total estimators
#' @export
estim_tot_m2 <- function(Yobs, modes, pi, p1, p2, phi = rep(1.0, length(Yobs)))
{
  maskSmr <- modes == "m2"

  sum((pi[maskSmr] * (1.0 - p1[maskSmr]) * p2[maskSmr])^-1L *
        phi[maskSmr] * Yobs[maskSmr])
}


#' Estimation of tphiy
#'
#' Estimation of the total of the weighted potential outcomes
#' of each mode of response \eqn{t_{y\phi}}.
#'
#' @inheritParams common_arguments
#' @param phi optional vector containing weights for the total.
#' The values must be between 0 and 1. The more the weight is close to one,
#' the more the first mode of reponse is important (numeric vector size N).
#' @return the estimation of the total \eqn{t_{y\phi}} (scalar).
#' @family total estimators
#' @export
estim_tot_phi_y <- function(Yobs, modes, pi, p1, p2, phi = rep(0.5, length(Yobs)))
{
  estim_tot_m1(Yobs, modes, pi, p1, phi) +
    estim_tot_m2(Yobs, modes, pi, p1, p2, 1.0 - phi)
}

#' Estimation of t0
#'
#' Estimation of the total \eqn{t_0} under an almost sure equality
#'  between the potential outcomes of each unit. The total can be estimated
#'  using the mode selection probabilities or the global answer probabilities.
#' @name estim_t0
#' @family total estimators
#' @return the estimation of the total \eqn{t_0} (scalar).
NULL

#' @describeIn estim_t0 Sum an estimator per mode. Do not use the first mode
#' selection probabilities.
#' @inheritParams common_arguments
#' @export
estim_tot_0_no_p1 <- function(Yobs, modes, pi, p2)
{
  p1 <- phi <- rep(1.0, length(Yobs))

  estim_tot_phi_y(Yobs, modes, pi, p1, p2, phi)
}


#' @describeIn estim_t0 Use the probability to answer by any mode instead
#' of by each mode.
#' @param responses vector indicating which unit answered.
#' If a unit answered then the corresponding value is "r"
#' (character vector or factor of size N).
#' @param pa vector containing the true or estimated global answer probability
#' (numeric vector of size N).
#' @export
estim_tot_a <- function(Yobs, responses, pi, pa)
{
  maskSa <- responses == "r"

  sum((pi[maskSa] * pa[maskSa])^-1L * Yobs[maskSa])
}
