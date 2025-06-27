#' Estimation of t_1
#'
#' @param Yobs vector of the observed outcomes. For the m1 non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The first mode of the protocol is defined as "m1"
#' (character vector or factor of size N).
#' @param pi vector containing the inclusion probabilities pi_k. Can be
#' equal to NA for the m1 non-respondents (numeric vector of size N).
#' @param p1 vector containing the true or estimated probabilities p_1k.
#' Can be equal to NA for the m1 non-respondents (numeric vector of size N).
#' @param phi optional vector containing weights for the total.
#' Must be known for the m1-respondents (numeric vector of [0,1]^N).
#' @return the estimation of the total t_1, possibly weighted (scalar).
#' @export
estim_tot_m1 <- function(Yobs, modes, pi, p1, phi = rep(1.0, length(Yobs)))
{
  maskSr <- modes == "m1"

  sum((pi[maskSr] * p1[maskSr])^-1L * phi[maskSr] * Yobs[maskSr])
}

#' Estimation of t_2
#' @param Yobs vector of the observed outcomes. For the m2 non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The second mode of the protocol is defined as "m2"
#' (character vector or factor of size N).
#' @param pi vector containing the inclusion probabilities pi_k. Can be
#' equal to NA for the m2 non-respondents (numeric vector of size N).
#' @param p1 vector containing the true or estimated probabilities p_1k.
#' Can be equal to NA for the m2 non-respondents (numeric vector of size N).
#' @param p2 vector containing the true or estimated probabilities p_2k.
#' Can be equal to NA for the m2 non-respondents (numeric vector of size N).
#' @param phi optional vector containing weights for the total. Must be
#' known for the m2-respondents (numeric vector of [0,1]^N).
#' @return the estimation of the total t_2, possibly weighted (scalar).
#' @export
estim_tot_m2 <- function(Yobs, modes, pi, p1, p2, phi = rep(1.0, length(Yobs)))
{
  maskSmr <- modes == "m2"

  sum((pi[maskSmr] * (1.0 - p1[maskSmr]) * p2[maskSmr])^-1L *
        phi[maskSmr] * Yobs[maskSmr])
}


#' Estimation of t_phiy
#' @param Yobs vector of the observed outcomes. For the non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The first mode of the protocol is defined as "m1" and the second as "m2"
#' (character vector or factor of size N).
#' @param pi vector containing the inclusion probabilities pi_k. Can be
#' equal to NA for the non-respondents (numeric vector of size N).
#' @param p1 vector containing the true or estimated probabilities p_1k.
#' Can be equal to NA for the non-respondents (numeric vector of size N).
#' @param p2 vector containing the true or estimated probabilities p_2k.
#' Can be equal to NA for the non-respondents (numeric vector of size N).
#' @param phi optional vector containing weights for the total. Must be
#' known for the respondents (numeric vector of [0,1]^N).
#' @return the estimation of the total t_phiy (scalar).
#' @export
estim_tot_m1_m2 <- function(Yobs, modes, pi, p1, p2, phi = rep(0.5, length(Yobs)))
{
  estim_tot_m1(Yobs, modes, pi, p1, phi) +
    estim_tot_m2(Yobs, modes, pi, p1, p2, 1.0 - phi)
}

#' Estimation of t_1
#'
#' @param Yobs vector of the observed outcomes. For the m1 non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The first mode of the protocol is defined as "m1"
#' (character vector or factor of size N).
#' @param pi vector containing the inclusion probabilities pi_k. Can be
#' equal to NA for the m1 non-respondents (numeric vector of size N).
#' @param p1 vector containing the true or estimated probabilities p_1k.
#' Can be equal to NA for the m1 non-respondents (numeric vector of size N).
#' @param phi optional vector containing weights for the total.
#' Must be known for the m1-respondents (numeric vector of [0,1]^N).
#' @return the estimation of the total t_1, possibly weighted (scalar).
#' @export
estim_tot_m1 <- function(Yobs, modes, pi, p1, phi = rep(1.0, length(Yobs)))
{
  maskSr <- modes == "m1"

  sum((pi[maskSr] * p1[maskSr])^-1L * phi[maskSr] * Yobs[maskSr])
}


#' Estimation of t_0
#'
#' @param Yobs vector of the observed outcomes. For the non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The first mode of the protocol is defined as "m1" and the second as "m2"
#' (character vector or factor of size N).
#' @param pi vector containing the inclusion probabilities pi_k. Can be
#' equal to NA for the non-respondents (numeric vector of size N).
#' @param p2 vector containing the true or estimated probabilities p_2k.
#' Can be equal to NA for the non-respondents by m2 (numeric vector of size N).
#' @return the estimation of the total t_0 (scalar).
#' @export
estim_tot_as_eq_no_p1 <- function(Yobs, modes, pi, p2)
{
  p1 <- phi <- rep(1.0, length(pi))

  estim_tot_m1_m2(Yobs, modes, pi, p1, p2, phi)
}

#' Estimation of t_ay
#' @param Yobs vector of the observed outcomes. For the non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param responses vector indicating which unit answered.
#' If a unit answered then the corresponding value is "r"
#' (character vector or factor of size N).
#' @param pi vector containing the inclusion probabilities pi_k. Can be
#' equal to NA for the non-respondents (numeric vector of size N).
#' @param pa vector containing the true or estimated probabilities p_ak.
#' Can be equal to NA for the non-respondents (numeric vector of size N).
#' @return the estimation of the total t_ay (scalar).
#' @export
estim_tot_a <- function(Yobs, responses, pi, pa)
{
  maskSa <- responses == "r"

  sum((pi[maskSa] * pa[maskSa])^-1L * Yobs[maskSa])
}
