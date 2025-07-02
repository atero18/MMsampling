#' Estimation of the variance of t_ephiy and the approximate variance
#' of t_pqphiy
#'
#'
#' HOMOSCEDASTICITY, INDEPENDENCE
#' @param Yobs vector of the observed outcomes. For the non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes vector of the selected mode of each unit.
#' The first mode of the protocol is defined as "m1" and the second as "m2"
#' (character vector or factor of size N).
#' @param pi_mat matrix containing the second-order inclusion probabilities
#' pi_kl. The values must be known on Sr^2 (except if `m2Only` = TRUE) and
#' Smr^2 (except if `m1Only` = TRUE) and can be equal to NA otherwise
#' (symmetric numeric matrix of order N).
#' @param p1 vector containing the true or estimated probabilities p_1k.
#' Can be equal to NA for the non-respondents or the units that didn't answer
#' by m1 if `m1Only` = TRUE (numeric vector of size N).
#' @param p2 vector containing the true or estimated probabilities p_2k.
#' Used if `m1Only` = FALSE. Can be equal to NA for the non-respondents
#' or the units that didn't answer by m2 if `m2Only` = TRUE (numeric vector of size N).
#' @param sd1 true or estimated value of the standard deviation of the m1 potential outcomes.
#' Used when `m2Only` = FALSE and `phi` != 0 (positive scalar).
#' @param sd2 true or estimated value of the standard deviation of the
#' m2 potential outcomes.  Used when `m1Only` = FALSE and `phi` != 1
#' (positive scalar).
#' @param cov12 the true or estimated covariance between the m1 and m2
#' potential outcomes. Default to `sd1` * `sd2` ; i.e. Y1 and Y2 are
#' conditionally collinear (scalar).
#' @param phi vector containing weights for the total. Must be known on the
#' set of respondents (numeric vector of [0,1]^N).
#' @param eqPotOut TRUE if we consider that y_1k = y_2k on the population
#' (boolean)
#' @param correcEW1 TRUE if we must apply a correction due to the
#' use of estimation of `p1` (boolean).
#' @param correcEW2 TRUE if we must apply a correction due to the
#' use of estimation of `p2` (boolean).
#' @param I vector containing the indicators I_k equal to TRUE if unit k
#' has been selected. Used only if `correcEW1` or `correcEW2` is true
#' (logical vector of size N).
#' @param Z design matrix. Used only when `correcEW1` or `correcEW2` is TRUE
#' (numeric matrix with N rows and q columns).
#' @param m1Only TRUE if we estimate the variance (resp. approximate variance)
#' of t_{e1} (resp. t_{pq1}) and only with the m1 respondents (boolean).
#' @param m2Only TRUE if we estimate the variance (resp. approximate variance)
#' of t_{e2} (resp. t_{pq2}) and only with the m2 respondents (boolean).
#' @param weightSrEstVp weight considering the impact of the sampling design
#' variance estimator made of the m1 respondents vs the m2 respondents. Used
#' only when `m1Only` and `m2Only` are FALSE (scalar in [0,1]).
#' @param weightSrEstVq1 similar to `weightSrEstVp` but for the m1 selection
#' variability (scalar in [0,1]).
#' @param ... arguments for the function `MatchIt::matchit`.
estim_AV_seq_m1_m2 <- function(Yobs,
                               modes,
                               pi_mat,
                               p1,
                               p2,
                               sd1 = 0.0,
                               sd2 = 0.0,
                               cov12 = sd1 * sd2,
                               phi = rep(0.5, length(Yobs)),
                               eqPotOut = FALSE,
                               correcEW1 = TRUE,
                               correcEW2 = TRUE,
                               I,
                               Z = matrix(1.0, nrow = length(Yobs), ncol = 1L),
                               m1Only = FALSE,
                               m2Only = FALSE,
                               weightSrEstVp = 0.5,
                               weightSrEstVq1 = weightSrEstVp,
                               ...)
{

  pi <- diag(pi_mat)
  N <- length(pi)

  maskSr <- modes == "m1"
  maskSmr <- modes == "m2"
  maskSa <- maskSr | maskSmr
  na <- sum(maskSa) # Number of respondents
  maskm1Sa <- modes[maskSa] == "m1" # mask of the m1 respondents within Sa


  m1Only <- isTRUE(m1Only)
  m2Only <- !m1Only && isTRUE(m2Only)

  # In the case we focus on the estimator t_{e1} or t_{pq1} and
  # we do not want to use the estimated potential outcomes of y_1k on
  # the set of the m2 respondents (i.e. use only the date (z_k, y_1k) in Sr).
  # The weights vector phi is not considered in fine but set to 1 to use
  # the equations that take phi into account.
  if (m1Only)
  {
    phi <- rep(1.0, N)
    phibar <- numeric(N)
  }

  # Same for t_{e2} or t_{pq2}, setting phi to 0
  else if (m2Only)
  {
    phi <- numeric(N)
    phibar <- rep(1.0, N)
  }
  else
  {
    phibar <- 1.0 - phi

    eqPotOut <- isTRUE(eqPotOut)

    if (eqPotOut)
    {
      sd2 <- sd1
      cov12 <- sd1^2L
    }
  }


  piSa <- pi[maskSa]
  p1Sa <- p1[maskSa]
  p1barSa <- 1.0 - p1Sa

  # If the data of the m1 respondents will be used
  if (!m2Only)
  {
    phiSr <- phi[maskSr]
    Y1Sr <- Yobs[maskSr]
    phiY1Sr <- phiSr * Y1Sr

    piSr <- pi[maskSr]
    p1Sr <- p1[maskSr]

    piSr_mat <- pi_mat[maskSr, maskSr, drop = FALSE]

    covarpSr <- pi2_to_covarInc(piSr_mat)
  }
  else
    sd1 <- 0.0


  # If the data of the m2 respondents will be used
  if (!m1Only)
  {
    phiSmr <- phi[maskSmr]
    Y2Smr <- Yobs[maskSmr]
    phibarY2Smr <- phibar[maskSmr] * Y2Smr

    piSmr <- pi[maskSmr]
    p1Smr <- p1[maskSmr]
    p1barSmr <- 1.0 - p1Smr
    p2Smr <- p2[maskSmr]
    p2Sa <- p2[maskSa]

    piSmr_mat <- pi_mat[maskSmr, maskSmr, drop = FALSE]
    covarpSmr <- pi2_to_covarInc(piSmr_mat)
  }
  else
  {
    sd2 <- 0.0
    p2 <- NULL # The m2 selection probabilities will not be used
  }

  # Sampling design variability (p, S)
  # There is no correction needed for probabilities estimation
  # If we use only the m1 outcomes
  # (we do not weight by phi)
  if (m1Only)
  {
    correctedY1Srp <- (piSr * p1Sr)^-1L * Y1Sr
    varpEst <-
      t(correctedY1Srp) %*%
      (covarpSr / piSr_mat) %*%
      correctedY1Srp %>%
      as.numeric()
    varpEst <- varpEst +
      sum((1.0 - piSr) * piSr^-2L * p1Sr^-1L * (1.0 - p1Sr^-1L) *
            Y1Sr^2L)

    # Second order inclusion and covariance matrices are
    # directly deleted due to their important size (#Sr^2).
    # Due to the conditional independence within the selection mechanisms
    # they are not used after.
    # rm(piSr_mat, covarpSr)
  }
  # If we use only the m2 outcomes
  # (we do not weight by 1 - phi)
  else if (m2Only)
  {
    Y2Smr <- Yobs[maskSmr]
    correctedY2Smrp <- (piSmr * p1barSmr * p2Smr)^-1L * Y2Smr
    varpEst <-
      t(correctedY2Smrp) %*%
      (covarpSmr / piSmr_mat) %*%
      correctedY2Smrp %>%
      as.numeric()
    varpEst <- varpEst +
      sum((1.0 - piSmr) * piSmr^-2L * (p1barSmr * p2Smr)^-1L *
            (1.0 - (p1barSmr * p2Smr)^-1L) * Y2Smr^2L)

    # Unused matrices of size #Smr^2 deleted
    rm(piSmr_mat, covarpSmr)
  }
  # If we use the m1 and m2 outcomes
  # We will need here have to estimate the counterfactuals
  # on Sr (y_2k) and Smr (y_2k) by single matching with Z
  else
  {
    # If Y1 = Y2 on U we can set Y1Est = Y2Est on Sa equal
    # to the vector of observed values Yobs
    if (eqPotOut)
      Y1Est <- Y2Est <- Yobs[maskSa]
    else
    {
      YEst <- estim_counterfactuals(Yobs, modes, Z, ...)
      Y1Est <- YEst[maskSa, "Y1Est"]
      Y2Est <- YEst[maskSa, "Y2Est"]
      rm(YEst)
    }


    phiY1SaEst <- phi[maskSa] * Y1Est
    phibarY2SaEst <- phibar[maskSa] * Y2Est

    piSa_mat <- pi_mat[maskSa, maskSa, drop = FALSE]
    covarpSa <- pi2_to_covarInc(piSa_mat)


    # The variance estimator is the weighted average of
    # the variance estimator based on Sr and the one based on Smr
    # Sr
    PhiYSr <- phiY1SaEst[maskm1Sa] + phibarY2SaEst[maskm1Sa]
    correctedPhiYSrp <-
      (piSr * p1Sr)^-1L * PhiYSr

    varpSrEst <-
      t(correctedPhiYSrp) %*%
      (covarpSr / piSr_mat) %*%
      correctedPhiYSrp %>%
      as.numeric()
    varpSrEst <- varpSrEst +
      sum((1.0 - piSr) * piSr^-2L * p1Sr^-1L * (1.0 - p1Sr^-1L) *
            PhiYSr^2L)

    # Smr
    PhiYSmr <- phiY1SaEst[!maskm1Sa] + phibarY2SaEst[!maskm1Sa]
    correctedPhiYSmrp <-
      (piSmr * p1barSmr * p2Smr)^-1L * PhiYSmr

    varpSmrEst <-
      t(correctedPhiYSmrp) %*%
      (covarpSmr / piSmr_mat) %*%
      correctedPhiYSmrp %>%
      as.numeric()
    varpSmrEst <- varpSmrEst +
      sum((1.0 - piSmr) * piSmr^-2L * (p1barSmr * p2Smr)^-1L *
            (1.0 - (p1barSmr * p2Smr)^-1L) * PhiYSmr^2L)


    varpEst <- weightSrEstVp * varpSrEst + (1.0 - weightSrEstVp) * varpSmrEst

    rm(piSr_mat, covarpSr, piSmr_mat, covarpSmr)
  }


  # m1 selection variability (q1, R1)
  if  (correcEW1)
  {
    ZSr <- Z[maskSr, , drop = FALSE]
    estLinCoefq1 <- 0.0
    if (!m2Only)
    {
      estVPhi11 <- -crossprod(ZSr,
                              (piSr * p1Sr)^-1L *
                                phiY1Sr * (1.0 - p1Sr))


      estbPhi11 <- solve(Fisher_Information_Matrix(p1, Z, I)) %*%
        estVPhi11

      estLinCoefq1 <- estLinCoefq1 + estbPhi11
    }

    if (!m1Only)
    {
      ZSmr <- Z[maskSmr, , drop = FALSE]
      estVPhibar21 <- crossprod(ZSmr,
                                (piSmr * (1.0 - p1Smr) * p2Smr)^-1L *
                                  phibarY2Smr * p1Smr)

      estbPhibar21 <- solve(Fisher_Information_Matrix(p1, Z, I)) %*%
        estVPhibar21

      estLinCoefq1 <- estLinCoefq1 + estbPhibar21
    }
  }

  # If we use only the m1 outcomes
  # (we do not weight by phi)
  if (m1Only)
  {

    correctedY1Srq1 <- correctedY1Srp

    if (correcEW1)
      correctedY1Srq1 <- correctedY1Srq1 + ZSr %*% estLinCoefq1

    varq1Est <- sum((1.0 - p1Sr) * correctedY1Srq1^2L)
  }
  # If we use only the m2 outcomes
  # (we do not weight by phi)
  else if (m2Only)
  {
    correctedY2Smrq1 <- (piSmr * p1barSmr)^-1L * Y2Smr

    if (correcEW1)
      correctedY2Smrq1 <- correctedY2Smrq1 - ZSmr %*% estLinCoefq1

    varq1Est <- sum(p1Smr * p2Smr^-1L * correctedY2Smrq1^2L)
  }
  # If we use the m1 and m2 outcomes, we will use the estimations
  # of the counterfactuals. As the estimation of the variability due
  # to the sampling design, we consider the weighted average of
  # the variance estimator based on Sr and the one based on Smr
  else
  {
    correctedPhiYSaq1 <- piSa^-1L *
      (p1Sa^-1L * phiY1SaEst - p1barSa^-1L * phibarY2SaEst)

    if (correcEW1)
    {
      correctedPhiYSaq1 <-
        correctedPhiYSaq1 + Z[maskSa, , drop = FALSE] %*% estLinCoefq1
    }

    weightsq1 <- numeric(na)
    weightsq1[maskm1Sa] <- weightSrEstVq1 * (1.0 - p1Sr)
    weightsq1[!maskm1Sa] <- (1.0 - weightSrEstVq1) * (p1Smr * p2Smr^-1L)

    varq1Est <- sum(weightsq1 * correctedPhiYSaq1^2L)
  }


  # m2 selection variability (q2, R2)
  # There is variability only if we consider m2. We do not need
  # to use counterfactuals estimators
  if (!m1Only)
  {
    correctedY2Smrq2 <- (piSmr * p1barSmr * p2Smr)^-1L * phibarY2Smr
    if (correcEW2)
    {
      estVPhibar22 <- -crossprod(ZSmr,
                                 (piSmr * (1.0 - p1Smr) * p2Smr)^-1L *
                                   (1.0 - p2Smr) * phibarY2Smr)

      maskSm <- I & !maskSr

      estbPhibar22 <-
        solve(Fisher_Information_Matrix(p2, Z, maskSm)) %*% estVPhibar22

      correctedY2Smrq2 <- correctedY2Smrq2 + ZSmr %*% estbPhibar22
    }

    varq2Est <- sum((1.0 - p2Smr) * correctedY2Smrq2^2L)
  }
  else
    varq2Est <- 0.0


  # Potential outcomes variability
  if (sd1 <= 0.0 || sd2 <= 0.0)
    cov12 <- 0.0

  # We use the mean of the phi_k in case some of them are unknown
  # (if they are all known it becomes ||phi||_2^2, <phi, 1 - phi>
  # and ||1 - phi||_2^2)).
  varphiYEst <-
    mean(phi^2L, na.rm = TRUE) * N * sd1^2L +
    2.0 * mean(phi * phibar, na.rm = TRUE) * N * cov12 +
    mean(phibar^2L, na.rm = TRUE) * N * sd2^2L


  varpEst + varq1Est + varq2Est + varphiYEst
}

#' Estimation of the variance of t_e1 and the approximate variance
#' of t_pq1
#'
#'
#' HOMOSCEDASTICITY
#' @param Yobs vector of the observed outcomes. For the non-respondents
#' the value is not considered and therefore can be equal to NA
#' (numeric vector of size N the size of the population).
#' @param modes
#' @param pi_mat matrix containing the second-order inclusion probabilities
#' pi_kl. The values must be known on Sr^2 and can be equal to NA otherwise
#' (symmetric numeric matrix of order N).
#' @param p1 vector containing the true or estimated probabilities p_1k.
#' Can be equal to NA for the units that didn't answer by m1
#' (numeric vector of size N).
#' @param sd1 true or estimated value (from a consistent estimator if possible)
#' of the standard deviation of the m1 potential outcomes (positive scalar).
#' @param correcEW1 TRUE if we must apply a correction due to the
#' use of estimations of `p1` (boolean).
#' @param I vector containing the indicators I_k equal to TRUE if unit k
#' has been selected. Used only if `correcEW1` = TRUE
#' (logical vector of size N).
#' @param Z design matrix. Used only when `correcEW1` = TRUE
#' (numeric matrix with N rows and q columns).
#' @param m1Only TRUE if we estimate the variance of t_{pq1} and only with
#' the m1 respondents (boolean).
#' @param m2Only TRUE if we estimate the variance of t_{pq2} and only with
#' the m2 respondents (boolean).
#' @param ... arguments for the function `MatchIt::matchit`.
estim_AV_seq_m1 <- function(Yobs, modes,
                            pi_mat,
                            p1, sd1,
                            correcEW1 = TRUE,
                            I,
                            Z = matrix(1.0,
                                       nrow = length(Yobs), ncol = 1L),
                            ...)
{
  estim_AV_seq_m1_m2(Yobs, modes,
                     pi_mat,
                     p1, p2 = NULL,
                     sd1,
                     sd2 = 0.0,
                     cov12 = 0.0,
                     phi = rep(1.0, length(Yobs)),
                     eqPotOut = NULL,
                     correcEW1,
                     correcEW2 = FALSE,
                     I,
                     Z,
                     m1Only = TRUE,
                     ...)
}


estim_AV_seq_m2 <- function(Yobs, modes,
                            pi_mat,
                            p1, p2, sd2,
                            correcEW1 = TRUE,
                            correcEW2 = TRUE,
                            I,
                            Z = matrix(1.0,
                                       nrow = length(Yobs), ncol = 1L),
                            ...)
{
  estim_AV_seq_m1_m2(Yobs, modes,
                     pi_mat,
                     p1, p2,
                     sd1 = 0.0,
                     sd2,
                     cov12 = 0.0,
                     phi = numeric(length(Yobs)),
                     eqPotOut = NULL,
                     correcEW1,
                     correcEW2,
                     I,
                     Z,
                     m2Only = TRUE,
                     ...)
}

#' Is assumed that Y1 = Y2
estim_AV_seq_a <- function(Yobs,
                           responses,
                           pi_mat,
                           pa,
                           sd = NULL,
                           correcEW = TRUE,
                           I,
                           Z = matrix(1.0, nrow = length(Yobs), ncol = 1L),
                           covarp = NULL,
                           n = 1L,
                           ...)
{
  responses[responses == "r"] <- "m1"

  estim_AV_seq_m1(Yobs, responses,
                  pi_mat, pa, sd,
                  correcEW,
                  I,
                  Z,
                  covarp = covarp,
                  n = n,
                  ...)
}

var_expansion_m1_m2 <- function(Y1exp, Y2exp,
                                pi_mat,
                                p1, p2,
                                sd1 = 0.0, sd2 = 0.0, cov12 = sd1 * sd2,
                                phi = rep(0.5, length(Y1exp)))
{
  pi <- diag(pi_mat)

  p1bar <- 1.0 - p1

  phibar <- 1.0 - phi

  covarPi <- pi2_to_covarInc(pi_mat)


  # Sampling design variability (p, S)
  correctedYp <- pi^-1L * (phi * Y1exp + phibar * Y2exp)
  varp <- as.numeric(t(correctedYp) %*% covarPi %*% correctedYp) +
    sum(pi^-1L * (1.0 - pi) *
          (phi^2L * sd1^2L + phibar^2L * sd2^2L +
             2.0 * phi * phibar * cov12))

  # m1 selection variability (q1, R1)
  # with the independence between p and q1
  # Variance of the term phik y1k / p1k - phikbar y2k / p1kbar
  varCorrectedDifference <-
    (p1^-2L * phi^2L * sd1^2L +
       p1bar^-2L * phibar^2L * sd2^2L -
       2.0 * (p1 * p1bar)^-1L * phi * phibar * cov12)
  varq1 <-
    sum(pi^-1L * p1 * p1bar *
          (varCorrectedDifference +
             (p1^-1L * phi * Y1exp - p1bar^-1L * phibar^2L * Y2exp)^2L))


  # m2 selection variability (q2, R2)
  # with the conditional independence between p and q1 ; and q1 and q2
  if (any(phibar > 0.0))
  {
    varPhibarY2 <- phibar^2L * (sd2^2L + Y2exp^2L) # Variance of each phikbar y_2k
    varq2 <- sum((pi * p1bar * p2)^-1L * (1.0 - p2) * varPhibarY2)
  }
  else
    varq2 <- 0.0

  # potential outcomes variability
  varPhiY <-
    sum(phi^2L) * sd1^2L +
    2.0 * sum(phi * phibar) * cov12 +
    sum(phibar^2L) * sd2^2L

  varp + varq1 + varq2 + varPhiY
}

var_expansion_m1 <- function(Y1exp,
                             pi_mat,
                             p1,
                             sd1 = 0.0,
                             phi = rep(1.0, length(Y1exp)))
{
  var_expansion_m1_m2(Y1exp, Y2exp = numeric(length(Y1exp)),
                      pi_mat,
                      p1, p2 = NULL,
                      sd1, sd2 = 0.0, cov12 = 0.0,
                      phi)
}

var_expansion_m2 <- function(Y2exp,
                             pi_mat,
                             p1, p2,
                             sd2 = 0.0,
                             phi = rep(1.0, length(Y2exp)))
{
  var_expansion_m1_m2(Y1exp = numeric(length(Y2exp)), Y2exp,
                      pi_mat,
                      p1, p2,
                      sd1 = 0.0, sd2, cov12 = 0.0,
                      1.0 - phi)
}

var_expansion_a <- function(Yexp,
                            pi_mat,
                            pa,
                            sd = 0.0)
{
  var_expansion_m1(Yexp, pi_mat, pa, sd)
}
