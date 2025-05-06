#' @export
estim_var_mean_phi1 <- function(Yobs,
                                modes, I,
                                piMat,
                                p1, Z,
                                phi = rep(1.0, length(Yobs)),
                                correcEstimWeights = TRUE)
{

  maskSr <- modes == "m1"

  weightsSr <- (diag(piMat)[maskSr] * p1[maskSr])^-1L
  Y1Sr <- Yobs[maskSr]
  # Estimator of the mean of y1 on the finite population U
  Hajekm1 <- sum(Y1Sr * weightsSr) / sum(weightsSr)


  errTerms <- rep(NA_real_, length(Yobs))
  errTerms[maskSr] <- Y1Sr - Hajekm1


  sumPhi <- sum(phi)

  # Returns the variance of Horvitz-Thompson total
  # (with or without estimated probabilities)of the variable
  # phi_k (y_1k - average on U of the y_1l) / sum of the phi_k on U
  estim_appr_var_seq_phi1(errTerms, modes, I, piMat,
                          p1, Z, phi,
                          estSD1 = 0.0,
                          correcEstimWeights) / sumPhi^2L
}

#' @export
estim_var_mean_phi2 <- function(Yobs,
                                modes, I,
                                piMat,
                                p1,
                                p2, Z,
                                phi = rep(1.0, length(Yobs)),
                                correcEstimWeights = TRUE)
{
  maskSmr <- modes == "m2"

  weightsSmr <-
    (diag(piMat)[maskSmr] * (1.0 - p1[maskSmr]) * p2[maskSmr])^-1L
  Y2mr <- Yobs[maskSmr]
  # Estimator of the mean of y2 on the finite population U
  Hajekm2 <- sum(Y2mr * weightsSmr) / sum(weightsSmr)

  errTerms <- rep(NA_real_, length(Yobs))
  errTerms[maskSmr] <- Y2mr - Hajekm2

  sumPhi <- sum(phi)

  # Returns the variance of Horvitz-Thompson total
  # (with or without estimated probabilities) of the variable
  # phi_k (y_2k - average on U of the y_2l) / sum of the phi_k on U
  estim_appr_var_seq_phi2(errTerms, modes, I, piMat,
                          p1, p2, Z, phi,
                          estSD2 = 0.0,
                          correcEstimWeights) / sumPhi^2L
}
