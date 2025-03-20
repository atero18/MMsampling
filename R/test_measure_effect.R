#' @export
estim_var_mean_phi1 <- function(Yobs,
                                modes, I,
                                piMat,
                                pq1Mat, Z,
                                phi = rep(1.0, length(Yobs)),
                                correcEstimWeights = TRUE,
                                independenceq1 = NULL,
                                independenceq2 = NULL,
                                ...)
{

  maskSr <- modes == "m1"

  args <- list(...)

  if ("p1" %in% names(args))
    p1 <- args[["p1"]]
  else
    p1 <- diag(pq1Mat)

  weightsWeb <- (diag(piMat)[maskSr] * p1[maskSr])^-1L

  # Estimator of the mean of y1 on the finite population U
  HajekWeb <- sum(Yobs[maskSr] * weightsWeb) / sum(weightsWeb)


  errTerms <- rep(NA_real_, length(Yobs))
  errTerms[maskSr] <- Yobs[maskSr] - HajekWeb


  sumPhi <- sum(phi)

  # Returns the variance of Horvitz-Thompson total
  # (with or without estimated probabilities)of the variable
  # phi_k (y_1k - average on U of the y_1l) / sum of the phi_k on U
  estim_appr_var_seq_phi1(errTerms, modes, I, piMat,
                          pq1Mat, Z, phi,
                          sd1 = 0.0,
                          correcEstimWeights,
                          independenceq1, independenceq2, ...) / sumPhi^2L
}

#' @export
estim_var_mean_phi2 <- function(Yobs,
                                modes, I,
                                piMat,
                                pq1Mat,
                                pq2Mat, Z,
                                phi = rep(1.0, length(Yobs)),
                                correcEstimWeights = TRUE,
                                independenceq1 = NULL,
                                independenceq2 = NULL,
                                ...)
{
  maskSmr <- modes == "m2"

  args <- list(...)

  if ("p1" %in% names(args))
    p1 <- args[["p1"]]
  else
    p1 <- diag(pq1Mat)

  if ("p2" %in% names(args))
    p2 <- args[["p2"]]
  else
    p2 <- diag(pq2Mat)

  weightsTel <-
    (diag(piMat)[maskSmr] * (1.0 - p1[maskSmr]) * p2[maskSmr])^-1L

  # Estimator of the mean of y2 on the finite population U
  HajekTel <- sum(Yobs[maskSmr] * weightsTel) / sum(weightsTel)

  errTerms <- rep(NA_real_, length(Yobs))
  errTerms[maskSmr] <- Yobs[maskSmr] - HajekTel

  sumPhi <- sum(phi)

  # Returns the variance of Horvitz-Thompson total
  # (with or without estimated probabilities) of the variable
  # phi_k (y_2k - average on U of the y_2l) / sum of the phi_k on U
  estim_appr_var_seq_phi2(errTerms, modes, I, piMat,
                          pq1Mat, pq2Mat, Z, phi,
                          sd2 = 0.0,
                          correcEstimWeights,
                          independenceq1, independenceq2, ...) / sumPhi^2L
}
