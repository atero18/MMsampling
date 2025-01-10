estim_var_mean_phi1 <- function(Yobs,
                                    modes, I,
                                    piMat,
                                    pq1Mat, Z,
                                    phi = rep(1.0, length(Yobs)),
                                    correcEstimWeights = TRUE)
{

  maskInt <- modes == "int"
  weightsWeb <- (diag(piMat)[maskInt] * diag(pq1Mat)[maskInt])^-1L

  # Estimator of the mean of y1 on the finite population U
  HajekWeb <- sum(Yobs[maskInt] * weightsWeb) / sum(weightsWeb)



  errTerms <- rep(NA_real_, length(Yobs))
  errTerms[maskInt] <- Yobs[maskInt] - HajekWeb

  sumPhi <- sum(phi)

  # Returns the variance of Horvitz-Thompson total
  # (with or without estimated probabilities)of the variable
  # phi_k (y_1k - average on U of the y_1l) / sum of the phi_k on U
  estim_appr_var_seq_phi1(errTerms, modes, I, piMat,
                          pq1Mat, Z, phi,
                          sd1 = 0.0,
                          correcEstimWeights) / sumPhi^2L
}

estim_var_mean_phi2 <- function(Yobs,
                                    modes, I,
                                    piMat,
                                    pq1Mat,
                                    pq2Mat, Z,
                                    phi = rep(1.0, length(Yobs)),
                                    correcEstimWeights = TRUE)
{
  maskSmr <- modes == "tel"
  weightsTel <- (diag(piMat)[maskSmr] *
                   (1.0 - diag(pq1Mat)[maskSmr]) *
                   diag(pq2Mat)[maskSmr])^-1L

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
                          correcEstimWeights) / sumPhi^2L
}
