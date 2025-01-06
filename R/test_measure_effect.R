estim_var_centered_phiweb <- function(Yobs,
                                      modes, I,
                                      piMat,
                                      pq1Mat, Z,
                                      phi = rep(1.0, length(Yobs)),
                                      correcEstimWeights = TRUE)
{

  maskInt <- modes == "int"
  weightsWeb <- (diag(piMat)[maskInt] * diag(pq1Mat)[maskInt])^-1L

  # Estimator of the mean of phi * y1 on the finite population U
  HajekWeb <- sum(Yobs[maskInt] * weightsWeb) / sum(weightsWeb)



  errTerms <- rep(NA_real_, length(Yobs))
  errTerms[maskInt] <- Yobs[maskInt] - HajekWeb

  sumPhi <- sum(phi)

  # Returns the variance of Hovitz-Thompson total
  # (with or without estimated probabilities)of the variable
  # phi_k (y_1k - average on U of the y_1l) / sum of the phi_k on U
  estim_appr_var_seq_phi1(errTerms / sumPhi, modes, I, piMat,
                          pq1Mat, Z, phi,
                          correcEstimWeights)
}
