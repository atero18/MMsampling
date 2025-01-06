estim_var_centered_phiweb_mean <- function(Yobs,
                                           modes, piMat,
                                           p1Mat,
                                           phi = rep(1.0, length(Yobs)),
                                           correcEstimWeights = TRUE)
{
  maskWeb <- modes == "web"

  sumPhi <- sum(phi)

  pi <- diag(piMat)
  weights <- (pi * diag(p1Mat))^-1L
  weightsWeb <- weights[maskWeb]

  # Estimator of the mean of phi * y1 on the finite population U
  HajekPhiWeb <- sum(phi[maskWeb] * Yobs[maskWeb] * weightsWeb) /
    sum(weightsWeb)

  errTermsWeb <- Yobs[maskWeb] - HajekPhiWeb

  covarPiWeb <- pi2_to_covarInc(piMat[maskWeb, maskWeb, drop = FALSE])
  varSEst <-
    t(errTermsWeb / pi[maskWeb]) %*%
    covarPiWeb /
    (piMat[maskWeb, maskWeb, drop = FALSE] *
       p1Mat[maskWeb, maskWeb, drop = FALSE])^-1L %*%
    (errTermsWeb / pi[maskWeb]) / sumPhi^2L


  covarq1Web <- pi2_to_covarInc(p1Mat[maskWeb, maskWeb, drop = FALSE])
  varq1Est <-
    t(errTermsWeb / weightsWeb) *
    piMat[maskWeb, maskWeb, drop = FALSE] *
    covarq1Web *
    (errTermsWeb / weightsWeb) / sumPhi^2L

  varSEst + varq1Est
}
