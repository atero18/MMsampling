#' Calculation of the inclusion covariance matrix using the
#' second order inclusion probabilities.
#'
#' @param pi_mat matrix containing the second-order inclusion
#' probabilities pi_kl (symmetric numeric matrix of order N).
#' @return the corresponding covariance matrix
#' (symmetric numeric matrix of order N).
#' @noRd
pi2_to_covarInc <- function(pi_mat)
{
  pi <- diag(pi_mat)
  pi_mat - outer(pi, pi)
}

## Useful ?
#' @importFrom tibble add_column
add_nr_prob <- function(probs, nrName = "nr")
{
  if (!nrName %in% colnames(probs))
  {
    probs <- cbind(probs, 1.0 - rowSums(probs))
    colnames(probs)[ncol(probs)] <- nrName
  }

  return(probs)
}
