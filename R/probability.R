#' Transformation de matrices d'interaction
#'
#' Les fonctions suivantes transforment des matrices contenant des
#' informations sur les interactions entre variables (ex : corrélation)
#' en d'autres matrices.
#' @keywords internal
#' @name mats_interaction
NULL

#' @describeIn mats_interaction Transformation d'une matrice de corrélation
#' en une matrice de variance.
#' @param cor Matrice de corrélation (matrice)
#' @param sigmas Vecteur des écarts-types. Doit être de même taille
#' que l'ordre de `cor`.
#' @importFrom checkmate testScalar
#' @export
cor2cov <- function(cor, sigmas)
{
  if (testRownamed(cor))
    names <- rownames(cor)
  else
    names <- names(sigmas)

  if (testScalar(cor))
  {
    assertSDs(sigmas, len = 1L)

    cov <- matrix(sigmas, nrow = 1L, ncol = 1L)
  }
  else
  {
    ## assertCorMat(cor) Vérifier les arguments
    p <- nrow(cor)
    assertSDs(sigmas, len = p)

    cov <- cor * (sigmas %*% t(sigmas))
  }

  colnames(cov) <- rownames(cov) <- names

  return(cov)
}

#' Calculation of the inclusion covariance matrix using the
#' second order inclusion probabilities.
#'
#' @param pi_mat matrix containing the second order inclusion
#' probabilities pi_kl (symmetric numeric matrix of order N).
#' @return the corresponding covariance matrix
#' (symmetric numeric matrix of order N).
pi2_to_covarInc <- function(pi_mat)
{
  pi <- diag(pi_mat)
  pi_mat - outer(pi, pi)
}

#' @describeIn mats_interaction Calcul of the second order probability matrix
#' in case of independency
#' @param pi First order probability vector
#' @export
pi_to_pi2 <- function(pi)
{
  # pi2 <- pi %*% t(pi)
  pi2 <- outer(pi, pi)

  # Replacement of the diagonal of pi2 by pi
  # (this affection seems faster than diag<-)
  N <- length(pi)
  # Linear indices for diagonal elements
  idx <- 1L + seq(from = 0, to = N - 1L) * (N + 1L)
  pi2[idx] <- pi
  # diag(pi2) <- pi

  pi2
}

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
