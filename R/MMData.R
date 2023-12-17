
# Data multimode (MMContext) ----------------------------------------------


#' @importFrom methods setRefClass
MMContext <-
  setRefClass("MMContext",
              fields = list(X = "data.frame", phi_tab = "matrix",
                            modesRef = "vector", Y_tab = "matrix",
                            M = "vector", U = "vector", K = "integer",
                            N = "integer", p = "integer"),
              methods = list(Xm = function(masqueMode)
                             {
                               .self$X[masqueMode, , drop = FALSE]
                             },
                             #' @importFrom checkmate assertFlag
                             phi = function(modes, answersOnly = TRUE)
                             {
                               assertFlag(answersOnly)
                               # Those who didn't answered have a phi equal
                               # to zero
                               phi <- numeric(.self$N)
                               maskTrueModes <- !is.na(modes) & modes != "nr"
                               phi[maskTrueModes] <-
                                 get_value_by_mode(
                                   .self$phi_tab[maskTrueModes, , drop = FALSE],
                                   modes[maskTrueModes])

                               if (answersOnly)
                                 phi[maskTrueModes]
                               else
                                 phi
                             })
  )

MMContext$lock("U", "K", "N", "p")

MMContext$methods(
  #' @importFrom purrr compose
  initialize = function(X, phi_tab, M = NULL,
                        modesRef = NULL, Y_tab = NULL, ...)
  {
    ## À terminer

    # Check if there is a constant column in X. In that case the column is
    # removed because the constant is used by default in each model used
    # here.
    checkConstVar <- apply(X, MARGIN = 2L, compose(length, unique)) == 1L
    if (any(checkConstVar))
      .self$X <- as.data.frame(X[, !checkConstVar, drop = FALSE])
    else
      .self$X <- as.data.frame(X)

    .self$N <- nrow(X)
    .self$p <- ncol(X)

    .self$phi_tab <- as.matrix(phi_tab)

    if (is.null(M))
      .self$M <- colnames(phi_tab)
    else
      .self$M <- unique(M[tolower(M) != "nr"])

    .self$K <- length(.self$M)

    if (is.null(modesRef))
      .self$modesRef <- .self$M[1L]
    else
      .self$modesRef <- unique(modesRef)

    if (is.null(Y_tab))
    {
      .self$Y_tab <- matrix(NA_real_, nrow = .self$N, ncol = .self$K)
      colnames(.self$Y_tab) <- .self$M
    }

    else
      .self$Y_tab <- as.matrix(Y_tab)

    if (testRownamed(X))
      names <- rownames(X)
    else if (testRownamed(phi_tab))
      names <- rownames(phi_tab)
    else if (testRownamed(Y_tab))
      names <- rownames(Y_tab)
    else
      names <- seq_len(.self$N)

    .self$U <- rownames(.self$phi_tab) <- rownames(.self$X) <-
      rownames(.self$phi_tab) <- names

    if (!testColnamed(X))
      colnames(.self$X) <- paste0("X", seq_len(.self$p))

  }
)

#' @importFrom checkmate checkSubset
validate_MMContext <- function(object)
{
  check <- checkM(object$M)

  if (!isTRUE(check))
    return(check)

  check <- checkSubset(object$modesRef, object$M)

  if (!isTRUE(check))
    return(check)

  if (object$K == 0L)
    return("There is no mode")

  if (length(object$modesRef) == 0L)
    return("There is no reference mode")

  if (length(object$modesRef) == object$K)
    return("There is no biased mode")

  check <- checkX(object$X, N = object$N, p = object$p)

  if (!isTRUE(check))
    return(check)

  checkPhi <- assertPhi_tab(object$phi_tab, N = object$N, M = object$M)

  if (!isTRUE(checkPhi))
    return(check)
}

#' @importFrom methods setValidity
setValidity("MMContext", validate_MMContext)

MMContext$methods(
  Xm = function(masqueMode) .self$X[masqueMode, , drop = FALSE],
  #' @importFrom checkmate assertFlag
  phi = function(modes, answersOnly = TRUE)
  {
    assertFlag(answersOnly)
    # Those who didn't answered have a phi equal
    # to zero
    phi <- numeric(.self$N)
    maskTrueModes <- !is.na(modes) & modes != "nr"
    phi[maskTrueModes] <-
      get_value_by_mode(
        .self$phi_tab[maskTrueModes, , drop = FALSE],
        modes[maskTrueModes])

    if (answersOnly)
      phi[maskTrueModes]
    else
      phi
  }
)

MMContext$methods(
  Ym = function(m) .self$Y_tab[, m],
  Yref = function() .self$Ym(.self$modesRef[1L]),
  yPhi = function() .self$Y_tab * .self$phi_tab,
  yPhim = function(mY, mPhi = mY) .self$Y_tab[, mY] * .self$phi_tab[, mPhi],
  #' @importFrom checkmate testChoice
  mesBias = function(m = NULL)
  {
    if (is.null(m))
    {
      .self$Y_tab -
        matrix(.self$Y,
               nrow = .self$N, ncol = .self$K,
               byrow = FALSE)
    }
    else if (m %in% .self$modesRef)
      numeric(.self$N)
    else
      .self$Ym(m) - .self$Yref()

  },
  totsY = function() colSums(.self$Y_tab),
  totYm = function(m) sum(.self$Y_tab[, m]),
  totYref = function() .self$totYm(.self$modesRef[1L])
)
