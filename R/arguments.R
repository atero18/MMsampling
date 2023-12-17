#' @include tools.R

checkX <- function(X, N = nrow(X), p = ncol(X))
{
  checkTable(X, nrows = N, ncols = p)
}

#' @importFrom checkmate makeAssertionFunction
assertX <- makeAssertionFunction(checkX)

checkY <- function(Y, N = length(Y))
{
  checkNumericVector(Y, finite = TRUE,
                     any.missing = FALSE, all.missing = FALSE,
                     min.len = 1L, len = N)
}

#' @importFrom checkmate makeAssertionFunction
assertY <- makeAssertionFunction(checkY)

#' @importFrom checkmate checkSetEqual
checkY_tab <- function(Y_tab, N = nrow(Y_tab), K = ncol(Y_tab), M, modesRef = NULL)
{
  check <-
    checkTable(Y_tab,
               nrows = N, ncols = K,
               mode = "numeric", col.names = "named")

  if (!isTRUE(check))
    return(check)

  assertM(M)

  checkCols <- checkSetEqual(colnames(Y_tab), M)

  if (!isTRUE(checkCols))
    return(checkCols)

  if (any(is.infinite(Y_tab)))
    return("At least one value for Y is infinite")

  if (!is.null(modesRef))
  {
    assertModesRef(modesRef, M)

    Kref <- length(modesRef)
    if (Kref > 1L)
    {
      mat <-
        as.matrix(Y_tab[, modesRef[-1L]]) -
        matrix(Y_tab[, modesRef[1L]], nrow = N, ncol = Kref - 1L, byrow = FALSE)

      if (any(mat != 0.0))
        return("Some reference modes gave different results.")

    }
  }
  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertY_tab <- makeAssertionFunction(checkY_tab)

checkPhi_tab <- function(phi_tab, N = nrow(phi_tab), M)
{
  assertM(M)

  checkProbTable(checkProbTable, N, col.names = M, addNR = FALSE)
}

#' @importFrom checkmate makeAssertionFunction
assertPhi_tab <- makeAssertionFunction(checkPhi_tab)

#' @importFrom checkmate assertInteger assertVector
set_I <- function(I, N)
{
  assertVector(I, any.missing = FALSE, all.missing = FALSE, min.len = 1L)

  ERRORMESSAGE <- "I is supposed to be a vector with logical or integer values"
  if (is.logical(I))
  {
    stopifnot(length(I) == N)
    return(I)
  }
  else if (is.integer(I))
  {
    if (all(I %in% c(0L, 1L)))
    {
      stopifnot(length(I) == N)
      return(I == 1L)
    }

    else
    {
      assertInteger(I, unique = TRUE, lower = 1L, upper = N)
      mask <- logical(length(I))
      mask[I] <- TRUE
      return(mask)
    }
  }
  else
    abort(ERRORMESSAGE)
}


#' @importFrom purrr partial
#' @importFrom checkmate checkVector
checkModes <-
  partial(checkmate::checkVector,
          any.missing = FALSE, all.missing = FALSE,
          min.len = 1L, null.ok = FALSE)

#' @importFrom purrr partial
#' @importFrom checkmate assertVector
assertModes <- partial(checkmate::assertVector,
                       any.missing = FALSE, all.missing = FALSE,
                       min.len = 1L, null.ok = FALSE)

#' @importFrom purrr partial
#' @importFrom checkmate checkVector
checkM <-
  partial(checkmate::checkVector,
          any.missing = FALSE, all.missing = FALSE,
          min.len = 2L, null.ok = FALSE, unique = TRUE)


#' @importFrom purrr partial
#' @importFrom checkmate assertVector
assertM <- partial(checkmate::assertVector,
                   any.missing = FALSE, all.missing = FALSE,
                   min.len = 2L, null.ok = FALSE, unique = TRUE)

#' @importFrom purrr partial
checkMode <- partial(checkModes, len = 1L)

#' @importFrom checkmate makeAssertionFunction
assertMode <- makeAssertionFunction(checkMode)

#' @importFrom checkmate checkVector checkSubset
checkModesRef <- function(modesRef, allModes = NULL)
{
  check <- checkVector(modesRef,
                       any.missing = FALSE, all.missing = FALSE,
                       min.len = 1L, unique = TRUE,
                       null.ok = FALSE)

  if (!isTRUE(check))
    return(check)

  if (!is.null(allModes))
    return(checkSubset(modesRef, allModes))
  else
    return(TRUE)
}

#' @importFrom checkmate makeAssertionFunction
assertModesRef <- makeAssertionFunction(checkModesRef)

#' @importFrom purrr partial
checkMeans <- partial(checkNumericVector, finite = TRUE,
                      any.missing = FALSE, all.missing = FALSE,
                      null.ok = FALSE)

#' @importFrom checkmate makeAssertionFunction
assertMeans <- partial(assertNumericVector, finite = TRUE,
                       any.missing = FALSE, all.missing = FALSE,
                       null.ok = FALSE)


#' @importFrom purrr partial
checkVariances <- checkSDs <-
  partial(checkNumericVector, lower = 0.0, finite = TRUE,
          any.missing = FALSE, all.missing = FALSE,
          null.ok = FALSE)

#' @importFrom purrr partial
assertVariances <- assertSDs <-
  partial(assertNumericVector, lower = 0.0, finite = TRUE,
          any.missing = FALSE, all.missing = FALSE,
          null.ok = FALSE)


#' @importFrom checkmate assertFlag
#' @keywords internal
checkProbabilityVec <- function(x,
                                len = NULL, min.len = 1L, max.len = NULL,
                                unique = FALSE,
                                striclyPos = FALSE, striclyUnsure = FALSE)
{
  assertFlag(striclyPos)
  assertFlag(striclyUnsure)

  checkProb <- checkNumericVector(x,
                                  min.len = min.len, len = len,
                                  max.len = max.len,
                                  lower = 0.0, upper = 1.0,
                                  finite = TRUE,
                                  unique = unique)

  if (!isTRUE(checkProb))
    return(checkProb)

  if (striclyPos && any(x == 0.0))
    return("probabilities are supposed to be stricly positive")

  if (striclyUnsure && any(x == 1.0))
    return("probabilities are supposed to be inferior to one")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction
#' @keywords internal
assertProbabilityVec <- makeAssertionFunction(checkProbabilityVec)

checkProbabilityScalar <- partial(checkProbabilityVec,
                                  len = 1L, max.len = 1L,
                                  unique = FALSE)

assertProbabilityScalar <- partial(assertProbabilityVec,
                                   len = 1L, max.len = 1L,
                                   unique = FALSE)

#' @importFrom checkmate checkMatrix
checkCovarMat <- function(mat, p = nrow(mat))
{
  check <- checkMatrix(mat, mode = "numeric",
                       any.missing = FALSE, all.missing = FALSE,
                       min.rows = 1L, nrows = p, ncols = p,
                       null.ok = FALSE)

  if (!isTRUE(check))
    return(check)


  if (any(diag(mat) < 0.0) || any(is.infinite(mat)))
    return("covariance matrix is supposed to have only finite positive values")

  weigths <- (1.0 / diag(mat)) %>% sqrt()
  covarStand <- weigths %*% mat %*% weigths

  if (any(abs(covarStand) > 1.0))
    return("Some covariances are higher than the variances")

  if (!isSymmetric(mat))
    return("covariance matrix is supposed to be symmetric")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertCovarMat <- makeAssertionFunction(checkCovarMat)

checkPi2 <- function(mat, N = nrow(mat), pi = diag(mat))
{
  check <- checkMatrix(mat, mode = "numeric",
                       any.missing = FALSE, all.missing = FALSE,
                       min.rows = 1L, nrows = N, ncols = N,
                       null.ok = FALSE)

  if (!isTRUE(check))
    return(check)

  if (any(mat < 0.0 | mat > 1.0))
    return("Values must be in [0,1]")


  if (!all(diag(mat) == pi))
    return("diagonal of pi2 must be equal to pi1")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertPi2 <- makeAssertionFunction(checkPi2)

checkCorMat <- function(mat, p = nrow(mat))
{
  checkCovar <- checkCovarMat(mat, p)

  if (!isTRUE(checkCovar))
    return(checkCovar)

  if (!all(diag(mat) == 1.0))
    return("Some elements on the diagonal are not equal to 0")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertCorMat <- makeAssertionFunction(checkCorMat)

#' @importFrom checkmate checkMatrix
checkBiasesMes <- function(biases, N, K = ncol(biases))
{

  if (is.matrix(biases))
  {
    checkMat <- checkMatrix(biases, mode = "numeric",
                            any.missing = FALSE, all.missing = FALSE,
                            min.rows = 1L, nrows = N,
                            min.cols = 1L, ncols = K,
                            col.names = "named",
                            null.ok = FALSE)

    if (!isTRUE(checkMat))
      return(checkMat)

    if (all(biases == 0.0))
      return("There is no biased mode")

    if (any(is.infinite(biases)))
      return("Biases must be finite")
  }
  else if (is.vector(biases))
  {
    checkVec <- checkNumericVector(biases, len = N, finite = TRUE,
                                   any.missing = FALSE, all.missing = FALSE)

    if (!isTRUE(checkVec))
      return(checkVec)
  }
  else
    return("Biases must be in a matrix or a numerical vector")


  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertBiasesMes <- makeAssertionFunction(checkBiasesMes)

#' @importFrom checkmate checkVector
#' @importFrom purrr partial
checkNames <- partial(checkmate::checkVector,
                      any.missing = FALSE, all.missing = FALSE,
                      unique = TRUE, min.len = 1L, null.ok = FALSE)

#' @importFrom checkmate assertVector
#' @importFrom purrr partial
assertNames <- partial(checkmate::assertVector,
                       any.missing = FALSE, all.missing = FALSE,
                       unique = TRUE, min.len = 1L, null.ok = FALSE)

checkRGH <- checkNames
assertRGH <- assertNames

set_RGH <- function(RGH, N = length(RGH))
{
  if (is.null(RGH))
    RGH <- rep(1L, N)
  else
    assertRGH(len = N)

  RGH
}
