#' @include tools.R

#' @title Definition of the mode list `M`
#' @param M list of expected modes. Must be non-empty, composed of strings
#' and has to be made of unique values. (character vector)
#' @keywords internal
#' @name param_M
NULL

#' Definition of the sample indicator vector `I`
#' @param I vector of boolean indicating for each unit if it has been selected
#' (TRUE) or not (FALSE) in the sample `S`.
#' @keywords internal
#' @name param_I
NULL

#' Definition of the number of units `N`
#' @param N number of statistical units. (stricly positive integer)
#' @keywords internal
#' @name param_N
NULL


#' Check covariates
#' @inheritParams param_N
#' @keywords internal
checkX <- function(X, N = nrow(X), p = ncol(X))
{
  checkTable(X, nrows = N, ncols = p)
}

#' @importFrom checkmate makeAssertionFunction
assertX <- makeAssertionFunction(checkX)

## À vérifier
#' Check outcomes
#' @keywords internal
checkY <- function(Y, N = length(Y))
{
  checkNumericVector(Y, finite = TRUE,
                     any.missing = FALSE, all.missing = FALSE,
                     min.len = 1L, len = N)
}

#' @importFrom checkmate makeAssertionFunction
assertY <- makeAssertionFunction(checkY)


#' Check the table of outcomes.
#'
#' One row is a unit and one column a mode
#' @param Y_tab table containing outcome values for the different modes
#' available. May contain missing values. Available values must be finite.
#' @inheritParams param_N
#' @inheritParams param_M
#' @keywords internal
#' @importFrom checkmate checkSetEqual
checkY_tab <- function(Y_tab, N = nrow(Y_tab),
                       K = ncol(Y_tab), M = colnames(Y_tab))
{
  # Assert the name of the modes are correct
  assertM(M)



  # Check that Y_tab is a table of numeric values with
  # the right number of rows and columns
  check <-
    checkTable(Y_tab,
               nrows = N, ncols = K,
               mode = "numeric", col.names = "named")

  if (!isTRUE(check))
    return(check)



  checkCols <- checkSetEqual(colnames(Y_tab), M)

  if (!isTRUE(checkCols))
    return(checkCols)

  # Values are supposed to be finite
  if (any(is.infinite(Y_tab)))
    return("At least one value for Y is infinite")


  TRUE
}

#' @importFrom checkmate makeAssertionFunction
assertY_tab <- makeAssertionFunction(checkY_tab)


#' Check if the ponderation for each unit is a convex weight
#' @inheritParams param_N
#' @inheritParams param_M
#' @keywords internal
checkPhi_tab <- function(phi_tab, N = nrow(phi_tab), M)
{
  assertM(M)

  checkProbTable(checkProbTable, N, col.names = M, addNR = FALSE)
}

#' @importFrom checkmate makeAssertionFunction
assertPhi_tab <- makeAssertionFunction(checkPhi_tab)

#' Assert the `I` vector
#' Assert if the indicator vector `I` is correct, and return it as a boolean
#' vector if it is made of integers
#' @importFrom checkmate assertInteger assertVector
#' @inheritParams param_I
#' @inheritParams param_N
#' @keywords internal
set_I <- function(I, N)
{
  assertVector(I, any.missing = FALSE, all.missing = FALSE, min.len = 1L)

  assertCount(N, positive = TRUE)

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

#' Check if a vector supposedly containing mode affectation for each
#' unit is correct.
#'
#' Can contain NA values. We must have at least one unit
#' which answered.
#' @keywords internal
#' @name check_modes
NULL

#' @describeIn check_modes Do a check
#' @importFrom purrr partial
#' @importFrom checkmate checkVector
#' @keywords internal
checkModes <-
  partial(checkmate::checkVector,
          all.missing = FALSE,
          min.len = 1L, null.ok = FALSE)

#' @describeIn check_modes Do an assertion
#' @importFrom purrr partial
#' @importFrom checkmate assertVector
#' @keywords internal
assertModes <- partial(checkmate::assertVector,
                       all.missing = FALSE,
                       min.len = 1L, null.ok = FALSE)

#' @describeIn check_modes Do a check for a unique mode
#' @keywords internal
#' @importFrom purrr partial
checkMode <- partial(checkModes, len = 1L)

#' @describeIn check_modes Do an assertion for a unique mode
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertMode <- makeAssertionFunction(checkMode)


#' Check if a vector of values is a vector of unique modes
#' @keywords internal
#' @name check_M
NULL

#' @describeIn check_M check the vector M
#' @keywords internal
#' @importFrom purrr partial
#' @importFrom checkmate checkVector
checkM <-
  partial(checkmate::checkVector,
          any.missing = FALSE, all.missing = FALSE,
          min.len = 2L, null.ok = FALSE, unique = TRUE)


#' @describeIn check_M assert the vector M
#' @keywords internal
#' @importFrom purrr partial
#' @importFrom checkmate assertVector
assertM <- partial(checkmate::assertVector,
                   any.missing = FALSE, all.missing = FALSE,
                   min.len = 2L, null.ok = FALSE, unique = TRUE)

#' Check if a vector is a vector of probability
#' @keywords internal
#' @name check_probability
NULL

#' @describeIn check_probability check a vector of probability
#' @importFrom checkmate assertFlag
#' @keywords internal
checkProbabilityVec <- function(x,
                                len = NULL, min.len = 1L, max.len = NULL,
                                unique = FALSE,
                                striclyPos = FALSE, striclyUnsure = FALSE,
                                any.missing = FALSE)
{
  assertFlag(striclyPos)
  assertFlag(striclyUnsure)

  checkProb <- checkNumericVector(x,
                                  min.len = min.len, len = len,
                                  max.len = max.len,
                                  lower = 0.0, upper = 1.0,
                                  finite = TRUE,
                                  unique = unique,
                                  any.missing = any.missing)

  if (!isTRUE(checkProb))
    return(checkProb)

  if (striclyPos && any(x == 0.0, na.rm = TRUE))
    return("probabilities are supposed to be stricly positive")

  if (striclyUnsure && any(x == 1.0, na.rm = TRUE))
    return("probabilities are supposed to be inferior to one")

  TRUE
}

#' @describeIn check_probability assert a vector of probability
#' @importFrom checkmate makeAssertionFunction
#' @keywords internal
assertProbabilityVec <- makeAssertionFunction(checkProbabilityVec)

#' @describeIn check_probability check a scalar of probability
#' @importFrom checkmate assertFlag
#' @keywords internal
checkProbabilityScalar <- partial(checkProbabilityVec,
                                  len = 1L, max.len = 1L,
                                  unique = FALSE)

#' @describeIn check_probability assert a scalar of probability
#' @importFrom checkmate assertFlag
#' @keywords internal
assertProbabilityScalar <- partial(assertProbabilityVec,
                                   len = 1L, max.len = 1L,
                                   unique = FALSE)

#' Check related to variances
#' @keywords internal
#' @name check_variances
NULL

#' @describeIn check_variances Check variances and standard deviations
#' @keywords internal
#' @importFrom purrr partial
checkVariances <- checkSDs <-
  partial(checkNumericVector, lower = 0.0, finite = TRUE,
          any.missing = FALSE, all.missing = FALSE,
          null.ok = FALSE)

#' @describeIn check_variances Assert variances and standard deviations
#' @keywords internal
#' @importFrom purrr partial
assertVariances <- assertSDs <-
  partial(assertNumericVector, lower = 0.0, finite = TRUE,
          any.missing = FALSE, all.missing = FALSE,
          null.ok = FALSE)

#' @describeIn check_variances Check a variance-covariance matrix
#' @keywords internal
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

#' @describeIn check_variances Assert a variance-covariance matrix
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCovarMat <- makeAssertionFunction(checkCovarMat)

#' @describeIn check_variances Check a correlation matrix
#' @keywords internal
checkCorMat <- function(mat, p = nrow(mat))
{
  checkCovar <- checkCovarMat(mat, p)

  if (!isTRUE(checkCovar))
    return(checkCovar)

  if (!all(diag(mat) == 1.0))
    return("Some elements on the diagonal are not equal to 0")

  TRUE
}

#' @describeIn check_variances Assert a correlation matrix
#' @keywords internal
#' @importFrom checkmate makeAssertionFunction
assertCorMat <- makeAssertionFunction(checkCorMat)

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

set_RGH <- function(RGH, N = length(RGH))
{

  if (is.null(RGH) || all(is.na(RGH)))
  {
    RGH <- rep(1L, N)
  }
  else if (anyNA(RGH))
    RGH[is.na(RGH)] <- max(RGH) + 1L

  return(RGH)
}

checkRGH <- checkNames
assertRGH <- assertNames
