#' @param seed graine aléatoire à utiliser pour la reproductibilité. Si vaut
#' `NULL` alors la graine n'est pas fixée.
#' @keywords internal
#' @importFrom checkmate assertCount
#' @name seed
fix_seed <- function(seed = NULL)
{
  assertCount(seed, null.ok = TRUE)

  if (!is.null(seed))
    set.seed(seed)
}

checkRownamed <- function(x)
{
  if (is.null(rownames(x)))
    return("argument doesn't have row names")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction makeTestFunction
assertRownamed <- makeAssertionFunction(checkRownamed)
testRownamed <- makeTestFunction(checkRownamed)

checkColnamed <- function(x)
{
  if (is.null(colnames(x)))
    return("argument doesn't have column names")

  TRUE
}

#' @importFrom checkmate makeAssertionFunction makeTestFunction
assertColnamed <- makeAssertionFunction(checkColnamed)
testColnamed <- makeTestFunction(checkColnamed)



# Creation of checks on numeric vectors -----------------------------------


#' @importFrom checkmate assertFlag checkVector checkNumeric
checkNumericVector <- function(x, lower = -Inf, upper = Inf, finite = FALSE,
                               any.missing = TRUE, all.missing = TRUE,
                               len = NULL, min.len = NULL, max.len = NULL,
                               unique = FALSE, sorted = FALSE,
                               names = NULL, typed.missing = FALSE,
                               null.ok = FALSE)
{
  assertFlag(null.ok)
  if (null.ok && is.null(x))
    return(TRUE)

  vectorCheck <- checkVector(x,
                             len = len,
                             min.len = min.len,
                             max.len = max.len,
                             null.ok = FALSE)

  if (!isTRUE(vectorCheck))
    return(vectorCheck)

  checkNumeric(x,
               any.missing = any.missing, all.missing = all.missing,
               lower = lower, upper = upper, finite = finite,
               unique = unique, sorted = sorted,
               names = names, typed.missing = typed.missing,
               null.ok = FALSE)

}

#' @importFrom checkmate makeAssertionFunction
assertNumericVector <- makeAssertionFunction(checkNumericVector)


#' @importFrom purrr partial
checkNumericScalar <- partial(checkNumericVector,
                              min.len = 1L, len = 1L, max.len = 1L,
                              any.missing = FALSE, all.missing = FALSE,
                              null.ok = FALSE, unique = FALSE, sorted = FALSE)

#' @importFrom purrr partial
assertNumericScalar <- partial(assertNumericVector,
                              min.len = 1L, len = 1L, max.len = 1L,
                              any.missing = FALSE, all.missing = FALSE,
                              null.ok = FALSE, unique = FALSE, sorted = FALSE)

#' @importFrom tibble is_tibble
#' @importFrom checkmate checkMatrix checkDataFrame checkTibble
checkTable <- function(table,
                       any.missing = FALSE, all.missing = FALSE,
                       min.rows = 1L, nrows = NULL, max.rows = NULL,
                       min.cols = 1L, ncols = NULL, max.cols = NULL,
                       null.ok = FALSE, mode = NULL, col.names = NULL)
{

  if (!is.null(col.names) && col.names != "named")
  {
    columnsNames <- col.names

    checkCols <- checkSetEqual(colnames(table), columnsNames)

    if (!isTRUE(checkCols))
      return(checkCols)

    col.names <- "named"
  }

  if (is_tibble(table))
    funCheck <- checkTibble

  else if (is.data.frame(table))
    funCheck <- checkDataFrame

  else if (is.matrix(table))
    funCheck <- checkMatrix

  else
    return("argument is supposed to be a Tibble, a data.frame or a matrix")

  if (is.null(mode))
    types <- character(0L) # nolint: object_usage_linter
  else
    types <- mode # nolint: object_usage_linter

  dataEnv <- as.list(environment())

  dataEnv <- dataEnv[intersect(names(dataEnv), names(formals(funCheck)))]

  formals(funCheck)[names(dataEnv)] <- dataEnv

  funCheck(table)
}

#' @importFrom checkmate makeAssertionFunction
assertTable <- makeAssertionFunction(checkTable)

#' @importFrom simstudy genCorMat
get_value_by_mode <- function(data, modes)
{

  if (is.vector(data))
    return(data)


  N <- nrow(data)
  res <- numeric(N)

  maskMissing <- is.na(modes) | modes == "nr"

  res[maskMissing] <- NA_real_

  if (any(maskMissing))
  {
    indexMatrix <-  cbind(which(!maskMissing, modes[!maskMissing]))

    res[!maskMissing] <-  data[indexMatrix]
  }

  return(res)
}

#' @importFrom checkmate assertFlag
data_proba_to_vec <- function(data, N = NULL, modes = NULL,
                              probVec = TRUE, addNR = TRUE,
                              striclyPos = FALSE, striclyUnsure = FALSE)
{

  if (is.matrix(data))
  {
    assertProbTable(data, probVec = probVec, addNR = addNR)
    assertModes(modes, len = N)
    assertSubset(modes, colnames(data))
    data <- get_value_by_mode(data, modes)
  }

  assertProbabilityVec(data, len = N, striclyPos = FALSE, striclyUnsure = FALSE)

  return(data)
}

#' @importFrom checkmate assertFlag assertCount
checkProbTable <- function(probsTable,
                           N = nrow(probsTable), p = ncol(probsTable),
                           probVec = TRUE, col.names = "named", addNR = TRUE)
{
  assertCount(N)
  assertFlag(probVec)

  checkType <- checkTable(probsTable,
                          mode = "numeric",
                          nrows = N, ncols = p,
                          col.names = col.names)

  if (!isTRUE(checkType))
    return(checkType)

  if (!testColnamed(probsTable))
    return("columns must be named")

  if (probVec)
  {
    assertFlag(addNR)
    if (addNR)
      probsTable <- add_nr_prob(probsTable, "nr")

    if (!all.equal(rowSums(probsTable) - 1.0, numeric(N), check.names = FALSE))
      return("Each row must be a probability")
  }


  if (!all(probsTable >= 0.0) || !all(probsTable <= 1.0))
    return("all values must be probabilities")

  TRUE

}

#' @importFrom checkmate makeAssertionFunction
assertProbTable <- makeAssertionFunction(checkProbTable)
