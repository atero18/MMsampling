#' @importFrom stats predict
estim_tot_imput <- function(X, model)
{
  Ypred <- predict(model, newdata = X)
  sum(Ypred)
}

#' @importFrom stats lm predict
MCO_Y <- function(sample)
{
  maskRespRef <- sample$respondents_ref()

  Ytilde <- sample$Yref() # nolint: object_usage_linter
  X <- sample$X %>% as.matrix() # nolint: object_usage_linter

  lm(Ytilde ~ X, subset = maskRespRef) %>%
    predict(newdata = sample$X)
}

#' @importFrom checkmate assertChoice
#' @importFrom stats coef
#' @export
estim_tot_MCO_Y <- function(plan, modesRef, X, Y, type = "pred")
{
  type <- tolower(type)
  assertChoice(type, c("pred", "beta"))
  modele <- MCO_Y(plan, modesRef, X, Y)

  if (type == "pred")
    estim_tot_imput(X, modele)

  else if (type == "beta")
  {
    betaSR <- coef(modele)
    (t(X) %*% X %*% betaSR)[1L]
  }
}

#' @importFrom checkmate assertVector
estim_tot_GH <- function(sample, groups)
{
  assertVector(groups,
               any.missing = FALSE, all.missing = FALSE,
               len = length(Y))

  if (length(unique(groups)) == 1L)
    return(estim_tot_1GH(plan))

  masqueRepRef <- sample$respondents_ref()

  groupSizes <- table(groups)
  groupNames <- names(groupSizes)
  groups <- groups[masque_repondants_reference]

  if (!all(groupNames %in% groups))
      stop("At least one group is not represented")

  crossprod(tapply(sample$Yref_resp(), groups, mean), groupSizes)
}

nearest_neighbor_reference <- function(sample, distance)
{
  N <- sample$N
  masqueRepRef <- sample$respondents_ref()
  masqueNonRep <- !masqueRepRef

  nn <- seq_len(N)

  XRR <- X[masqueRepRef, ]
  XnRR <- X[masqueNonRep, ]

  distances <- apply(XnRR, MARGIN = 1L,
                     function(xnRR) apply(XRR, MARGIN = 1L,
                                          function(xRR) distance(xnRR, XRR)))
  nn[masqueNonRep] <-
    which(masqueRepRef)[apply(distances, MARGIN = 1L, which.min)]

  nn
}

Y_matching <- function(sample, method)
{

  if (!requireNamespace("MatchIt", quietly = TRUE))
    abort("Package MatchIt is needed")

  masqueRepsRef <- sample$respondents_ref()
  class <- !masqueRepsRef %>% as.integer() # nolint: object_usage_linter

  X <- sample$X
  data <- cbind(class, X)
  resImp <- MatchIt::matchit(class ~ ., data = data, replace = TRUE)$match.matrix[, 1L]

  Ycf <- sample$Ytilde()
  Ycf[!masqueRepsRef] <- sample$Yref()[resImp]

  return(Ycf)
}
