#' Estimation of mode choising probabilities
#' @inheritParams sample
#' @keywords internal
#' @name proba_mode_estim
NULL

#' @describeIn proba_mode_estim Estimation by a multimodal point of view
#' @importFrom checkmate assertFlag
#' @importFrom VGAM vglm multinomial predictvglm
estim_response_prob_global <- function(I, mode, X,
                                       RGH = NULL, constRGH = TRUE,
                                       chosenOnly = FALSE)
{

  ## Ajouter vérification arguments ?
  modes <- as.factor(modes)
  M <- levels(modes)
  N <- length(mode)
  I <- set_I(I, N)

  assertX(X, N = N)

  assertFlag(constRGH)
  assertFlag(chosenOnly)

  RGH <- set_RGH(RGH)
  RGHNames <- unique(RGH)

  data <- cbind(mode = as.factor(modes), X)
  fittedProbs <- numeric(N)


  # For each RGH we do a multinomial regression on the whole set of modes
  for (group in RGHNames)
  {
    maskGroup <- RGH == group

    modelGroup <- vglm(mode ~ X, data = data,
                       family = multinomial, subset = I & maskGroup)

    ## À tester pour savoir quel résultat est rendu ; va bloquer car
    ## rendra sans doute une matrice. Est-ce que l'on obtient bien tous
    ## les modes pour tous les groupes ?
    fittedProbsGroup <-
      predictvglm(modelGroup,
                  newdata = X[maskGroup, , drop = FALSE],
                  type = "response")

    if (constRGH)
      fittedProbs[maskGroup] <- mean(fittedProbsGroup)
    else
      fittedProbs[maskGroup] <- fittedProbsGroup
  }

  ## Ça fait quoi ça ?
  if (chosenOnly)
    get_value_by_mode(fittedProbs, modes)
  else
    fittedProbs
}

#' @describeIn proba_mode_estim Estimation sequentialy with conditional
#' probabilities (interesting particularly with sequential mixed-modes protocoles)
#' @param orderModes the order modes are considered. The first mode will have
#' an absolute probability of selection, the second a probability a selection
#' conditionaly to the non-selection of the first one, etc.
#' @importFrom stats glm binomial predict.glm predict
#' @importFrom checkmate assertFlag assertVector
estim_response_prob_sequential <- function(I, X, modes, orderModes,
                                           RGH = NULL, constRGH = TRUE,
                                           link = "logit",
                                           chosenOnly = FALSE)
{

  ## Ajouter vérification arguments ?
  modes <- as.factor(modes)
  M <- levels(modes)
  N <- length(mode)
  I <- set_I(I, N)

  assertX(X, N = N)

  assertFlag(constRGH)
  assertFlag(chosenOnly)

  RGH <- set_RGH(RGH)

  RGHNames <- unique(RGH)

  # If the package `fastglm` is installed we will prefer to use it
  usefastglm <- requireNamespace("fastglm", quietly = TRUE)

  orderModes <- orderModes[tolower(orderModes) != "nr"]
  assertVector(orderModes,
               any.missing = FALSE, all.missing  = FALSE,
               min.len = 1L, unique = TRUE)

  assertSubset(orderModes, M)


  # Will contain the conditional probs for each mode
  conditionalProbs <- matrix(0.0,
                             nrow = N,
                             ncol = length(orderModes) + 1L)

  rownames(conditionalProbs) <- names(I)
  colnames(conditionalProbs) <- M

  # Will contain the unconditional probs for each mode
  unconditionalProbs <- conditionalProbs
  unconditionalVec <- rep(1.0, N)

  subset <- rep(TRUE, N) # TRUE if the individual is used for the regression
  # TRUE if the individual answered by the considered mode
  response <- rep(FALSE, N) # TRUE is the unit answered by the mode or already
  # answered

  if (usefastglm)
    predictFun <- predict
  else
    predictFun <- predict.glm

  # For each mode (in the correct order) we will do the regression, considering
  # only the people which didn't answer yet
  for (i in seq_len(orderModes))
  {
    mode <- orderModes[i]
    # We set to TRUE the unit in the subset that used this mode to answer
    response[subset] <- sample$mode[subset] == mode

    # One evaluation is made per mode, per RGH
    for (group in RGHNames)
    {
      maskGroup <- RGH == group
      nbElemsGroup <- sum(maskGroup)

      if (!any(subset & maskGroup))
      {
        glue("In group {group} nobody answered with mode {mode}. \\
             Probability fixed to zero") %>% inform()

        conditionalProbs[maskGroup, orderModes[i]] <- 0.0
        unconditionalProbs[maskGroup, orderModes[i]] <- 0.0
      }
      else if (all(maskGroup[subset]))
      {
        group("Everyone in {group} answered with mode {mode}. \\
              Probability fixed to one") %>% inform()

        conditionalProbs[maskGroup, orderModes[i]] <- 1.0
        unconditionalProbs[maskGroup, orderModes[i]] <-
          unconditionalProbs[maskGroup, i - 1L]

        unconditionalVec[maskGroup] <- 0.0
      }
      else
      {

        if (usefastglm)
        {
          subResponse <- response[subset & maskGroup] # nolint: object_usage_linter
          subX <- X[subset & maskGroup, , drop = FALSE] # nolint: object_usage_linter
          modelGroup <-
            fastglm::fastglm(subResponse ~ subX,
                             family = binomial, link = link)

          fittedProbsGroup <- modelGroup$fitted.values
        }
        else
        {
          modelGroup <- glm(response ~ X, family = binomial,
                            subset = subset & maskGroup, link = link)

          fittedProbsGroup <- predict.glm(modelGroup,
                                      newdata = X[maskGroup, , drop = FALSE],
                                      type = "response")
        }

        if (constRGH)
          fittedProbsGroup <- rep(mean(fittedProbsGroup), nbElemsGroup)

        conditionalProbs[maskGroup, orderModes[i]] <- fittedProbsGroup

        unconditionalProbs[maskGroup, mode] <-
          unconditionalVec[maskGroup] * fittedProbsGroup

        unconditionalVec[maskGroup] <-
          unconditionalVec[maskGroup] * (1.0 - fittedProbsGroup)

      }
    }

    subset[response] <- FALSE
  }

  # Non-response probabilities are calculated
  conditionalProbs[, "nr"] <-
    1.0 - conditionalProbs[, mode]

  unconditionalProbs$nr <- 1.0 - apply(unconditionalProbs, MARGIN = 1L, sum)

  if (chosenOnly)
  {
    unconditionalProbs <- get_value_by_mode(unconditionalProbs, sample$mode)
    conditionalProbs <- get_value_by_mode(conditionalProbs, sample$mode)
  }

  list(unconditional = unconditionalProbs,
       conditional = conditionalProbs)
}
