#' Estimation of mode choising probabilities
#' @inheritParams sample
#' @name proba_mode_estim
NULL

#' @describeIn proba_mode_estim Estimation by a multimodal point of view
#' @importFrom checkmate assertFlag
#' @importFrom VGAM vglm multinomial predictvglm
estim_response_prob_global <- function(sample,
                                       RGH = NULL, constRGH = TRUE,
                                       chosenOnly = FALSE)
{

  N <- sample$N
  masqueEchant <- sample$I
  modes <- sample$mode

  assertFlag(constRGH)
  assertFlag(chosenOnly)

  set_RGH(RGH)

  X <- sample$X
  data <- cbind(modes, X)
  fittedProbs <- numeric(N)
  RGHNames <- unique(RGH)
  for (group in RGHNames)
  {
    maskGroup <- RGH == group

    modelGroup <- vglm(modes ~ X, data = data,
                       family = multinomial, subset = masqueEchant & maskGroup)

    fittedProbsGroup <-
      predictvglm(modelGroup, newdata = X[maskGroup, ], type = "response")

    if (constRGH)
      fittedProbs[maskGroup] <- mean(fittedProbsGroup)
    else
      fittedProbs[maskGroup] <- fittedProbsGroup
  }

  if (chosenOnly)
    get_value_by_mode(fittedProbs, modes)
  else
    fittedProbs
}

#' @describeIn proba_mode_estim Estimation sequentialy with conditional
#' probabilities (interesting particularly with sequential multimode protocoles)
#' @param orderModes the order modes are considered. The first mode will have
#' an absolute probability of selection, the second a probability a selection
#' conditionaly to the non-selection of the first one, etc.
#' @importFrom stats glm binomial predict.glm predict
#' @importFrom checkmate assertFlag assertVector
estim_response_prob_sequential <- function(sample, orderModes,
                                           RGH = NULL, constRGH = TRUE,
                                           link = "logit",
                                           chosenOnly = FALSE)
{
  assertFlag(constRGH)
  assertFlag(chosenOnly)

  set_RGH(RGH)

  usefastglm <- requireNamespace("fastglm", quietly = TRUE)

  orderModes <- orderModes[tolower(orderModes) != "nr"]
  assertVector(orderModes,
               any.missing = FALSE, all.missing  = FALSE,
               min.len = 1L, unique = TRUE)

  N <- sample$N

  conditionalProbs <- matrix(0.0,
                             nrow = N,
                             ncol = length(orderModes) + 1L)

  rownames(conditionalProbs) <- sample$U
  colnames(conditionalProbs) <- c(orderModes, "nr")

  unconditionalProbs <- conditionalProbs

  X <- sample$X

  unconditionalVec <- rep(1.0, N)

  subset <- rep(TRUE, N)

  if (usefastglm)
    predictFun <- predict
  else
    predictFun <- predict.glm

  for (i in seq_len(orderModes))
  {
    mode <- orderModes[i]
    response[subset] <- sample$mode[subset] == mode

    for (group in RGH)
    {
      maskGroup <- RGH == group
      nbElemsGroup <- sum(maskGroup)

      if (!any(subset & maskGroup))
      {
        glue("In group {group} nobody answered with mode {mode}. \\
             Probability fixed to zero") %>% inform()
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
        }
        else
        {
          modelGroup <- glm(response ~ X, family = binomial,
                            subset = subset & maskGroup, link = link)
        }

        fittedProbsGroup <- predictFun(modelGroup,
                                       newdata = X[maskGroup, , drop = FALSE],
                                       type = "response")

        if (constRGH)
          fittedProbsGroup <- rep(mean(fittedProbsGroup), nbElemsGroup)

        conditionalProbs[maskGroup, orderModes[i]] <- fittedProbsGroup

        unconditionalProbs[maskGroup, mode] <-
          unconditionalVec[maskGroup] * fittedProbsGroup

        unconditionalVec[maskGroup] <-
          unconditionalVec[maskGroup] * (1.0 - fittedProbsGroup)

      }
    }

    subset <- response == 0.0
  }

  conditionalProbs[, "nr"] <-
    1.0 - unconditionalProbs[, ncol(unconditionalProbs)]

  if (chosenOnly)
  {
    unconditionalProbs <- get_value_by_mode(unconditionalProbs, sample$mode)
    conditionalProbs <- get_value_by_mode(conditionalProbs, sample$mode)
  }

  list(unconditional = unconditionalProbs,
       conditional = conditionalProbs)
}
