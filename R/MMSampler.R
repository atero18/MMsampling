#' @include MMData.R

# Data for samplers (MMSamplerData) ---------------------------------------

#' @importFrom methods setRefClass
MMSamplerData <-
  setRefClass("MMSamplerData",
              contains = "MMContext",
              fields = list(pi = "vector",
                            pi2 = "matrix",
                            incPlan = "function",
                            fixedSize = "logical",
                            probaModes = "matrix")
  )

MMSamplerData$methods(
  initialize = function(pi, pi2, incPlan, fixedSize = FALSE, probaModes, ...)
  {
    callSuper(...)
    .self$pi <- pi
    .self$pi2 <- pi2
    .self$incPlan <- incPlan
    .self$fixedSize <- fixedSize
    .self$probaModes <- probaModes
  }
)

#' @importFrom checkmate checkFlag checkSubset
valid_MMSamplerData <- function(object)
{
  checkPi <- checkProbabilityVec(object$pi, len = object$N)

  if (!isTRUE(checkPi))
    return(checkPi)

  if (!any(object$pi > 0.0))
    return("No unit can be select (all pi are equal to zero)")

  checkPi2 <- checkPi2(object$pi2, object$N, object$pi)

  if (!isTRUE(checkPi2))
    return(checkPi2)

  if (!is.function(object$incPlan))
    return("incPlan must be a function")


  checkFixedSize <- checkFlag(object$fixedSize)

  if (!isTRUE(checkFixedSize))
    return(checkFixedSize)

  checkProbaModes <- checkProbTable(object$probaModes, N = object$N)

  if (!isTRUE(checkProbaModes))
    return(checkProbaModes)

  checkProbaModes <- checkSubset(colnames(object$probaModes), c(object$M, "nr"))

  if (!isTRUE(checkProbaModes))
    return(checkProbaModes)

  TRUE
}

#' @importFrom methods setValidity
setValidity("MMSamplerData", valid_MMSamplerData)

MMSamplerData$methods(
  esp_n = function() sum(.self$pi),
  var_n = function()
  {
    if (fixedSize)
      0.0
    else if (is.null(pi2))
      NA_real_
    else
      pi2_to_covarInc(.self$pi2) %>% sum()
  })

MMSamplerData$methods(
  is_HT_unbiased = function() all(pi %*% .self$pModes(withNR = FALSE) > 0.0 |
                                    .self$phi == 0.0) ## À vérifier,

)

MMSamplerData$methods(
  covarInc = function() pi2_to_covarInc(.self$pi2),
  isVarUnbiased = function() .self$HTunbiased() && all(.self$pi2 > 0.0) ## À vérifier
)

MMSamplerData$methods(
  proba_resp_ref = function() .self$probaModes[, .self$modesRef, drop = TRUE]
)

MMSamplerData$methods(
  weights =
  function(modes, inv = TRUE)
  {
    HT_mm_weights(pi = .self$pi, probaMode = .self$pModes(modes), inv = inv)
  },
  weights_ref = function(inv = TRUE)
  {
    weights <- .self$pi *
      rowSums(.self$probaModes[, .self$modesRef, drop = FALSE])

    if (inv)
      weights
    else
      1.0 / weights

  },
  #' @importFrom checkmate assertSubset
  completeWeights =
  function(modes = NULL)
  {
    assertSubset(modes, .self$modes, empty.ok = TRUE)
    if (is.null(modes))
    {
      inclusionWeights <-
        matrix(1.0 / pi, nrow = .self$N, ncol = .self$K, byrow = FALSE)

      modeAndNRWeights <- .self$probaModes
    }
    else
    {
      stopifnot(length(modes) == .self$N)
      respondants <- modes != "nr"
      inclusionWeights <- numeric(.self$N)
      inclusionWeights[respondants] <- 1.0 / .self$pi[respondants]

      modeAndNRWeights <- numeric(.self$N)
      modeAndNRWeights[respondants] <-
        1.0 /
        get_value_by_mode(.self$probaModes[respondants, ],
                          modes[respondants])

    }

    inclusionWeights * modeAndNRWeights
  }
)

# Variance calculation ----------------------------------------------------

MMSamplerData$methods(
  var_plan = function(Y)
  {
    YdivPi <- Y / .self$pi
    varPlan <- (YdivPi %*% t(YdivPi)) * .self$covarInc

    sum(varPlan)
  },
  var_mode = function(Y)
  {
    varMode <-
      (Y^2L / .self$pi) %*%
      (.self$phi^2L / .self$pModes(withNR = FALSE) - 1.0)

    sum(varMode)
  },
  var_HT = function(Y)
  {
    assertY(Y, N = .self$N())
    .self$var_plan(Y) + .self$var_mode(Y)
  }
)


# Generation of samples (MMSampler) ---------------------------------------

#' @importFrom methods setRefClass
MMSampler <-
  setRefClass("MMSampler", contains = "MMSamplerData",
              methods = list(
                initialize = function(...) callSuper(...))
)

MMSampler$methods(
  #' @inheritParams seed
  samplingInclusion = function(seed = NULL)
  {
    fix_seed(seed)

    I <- .self$incPlan(pi = .self$pi, X = .self$X)

    if (all(I %in% c(0L, 1L)))
      I <- I == 1L

    set_I(I, N = .self$N)

    ##if (.self$fixedSize)
    ##  stopifnot(sum(I) == as.integer(.self$n))

    return(I)
  },
  #' @inheritParams seed
  samplingMode = function(I = rep(TRUE, .self$N),
                          C = logical(.self$N), seed = NULL)
  {
    fix_seed(seed)

    plan <- data.frame(I = I, pi = .self$pi, C = C)
    rownames(plan) <- .self$U

    sim_choix_mode_Poisson(plan = plan, p = .self$probaModes)
  },
  #' @inheritParams seed
  #' @importFrom rlang inject
  sampling = function(seed = NULL)
  {
    fix_seed(seed)

    I <- .self$samplingInclusion(seed = NULL)

    samplingModes <- .self$samplingMode(I, seed = NULL)

    Ytilde <- get_value_by_mode(.self$Y_tab, samplingModes$mode)

    inject(MMSample$new(I = I, R = samplingModes$R,
                        mode = samplingModes$mode,
                        Ytilde = Ytilde,
                        !!!as.list.environment(.self)))
  },
  #' @importFrom checkmate assertCount
  #' @importFrom parallel mclapply detectCores
  #' @inheritParams seed
  samplings = function(B = 1000L, seed = NULL, nbCores = detectCores() - 1L)
  {
    assertCount(B, positive = TRUE)
    assertCount(nbCores)
    fix_seed(seed)

    onWindows <- .Platform$OS.type == "windows"

    if (!onWindows)
      assertCount(nbCores, positive = TRUE)

    if (onWindows || nbCores <= 1L)
      plans <- replicate(B, .self$sampling(seed = NULL), simplify = FALSE)

    else
    {
      oldRandom <- .Random.seed

      RNGkind(kind = "L'Ecuyer-CMRG")

      plans <-
        mclapply(seq_len(B),
                 .self$sampling,
                 mc.cores = nbCores,
                 mc.set.seed = TRUE)

      # We set as random generator the old one
      # and ask for B random numbers to ensure
      # reproductibility and remove the deterministic
      # characteristic.
      .Random.seed <- oldRandom
      runif(B)
    }

    return(plans)
  }
)



# SRS sampling case (MMSRS) -----------------------------------------------


#' @importFrom methods setRefClass
MMSRS <-
  setRefClass("MMSRS",
              contains = "MMSampler",
              fields = list(n = "integer"),
              methods = list(
                #' @importFrom checkmate assertCount assertInt
                initialize = function(n, N, ...)
                {
                  assertCount(N, positive = TRUE)
                  assertInt(n, lower = 1L, upper = N)
                  .self$n <- as.integer(n)
                  piTemp <- rep(n / N, N)
                  pi2Temp <- pi2_SRS(n, N)
                  incPlanTemp <- function(...)
                  {
                    s <- logical(N)
                    s[sample.int(N, n)] <- TRUE
                    return(s)
                  }
                  callSuper(fixedSize = TRUE, pi = piTemp, pi2 = pi2Temp,
                                 incPlan = incPlanTemp, ...)
                },
                esp_n = function() .self$n,
                var_n = function() 0.0)
              )

#' @importFrom checkmate assertCount assertInt
#' @export
in_out_probs_SRS <- function(n, N, nbIncluded = 1L,
                            nbExcluded = 0L)
{
  assertCount(N, positive = TRUE)
  assertInt(n, lower = 1L, upper = N)
  assertInt(nbIncluded, lower = 0L, upper = n - 1L)
  assertInt(nbExcluded, lower = 0L, upper = N - n)

  if (n == N)
    return(1.0)

  if (nbExcluded == 0L)
    prod(seq(n - nbIncluded, n)) / prod(seq(N - nbIncluded, N))

  else
  {
    factorial(N - nbIncluded - nbExcluded) /
      (factorial(n - nbIncluded) * factorial(N - n - nbExcluded)) *
      (factorial(n) * factorial(N - n)) / factorial(N)
  }

}

inclusion_probs_SRS <-
  function(n, N, k) in_out_probs_SRS(n, N, k, nbExcluded = 0L)


#' @export
pi2_SRS <- function(n, N)
{
  pikl <- inclusion_probs_SRS(n, N, k = 2L)
  mat <- matrix(pikl, nrow = N, ncol = N)
  diag(mat) <- n / N
  return(mat)
}



# Poisson sampling case (MMPoisson) ---------------------------------------


#' @importFrom methods setRefClass
MMPoisson <-
  setRefClass("MMPoisson",
              contains = "MMSampler",
              methods = list(
                initialize = function(...)
                {
                  funPoisson <- function(...) runif(.self$N) <= .self$pi

                  callSuper(fixedSize = FALSE,
                            incPlan = funPoisson, ...)
                })
              )

#' @export
pi2_Poisson <- function(pi)
{
  assertProbabilityVec(pi)

  mat <- pi %*% t(pi)
  diag(mat) <- pi

  return(mat)
}
