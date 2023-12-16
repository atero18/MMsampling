#' @importFrom checkmate assertCount assertVector
#' @importFrom simstudy genCorMat
#' @inheritParams seed
#' @export
genVarMatrix <- function(p, seed = NULL, names = NULL,
                         sigmas = NULL, maxSigma = 10.0)
{
  assertCount(p, positive = TRUE)
  fix_seed(seed)

  if (is.null(sigmas))
  {
    assertNumericScalar(maxSigma, lower = 0.0, finite = TRUE)
    sigmas <- runif(p, 0.0, maxSigma)
  }
  else
    assertSDs(sigmas, len = p)

  if (!is.null(names))
  {
    assertVector(names,
                 len = p, unique = TRUE,
                 all.missing = FALSE, any.missing = FALSE)

    names(sigmas) <- names
  }

  if (p == 1L)
    corMat <- 1.0
  else
    corMat <- genCorMat(p)

  cor2cov(corMat, sigmas)
}

#' @importFrom checkmate assertCount testScalar
#' @importFrom simstudy genCorData
#' @importFrom stats cov2cor
#' @inheritParams seed
#' @export
lin_sim <- function(N = 1000L,
                    meanX = NULL, covarX = NULL,
                    betasY = NULL, varY = NULL,
                    Kbar = 1L, covarYMbar = NULL,
                    betasYMbar = NULL,
                    meanZ = NULL, covarZ = NULL,
                    gammaYMbar = NULL,
                    pX = nrow(covarX), pZ = length(meanZ),
                    seed = NULL)
{
  assertCount(N, positive = TRUE)
  assertCount(Kbar, positive = TRUE)
  assertCount(pX, positive = TRUE)
  assertCount(pZ)

  fix_seed(seed)

  if (is.null(meanX))
    meanX <- runif(pX, -10.0, +10.0)
  else
    assertMeans(meanX, len = pX)

  if (is.null(covarX))
    covarX <- genVarMatrix(pX)
  else
  {
    covarX <- as.matrix(covarX)
    assertCovarMat(covarX, pX)
  }

  sigmaX <- diag(covarX) %>% sqrt()

  nomsVarX <- paste0("X", seq_along(meanX))


  if (is.null(betasY))
    betasY <- runif(pX + 1L, -10.0, +10.0)
  else
  {
    assertNumericVector(betasY, finite = FALSE, len = pX + 1L,
                        any.missing = FALSE, all.missing = FALSE,
                        null.ok = FALSE)
  }

  names(betasY) <- nomsBetas <- c("cst", nomsVarX)



  if (is.null(varY))
    varY <- runif(1L, 0.0, 10.0)
  else
    assertVariances(varY, len = 1L)

  sigmaY <- sqrt(varY)

  nomsVarYMbar <- paste0("Y", seq_len(Kbar))
  nomsEpsYMbar <- paste0("eps", seq_len(Kbar))

  if (is.null(betasYMbar))
  {
    betasYMbar <-
      runif(Kbar * (pX + 1L), -10.0, +10.0) %>%
      matrix(nrow = pX + 1L, ncol = Kbar)
  }
  else if (Kbar == 1L && is.vector(betasYMbar))
  {
    assertNumericVector(betasYMbar, finite = TRUE,
                        any.missing = FALSE, all.missing = FALSE,
                        len = Kbar + 1L)
    betasYMbar <- matrix(betasYMbar, nrow = pX + 1L, ncol = 1L)
  }
  else
    assertTable(betasYMbar, ncols = Kbar, nrows = pX + 1L)

  rownames(betasYMbar) <- nomsBetas
  colnames(betasYMbar) <- nomsVarYMbar

  if (is.null(covarYMbar))
    covarYMbar <- genVarMatrix(Kbar)
  else
  {
    covarYMbar <- as.matrix(covarYMbar)
    assertCovarMat(covarX, Kbar)
  }

  sigmaYMbar <- diag(covarYMbar) %>% sqrt()

  if (pZ > 0L)
  {
    if (is.null(meanZ))
      meanZ <- runif(pZ, -10.0, +10.0)

    else
      assertMeans(meanZ, len = pZ)

    if (is.null(gammaYMbar))
    {
      gammaYMbar <-
        runif(Kbar * pZ, -10.0, +10.0) %>% matrix(nrow = pZ, ncol = Kbar)
    }

    nomsVarZ <- nomsGamma <- paste0("Z", seq_len(pZ))

    rownames(gammaYMbar) <- nomsVarZ
    colnames(gammaYMbar) <- nomsVarYMbar

    if (is.null(covarZ))
      covarZ <- genVarMatrix(pZ)
    else
      assertCovarMat(covarZ, pZ)

    sigmaZ <- diag(covarZ) %>% sqrt()

  }
  else
  {
    nomsVarZ <- nomsGamma <- character(0L)
    sigmaZ <- numeric(0L)
  }

  covar <- matrix(0.0,
                  ncol = pX + 1L + Kbar + pZ,
                  nrow = pX + 1L + Kbar + pZ)
  rownames(covar) <-
    colnames(covar) <- c(nomsVarX, "eps", nomsVarZ, nomsEpsYMbar)

  covar[nomsVarX, nomsVarX] <- covarX
  covar["eps", "eps"] <- varY
  covar[nomsEpsYMbar, nomsEpsYMbar] <- covarYMbar

  if (pZ > 0L)
    covar[nomsVarZ, nomsVarZ] <- covarZ

  sigmas <- c(sigmaX, sigmaY, sigmaYMbar, sigmaZ)

  corrMat <- cov2cor(covar)

  # We delete the potentials deltas making
  # the correlation matrix not symmetric
  if (!isSymmetric(corrMat))
    corrMat <- (corrMat +  t(corrMat)) / 2.0


  data <- genCorData(N,
                     mu = c(meanX, numeric(1L + Kbar), meanZ),
                     sigma = sigmas, corMatrix = corrMat, rho = 0.0,
                     cnames = c(nomsVarX, "eps",
                                nomsEpsYMbar, nomsGamma))[, -1L]
  data <- as.matrix(data)

  Y <- betasY[1L] + data[, nomsVarX, drop = FALSE] %*% betasY[-1L] + data[, "eps"]
  colnames(Y) <- "Y"

  data <- cbind(data, Y = Y)

  YMbar <-
    betasYMbar["cst", , drop = TRUE] +
    data[, nomsVarX, drop = FALSE] %*%
    betasYMbar[-1L, , drop = FALSE] +
    data[, nomsEpsYMbar]

  if (pZ > 0L)
    YMbar <- YMbar + data[, nomsVarZ, drop = FALSE] %*% as.matrix(gammaYMbar)

  colnames(YMbar) <- nomsVarYMbar

  data <- cbind(data, YMbar)

  coefsReg <- rbind(betasYMbar, gammaYMbar)
  coefsReg <- cbind(Y = c(betasY, numeric(pZ)), coefsReg)
  rownames(coefsReg) <- c(nomsBetas, nomsVarZ)


  probsM <- runif(N * (Kbar + 2L)) %>% matrix(nrow = N, ncol = Kbar + 2L)
  probsM <-
    apply(probsM, MARGIN = 1L, function(probs) probs / sum(probs)) %>% t()
  colnames(probsM) <- c("Y", nomsVarYMbar, "nr")

  phi_tab <- matrix(1.0 / (Kbar + 1L), ncol = Kbar + 1L, nrow = N)
  colnames(phi_tab) <- c("Y", nomsVarYMbar)

  list(covar = covar, data = data,
       phi_tab = phi_tab,
       coefs = coefsReg,
       probsM = probsM)
}

#' @importFrom purrr partial
sim_2_modes_1_X_1_Z <- partial(lin_sim, covarX = NULL,
                               Kbar = 1L, covarYMbar = 1L,
                               covarZ = NULL, pX = 1L, pZ = 1L)

#' @importFrom purrr compose
sim_internet_telephone <-
  compose(
    function(df)
    {
      colnames(df)[colnames(df) == c("Y", "Y1")] <- c("Int", "Tel")
      df
    },
    sim_2_modes_1_X_1_Z
  )

#' @importFrom checkmate assertString
#' @keywords internal
extract_data_from_sim <- function(data, toExtract = "X")
{

  if (is.list(data))
    data <- data$data

  assertString(toExtract, min.chars = 1L)
  toExtract <- toupper(toExtract)
  if (toExtract %in% c("X", "Y_TAB"))
  {
    pattern <- switch(toExtract,
                      X = "X",
                      Y_TAB = "Y")
    pattern <- paste0("^", pattern)
    correspondingColumns <- grep(pattern, colnames(data))
    return(data[, correspondingColumns, drop = FALSE])
  }


  NULL
}


#' @importFrom purrr partial
#' @keywords internal
extract_X_from_sim <-
  partial(extract_data_from_sim, toExtract = "X")

#' @importFrom purrr partial
#' @keywords internal
extract_Y_tab_from_sim <-
  partial(extract_data_from_sim, toExtract = "Y_tab")


#' @importFrom checkmate assertCount
#' @importFrom tibble tibble
#' @keywords internal
MC_mm <- function(sampler, B = 1000L,
                  samples = NULL, seed = NULL, ...)
{
  assertCount(B, positive = TRUE)

  fix_seed(seed)

  if (is.null(samples))
  {
    samples <- sampler$samplings(B)

    dataSamples <- lapply(samples, function(sample) sample$show())

    dataSamples <- do.call("rbind", dataSamples)
  }
  else
    dataSamples <- NULL


  estimFun <- function(sample) HT_mm_with_fitting(sample, ...)
  estimsTot <- vapply(samples, estimFun, numeric(1L))

  res <- tibble(itMC = seq_len(B),
                HT = estimsTot, totY = problem$totY,
                bias = estimsTot - problem$totY)

  if (!is.null(dataSamples))
    res <- res %>% cbind(dataSamples)

  return(res)
}

#' Évaluation des estimations via Monte-Carlo
#' @param estims vecteur numérique des estimations
#' @param target valeur cible / vraie valeur (réelle)
#' @name estim_MC
NULL

#' @describeIn estim_MC Estimation du biais
#' @export
bias_MC <- function(estims, target) mean(estims) - target

#' @describeIn estim_MC Variance empirique corrigée des estimations
#' @importFrom stats var
#' @export
var_MC <- function(estims) var(estims)

#' @describeIn estim_MC Écart-type empirique corrigé des estimations
#' @importFrom purrr compose
#' @export
sd_MC <- compose(sqrt, var_MC)

#' @describeIn estim_MC Erreur quadratique moyenne par rapport à l'estimateur
MSE_MC <- function(estims, target)
{
  bias_MC(estims, target)^2L + mean((estims - target)^2L)
}

#' @describeIn estim_MC Coefficient de variation (par défaut estimé via la
#' moyenne des estimations)
#' @export
CV_MC <- function(estims, target = mean(estims)) sd_MC(estims) / target

#' @importFrom tibble tibble
#' @importFrom tidyr expand_grid
#' @importFrom dplyr filter mutate distinct full_join select join_by inner_join
#' @importFrom dplyr group_by cur_group_id ungroup
#' @inheritParams seed
#' @export
make_grid_sim <- function(N = 10000L, seed = NULL)
{
  fix_seed(seed)

  phiList <- c(0.0, NA_real_)
  pXList <- seq_len(3L)
  pZList <- seq(0L, 2L)
  KRefBarList <- seq_len(2L)
  KRefList <- 1L ##KRefList <- seq_len(3L)
  propList <- c(5.0, 10.0, 20.0, 40.0) / 100.0
  samplerTypeList <- "SRS" ## c("SRS", "Poisson")
  checkEqualityList <- c("no", "MCO", "MCO_agreg")
  imputationList <- c("true_values", "nearest", "optimal", "genetic", "MCO")
  checkNullityBiasList <- c("no", "MCO_Ym", "MCO_CF")
  deltaEstimList <- c("CF", "MCO", "MCO_Ym", "G-COMP")
  deltaGHBList <- c("no", "k-means", "k-means_agreg")
  pmEstimList <- c("true_values", "multinomial")
  pmRGHList <- c("no", "k-means", "k-means_agreg")

  nbPis <- 2L


  probasGrid <- rep(N, nbPis) %>% sort()

  listPis <- list()

  for (i in seq_along(probasGrid))
    listPis[[i]] <- runif(probasGrid[i])

  probasGrid <- tibble(N = probasGrid, pi = listPis)


  phi <- N <- KRef <- Kbar <- pX <- pZ <- samplerType <- prop <-
    imputation <- deltaEstim <- n <- NULL
  # Creation of the hyperparameters grid
  grid <-
    expand_grid(phi = phiList,
                N = N,
                KRef = KRefList,
                Kbar = KRefBarList,
                pX = pXList,
                pZ = pZList,
                samplerType = samplerTypeList,
                prop = propList,
                imputation = imputationList,
                checkEquality = checkEqualityList,
                deltaEstim = deltaEstimList,
                checkNullityBias = checkNullityBiasList,
                pmEstim = pmEstimList)

  grid <- inner_join(grid, probasGrid, by = join_by(N))

  # Remove redundancies
  grid <- grid %>% filter(phi == 0.0 | KRef > 1L)
  grid <- grid %>% mutate(pi = ifelse(samplerType != "SRS", pi, NA))

  grid <- grid %>%
    mutate(n = ifelse(samplerType == "SRS",
                      ceiling(N * prop), NA_integer_)) %>%
    select(-prop)

  # With estimation for measure bias with G-COMP
  # There is no method of imputation used
  grid <-
    grid %>%
    mutate(imputation = ifelse(deltaEstim != "G-COMP",
                               imputation, NA_character_))

  grid <- distinct(grid)

  grid <-
    grid %>%
    group_by(N, KRef, Kbar, pX, pZ, phi) %>%
    mutate(idProblem = cur_group_id()) %>%
    group_by(n, .add = TRUE) %>%
    mutate(idSample = cur_group_id()) %>%
    ungroup()

  return(grid)

}

#' @importFrom purrr partial
make_grid_int_tel <- function(N = 10000L, seed = NULL)
{
  make_grid_sim(seed) %>%
    filter(pZ <= 1L, pX <= 2L, KRef == 1L, Kbar == 1L) # nolint: object_usage_linter
}

## À supprimer
sim_CH <- function(N = 10L, delta = 0.25, ratioVar = 1.0, seed = 123L)
{
    fix_seed(seed)

    assertNumericScalar(delta, finite = TRUE)
    assertNumericScalar(ratioVar, lower = 1L, finite = TRUE)

    meanX <- 1.0
    varX <- 1.5
    betasYtel <- c(5.0, 2.0) # Intercept + fonction linéaire X pour réponse téléphone
    betasYint <- c(5.0, 2.0 + delta) # intercept + fonction linéaire X pour réponse internet
    varYtel <- 1.0 # Variance des résidus pour le mode téléphone
    varYint <- varYtel * ratioVar # Variance des résidus pour le mode internet
    simulation <- lin_sim(N, pX = 1L, pZ = 0L, Kbar = 1L,
                          meanX = meanX, covarX = varX,
                          betasY = betasYtel, varY = varYtel,
                          betasYMbar = betasYint, covarYMbar = varYint,
                          seed = NULL)

    NOMSMODES <- c("tel", "int")
    colnames(simulation$data)[ncol(simulation$data) + -1:0] <- NOMSMODES
    colnames(simulation$phi_tab) <- NOMSMODES
    colnames(simulation$coefs) <- NOMSMODES
    colnames(simulation$probsM) <- c(NOMSMODES, "nr")

    return(simulation)

}


#' @importFrom dplyr select
#' @importFrom tidyr expand_grid
#' @inheritParams seed
#' @export
grid_sim <- function(B = 100L, grid = NULL, seed = NULL)
{
  fix_seed(seed)

  if (is.null(grid))
    grid <- make_grid_sim()

  sizeGrid <- nrow(grid)
  grid$idModel <- seq_len(sizeGrid)

  nbParameters <- ncol(grid)

  nbProblems <- max(grid$idProblem)

  dataProblems <- grid %>% select(idProblem, N, pX, pZ, phi)
  dataProblems$totY <- 0.0

  dataSamples <- grid %>% select(idSample, samplerType, n)

  nbSamples <- max(dataSamples$idSample)

  dataModes <- grid %>% select(idProblem)

  dataModels <- grid %>% select(idProblem, idSample, idModel)

  dataMC <-
    expand_grid(idModel = seq_len(nbSamples),
                itMC = seq_len(B))

  dataMCbis <- NULL

  dataMC$HT <- 0.0

  idModel <- 0.0 ##

  useProgressBar <- requireNamespace("progress", quietly = TRUE)

  if (useProgressBar)
  {
    bar <-
      progress::progress_bar$new(
        total = sizeGrid,
        format = paste0(":current / :total,:percent - ",
                        "elapsed: :elapsedfull | eta : :eta"))
    bar$tick(0L)
  }

  for (idProblem in unique(grid$idProblem))
  {
    # Select hyperparameters of the problem
    paramsProblem <- grid[which(grid$idProblem == idProblem)[1L], ]

    N <- paramsProblem$N
    pX <- paramsProblem$pX
    pZ <- paramsProblem$pZ
    phi <- paramsProblem$phi

    if (is.na(phi))
      phi <- NULL

    simulation <- lin_sim(N = N, pX = pX, pZ = pZ)

    X <- extract_data_from_sim(simulation, "X")

    Y_tab <- extract_data_from_sim(simulation, "Y")

    modes <- colnames(simulation$probsM)
    modesRef <- "Y"

    probaModes <- simulation$probsM

    problem <-
      MMProblem$new(X = X, Y_tab = Y_tab,
                    probaModes = simulation$probsM, modesRef = modesRef)

    dataProblems$totY[idProblem] <- problem$totY

    idSamples <-
      grid$idSample[grid$idProblem == idProblem] %>% unique()

    for (idSample in idSamples)
    {

      paramsSampling <-
        grid[grid$idSample == idSample, ][1L, ]

      n <- paramsSampling$n
      samplerType <- paramsSampling$samplerType
      pi <- paramsSampling$pi[[1L]]

      if (samplerType == "SRS")
        sampler <- MMSRS$new(X = X, n = n, N = N, probaModes = probaModes)

      else if (samplerType == "Poisson")
        sampler <- MMPoisson$new(X = X, N = N, pi = pi, probaModes = probaModes)


      plans <- sampler$samplings(B)

      infosPlans <-
        lapply(plans, about_plan, modes = modes, modesRef = modesRef)

      dataMCbis <- rbind(dataMCbis, do.call("rbind", infosPlans))

      idModels <-
        which(grid$idProblem == idProblem & grid$idSample == idSample)

      for (i in idModels)
      {
        paramsModel <- grid[i, ]

        imputation <- paramsModel$imputation
        deltaEstim <- paramsModel$deltaEstim
        pmEstim <- paramsModel$pmEstim

        if (pmEstim == "true_values")
          pmEstim <- NULL

        resMM <- MC_mm(problem = problem, B = B, plans = plans)
        dataMC$HT[seq((i - 1L) * B + 1L, i * B)] <- resMM$HT

        if (useProgressBar)
          bar$tick(1L)
      }
    }
  }

  list(problems = dataProblems,
       samples = dataSamples,
       models = dataModels,
       MC = cbind(dataMC, dataMCbis))
}
