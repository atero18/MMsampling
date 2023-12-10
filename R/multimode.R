

#' @importFrom checkmate assertScalar assertSubset assertFlag
#' @importFrom dplyr filter
#' @importFrom stats lm coef vcov qchisq
check_modes_equality_MCO <- function(sample, modes = NULL,
                                     alpha = 0.05, mergeEquivalents = TRUE)
{
  assertProbabilityVec(alpha, len = 1L, striclyPos = TRUE, striclyUnsure = TRUE)
  assertFlag(mergeEquivalents)

  masqueRepondants <- sample$R

  plan <- plan[masqueRepondants, ]
  X <- plan$X[masqueRepondants, , drop = FALSE]
  Ytilde <- plan$Ytilde()[masqueRepondants]

  if (!is.null(problem$modesRef))
  {
    masqueRepondantsRef <- masque_repondants_reference(plan, problem$modesRef)

    if (any(masqueRepondantsRef) && length(problem$modesRef) > 1L)
      plan$modes[plan$modes %in% problem$modesRef] <- problem$modesRef[1L]
  }

  usedModes <- modes_utilises(plan)
  if (!is.null(modes))
  {
    assertSubset(modes, usedModes)
    usedModes <- modes
  }

  Kused <- length(usedModes)

  p <- ncol(X)

  betas <- matrix(NA_real_, nrow = Kused, ncol = p + 1L)

  rownames(betas) <- usedModes

  if (!testColnamed(X))
    colnames(X) <- seq_len(p)

  colnames(betas) <- c("cst", colnames(X))

  varsBetas <- replicate(Kused, NULL, simplify = FALSE)
  names(varsBetas) <- usedModes

  uChi2 <- qchisq(1.0 - alpha, p + 1L)

  modifs <- character(0L)

  i <- 1L
  while (i < nrow(betas))
  {
    m1 <- rownames(betas)[i]
    masquem1 <- masque_repondants_mode(plan, m1)

    if (anyNA(betas[m1, ]))
    {
      model <- lm(Ytilde ~ X, subset = masquem1)
      betas[m1, ] <- coef(model)
      varsBetas[[m1]] <- vcov(model)
    }

    j <- i + 1L
    while (j <= nrow(betas))
    {
      m2 <- rownames(betas)[j]
      masquem2 <- masque_repondants_mode(plan, m2)

      if (anyNA(betas[m2, ]))
      {
        model <- lm(Ytilde ~ X, subset = masquem2)
        betas[m2, ] <- coef(model)
        varsBetas[[m2]] <- vcov(model)
      }

      deltasCoefs <- betas[m1, , drop = TRUE] - betas[m2, , drop = TRUE]
      cumulatedVar <- varsBetas[[m1]] + varsBetas[[m2]]

      chi2Value <- t(deltasCoefs) %*% solve(cumulatedVar) %*% deltasCoefs
      chi2Value <- as.vector(chi2Value)

      if (chi2Value <= uChi2)
      {
        plan$mode[masquem2] <- m1
        betas <- betas[-j, , drop = FALSE]
        problem$phi[, m1] <- problem$phi[, m1] + problem$phi[, m2]
        problem$phi[, m2] <- 0.0
        problem$Y_tab[masquem2, m1] <- Ytilde[masquem2]
        problem$Y_tab <- problem$Y_tab[, -j, drop = FALSE]
        problem$probaModes[, m1] <-
          problem$probaModes[, m1] + problem$probaModes[, m2]

        modifs <- c(modifs, paste0(m2, "->", m1))

        if (mergeEquivalents)
        {
          betas[m1, ] <- NA_real_
          i <- 0L
          j <- nrow(betas) + 1L
        }
      }
      else
        j <- j + 1L
    }
    i <- i + 1L
  }

  list(changes = modifs, remainingModes = rownames(betas), plan = plan)
}
