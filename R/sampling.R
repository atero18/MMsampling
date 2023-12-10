#' @inheritParams seed
ind_controle <- function(plan, nu, seed = NULL)
{
  N <- taille_pop(plan)
  n <- taille_echant(plan)

  if (n == 0L)
  {
    plan$mode <- logical(N)
    return(plan)
  }

  assertProbTable(nu, N)
  nu <- add_nr_prob(nu, "nC")

  modesC <- colnames(nu)

  fix_seed(seed)

  C <- logical(N)
  modeC <- apply(nu[I, ], 1L,
                 function(nus) sample(modesC, size = 1L, prob = nus))

  C[I] <- modeC != "nC"

  mode <- rep(NA, N)
  mode[C] <- modeC[mode != "nC"]

  plan$mode <- mode

  return(plan)
}

#' Simulation du choix des modes sous un modèle de Poisson
#' @name sampling_mode
NULL

#' @describeIn sampling_mode Simulation d'un choix de mode pour chaque
#' individu de manière indépendante
#' @inheritParams seed
#' @export
sim_ind_MM <- function(p, seed = NULL)
{
  p <- add_nr_prob(p)
  assertProbTable(p)

  N <- nrow(p)
  M <- colnames(p)

  # Special case when there is only two choices : respond or not respond
  # (Use of Poisson sampling)
  if (length(M) == 2L)
  {
    chosenModes <- M[runif(N) <= p[, 2L] + 1L]
  }
  else
  {
    chosenModes <- apply(p, MARGIN = 1L,
                         function(prob) sample(M, size = 1L, prob = prob))
  }

  factor(chosenModes, levels = c(M, NA))
}

#' @describeIn sampling_mode Simulation de sélection avec probabilités égales
#' pour tous les modes et tous les individus
#' @param N Taille de la population
#' @param modes Liste des différents modes disponibles (vecteur). Possibilité
#' que ce soit un entier, auquel cas on considérera autant de modes, plus
#' éventuellement de la non-réponse.
#' @importFrom checkmate assertCount testInt assertVector assertFlag
#' @export
sim_ind_MM_eq_prob <- function(N, modes, addNR = TRUE)
{
  assertCount(N, positive = TRUE)
  assertFlag(addNR)

  if (testInt(modes, lower = 2L - addNR))
    modes <- seq_len(modes) %>% as.character()
  else
  {
    assertVector(modes,
                 any.missing = FALSE, all.missing = FALSE,
                 min.len = 2L - addNR, unique = TRUE)
  }

  if (addNR && !"nr" %in% modes)
    modes <- c(modes, "nr")

  res <- sample(modes, size = N, replace = TRUE)
  factor(res, levels = modes)
}

#' @importFrom stats runif
sim_choix_mode_Poisson <- function(plan, p, q = NULL, seed = NULL)
{
  N <- taille_pop(plan)

  fix_seed(seed)

  plan$mode <- factor(NA,
                      levels = c(colnames(p),
                                 colnames(q), "nr", NA) %>% unique())

  if (!"C" %in% colnames(plan))
    plan$C <- FALSE

  echantControle <- masque_controle(plan)
  echantNonControle <- masque_echant(plan) & !echantControle

  R <- logical(N)

  if (any(echantControle))
  {
    assertProbTable(q, N)
    q <- add_nr_prob(q, "nr")

    R[echantControle] <-
    vapply(which(echantControle),
           function(i) runif(1L) <= q[i, plan$mode[i]], logical(1L))
  }

  if (any(echantNonControle))
  {

    checkProbTable(p, N)
    ##assertProbTable(p, N)
    p <- add_nr_prob(p, "nr")

    plan$mode[echantNonControle] <-
      sim_ind_MM(p[echantNonControle, , drop = FALSE])

    R[echantNonControle] <- plan$mode[echantNonControle] != "nr"
  }

  plan$R <- R

  return(plan)
}
