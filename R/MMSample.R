#' @include MMSampler.R


#' @param sample object containing data about a sample and its context
#' @name sample
NULL

# Result of a sampler (MMSample) ----------------------------------------

#' @importFrom methods setRefClass
MMSample <-
  setRefClass(
    "MMSample",
    contains = "MMSamplerData",
    fields = list(I = "logical", R = "logical",
                  mode = "factor",
                  n = "integer", nbRespondents = "integer",
                  nbRef = "integer")
  )

MMSample$lock("n", "nbRespondents", "nbRef")

MMSample$methods(
  initialize = function(I, R, mode, ...)
  {
    callSuper(...)
    .self$I <- set_I(I, N = .self$N)
    .self$n <- sum(.self$I)
    .self$R <- R
    .self$nbRespondents <- sum(.self$R)
    .self$mode <- mode
    levels(.self$mode) <- c(.self$M, "nr", NA)
    .self$nbRef <- sum(.self$respondents_ref())
  }
)

valid_MMSample <- function(object)
{
  if (length(object$I) != object$N)
    return("I doesn't have the correct length")

  object$I[is.na(object$I) & object$pi == 0.0] <- FALSE
  object$I[is.na(object$I) & object$pi == 1.0] <- TRUE

  if (anyNA(object$I))
    return("I has missing value(s)")

  if (any(object$pi == 0.0 & object$I))
    return("At least one element is selected with a probability of zero")

  if (any(object$pi == 1.0 & !object$I))
    return("At least one element isn't selected with a probability of one")

  if (length(object$R) != object$N)
    return("R doesn't have the correct length")

  if (any(is.na(object$R) & !object$I))
    object$R[is.na(object$R) & !object$I] <- FALSE

  if (anyNA(object$R))
    return("R has missing value(s)")

  if (any(object$R > object$I))
    return("At least one unit answered without being in the sample")

  if (length(object$mode) != object$N)
    return("mode doesn't have the correct length")

  if (any(object$R & !object$mode %in% object$M))
    return("At least one element answered but its mode is 'unanswered' or unknown")

  object$mode[object$I & !object$R] <- "nr"
  object$mode[!object$I] <- NA

  if (!testSubset(levels(object$mode), c(object$M, "nr", NA)))
    return("Invalid values for mode")

  TRUE
}

#' @importFrom methods setValidity
setValidity("MMSample", valid_MMSample)

MMSample$methods(
  add_answers = function(Ytilde)
  {
    if (length(Ytilde) == .self$N)
      Ytilde <- Ytilde[.self$R]

    assertY(Ytilde, .self$n)

    modes <- .self$mode[.self$R]

    matrixPos <- cbind(which(.self$R), modes)
    .self$Y_tab[matrixPos] <- Ytilde
  }
)

#' @importFrom purrr compose
MMSample$methods(
  respondents_mode = function(mode) .self$R & .self$mode == mode,
  respondents_modes = function(modes) .self$R & .self$mode %in% modes,
  respondents_ref = function() .self$respondents_modes(.self$modesRef),
  respondents_with_biais = function() .self$R & !.self$mode %in% .self$modesRef,
  non_respondents = function() .self$I & !.self$R,
  used_modes = function() .self$mode[!.self$mode %in% c(NA, "nr")] %>% unique(),
  used_biased_modes = function() setdiff(.self$used_modes(), .self$modesRef)
)

MMSample$methods(
  Xm = function(mode) .self$X[.self$respondents_mode(mode), , drop = FALSE],
  Xref = function() .self$X[.self$respondents_ref(), , drop = FALSE],
  #' @importFrom checkmate assertFlag
  Ytilde = function(answersOnly = FALSE)
  {
    assertFlag(answersOnly)
    Ytilde <- get_value_by_mode(.self$Y_tab, .self$mode)

    if (answersOnly)
      Ytilde[!is.na(Ytilde)]
    else
      Ytilde
  },

  #' @importFrom checkmate assertFlag
  extract_Y = function(mask, answersOnly = TRUE)
  {
    assertFlag(answersOnly)

    if (answersOnly)
      .self$Ytilde()[mask]
    else
      ifelse(mask, .self$Ytilde(), NA_real_)
  },
  Ym_resp = function(mode, answersOnly = TRUE)
  {
    .self$extract_Y(.self$respondents_mode(mode), answersOnly)
  },
  Yref_resp = function(answersOnly = TRUE)
  {
    .self$extract_Y(.self$respondents_ref(), answersOnly)
  },
  HT_Ym = function(mode)
  {
    HT_Ym(self$pi, .self$Ytilde(), mode, .self$phi, .self$probaModes, .self$mode)
  }
)


# Diplay results of sampling ----------------------------------------------


#' @importFrom purrr compose
MMSample$methods(
  show = function()
  {
    infos <- data.frame(N = .self$N,
                        n = .self$n,
                        nEsp = .self$esp_n(),
                        nVar = .self$var_n(),
                        nbResp = .self$nbRespondents,
                        nbNonResp = .self$n - .self$nbRespondents,
                        nbRef = .self$nbRef,
                        YREsp_mean = mean(.self$Ytilde(), na.rm = TRUE))

    nbModes <- vapply(.self$M, compose(sum, .self$respondents_mode), integer(1L))
    names(nbModes) <- paste0("nb_", .self$M)

    modesResp <- .self$mode[.self$R] %>% droplevels("nr")
    meansYm <- tapply(.self$Ytilde()[.self$R], modesResp, mean)
    names(meansYm) <- paste0("Ymean_", names(meansYm))

    cbind(infos, t(nbModes), t(meansYm))
  }
)

#' @export
plot.MMSample <- function(x, ...)
{
  loadNamespace("igraph")

  nodes <- c(glue("U (N = {x$N})"), "S", x$M, "nr")
  colors <- c("purple", "red")

  nbPerMode <- table(x$mode)
  nbPerMode <- nbPerMode[c(x$M, "nr")]
  sizes <- c(x$N, x$n, nbPerMode)

  useLaTeX <- requireNamespace("latex2exp", quietly = TRUE)

  if (useLaTeX)
    symbolY <- "$\\bar{Y}$" # nolint: nonportable_parth_linter
  else
    symbolY <- "Ymean"

  meanYGlobal <- mean(x$Ytilde(), na.rm = TRUE) %>% ceiling()

  edgeLabels <- glue("{x$n} ({symbolY} = {meanYGlobal})")

  for (mode in x$M)
  {
    if (nbPerMode[mode] == 0L)
      edgeLabels <- c(edgeLabels, "0")
    else
    {
      meanYMode <- x$Ym_resp(mode) %>% mean() %>% ceiling()
      edgeLabels <- c(edgeLabels,
                      glue("{nbPerMode[mode]} ({symbolY} = {meanYMode})"))
    }
  }

  edgeLabels <- c(edgeLabels, as.character(sum(x$non_respondents())))

  colors <- c("purple", "red",
              ifelse(x$M %in% x$modesRef, "blue", "green"),
              "grey")

  edges <- c(1L, 2L)

  edges <- rep(2L, 2L * (x$K + 1L) + 2L)
  edges[1L] <- 1L
  edges[2L + 2L * seq_len(x$K + 1L)] <- seq_len(x$K + 1L) + 2L

  layout <- matrix(c(0.0, 1.0, # U
                     0.0, 0.5), # S
                   ncol = 2L, byrow = TRUE)

  xModes <- -1.0 + 2.0 * (0L:x$K) / x$K

  layout <- rbind(layout,
                  cbind(xModes, 0.0))

  graph <- igraph::make_graph(edges, directed = TRUE)
  igraph::V(graph)$name <- nodes
  #V(graph)$size <- log(sizes)
  igraph::V(graph)$color <- colors
  igraph::E(graph)$label <- edgeLabels
  ##
  # https://igraph.org/r/doc/plot.common.html

  igraph::plot.igraph(graph, layout = layout, label.dist = 1.0)
}
