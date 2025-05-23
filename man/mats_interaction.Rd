% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probability.R
\name{mats_interaction}
\alias{mats_interaction}
\alias{cor2cov}
\alias{pi2_to_covarInc}
\alias{pi_to_pi2}
\title{Transformation de matrices d'interaction}
\usage{
cor2cov(cor, sigmas)

pi2_to_covarInc(pi2)

pi_to_pi2(pi)
}
\arguments{
\item{cor}{Matrice de corrélation (matrice)}

\item{sigmas}{Vecteur des écarts-types. Doit être de même taille
que l'ordre de \code{cor}.}

\item{pi2}{Matrice d'inclusion du second ordre}

\item{pi}{First order probability vector}
}
\description{
Les fonctions suivantes transforment des matrices contenant des
informations sur les interactions entre variables (ex : corrélation)
en d'autres matrices.
}
\section{Functions}{
\itemize{
\item \code{cor2cov()}: Transformation d'une matrice de corrélation
en une matrice de variance.

\item \code{pi2_to_covarInc()}: Calcul de la matrice de covariance
d'inclusion d'individus à partir des probabilités d'inclusion du
second ordre.

\item \code{pi_to_pi2()}: Calcul of the second order probability matrix
in case of independency

}}
\keyword{internal}
