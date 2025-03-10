% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arguments.R
\name{check_modes}
\alias{check_modes}
\alias{checkModes}
\alias{assertModes}
\alias{checkMode}
\alias{assertMode}
\title{Check if a vector supposedly containing mode affectation for each
unit is correct.}
\usage{
checkModes(...)

assertModes(...)

checkMode(...)

assertMode(..., .var.name = checkmate::vname(...), add = NULL)
}
\description{
Can contain NA values. We must have at least one unit
which answered.
}
\section{Functions}{
\itemize{
\item \code{checkModes()}: Do a check

\item \code{assertModes()}: Do an assertion

\item \code{checkMode()}: Do a check for a unique mode

\item \code{assertMode()}: Do an assertion for a unique mode

}}
\keyword{internal}
