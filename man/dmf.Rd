% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmf.R
\name{dmf}
\alias{dmf}
\title{Perform deviance matrix factorization.}
\usage{
dmf(
  x,
  family = gaussian(),
  rank = ncol(x),
  weights = 1,
  offset = zeros(x),
  control = glm.control(epsilon = 1e-06, maxit = 100)
)
}
\arguments{
\item{x}{Input matrix to be factorized.}

\item{family}{Family object to specify deviance loss.}

\item{rank}{Decomposition rank.}

\item{weights}{Entrywise weight.}

\item{offset}{Entrywise offset.}

\item{control}{Algorithm control parameters (see \code{glm.control}).}
}
\value{
DMF structure.
}
\description{
Perform deviance matrix factorization.
}