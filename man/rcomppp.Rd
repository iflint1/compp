% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcompp.R
\name{rcomppp}
\alias{rcomppp}
\title{Homogeneous COM-Poisson point process simulation}
\usage{
rcomppp(n = 1, lambda = 1, nu = 1, window = owin(), drop = TRUE)
}
\arguments{
\item{n}{Number of replications.}

\item{lambda}{Quasi-intensity.}

\item{nu}{Dispersion parameter}

\item{window}{Observation window.}

\item{drop}{Logical. If n = 1 and drop = TRUE (the default), the result will be a point pattern,
rather than a list containing a single point pattern.}
}
\value{
List of \code{spatstat::ppp} samples of a homogeneous COM-Poisson point process.
}
\description{
Homogeneous COM-Poisson point process simulation
}
