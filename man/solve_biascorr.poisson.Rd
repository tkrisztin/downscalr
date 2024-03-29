% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_biascorr.R
\name{solve_biascorr.poisson}
\alias{solve_biascorr.poisson}
\title{Bias correction solver for Poisson type problems}
\usage{
solve_biascorr.poisson(targets, xmat, betas, options = downscale_control_pop())
}
\arguments{
\item{targets}{A dataframe with columns pop.type and value (all targets >= 0)}

\item{xmat}{A dataframe of explanatory variables with columns ks and value.}

\item{betas}{A dataframe of coefficients with columns ks, lu.from, lu.to & value}

\item{options}{A list with solver options. Call \code{\link{downscale_control_pop}} for default options and for more detail.}
}
\value{
A list containing
\itemize{
\item \code{out.res} A \code{n x p} matrix of area allocations
\item \code{out.solver} A list of the solver output
}
}
\description{
Bias correction solver for Poisson type problems
}
\details{
Given \code{J} targets matches  the projections from a Poisson-type model.

You should not call this functions directly, call \code{\link{downscale_pop}} instead.

Targets should be specified in a dataframe with  a pop.type and a value column.
All targets must be larger or equal to zero.
}
\examples{
## A basic example
}
