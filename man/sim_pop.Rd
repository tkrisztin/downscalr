% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{sim_pop}
\alias{sim_pop}
\title{Simulating population from a data generating process}
\usage{
sim_pop(
  n,
  J = 2,
  tt = 2,
  k = 2,
  alphas = matrix(stats::runif(tt * (J - 1)), tt, J) * 2,
  xmat = matrix(stats::rnorm(n * k), n, k),
  betas = matrix(stats::rnorm(k * J), k, J)
)
}
\arguments{
\item{n}{Number of spatial observations \eqn{n}.}

\item{J}{Number of population types \eqn{J}. Defaults to two}

\item{tt}{Number of time steps to simulate \eqn{T}. Defaults to two.}

\item{k}{Number of covariates to simulate \eqn{k}. Defaults to two.}

\item{alphas}{Matrix of dimensions \eqn{T \times J}. Provides the true values
for all \eqn{\alpha_{jt}}. Defaults to a random uniform matrix.}

\item{xmat}{Matrix of dimensions \eqn{n \times k}. Provides the values for the
matrix \eqn{X}. Defaults to a matrix of random, normally distributed
values (with zero mean and uniform variance).}

\item{betas}{Matrix of dimensions \eqn{k \times J}. Provides the values for \eqn{\beta_1,...,\beta_J}.
Defaults to a matrix of random, normally distributed
values (with zero mean and uniform variance).}
}
\value{
A list with the generated \eqn{\alpha} (\code{alphas}), \eqn{X}  (\code{xmat}),
\eqn{\beta}  (\code{betas}), \eqn{Y}  (\code{Y}),
and \eqn{T} (\code{targets}). The returned values are in long format as required by
\code{\link{downscale_pop}}.
}
\description{
This function can be used to generate population data from a data generating process for testing downscaling.
}
\details{
The generated Poisson model over \eqn{J} distinct population tyes
(with \eqn{j',j = 1,...,J}) takes the form:

\deqn{
   y_{ijt} = \exp{\exp(\alpha_{jt} + X_i \beta_{j})}
}

with \eqn{X_i} being the \eqn{i}-th row of the \eqn{n \times k} matrix \eqn{X}. \eqn{\beta_j} is the
\eqn{k \times J} matrix of population type specific slope parameters. \eqn{\alpha_{jt}} are population type
specific intercepts, which are typically optimised using the bias correction method. These are assumed
to vary over time.

The function generates the \eqn{NT \ times J} matrix \eqn{Y}. Based on this, the function
generates a set of targets \eqn{T_{jt}} for downscaling. These are calculated by:

\deqn{
 T_{jt} = \sum^n_{i=1} y_{ijt}
}

The output of this function can be directly used in \code{\link{downscale_pop}}.
}
\examples{
dgp1 = sim_pop(n = 100)
}
