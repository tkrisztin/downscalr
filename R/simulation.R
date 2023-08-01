


#' Simulating land-use from a data generating process
#'
#' This function can be used to generate land use data from a data generating process
#'
#' The generated multinomial logit land use model over \eqn{J} distinct land uses
#' (with \eqn{j',j = 1,...,J}) takes the form:
#'
#' \deqn{
#'    y_{ijt} = \frac{\exp(\alpha_{jt} + X_i \beta_j)}{ \sum^{J}_{j' = 1} \exp(\alpha_{jt} + X_i \beta_{j'})}
#' }
#'
#' with \eqn{X_i} being the \eqn{i}-th row of the \eqn{n \times k} matrix \eqn{X}. \eqn{\beta_j} is the
#' \eqn{k \times J} matrix of land use class specific slope parameters. \eqn{\alpha_{jt}} are class
#' specific intercepts, which are typically optimised using the bias correction method. These are assumed
#' to vary over time.
#'
#' For the purpose of identification the slope parameters \eqn{\beta_J} and intercepts
#' \eqn{\alpha_J} associated with the \eqn{J}-th class have to be zero.
#'
#' The function generates the \eqn{NT \ times J} matrix \eqn{Y}. Based on this, the function
#' generates a set of targets \eqn{T_{jt}} for downscaling. These are calculated by:
#'
#' \deqn{
#'  T_{jt} = \sum^n_{i=1} y_{ijt} a_i
#' }
#'
#' The output of this function can be directly used in \code{\link{downscale}}.
#'
#'
#' @param n Number of spatial observations \eqn{n}.
#' @param J Number of land use classes \eqn{J}. Defaults to three.
#' @param tt Number of time steps to simulate \eqn{T}. Defaults to two.
#' @param k Number of covariates to simulate \eqn{k}. Defaults to one.
#' @param areas Vector of dimensions \eqn{n \times 1}. Provides the spatial observation
#'              areas \eqn{a_i} for calculating targets \eqn{T_1,...,T_J}.
#' @param alphas Matrix of dimensions \eqn{T \times J}. Provides the true values
#'              for all \eqn{\alpha_{jt}}. The \eqn{J}-th column of the matrix
#'              has to be equal to zero. Defaults to a random uniform matrix.
#' @param xmat Matrix of dimensions \eqn{n \times k}. Provides the values for the
#'              matrix \eqn{X}. Defaults to a matrix of random, normally distributed
#'              values (with zero mean and uniform variance).
#' @param betas Matrix of dimensions \eqn{k \times J}. Provides the values for \eqn{\beta_1,...,\beta_J}.
#'              The \eqn{J}-th column of the matrix
#'              has to be equal to zero. Defaults to a matrix of random, normally distributed
#'              values (with zero mean and uniform variance).
#'
#' @return A list with the generated \eqn{\alpha} (\code{alphas}), \eqn{X}  (\code{xmat}),
#'          \eqn{\beta}  (\code{betas}), \eqn{Y}  (\code{Y}),
#'          and \eqn{T} (\code{targets}). The returned values are in long format as required by
#'          \code{\link{downscale}}.
#'
#' @export sim_lu
#' @import tidyr
#' @import dplyr
#' @importFrom stats rnorm runif
#'
#' @examples
#' dgp1 = sim_lu(n = 100)
sim_lu = function(n,J=3,tt=2,k = 1,
                   areas = stats::runif(n),
                   alphas = cbind(matrix(stats::runif(tt * (J-1)),tt,J-1),0),
                   xmat = matrix(stats::rnorm(n*k),n,k),
                   betas = cbind(matrix(stats::rnorm(k*(J-1)),k,J-1),0)) {

  ns = as.character(1:n)
  lu.to = as.character(paste0("lu",1:J))
  ks = as.character(paste0("k",1:k))
  times = as.character(1:tt)
  areas2 = data.frame(ns = ns,value = areas)

  gdp_luc = data.frame()
  for (t in 1:tt) {
    mu = exp( matrix(alphas[t,],n,J,byrow = TRUE) + xmat %*% betas)
    yhat = mu / rowSums(mu) * areas
    colnames(yhat) = lu.to
    gdp_luc = gdp_luc %>% bind_rows(
      data.frame(times = times[t],ns = ns,yhat) %>%
        pivot_longer(cols = -c(1:2),names_to = "lu.from")
    )
  }
  lu.from = value = NULL
  targets = gdp_luc %>%
    group_by(times,lu.from) %>%
    summarise(value = sum(value))
  return(list(n = n, J = J, tt = tt,k=1,
              ns = ns,
              lu.to = lu.to,
              ks = ks,
              times = times,
              areas = areas2,
              alphas = alphas,
              xmat = xmat,
              betas = betas,
              Y = gdp_luc,
              targets = targets))
}
