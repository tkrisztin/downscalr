
#' Simulating land-use change from a data generating process
#'
#' This function can be used to generate land use change data from a data generating process
#'
#' The function simulates a set of \eqn{J} different multinomial logit land use models over \eqn{J} distinct land uses
#' (with \eqn{j',j = 1,...,J}). Each of the models takes the form:
#'
#' \deqn{
#'    y_{iljt} = \frac{\exp(\alpha_{ljt} + X_i \beta_{lj})}{ \sum^{J}_{j' = 1} \exp(\alpha_{ljt} + X_i \beta_{lj'})}
#' }
#'
#' where \eqn{l} denotes the individual multinomial logit models. \eqn{X_i} is the \eqn{i}-th
#' row of the \eqn{n \times k} matrix \eqn{X}. \eqn{\beta_j} is the
#' \eqn{k \times J} matrix of land use class specific slope parameters. \eqn{\alpha_{jt}} are class
#' specific intercepts, which are typically optimised using the bias correction method. These are assumed
#' to vary over time.
#'
#' For the purpose of identification the slope parameters \eqn{\beta_J} and intercepts
#' \eqn{\alpha_J} associated with the \eqn{J}-th class have to be zero.
#'
#' Each of the \eqn{J} different models captures land use being converted from the \eqn{l}-th land use class
#' to \eqn{J} other land uses.
#'
#' The function generate \eqn{J} different  \eqn{NT \ times J} matrices \eqn{Y_l}. Based on these, the function
#' generates a set of targets \eqn{T_{ljt}} for downscaling. These are calculated by:
#'
#' \deqn{
#'  T_{ljt} = \sum^n_{i=1} y_{iljt} a_{ilt}.
#' }
#'
#' The initial values for \eqn{a_il1} are randomly generated. After this the values of \eqn{a_{ilt}} are dynamicall
#' updated through:
#'
#' \deqn{
#' a_{ilt+1} = \sum^n_{j=1} y_{iljt} a_{ilt}
#' }
#'
#' The output of this function can be directly used in \code{\link{downscale}}.
#'
#'
#' @inheritParams sim_lu
#' @param areas A matrix of dimensions \eqn{n \times J}.
#'              Provides the spatial observation of
#'              areas \eqn{a_{ij1}} in \eqn{t=1} for calculating targets \eqn{T_1,...,T_J}. Defaults to a
#'              random uniform matrix.
#' @param alphas A list of length \eqn{J}. Each list item contains a matrix of dimensions \eqn{T \times J}.
#'              These provide the true values for all \eqn{\alpha_{ljt}}. The \eqn{J}-th column of the matrix
#'              has to be equal to zero. Defaults to a list containing random uniform matrices.
#' @param betas A list of length \eqn{J}. Each list item contains a matrix of dimensions
#'              \eqn{k \times J}. Provides the values for each \eqn{\beta_{1l},...,\beta_{jl}}.
#'              The \eqn{J}-th column of the matrix
#'              has to be equal to zero. Defaults to a list containing \eqn{J} matrices of
#'              random, normally distributed values (with zero mean and uniform variance).
#'
#' @return A list with the generated \eqn{\alpha} (\code{alphas}), \eqn{X}  (\code{xmat}),
#'          \eqn{\beta}  (\code{betas}), \eqn{Y}  (\code{Y}),
#'          and \eqn{T} (\code{targets}). The returned values are in long format as required by
#'          \code{\link{downscale}}.
#' @export
#'
#' @examples
#' dgp1 = sim_luc(n = 100)
sim_luc = function(n,J=3,tt=2,k = 1,
                  areas = matrix(runif(n*J),n,J),
                  alphas = replicate(J,
                                     cbind(matrix(stats::runif(tt * (J-1)),tt,J-1),0),simplify = FALSE),
                  xmat = matrix(stats::rnorm(n*k),n,k),
                  betas = replicate(J,
                                    cbind(matrix(stats::rnorm(k*(J-1)),k,J-1),0),simplify = FALSE)) {

  value = lu.to = NULL
  ns = as.character(1:n)
  ks = as.character(paste0("k",1:k))
  lu.from = as.character(paste0("lu",1:J))
  lu.to = as.character(paste0("lu",1:J))
  times = 1:tt

  colnames(areas) = lu.from
  start.areas = data.frame(ns = ns,areas) %>%
    pivot_longer(cols = -c(1), names_to = "lu.from")

  dgp_luc = data.frame()
  for (t in 1:tt) {
    for (jj in 1:J) {
      curr_res = sim_lu(n = n,
                        J = J, tt = 1, k = k,
                        alphas = alphas[[jj]][t,,drop=FALSE],
                        areas = areas[,jj],
                        xmat = xmat,
                        betas = betas[[jj]])
      dgp_luc = dgp_luc %>% bind_rows(
        data.frame(times = as.character(times[t]),lu.from = lu.from[jj],
                   curr_res$Y %>% select(-times))
      )
    }
    areas = dgp_luc %>%
      filter(times == as.character(times[t])) %>%
      group_by(ns,lu.to) %>%
      summarise(value = sum(value)) %>%
      ungroup() %>%
      pivot_wider(names_from = lu.to) %>%
      select(-ns) %>% as.matrix()
  }
  targets = dgp_luc %>%
    group_by(times,lu.to,lu.from) %>%
    summarise(value = sum(value)) %>%
    dplyr::ungroup()
  all_betas = data.frame()
  for (jj in 1:J) {
    curr.betas = betas[[jj]]
    colnames(curr.betas) = lu.to
    curr.betas = data.frame(ks = ks,curr.betas) %>%
      pivot_longer(cols = -c(1),names_to = "lu.to")
    all_betas = all_betas %>% bind_rows(
      data.frame(lu.from = lu.from[jj],curr.betas)
    )
  }
  colnames(xmat) = ks
  xmat = data.frame(ns = ns,xmat) %>%
    pivot_longer(cols = -c(1),names_to = "ks")
  return(list(n = n, J = J, tt = tt,k=k,
              ns = ns,
              lu.to = lu.to,
              ks = ks,
              times = times,
              start.areas = start.areas,
              alphas = alphas,
              xmat = xmat,
              betas = all_betas,
              Y = dgp_luc,
              targets = targets))
}

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
#'              areas \eqn{a_i} for calculating targets \eqn{T_1,...,T_J}. Defaults to a
#'              random uniform matrix.
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

  dgp_luc = data.frame()
  for (t in 1:tt) {
    mu = exp( matrix(alphas[t,],n,J,byrow = TRUE) + xmat %*% betas)
    yhat = mu / rowSums(mu) * areas
    colnames(yhat) = lu.to
    dgp_luc = dgp_luc %>% bind_rows(
      data.frame(times = times[t],ns = ns,yhat) %>%
        pivot_longer(cols = -c(1:2),names_to = "lu.to")
    )
  }
  value = NULL
  targets = dgp_luc %>%
    group_by(times,lu.to) %>%
    summarise(value = sum(value)) %>%
    dplyr::ungroup()
  colnames(betas) = lu.to
  betas = data.frame(ks = ks,betas) %>%
    pivot_longer(cols = -c(1),names_to = "lu.to")
  colnames(xmat) = ks
  xmat = data.frame(ns = ns,xmat) %>%
    pivot_longer(cols = -c(1),names_to = "ks")
  return(list(n = n, J = J, tt = tt,k=k,
              ns = ns,
              lu.to = lu.to,
              ks = ks,
              times = times,
              areas = areas2,
              alphas = alphas,
              xmat = xmat,
              betas = betas,
              Y = dgp_luc,
              targets = targets))
}

#' Simulating population from a data generating process
#'
#' This function can be used to generate population data from a data generating process for testing downscaling.
#'
#' The generated Poisson model over \eqn{J} distinct population tyes
#' (with \eqn{j',j = 1,...,J}) takes the form:
#'
#' \deqn{
#'    y_{ijt} = \exp{\exp(\alpha_{jt} + X_i \beta_{j})}
#' }
#'
#' with \eqn{X_i} being the \eqn{i}-th row of the \eqn{n \times k} matrix \eqn{X}. \eqn{\beta_j} is the
#' \eqn{k \times J} matrix of population type specific slope parameters. \eqn{\alpha_{jt}} are population type
#' specific intercepts, which are typically optimised using the bias correction method. These are assumed
#' to vary over time.
#'
#' The function generates the \eqn{NT \ times J} matrix \eqn{Y}. Based on this, the function
#' generates a set of targets \eqn{T_{jt}} for downscaling. These are calculated by:
#'
#' \deqn{
#'  T_{jt} = \sum^n_{i=1} y_{ijt}
#' }
#'
#' The output of this function can be directly used in \code{\link{downscale_pop}}.
#'
#'
#' @param n Number of spatial observations \eqn{n}.
#' @param J Number of population types \eqn{J}. Defaults to two
#' @param tt Number of time steps to simulate \eqn{T}. Defaults to two.
#' @param k Number of covariates to simulate \eqn{k}. Defaults to two.
#' @param alphas Matrix of dimensions \eqn{T \times J}. Provides the true values
#'              for all \eqn{\alpha_{jt}}. Defaults to a random uniform matrix.
#' @param xmat Matrix of dimensions \eqn{n \times k}. Provides the values for the
#'              matrix \eqn{X}. Defaults to a matrix of random, normally distributed
#'              values (with zero mean and uniform variance).
#' @param betas Matrix of dimensions \eqn{k \times J}. Provides the values for \eqn{\beta_1,...,\beta_J}.
#'              Defaults to a matrix of random, normally distributed
#'              values (with zero mean and uniform variance).
#'
#' @return A list with the generated \eqn{\alpha} (\code{alphas}), \eqn{X}  (\code{xmat}),
#'          \eqn{\beta}  (\code{betas}), \eqn{Y}  (\code{Y}),
#'          and \eqn{T} (\code{targets}). The returned values are in long format as required by
#'          \code{\link{downscale_pop}}.
#'
#' @export sim_pop
#' @import tidyr
#' @import dplyr
#' @importFrom stats rnorm runif
#'
#' @examples
#' dgp1 = sim_pop(n = 100)
sim_pop = function(n,J=2,tt=2,k = 2,
                  alphas = matrix(stats::runif(tt * (J-1)),tt,J) * 2,
                  xmat = matrix(stats::rnorm(n*k),n,k),
                  betas = matrix(stats::rnorm(k*J),k,J)) {

  ns = as.character(1:n)
  pop.types = as.character(paste0("pop",1:J))
  ks = as.character(paste0("k",1:k))
  times = as.character(1:tt)

  dgp_pop = data.frame()
  for (t in 1:tt) {
    yhat = exp( matrix(alphas[t,],n,J,byrow = TRUE) + xmat %*% betas)
    colnames(yhat) = pop.types
    dgp_pop = dgp_pop %>% bind_rows(
      data.frame(times = times[t],ns = ns,yhat) %>%
        pivot_longer(cols = -c(1:2),names_to = "pop.type")
    )
  }
  pop.type = value = NULL
  targets = dgp_pop %>%
    group_by(times,pop.type) %>%
    summarise(value = sum(value)) %>%
    dplyr::ungroup()
  colnames(betas) = pop.types
  betas = data.frame(ks = ks,betas) %>%
    pivot_longer(cols = -c(1),names_to = "pop.type")
  colnames(xmat) = ks
  xmat = data.frame(ns = ns,xmat) %>%
          pivot_longer(cols = -c(1),names_to = "ks")
  return(list(n = n, J = J, tt = tt,k=k,
              ns = ns,
              pop.type = pop.type,
              ks = ks,
              times = times,
              alphas = alphas,
              xmat = xmat,
              betas = betas,
              Y = dgp_pop,
              targets = targets))
}
