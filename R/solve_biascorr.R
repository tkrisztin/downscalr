#' Bias correction solver for multinomial logit type problems
#'
#' @param targets A p x 1  vector of downscaling targets (with at least one targets >= 0)
#' @param areas An n x 1  vector of grid areas, with all areas >= 0 and with sum(areas) >= sum(targets)
#' @param xmat An n x k  matrix of explanatory variables
#' @param betas A k x p1 matrix (with p = p1 + p2)
#' @param priors A n x p2 matrix of priors (with priors >= 0)
#' @param restrictions A n x p binary matrix of restrictions; if restrictions are one, the MNL function is set to zero
#' @param options A list with solver options. Call \code{\link{get_opts.solve_biascorr}} for default options and for more detail.
#'
#' @details Given p targets matches either the projections from an MNL-type model or exogeneous priors.
#'
#' min \deqn{\sum  (  z ij areas_i  - targets_j )^2}
#' s.t. \deqn{ z_ij = \mu_ij}
#' \deqn{\mu_ij = \lambda_ij / (1 + \sum \lambda_ij )}
#' \deqn{\lambda_ij = x_j + \exp xmat_i betas_j} or \eqn{\lambda_ij = x_j + priors_ij}
#' \deqn{x_j >= 0}
#' with \eqn{i = 1,...,n} and  \eqn{j = 1,...,n}. #' For each target either betas and xmats or priors have to be supplied. Priors have to be strictly larger or equal to zero.
#'
#' When \code{cutoff} is specified, \eqn{z_ij} is defined as above if \eqn{mu_ij > cutoff}. If \eqn{mu_ij <= cutoff} then \eqn{z_ij = 0}. Per default \code{cutoff} is set to zero.
#'
#' Restrictions are binary and optional. If restrictions are supplied, in case \eqn{restrictions_ij = 1} then  \eqn{z_ij = 0}.
#'
#' @return A list containing
#' * \code{out.res} A n x p matrix of area allocations
#' * \code{out.solver} A list of the solver output
#'
#' @export solve_biascorr.mnl
#' @import nloptr
#'
#' @examples
#' ## A basic example
#' lu.to = c("crop","grass","forest","other")
#' n <- 100; k = 3
#' areas = 5+runif(n)
#' targets = runif(length(lu.to),1,sum(areas) / length(lu.to))
#' names(targets) = lu.to
#' betas = array(sample(-2:2,k*length(lu.to),replace = TRUE),c(k,length(lu.to)))
#' row.names(betas) = 1:k; colnames(betas) = lu.to
#' xmat = matrix(rnorm(n*k),n,k)
#' row.names(xmat) = 1:n; colnames(xmat) = 1:k
#' res1 = solve_biascorr.mnl(targets = targets,areas = areas,xmat = xmat,betas = betas)
#'
#' ## An example using priors
#' betas2 = betas[,-1]
#' priors = matrix(runif(n),n,1); colnames(priors) = lu.to[1]
#' res2 = solve_biascorr.mnl(targets = targets,areas = areas,xmat = xmat,
#'                          betas = betas2,priors = priors)
#'
#' ## An example with restrictions
#' restrictions = matrix(0,n,length(lu.to));
#' row.names(restrictions) = 1:n; colnames(restrictions) = lu.to
#' restrictions[,2] = sample(c(0,1),replace = TRUE,size = n)
#' res3 = solve_biascorr.mnl(targets = targets,areas = areas,xmat = xmat,
#'                          betas = betas2,priors = priors,
#'                          restrictions = restrictions)
solve_biascorr.mnl = function(targets,areas,xmat,betas,priors = NULL,restrictions=NULL,
                              options = get_opts.solve_biascorr()) {
  err.txt = options$err.txt;
  if (!all(targets >=0)) {
    targets[targets<0] = 0
    warning(paste0(err.txt,"Set negative targets to 0."))
  }
  if (!all(areas >= 0)) {stop(paste0(err.txt,"All areas must be larger or equal to zero."))}
  if (!sum(areas) >= sum(targets)) {
    targets = targets / sum(targets) * sum(areas)
    warning(paste0(err.txt,"Sum of areas larger than sum of targets: lowered all targets to equal sum(areas)."))
  }

  p = length(targets)
  n = length(areas)
  p1 = ncol(betas)
  k = nrow(betas)
  if (ncol(xmat)!=k || nrow(xmat)!=n) {stop(paste0(err.txt,"Dimensions of xmat, ares and betas do not match."))}
  if (!is.null(priors)) {
    p2 = ncol(priors)
    if (any(priors<0)) {stop(paste0(err.txt,"Priors must be strictly non-negative."))}
  } else {p2 = 0}
  if (!p == p1 + p2) {stop(paste0(err.txt,"Dimensions of betas, targets and priors does not match."))}
  if (is.null(colnames(betas)) || is.null(names(targets)) || (p2>0 && is.null(colnames(priors)))) {
    stop(paste0(err.txt,"Betas, targets and priors need colnames."))}
  if (!all(colnames(betas) %in% names(targets))) {stop(paste0(err.txt,"Names of betas and targets do not match."))}

  # check restrictions for consistency
  if (!is.null(restrictions)) {
    if (!all(colnames(restrictions) %in% names(targets))) {stop(paste0(err.txt,"Names of restrictions and targets do not match."))}
    if (any(!restrictions %in% c(0,1))) {stop(paste0(err.txt,"Restrictions must be binary only!"))}
    if (nrow(restrictions)!=n) {stop(paste0(err.txt,"Dimension of restrictions does not match."))}
    restr.mat = matrix(0,n,p); colnames(restr.mat) = names(targets)
    restr.mat[,colnames(restrictions)] = restrictions
  } else {restr.mat = NULL}

  # out.res contains downscaled estimates; priors.mu econometric & other priors for estimation
  out.res = priors.mu = matrix(0,n,p)
  colnames(out.res) = colnames(priors.mu) = names(targets)
  # match econometric priors
  priors.mu[,colnames(betas)] = xmat %*% betas
  # make sure the priors are numerically well behaved
  priors.mu[priors.mu > options$MAX_EXP] = options$MAX_EXP
  priors.mu[priors.mu < -options$MAX_EXP] = -options$MAX_EXP
  priors.mu[,colnames(betas)] = exp(priors.mu[,colnames(betas)])
  # match other priors (if they exist)
  if (p2 > 0) {
    priors.mu[,colnames(priors)] = priors
  }
  # remove targets that are all zero
  not.zero = (targets != 0)
  if (!all(not.zero)) {
    if (all(targets == 0)) {return(list(out.res = out.res, out.solver = NULL))}
    targets = targets[not.zero]
    priors.mu = priors.mu[,not.zero,drop = FALSE]
    if (!is.null(restrictions)) {restr.mat = restr.mat[,not.zero,drop = FALSE]}
  }
  x0 = targets / sum(targets + 1)
  opts <- list(algorithm = options$algorithm,
               xtol_rel = options$xtol_rel,
               xtol_abs = options$xtol_abs,
               maxeval = options$maxeval
  )
  redo = TRUE; countr = 1;
  while (redo) {
    res.x = nloptr::nloptr(x0,sqr_diff.mnl,
                   lb = rep(exp(-options$MAX_EXP),length(x0)),
                   ub = rep(exp(options$MAX_EXP),length(x0)),
                   opts=opts,
                   mu = priors.mu,areas = areas,targets = targets,restrictions = restr.mat,cutoff = options$cutoff)
    if (res.x$objective < 10^-8 || countr > 3) {
      redo =FALSE
      res.x$par = res.x$solution
    } else {
      countr = countr + 1
      x0 = res.x$solution
    }
  }
  out.mu = mu.mnl(res.x$solution[1:length(targets)],priors.mu,areas,restr.mat,options$cutoff)
  if (all(not.zero)) {out.res = out.mu
  } else {out.res[,not.zero] = out.mu}
  return(list(out.res = out.res, out.solver = res.x))
}



#' Multinomial logit mean calculation with restrictions and cutoff
#'
#' @param x Bias correction parameters
#' @param mu Means of each class
#' @param areas Vector of areas
#' @param restrictions Optional restrictions. Binary, if 1, the log-odds are set to zero
#' @param cutoff Optional, defaults to zero. If set, all log-odds below cutoff are set to zero
#'
#' @export mu.mnl
#'
#' @return Returns the vector of log-odds times the areas
mu.mnl = function(x,mu,areas,restrictions = NULL,cutoff = 0) {
  mu = mu * matrix(x,nrow(mu),ncol(mu),byrow = TRUE)
  mu = mu / rowSums(cbind(1,mu)) * areas
  if (!is.null(restrictions)) {
    mu[restrictions == 1] = 0
  }
  mu[mu<cutoff] = 0
  return(mu)
}

#' Squared differences optimization function with multinomial logit
#'
#' @param x Bias correction parameters
#' @param mu Means of each class
#' @param areas Vector of areas
#' @param targets Targets to match
#' @param restrictions Optional restrictions. Binary, if 1, the log-odds are set to zero
#' @param cutoff Optional, defaults to zero. If set, all log-odds below cutoff are set to zero
#'
#' @export sqr_diff.mnl
#'
#' @return Squared differences of optimized values and targets
sqr_diff.mnl = function(x,mu,areas,targets,restrictions,cutoff = 0) {
  x.mu = x[1:length(targets)]
  mu = mu.mnl(x.mu,mu,areas,restrictions,cutoff)
  return( sum((colSums(mu) - targets)^2) )
}

#' Get default options for bias corrections solver
#'
#' @return List with default options for bias correction solver
#'
#' @details Call this function if you want to change default options for the bias corrections solver.
#'
#' The default algorithm is NLOPT_LB_SBPLX (see the \code{\link[nloptr]{nloptr}} package for documentation and more detail). Relative and absolute tolerance are set to 1.0e-20 and \code{maxeval} to 1600.
#'
#' @export get_opts.solve_biascorr
#'
#' @examples
#' opts1 = get_opts.solve_biascorr()
get_opts.solve_biascorr = function() {
  return(list(
    "algorithm" = "NLOPT_LN_SBPLX",
    "xtol_rel" = 1.0e-20,
    "xtol_abs" = 1.0e-20,
    "maxeval" = 1600,
    MAX_EXP = log(.Machine$double.xmax),
    cutoff = 0,
    err.txt = ""
  ))
}

