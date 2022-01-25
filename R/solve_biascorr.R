#' Bias correction solver for multinomial logit type problems
#'
#' @param targets A dataframe with columns lu.from (optional), lu.to and value (all targets >= 0)
#' @param areas A dataframe of areas with columns lu.from (optional), ns and value, with all areas >= 0 and with sum(areas) >= sum(targets)
#' @param xmat A dataframe of explanatory variables with columns ks and value
#' @param betas A dataframe of coefficients with columns ks, lu.from (optional), lu.to & value
#' @param priors A dataframe of priors (if no \code{betas} were supplied) with columns ns, lu.from (optional), lu.to (with priors >= 0)
#' @param restrictions A dataframe with columns ns, lu.from (optional), lu.to and value. Values must be zero or one. If restrictions are one, the MNL function is set to zero
#' @param options A list with solver options. Call \code{\link{solve_biascorr_control}} for default options and for more detail.
#'
#' @details Given \code{p} targets matches either the projections from an MNL-type model or exogeneous priors.
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
#' Targets should be specified in a dataframe with at least a lu.to and a value column. All targets must be larger or equal to zero. If an lu.from column is supplied it has to be specified in all arguments. In this case the bias correction is performed for all lu.from classes.
#'
#' Areas correspond to either an areas per pixel (ns), with value or optionally the are of lu.from in a pixel. All areas must be larger The function expects lu.from
#'
#' Restrictions are binary and optional. If restrictions are supplied, in case \eqn{restrictions_ij = 1} then  \eqn{z_ij = 0}.
#'
#' @return A list containing
#' * \code{out.res} A \code{n x p} matrix of area allocations
#' * \code{out.solver} A list of the solver output
#'
#' @export solve_biascorr.mnl
#' @import nloptr
#' @import tidyr
#' @import dplyr
#' @import tibble
#'
#' @examples
#' ## A basic example
solve_biascorr.mnl = function(targets,areas,xmat,betas,priors = NULL,restrictions=NULL,
                              options = solve_biascorr_control()) {
  err.txt = options$err.txt;
  if (!"lu.from" %in% colnames(targets)) {
    targets = cbind(lu.from = "1",targets)
    areas = cbind(lu.from = "1",areas)
    betas = cbind(lu.from = "1",betas)
    if (!is.null(priors)) {priors = cbind(lu.from = "1",priors)}
    if (!is.null(restrictions)) {restrictions = cbind(lu.from = "1",restrictions)}
    null_lu.from = TRUE
  } else {null_lu.from = FALSE}

  if (!all(areas$value >= 0)) {stop(paste0(err.txt,"All areas must be larger or equal to zero."))}
  if (!all(targets >=0)) {
    targets[targets<0] = 0
    warning(paste0(err.txt,"Set negative targets to 0."))
  }
  if (!sum(areas$value) >= sum(targets$value)) {
    targets$value = targets$value / sum(targets$value) * sum(areas$value)
    warning(paste0(err.txt,"Sum of areas larger than sum of targets: lowered all targets to equal sum(areas)."))
  }

  lu.from <- unique(targets$lu.from)
  lu.to <- unique(targets$lu.to)
  ks = unique(betas$ks)
  # check if all targets are covered as either betas or priors
  chck.names = targets  %>% left_join(
    betas %>% group_by(lu.from,lu.to) %>% summarize(n = n(),.groups = "keep"),by = c("lu.from", "lu.to"))
  chck.names$n[is.na(chck.names$n)] = 0
  if (!is.null(priors)) {
    if (any(paste0(priors$lu.from) == paste0(priors$lu.to))) {stop(paste0(err.txt,"Priors lu.from must be unequal to lu.to."))}
    chck.names = chck.names %>%
      left_join(
        priors %>% group_by(lu.from,lu.to) %>% summarize(n2 = n(),.groups = "keep"),by = c("lu.from", "lu.to"))
    chck.names$n2[is.na(chck.names$n2)] = 0
    chck.names$n = chck.names$n + chck.names$n2
  }
  if (any(chck.names$n == 0)) {stop(paste0(err.txt,"Missing betas or priors for targets!"))}

  out.solver <- list()
  curr.lu.from <- lu.from[1]
  for(curr.lu.from in lu.from){
    # Extract targets
    curr.targets = dplyr::filter(targets,lu.from == curr.lu.from)$value
    names(curr.targets) <- targets$lu.to[targets$lu.from == curr.lu.from]

    # Extract areas
    curr.areas = dplyr::filter(areas,lu.from == curr.lu.from)$value
    names(curr.areas) <- areas$ns[areas$lu.from == curr.lu.from]

    # Extract betas
    curr.betas = dplyr::filter(betas,lu.from == curr.lu.from) %>%
      tidyr::pivot_wider(names_from = "lu.to",values_from = "value",id_cols = "ks") %>%
      tibble::column_to_rownames(var = "ks")
    curr.betas = as.matrix(curr.betas)

    # Extract xmat
    curr.xmat = dplyr::filter(xmat,ks %in% row.names(curr.betas)) %>%
      tidyr::pivot_wider(names_from = "ks",values_from = "value",id_cols = "ns")  %>%
      tibble::column_to_rownames(var = "ns")
    curr.xmat = as.matrix(curr.xmat)

    # Extract priors
    if (!is.null(priors) && any(priors$lu.from == curr.lu.from)) {
      curr.priors = dplyr::filter(priors,lu.from == curr.lu.from) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.priors = curr.priors[match(names(curr.areas),row.names(curr.priors)),,drop = FALSE]
      # check if betas have been provided for priors already
      if (any(colnames(curr.betas) %in% colnames(curr.priors))) {
        warning(paste0(err.txt,
                       "Priors provided for lu.from/lu.to combinations for which betas exist.\n These will be overwritten."))
        curr.betas = curr.betas[,-which(colnames(curr.betas) %in% colnames(curr.priors)),drop = FALSE]
      }
      curr.priors = as.matrix(curr.priors)
    } else {curr.priors = NULL}

    # Extract restrictions
    if (!is.null(restrictions) && any(restrictions$lu.from == curr.lu.from)) {
      curr.restrictions = dplyr::filter(restrictions,lu.from == curr.lu.from) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.restrictions = curr.restrictions[match(names(curr.areas),row.names(curr.restrictions)),,drop = FALSE]
      curr.restrictions = as.matrix(curr.restrictions)
    } else {curr.restrictions = NULL}

    p = length(curr.targets)
    n = length(curr.areas)
    p1 = ncol(curr.betas)
    k = nrow(curr.betas)

    if (ncol(curr.xmat)!=k || nrow(curr.xmat)!=n) {stop(paste0(err.txt,"Dimensions of xmat, ares and betas do not match."))}
    if (!is.null(priors)) {
      p2 = ncol(curr.priors)
      if (any(curr.priors<0)) {stop(paste0(err.txt,"Priors must be strictly non-negative."))}
    } else {p2 = 0}
    if (!p == p1 + p2) {stop(paste0(err.txt,"Dimensions of betas, targets and priors does not match."))}

    # check restrictions for consistency
    if (!is.null(restrictions)) {
      if (!all(colnames(curr.restrictions) %in% names(curr.targets))) {stop(paste0(err.txt,"Names of restrictions and targets do not match."))}
      if (any(!curr.restrictions %in% c(0,1))) {stop(paste0(err.txt,"Restrictions must be binary only!"))}
      if (nrow(curr.restrictions)!=n) {stop(paste0(err.txt,"Dimension of restrictions does not match."))}
      restr.mat = matrix(0,n,p); colnames(restr.mat) = names(curr.targets)
      restr.mat[,colnames(curr.restrictions)] = curr.restrictions
    } else {restr.mat = NULL}

    # out.res contains downscaled estimates; priors.mu econometric & other priors for estimation
    out.res = priors.mu = matrix(0,n,p)
    colnames(out.res) = colnames(priors.mu) = names(curr.targets)
    # match econometric priors
    priors.mu[,colnames(curr.betas)] = curr.xmat %*% curr.betas
    # make sure the priors are numerically well behaved
    priors.mu[priors.mu > options$MAX_EXP] = options$MAX_EXP
    priors.mu[priors.mu < -options$MAX_EXP] = -options$MAX_EXP
    priors.mu[,colnames(curr.betas)] = exp(priors.mu[,colnames(curr.betas)])
    # match other priors (if they exist)
    if (p2 > 0) {
      priors.mu[,colnames(curr.priors)] = curr.priors
    }
    # remove targets that are all zero
    not.zero = (curr.targets != 0)
    if (!all(not.zero)) {
      if (all(curr.targets == 0)) {return(list(out.res = out.res, out.solver = NULL))}
      curr.targets = curr.targets[not.zero]
      priors.mu = priors.mu[,not.zero,drop = FALSE]
      if (!is.null(curr.restrictions)) {restr.mat = restr.mat[,not.zero,drop = FALSE]}
    }
    x0 = curr.targets / sum(curr.targets + 1)
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
                             mu = priors.mu,areas = curr.areas,targets = curr.targets,
                             restrictions = restr.mat,cutoff = options$cutoff)
      if (res.x$objective < options$max_diff || countr > options$redo) {
        redo =FALSE
        res.x$par = res.x$solution
      } else {
        countr = countr + 1
        x0 = res.x$solution
      }
    }
    out.mu = mu.mnl(res.x$solution[1:length(curr.targets)],priors.mu,curr.areas,restr.mat,options$cutoff)
    if (all(not.zero)) {out.res = out.mu
    } else {out.res[,not.zero] = out.mu}
    out.solver[[curr.lu.from]] = res$out.solver

    # add residual own flows in output
    res$out.res2 = data.frame(ns = names(curr.areas),
                              curr.areas - rowSums(res$out.res),res$out.res)
    colnames(res$out.res2)[2] = paste0(curr.lu.from)
    # pivot into long format
    res.agg <- data.frame(
      lu.from=curr.lu.from,res$out.res2 %>%
        pivot_longer(cols = -c("ns"),names_to = "lu.to"))

    # aggregate results over dataframes
    if(curr.lu.from==lu.from[1]){out.res <- res.agg
    } else {
      out.res = bind_rows(out.res,res.agg)
    }
  }
  if (null_lu.from) {
    res.agg = res.agg %>% dplyr::select(-lu.from)
  }
  return(list(out.res = out.res, out.solver = out.solver))
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
#' @param algorithm Solver algorithm (see the \code{\link[nloptr]{nloptr}} package for documentation and more detail).
#' @param xtol_rel Relative tolerance of solver.
#' @param xtol_abs Absolute solver tolerance.
#' @param maxeval Maximum evaluation of solver.
#' @param MAX_EXP Numerical cutoff for MNL function.
#' @param cutoff Optional cutoff to avoid MNL values close to zero.
#' @param err.txt Error text for caller identification (used for debugging)
#' @param max_diff If difference to targets is larger, redo the estimation (helps to avoid convergence errors)
#' @param redo Maximum number of repeats
#'
#' @return List with default options for bias correction solver
#'
#' @details Call this function if you want to change default options for the bias corrections solver.
#'
#' @export solve_biascorr_control
#'
#' @examples
#' opts1 = solve_biascorr_control()
solve_biascorr_control = function(algorithm = "NLOPT_LN_SBPLX",
                                   xtol_rel = 1.0e-20,xtol_abs = 1.0e-20,maxeval = 1600,
                                   MAX_EXP = log(.Machine$double.xmax),cutoff = 0,
                                   redo = 2,max_diff = 10^-8,err.txt = "") {
  return(list(algorithm = algorithm,xtol_rel = xtol_rel,xtol_abs = xtol_abs,
              maxeval = maxeval,MAX_EXP = MAX_EXP,cutoff = cutoff,redo = redo,
              max_diff = max_diff,err.txt = err.txt
  ))
}

