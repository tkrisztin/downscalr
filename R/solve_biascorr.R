#' Bias correction solver for multinomial logit type problems
#'
#' @param targets A dataframe with columns lu.from, lu.to and value (all targets >= 0)
#' @param areas A dataframe of areas with columns lu.from, ns and value, with all areas >= 0
#'   and with sum(areas) >= sum(targets)
#' @param xmat A dataframe of explanatory variables with columns ks and value.
#' @param betas A dataframe of coefficients with columns ks, lu.from, lu.to & value
#' @param priors A dataframe of priors (if no \code{betas} were supplied) with columns ns, lu.from, lu.to (with priors >= 0)
#' @param restrictions A dataframe with columns ns, lu.from, lu.to and value. Values must be zero or one. If restrictions are one, the MNL function is set to zero
#' @param options A list with solver options. Call \code{\link{downscale_control}} for default options and for more detail.
#'
#' @details Given \code{p} targets matches either the projections from an MNL-type model or exogeneous priors.
#'
#' You should not call this functions directly, call \code{\link{downscale}} instead.
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
                              options = downscale_control()) {
  lu.from <- unique(targets$lu.from)
  lu.to <- unique(targets$lu.to)
  ks = unique(betas$ks)

  out.solver <- list()
  curr.lu.from <- lu.from[1]
  full.out.res = NULL
  for(curr.lu.from in lu.from){
    err.txt = paste0(curr.lu.from," ",options$err.txt)

    # Extract targets
    curr.targets = dplyr::filter(targets,lu.from == curr.lu.from)$value
    names(curr.targets) <- targets$lu.to[targets$lu.from == curr.lu.from]
    curr.lu.to = names(curr.targets)

    # Extract betas
    curr.betas = dplyr::filter(betas,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
      tidyr::pivot_wider(names_from = "lu.to",values_from = "value",id_cols = "ks") %>%
      tibble::column_to_rownames(var = "ks")
    curr.betas = as.matrix(curr.betas)

    # Extract xmat
    ## IMPORTANT CHECK ORDER OF VARIABLES FIXED
    curr.xmat = dplyr::filter(xmat,ks %in% row.names(curr.betas)) %>%
      tidyr::pivot_wider(names_from = "ks",values_from = "value",id_cols = "ns")  %>%
      tibble::column_to_rownames(var = "ns") %>%
      dplyr::select(rownames(curr.betas))
    curr.xmat = as.matrix(curr.xmat)

    # Extract areas
    curr.areas = dplyr::filter(areas,lu.from == curr.lu.from)$value
    names(curr.areas) <- areas$ns[areas$lu.from == curr.lu.from]
    # BUGFIX: MW, order of curr.areas wrong, need to re-arrange based on xmat
    if (nrow(curr.xmat) > 0) {
      curr.areas = curr.areas[match(rownames(curr.xmat),names(curr.areas))]
    }

    # Extract priors
    ## IMPORTANT CHECK ORDER OF VARIABLES SIMILARLY TO XMAT
    if (!is.null(priors) && any(priors$lu.from == curr.lu.from)) {
      curr.priors = dplyr::filter(priors,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "value",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      curr.prior_weights = dplyr::filter(priors,lu.from == curr.lu.from & lu.to %in% curr.lu.to) %>%
        tidyr::pivot_wider(names_from = lu.to,values_from = "weight",id_cols = "ns") %>%
        tibble::column_to_rownames(var = "ns")
      name_match = match(names(curr.areas),row.names(curr.priors))
      curr.priors = curr.priors[name_match,,drop = FALSE]
      curr.prior_weights = curr.prior_weights[name_match,,drop = FALSE]
      # check if betas have been provided for priors already
      mixed_priors =c()
      nonmixed_priors = colnames(curr.priors)
      if (any(colnames(curr.betas) %in% colnames(curr.priors))) {
        mixed_priors = colnames(curr.betas)[which(colnames(curr.betas) %in% colnames(curr.priors))]
        nonmixed_priors = nonmixed_priors[!nonmixed_priors %in% mixed_priors]
        #warning(paste0(err.txt,
        #               "Priors provided for lu.from/lu.to combinations for which betas exist.\n These will be overwritten."))
        #curr.betas = curr.betas[,-which(colnames(curr.betas) %in% colnames(curr.priors)),drop = FALSE]
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

    if (p1 > 0 ) {
      if (ncol(curr.xmat)!=k || nrow(curr.xmat)!=n) {
        stop(paste0(err.txt,"Dimensions of xmat, areas and betas do not match."))
      }
    }
    if (!is.null(curr.priors)) {
      #p2 = ncol(curr.priors)
      p2 = length(nonmixed_priors)
      p2_mixed = length(mixed_priors)
      if (any(curr.priors<0)) {stop(paste0(err.txt,"Priors must be strictly non-negative."))}
    } else {p2 = 0;p2_mixed = 0}

    # check restrictions for consistency
    if (!is.null(curr.restrictions) & any(colnames(curr.restrictions) %in% names(curr.targets))) {
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
      priors.mu[,nonmixed_priors] = curr.priors[,nonmixed_priors]
    }
    if (p2_mixed > 0) {
      w1 = curr.prior_weights[,mixed_priors,drop = FALSE]
      #re-scale exogeneous prior to priors.mu
      eco.priors_sum = apply(priors.mu[,mixed_priors,drop = FALSE],c(2),sum)
      exo.priors = curr.priors[,mixed_priors,drop = FALSE]
      exo.priors_sum = apply(exo.priors, 2, sum)
      exo.priors = t(
        (t(exo.priors) / exo.priors_sum) * eco.priors_sum   )
      priors.mu[,mixed_priors] = as.matrix((1-w1)*priors.mu[,mixed_priors] + w1*exo.priors)
    }
    # remove targets that are all zero
    not.zero = (curr.targets != 0)
    if (all(curr.targets == 0)) {

      #catch case if all targets are equal zero
      out.solver[[curr.lu.from]] = NULL
    } else {

      #cut out zero targets from targets and priors
      if (any(curr.targets == 0) && !all(curr.targets == 0)) {
        curr.targets = curr.targets[not.zero]
        priors.mu = priors.mu[,not.zero,drop = FALSE]
        if (!is.null(curr.restrictions)) {restr.mat = restr.mat[,not.zero,drop = FALSE]}
      }

      #proceed with bias correction
      x0 = curr.targets / sum(curr.targets + 1)
      opts <- list(algorithm = options$algorithm,
                   xtol_rel = options$xtol_rel,
                   xtol_abs = options$xtol_abs,
                   maxeval = options$maxeval
      )
      # check if the optimiser uses gradient or not
      if (grepl("_LD_",opts$algorithm)) {
        eval_grad_f = grad_sqr_diff.mnl
      } else {
        eval_grad_f = NULL
      }
      res.x = nloptr::nloptr(x0 = x0,
                             eval_f = sqr_diff.mnl,
                             eval_grad_f = eval_grad_f,
                             lb = rep(exp(-options$MAX_EXP),length(x0)),
                             ub = rep(exp(options$MAX_EXP),length(x0)),
                             opts=opts,
                             mu = priors.mu,areas = curr.areas,targets = curr.targets,
                             restrictions = restr.mat,cutoff = options$cutoff)
      res.x$par = res.x$solution

      out.mu = mu.mnl(res.x$solution[1:length(curr.targets)],priors.mu,curr.areas,restr.mat,options$cutoff)

      # If differences are still too large, do individual logit models boosted by grid search
      if (res.x$objective > options$max_diff) {
        # Find targets where out.mu deviates more than max_diff from the target
        error_ind = (curr.targets - colSums(out.mu))^2 > options$max_diff
        error_targets = curr.targets[error_ind]
        error_restrictions = restr.mat[,error_ind]

        # Calculate residual areas -  res_areas
        res_areas = curr.areas
        if (any(!error_ind)) {
          res_areas = res_areas - rowSums(out.mu[,!error_ind])
        }

        # Loop over remaining targets
        for (ccc in 1:length(error_targets)) {
          curr_error_target = error_targets[ccc]
          curr_error_restrictions = error_restrictions[,ccc]
          curr_error_mu = priors.mu[,names(curr_error_target),drop=FALSE]

          # Do an iterated grid search to find correct scaling coefficient for priors
          curr_scaling = iterated_grid_search(min_param = -options$MAX_EXP,
                                              max_param = options$MAX_EXP,
                                              func = sqr_diff.mnl,
                                              max_iterations = 10,
                                              precision_threshold = 1e-3,
                                              exp_transform = TRUE,
                                              mu = curr_error_mu,areas = res_areas,
                                              targets = curr_error_target,
                                              restrictions = curr_error_restrictions,
                                              cutoff = options$cutoff)

          # Re-scale prior
          curr_error_mu = curr_error_mu * curr_scaling$best_param

          # Optimize with scaled priors
          res.x = nloptr::nloptr(x0 = 1,
                                 eval_f = sqr_diff.mnl,
                                 eval_grad_f = eval_grad_f,
                                 lb = exp(-options$MAX_EXP),
                                 ub = exp(options$MAX_EXP),
                                 opts=opts,
                                 mu = curr_error_mu,areas = res_areas,targets = curr_error_target,
                                 restrictions = curr_error_restrictions,cutoff = options$cutoff)

          # Calculate areas of current target with mu.mnl
          curr_error_out.mu =
            mu.mnl(res.x$solution,
                   curr_error_mu,res_areas,curr_error_restrictions,
                   options$cutoff)

          # Add calculated areas to out.mu
          out.mu[,names(curr_error_target)] = curr_error_out.mu

          # Substract areas from res.areas
          res_areas = res_areas - curr_error_out.mu

          # Add note to res.x
          res.x$message = "INDIVIDUAL LOGIT BOOSTED BY GRID SEARCH: Standard optimization failed to converge"
        }
      }

      if (all(not.zero)) {out.res = out.mu
      } else {out.res[,not.zero] = out.mu}
      out.solver[[curr.lu.from]] = res.x
    }

    # add residual own flows in output
    out.res2 = data.frame(ns = names(curr.areas),
                              curr.areas - rowSums(out.res),out.res)
    colnames(out.res2)[2] = paste0(curr.lu.from)
    # pivot into long format
    res.agg <- out.res2 %>%
        pivot_longer(cols = -c("ns"),names_to = "lu.to") %>%
      bind_cols(lu.from = curr.lu.from)

    # aggregate results over dataframes
    if(curr.lu.from==lu.from[1]){
      full.out.res <- res.agg
    } else {
      full.out.res = bind_rows(full.out.res,res.agg)
    }
  }
  return(list(out.res = full.out.res, out.solver = out.solver))
}


#' Bias correction solver for Poisson type problems
#'
#' @param targets A dataframe with columns pop.type and value (all targets >= 0)
#' @param xmat A dataframe of explanatory variables with columns ks and value.
#' @param betas A dataframe of coefficients with columns ks, lu.from, lu.to & value
#' @param options A list with solver options. Call \code{\link{downscale_control_pop}} for default options and for more detail.
#'
#' @details Given \code{J} targets matches  the projections from a Poisson-type model.
#'
#' You should not call this functions directly, call \code{\link{downscale_pop}} instead.
#'
#' Targets should be specified in a dataframe with  a pop.type and a value column.
#' All targets must be larger or equal to zero.
#'
#' @return A list containing
#' * \code{out.res} A \code{n x p} matrix of area allocations
#' * \code{out.solver} A list of the solver output
#'
#' @export solve_biascorr.poisson
#' @import nloptr
#' @import tidyr
#' @import dplyr
#' @import tibble
#'
#' @examples
#' ## A basic example
solve_biascorr.poisson = function(targets,xmat,betas,
                              options = downscale_control_pop()) {
  pop.type <- unique(targets$pop.type)
  ks = unique(betas$ks)
  value = NULL

  out.solver <- list()
  curr.pop.type <- pop.type[1]
  full.out.res = data.frame()
  for(curr.pop.type in pop.type){
    err.txt = paste0(curr.pop.type," ",options$err.txt)

    # Extract targets
    curr.targets = dplyr::filter(targets,pop.type == curr.pop.type)$value

    # Extract betas
    curr.betas = dplyr::filter(betas,pop.type == curr.pop.type) %>%
      tibble::column_to_rownames(var = "ks") %>%
      dplyr::select(value)
    curr.betas = as.matrix(curr.betas)

    # Extract xmat
    ## IMPORTANT CHECK ORDER OF VARIABLES FIXED
    curr.xmat = dplyr::filter(xmat,ks %in% row.names(curr.betas)) %>%
      tidyr::pivot_wider(names_from = "ks",values_from = "value",id_cols = "ns")  %>%
      tibble::column_to_rownames(var = "ns") %>%
      dplyr::select(rownames(curr.betas))
    curr.xmat = as.matrix(curr.xmat)

    n = nrow(curr.xmat)
    k = nrow(curr.betas)


    if (ncol(curr.xmat)!=k || nrow(curr.xmat)!=n) {
      stop(paste0(err.txt,"Dimensions of xmat, areas and betas do not match."))
    }

    # out.res contains downscaled estimates; priors.mu econometric & other priors for estimation
    out.res = priors.mu = matrix(0,n,1)
    colnames(out.res) = colnames(priors.mu) = names(curr.pop.type)
    # match econometric priors
    priors.mu = curr.xmat %*% curr.betas
    # make sure the priors are numerically well behaved
    priors.mu[priors.mu > options$MAX_EXP] = options$MAX_EXP
    priors.mu[priors.mu < -options$MAX_EXP] = -options$MAX_EXP
    priors.mu = exp(priors.mu)

    res.x = list()
    if (curr.targets == 0) {
      #catch case if all targets are equal zero
      res.x$par = -Inf
      out.solver[[curr.pop.type]] = res.x
      out.res = priors.mu * 0
    } else {
      res.x$par = log(curr.targets / sum(priors.mu))

      out.res = exp(res.x$par) * priors.mu
      out.solver[[curr.pop.type]] = res.x
    }

    # add residual own flows in output
    res.agg = data.frame(ns = rownames(curr.xmat),
                          pop.type = curr.pop.type, value = out.res)

    # aggregate results over dataframes
    full.out.res = bind_rows(full.out.res,res.agg)
  }
  return(list(out.res = full.out.res, out.solver = out.solver))
}




