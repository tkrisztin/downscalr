#' Downscaling of land-use (change) data
#'
#' @param targets A dataframe with columns times, lu.from (optional), lu.to and value (all targets >= 0)
#' @param start.areas  A dataframe of areas with columns lu.from (optional), ns and value, with all areas >= 0 and with sum(areas) >= sum(targets)
#' @param xmat A dataframe of explanatory variables with columns ns, ks and value
#' @param betas A dataframe of coefficients with columns ks, lu.from (optional), lu.to & value
#' @param areas.update.fun function providing update for dynamic xmat columns, must take as arguments res, curr.areas, priors, xmat.proj, must return dataframe with columns ns, ks & value defaults to areas.sum_to() which sums over lu.to
#' @param xmat.coltypes ks vector, each can be either "static", "dynamic", or "projected"
#' @param xmat.proj dataframe with columns times, ns, ks, must be present for each xmat.coltype specified as projected
#' @param xmat.dyn.fun function providing update for dynamic xmat columns, must take as arguments res, curr.areas, priors, xmat.proj must return ns x ks(dynamic) columns
#' @param priors A dataframe of priors (if no \code{betas} were supplied) with columns ns, lu.from (optional), lu.to (with priors >= 0)
#' @param restrictions A dataframe with columns ns, lu.from (optional), lu.to and value. Values must be zero or one. If restrictions are one, the MNL function is set to zero
#' @param options A list with solver options. Call \code{\link{downscale_control}} for default options and for more detail.
#'
#' @details Given \code{p} targets matches either the projections from an MNL-type model or exogenous priors.
#'
#' @return A list containing
#' * \code{out.res} Dataframe with columns times, ns, lu.from, lu.to & value (area allocation)
#' * \code{out.solver} A list of the solver output
#'
#' @export downscale
#' @import nloptr
#' @import tidyr
#' @import dplyr
#' @import tibble
#'
#' @examples
#' ## A basic example
downscale = function(targets,start.areas,xmat,betas,
                     areas.update.fun = areas.sum_to,
                     xmat.coltypes = NULL,
                     xmat.proj = NULL,xmat.dyn.fun = xmat.sum_to,
                     priors = NULL,restrictions=NULL,
                     options = downscale_control()) {
  # Handle input checking
  err.txt = options$err.txt
  targets = complete_targets(targets)
  start.areas = complete_areas(start.areas)
  xmat = complete_xmat(xmat)
  betas = complete_betas(betas)
  complete_xmat.coltypes = complete_xmat.coltypes(xmat.coltypes,xmat)
  if (!is.null(priors)) {priors = complete_priors(priors,xmat)}
  if (!is.null(restrictions)) {restrictions = complete_restrictions(restrictions,xmat)}
  if (!is.null(xmat.proj )) {xmat.proj = complete_xmat.proj(xmat.proj)}
  err_check_inputs(targets,start.areas,xmat,betas,
                   areas.update.fun,xmat.coltypes,
                   xmat.proj,xmat.dyn.fun,
                   priors,restrictions,err.txt)

  # save column types of xmat
  if (any(xmat.coltypes$value == "projected")) {proj.colnames = filter(xmat.coltypes,.data$value == "projected")$ks}
  if (any(xmat.coltypes$value == "dynamic")) {dyn.colnames = filter(xmat.coltypes,.data$value == "dynamic")$ks}

  # Set starting values
  curr.areas = start.areas
  curr.xmat = xmat
  curr.priors = priors
  curr.restrictions = restrictions

  times = unique(targets$times)
  out.solver <- list()
  for (curr.time in times) {
    # Extract targets
    curr.targets = filter(targets,times == curr.time) %>% select(-times)

    if (options$solve_fun == "solve_biascorr") {
      curr.options = options
      curr.options$err.txt = paste0(curr.time," ",curr.options$err.txt)
      res = solve_biascorr.mnl(targets = curr.targets,
                               areas = curr.areas,
                               xmat = curr.xmat,
                               betas = betas,
                               priors = curr.priors,
                               restrictions=curr.restrictions,options = curr.options)
      out.solver[[as.character(curr.time)]] = res$out.solver
    }

    # update curr.area
    curr.areas = areas.update.fun(res, curr.areas, priors, xmat.proj)

    # update projected xmats
    if (any(xmat.coltypes$value == "projected")) {
      xmat = xmat %>%
        left_join(dplyr::filter(xmat.proj,times == curr.time) %>%
                    dplyr::select(-times) %>% rename("proj" = "value") %>%
                    mutate(ns = as.character(.data$ns)),by = c("ks","ns"))
      xmat = xmat %>% mutate(value = ifelse(!is.na(.data$proj),.data$proj,value)) %>%
        dplyr::select(-proj)
    }
    # update dynamic xmats
    if (any(xmat.coltypes$value == "dynamic")) {
      tmp.proj = xmat.dyn.fun(res, curr.areas, priors, xmat, xmat.proj)
      xmat = xmat %>%
        left_join(tmp.proj %>% rename("dyn" = "value") %>%
                    filter(ks %in% xmat.coltypes$ks[xmat.coltypes$value == "dynamic"]),by = c("ks","ns"))
      xmat = xmat %>% mutate(value = ifelse(!is.na(.data$dyn),.data$dyn,value)) %>%
        dplyr::select(-dyn)
    }
    # aggregate results over dataframes
    res.agg = data.frame(times = curr.time,res$out.res)
    if(curr.time==times[1]){out.res <- res.agg
    } else {
      out.res = bind_rows(out.res,res.agg)
    }
  }

  # Add back default columns
  # TODO

  return(list(out.res = out.res,out.solver = out.solver))
}
