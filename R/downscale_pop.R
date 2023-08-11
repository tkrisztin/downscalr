

#' Downscaling of population data
#'
#' @inheritParams downscale
#' @param targets A dataframe with columns times, pop.type (optional), and value (all targets >= 0)
#' @param xmat A dataframe of explanatory variables with columns ns, ks and value.   Defaults to NULL.
#' \code{xmat} and \code{betas} have to be provided for each
#' \code{pop.type} in \code{targets}.
#' @param betas A dataframe of coefficients with columns ks, pop.type (optional), and value. Defaults to NULL.
#' \code{xmat} and \code{betas} have to be provided for each
#' \code{pop.type} in \code{targets}.
#' @param options A list with solver options. Call \code{\link{downscale_control_pop}} for default options
#' and for more detail.
#'
#' @return A list containing
#' * \code{out.res} Dataframe with columns times, ns, pop.type & value (population allocation)
#' * \code{out.solver} A list of the solver output
#' * \code{ds.inputs} A list documenting all the downscale function inputs
#' @export
#'
#' @examples
#' dgp1 = sim_pop(1000)
#' res1 = downscale_pop(targets = dgp1$targets,xmat = dgp1$xmat,betas = dgp1$betas)
downscale_pop = function(targets,
                     times = NULL,
                     xmat,
                     betas,
                     xmat.coltypes = NULL,
                     xmat.proj = NULL,
                     xmat.dyn.fun = xmat.identity,
                     options = downscale_control_pop()) {
  pop.type = value = proj = ks = dyn = lu.to = NULL
  # Handle input checking
  err.txt = options$err.txt
  targets = complete_targets_pop(targets)
  xmat = complete_xmat(xmat)
  betas = complete_betas_pop(betas)

  complete_xmat.coltypes = complete_xmat.coltypes(xmat.coltypes, xmat)
  if (!is.null(xmat.proj)) {
    xmat.proj = complete_xmat.proj(xmat.proj)
  }
  err_check_inputs_pop(
    targets,
    xmat,
    betas,
    xmat.coltypes,
    xmat.proj,
    xmat.dyn.fun,
    err.txt
  )

  # save column types of xmat
  if (any(xmat.coltypes$value == "projected")) {
    proj.colnames = filter(xmat.coltypes, .data$value == "projected")$ks
  }
  if (any(xmat.coltypes$value == "dynamic")) {
    dyn.colnames = filter(xmat.coltypes, .data$value == "dynamic")$ks
  }

  # Set starting values
  curr.xmat = xmat

  if (is.null(times)) {
    times = unique(targets$times)
  }
  out.solver <- list()
  for (curr.time in times) {
    # Extract targets
    curr.targets = filter(targets, times == curr.time) %>% dplyr::select(-times)

    if (options$solve_fun == "solve_biascorr") {
      curr.options = options
      curr.options$err.txt = paste0(curr.time, " ", curr.options$err.txt)
      res = solve_biascorr.poisson(
        targets = curr.targets,
        xmat = curr.xmat,
        betas = betas,
        options = curr.options
      )
      out.solver[[as.character(curr.time)]] = res$out.solver
    }

    # update projected xmats
    if (any(xmat.coltypes$value == "projected")) {
      xmat = xmat %>%
        left_join(
          dplyr::filter(xmat.proj, times == curr.time) %>%
            dplyr::select(-times) %>% rename("proj" = "value") %>%
            mutate(ns = as.character(.data$ns)),
          by = c("ks", "ns")
        )
      xmat = xmat %>% mutate(value = ifelse(!is.na(.data$proj), .data$proj, value)) %>%
        dplyr::select(-proj)
    }
    # update dynamic xmats
    if (any(xmat.coltypes$value == "dynamic")) {
      tmp.proj = xmat.dyn.fun(res, NULL, NULL, xmat, xmat.proj)
      xmat = xmat %>%
        left_join(
          tmp.proj %>% rename("dyn" = "value") %>%  mutate(ns = as.character(.data$ns)) %>%
            filter(ks %in% xmat.coltypes$ks[xmat.coltypes$value == "dynamic"]),
          by = c("ks", "ns")
        )
      xmat = xmat %>% mutate(value = ifelse(!is.na(.data$dyn), .data$dyn, value)) %>%
        dplyr::select(-dyn)
    }
    # aggregate results over dataframes
    res.agg = data.frame(times = curr.time, res$out.res)
    if (curr.time == times[1]) {
      out.res <- res.agg %>%  mutate(ns = as.character(.data$ns))
    } else {
      out.res = bind_rows(out.res, res.agg %>%  mutate(ns = as.character(.data$ns)))
    }
  }

  # Remove default values columns
  if (any(out.res$times == PLCHOLD_T)) {
    out.res = out.res %>% dplyr::select(-times)
  }
  if (any(out.res$pop.type == PLCHOLD_POPT)) {
    out.res = out.res %>% dplyr::select(-pop.type)
  }

  ret <- list(out.res = out.res,
              out.solver = out.solver,
              ds.inputs = list(
                targets = targets,
                xmat = xmat,
                betas = betas,
                xmat.coltypes = xmat.coltypes,
                xmat.proj = xmat.proj,
                xmat.dyn.fun = xmat.dyn.fun,
                options = options))

  class(ret) = "downscalr"

  return(ret)
}
