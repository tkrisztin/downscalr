#' Downscaling of Land-Use (Change) Data
#'
#' Performs downscaling of land-use data over specified time steps using a range of inputs, including targets, areas, explanatory variables, and priors. It supports both bias correction and non-targeted downscaling methods.
#'
#' @param targets A dataframe with mandatory columns: `times` (character), `lu.to` (character), and `value` (numeric, all values >= 0). Optional column: `lu.from` (character). Represents the downscaling targets for each time step and land-use change.
#' @param start.areas A dataframe with starting areas. Includes mandatory columns: `ns` (character, representing grid IDs for downscaling) and `value` (numeric, all areas >= 0 and sum of areas >= sum of targets). Optional column: `lu.from` (character).
#' @param times A character vector of time steps for downscaling. The first time step must be present in `targets`. If `NULL`, times are derived from unique values in `targets`. Default is `NULL`.
#' @param xmat A dataframe with explanatory variables for econometric priors. Includes columns: `ns` (character), `ks` (character), and `value` (numeric). If `NULL`, a placeholder is used. Default is `NULL`.
#' @param betas A dataframe of coefficients for econometric priors. Includes columns: `ks` (character), `lu.to` (character), and `value` (numeric). Optional column: `lu.from` (character). If `NULL`, a placeholder is used. Default is `NULL`.
#' @param areas.update.fun A function providing an update for dynamic xmat columns. Takes as arguments `res`, `curr.areas`, `priors`, `xmat.proj` and must return a dataframe with columns `ns`, `ks`, and `value`. Defaults to `areas.sum_to()`, which sums over `lu.to`.
#' @param xmat.coltypes A vector `ks`, with each element being either "static", "dynamic", or "projected". Determines how different columns in `xmat` are treated during the downscaling process.
#' @param xmat.proj A dataframe with projections. Includes columns: `times` (character), `ns` (character), `ks` (character), and `value` (numeric). Required for each `xmat.coltype` specified as projected.
#' @param xmat.dyn.fun A function providing updates for dynamic xmat columns. Takes as arguments `res`, `curr.areas`, `priors`, `xmat.proj` and must return a dataframe with `ns x ks(dynamic)` columns.
#' @param priors A dataframe with exogenous priors. Includes columns: `times` (character, optional), `ns` (character), `lu.from` (character, optional), `lu.to` (character), and `value` (numeric, >= 0). An optional `weight` column (numeric, 0 <= weight <= 1) can be supplied to adjust the influence of exogenous priors.
#' @param restrictions A dataframe with restrictions. Includes columns: `ns` (character), `lu.from` (character, optional), `lu.to` (character), and `value` (numeric). Values must be either zero or one, indicating whether the MNL function should be set to zero for certain combinations.
#' @param options A list of solver options. Use `\link{downscale_control}` to obtain default options and for more detailed information.
#'
#' @details The function integrates various data inputs to match `p` targets using either projections from an MNL-type model or exogenous priors. Appropriate input validation and preprocessing are performed before downscaling.
#'
#' @return A list containing three elements:
#' * `out.res`: A dataframe with columns `times`, `ns`, `lu.from`, `lu.to`, and `value` (area allocation).
#' * `out.solver`: A list detailing the solver output.
#' * `ds.inputs`: A list documenting all the inputs used in the downscaling function.
#'
#' @export downscale
#' @import nloptr
#' @import tidyr
#' @import dplyr
#' @import tibble
#'
#' @examples
#' require(dplyr)
#' require(tidyr)
#' require(tibble)
#' betas = NULL
#' for (jj in unique(argentina_luc$lu.from)) {
#'  Y = dplyr::filter(argentina_luc,lu.from == jj & Ts == 2000) %>%
#'    pivot_wider(names_from = lu.to)
#'  X = argentina_df$xmat %>% tidyr::pivot_wider(names_from = "ks") %>%
#'    dplyr::arrange(match(ns,Y$ns))
#'  Y = Y %>% dplyr::select(-c(lu.from,Ts,ns))
#'  X = X %>% dplyr::select(-c(ns))
#'  res1 <- mnlogit(as.matrix(X), as.matrix(Y),baseline = which(colnames(Y) == jj),
#'           niter = 3,nburn = 2)
#'  betas = betas %>% dplyr::bind_rows(
#'   apply(res1$postb, c(1, 2), mean) %>%
#'   as.data.frame() %>% tibble::rownames_to_column("ks") %>%
#'   pivot_longer(cols = -c(1),names_to = "lu.to") %>%
#'   dplyr::mutate(lu.from = jj,.before="lu.to")
#'  )
#' }
#' ns = unique(argentina_df$lu_levels$ns)
#' priors = data.frame(ns = as.character(ns),lu.from="Cropland",
#'             lu.to="Forest",value = as.numeric(runif(length(ns))))
#' res1 = downscale(targets = argentina_FABLE %>% dplyr::filter(times == "2010"),
#'          start.areas = argentina_df$lu_levels,
#'          xmat = argentina_df$xmat,
#'          betas = betas %>% dplyr::filter(lu.from!="Cropland" | lu.to!="Forest"),
#'          priors = priors)
#'
#'  dgp1 = sim_luc(1000,tt = 3)
#'  res1 = downscale(targets = dgp1$targets,start.area = dgp1$start.areas,
#'         xmat = dgp1$xmat,betas = dgp1$betas,times = c(1:3))
downscale = function(targets,
                     start.areas,
                     times = NULL,
                     xmat = NULL,
                     betas = NULL,
                     areas.update.fun = areas.sum_to,
                     xmat.coltypes = NULL,
                     xmat.proj = NULL,
                     xmat.dyn.fun = xmat.sum_to,
                     priors = NULL,
                     restrictions = NULL,
                     options = downscale_control()) {
  lu.from = value = proj = ks = dyn = lu.to = NULL
  # Handle input checking
  err.txt = options$err.txt
  targets = complete_targets(targets)
  start.areas = complete_areas(start.areas)
  target_area_check(targets, start.areas)
  if (is.null(xmat)) {
    xmat = data.frame(ns = unique(start.areas$ns),
                      ks = PLCHOLD_K,
                      value = 0)
  }
  xmat = complete_xmat(xmat)
  if (is.null(betas)) {
    betas = data.frame(
      lu.from = paste0(PLCHOLD_LU, 1),
      lu.to = paste0(PLCHOLD_LU, 2),
      ks = unique(xmat$ks),
      value = 0
    )
  }
  betas = complete_betas(betas)
  xmat.coltypes = complete_xmat.coltypes(xmat.coltypes, xmat)
  if (!is.null(priors)) {
    priors = complete_priors(priors, xmat, targets)
  }
  if (!is.null(restrictions)) {
    restrictions = complete_restrictions(restrictions, xmat)
  }
  if (!is.null(xmat.proj)) {
    xmat.proj = complete_xmat.proj(xmat.proj)
  }
  err_check_inputs(
    targets,
    start.areas,
    xmat,
    betas,
    areas.update.fun,
    xmat.coltypes,
    xmat.proj,
    xmat.dyn.fun,
    priors,
    restrictions,
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
  curr.areas = start.areas
  curr.xmat = xmat
  curr.restrictions = restrictions

  if (is.null(times)) {
    times = unique(targets$times)
  }
  out.solver <- list()
  for (curr.time in times) {
    # Extract targets
    curr.targets = filter(targets, times == curr.time) %>% dplyr::select(-times)

    # Extract priors
    if (any(priors$times == curr.time)) {
      curr.priors = filter(priors, times == curr.time) %>% dplyr::select(-times)
    } else {
      curr.priors = NULL
    }

    # BUGFIX TK: added check for missing targets for every time iteration
    #   This was before only done for the initial period, resulting in missed
    #     land cover.
    # check if all lu in curr.areas in curr.targets
    if (!all(unique(curr.areas$lu.from) %in% unique(curr.targets$lu.from))) {
      # if not save them to a table and add them back manually
      missing_luc = unique(curr.areas$lu.from)[!unique(curr.areas$lu.from) %in% unique(curr.targets$lu.from)]
      missing_luc = dplyr::filter(curr.areas, lu.from %in% missing_luc) %>%
        mutate(lu.to = lu.from)
    } else {
      missing_luc = NULL
    }

    if (options$solve_fun == "solve_biascorr") {
      curr.options = options
      curr.options$err.txt = paste0(curr.time, " ", curr.options$err.txt)
      res = solve_biascorr.mnl(
        targets = curr.targets,
        areas = curr.areas,
        xmat = curr.xmat,
        betas = betas,
        priors = curr.priors,
        restrictions = curr.restrictions,
        options = curr.options
      )
      out.solver[[as.character(curr.time)]] = res$out.solver
    } else if (options$solve_fun == "solve_notarget") {
      curr.options = options
      curr.options$err.txt = paste0(curr.time, " ", curr.options$err.txt)
      res = solve_notarget.mnl(
        targets = curr.targets,
        areas = curr.areas,
        xmat = curr.xmat,
        betas = betas,
        restrictions = curr.restrictions,
        options = curr.options
      )
      out.solver[[as.character(curr.time)]] = res$out.solver
    }

    # add not covered land-uses (because no targets exist)
    if (!is.null(missing_luc)) {
      res$out.res = res$out.res %>% bind_rows(missing_luc)
    }

    # update curr.area
    curr.areas = areas.update.fun(res, curr.areas, priors, xmat.proj)

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
      tmp.proj = xmat.dyn.fun(res, curr.areas, priors, xmat, xmat.proj)
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
  if (any(out.res$lu.from == PLCHOLD_LU)) {
    out.res = out.res %>% dplyr::select(-lu.from)
  }
  if (any(out.res$lu.to == PLCHOLD_LU)) {
    out.res = out.res %>% dplyr::filter(lu.to != PLCHOLD_LU)
  }

  ret <- list(
    out.res = out.res,
    out.solver = out.solver,
    ds.inputs = list(
      targets = targets,
      start.areas = start.areas,
      xmat = xmat,
      betas = betas,
      areas.update.fun = areas.update.fun,
      xmat.coltypes = xmat.coltypes,
      xmat.proj = xmat.proj,
      xmat.dyn.fun = xmat.dyn.fun,
      priors = priors,
      restrictions = restrictions,
      options = options
    )
  )

  class(ret) = "downscalr"

  return(ret)
}

