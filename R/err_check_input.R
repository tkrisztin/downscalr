
PLCHOLD_REGION = "NA_REGION"
PLCHOLD_LU = "NA_LU"
PLCHOLD_T = "NA_TIME"

#' Error check inputs
#'
#' @param targets A dataframe with columns times, lu.from (optional), lu.to and value (all targets >= 0)
#' @param areas  A dataframe of areas with columns lu.from (optional), ns and value, with all areas >= 0 and with sum(areas) >= sum(targets)
#' @param xmat A dataframe of explanatory variables with columns ks and value
#' @param betas A dataframe of coefficients with columns ks, lu.from (optional), lu.to & value
#' @param areas.update.fun function providing update for dynamic xmat columns, must take as arguments res, curr.areas, priors, xmat.proj, must dataframe with columns ns, lu.from & value defaults to areas.sum_to() which sums over lu.to
#' @param xmat.coltypes A dataframce with columns ks and string value, can be either "static", "dynamic", or "projected"
#' @param xmat.proj dataframe with columns times, ns, ks, must be present for each xmat.coltype specified as projected
#' @param xmat.dyn.fun function providing update for dynamic xmat columns, must take as arguments res, curr.areas, priors, xmat.proj must return ns x ks(dynamic) columns
#' @param priors A dataframe of priors (if no \code{betas} were supplied) with columns ns, lu.from (optional), lu.to (with priors >= 0)
#' @param restrictions A dataframe with columns ns, lu.from (optional), lu.to and value. Values must be zero or one. If restrictions are one, the MNL function is set to zero
#'
#' Internal function. Must throw errors, no return value if inputs do not match the specification. Handle all error checking here
#' Use this for all error checking of inputs.
#'
#' @keywords internal
err_check_inputs = function(targets,areas,xmat,betas,
                            areas.update.fun,xmat.coltypes,
                            xmat.proj,xmat.dyn.fun,
                            priors,restrictions,err.txt) {
  # check NA
  if (any(is.na(targets)) ||
      any(is.na(areas)) ||
      any(is.na(xmat)) ||
      any(is.na(betas))) {stop(paste0(err.txt,"Input contains NA values"))}
  if (!is.null(priors) && any(is.na(priors))) {stop(paste0(err.txt,"Input contains NA values"))}
  if (!is.null(restrictions) && any(is.na(restrictions))) {stop(paste0(err.txt,"Input contains NA values"))}
  if (!is.null(xmat.coltypes) && any(is.na(xmat.coltypes))) {stop(paste0(err.txt,"Input contains NA values"))}
  if (!is.null(xmat.proj) && any(is.na(xmat.proj))) {stop(paste0(err.txt,"Input contains NA values"))}

  # check rows
  if (nrow(targets) < 1) {stop(paste0(err.txt,"No observations in targets!"))}
  if (nrow(areas) < 1) {stop(paste0(err.txt,"No observations in areas!"))}
  if (nrow(xmat) < 1) {stop(paste0(err.txt,"No observations in xmat!"))}
  if (nrow(betas) < 1) {stop(paste0(err.txt,"No observations in betas!"))}
  if (!is.null(priors)) {
    if (nrow(priors) < 1) {stop(paste0(err.txt,"No observations in priors!"))}
  }
  if (!is.null(restrictions)) {
    if (nrow(restrictions) < 1) {stop(paste0(err.txt,"No observations in restrictions!"))}
  }
  if (!is.null(xmat.proj)) {
    if (nrow(xmat.proj) < 1) {stop(paste0(err.txt,"No observations in xmat.proj!"))}
  }
  if (!is.null(xmat.coltypes)) {
    if (nrow(xmat.coltypes) < 1) {stop(paste0(err.txt,"No observations in xmat.coltypes!"))}
  }

  # check correct names
  check_names = all(tibble::has_name(targets, c("lu.from","lu.to","times","value")))
  if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in targets."))}
  check_names = all(tibble::has_name(xmat, c("ks","ns","value")))
  if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in xmat"))}
  check_names = all(tibble::has_name(areas, c("lu.from","ns","value")))
  if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in areas"))}
  check_names = all(tibble::has_name(betas, c("ks","lu.from","lu.to","value")))
  if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in betas"))}
  if (!is.null(priors)) {
    check_names = all(tibble::has_name(priors, c("ns","lu.from","lu.to","value")))
    if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in priors"))}
  }
  if (!is.null(restrictions)) {
    check_names = all(tibble::has_name(restrictions, c("ns","lu.from","lu.to","value")))
    if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in restrictions"))}
  }
  if (!is.null(xmat.proj)) {
    check_names = all(tibble::has_name(xmat.proj, c("ns","ks","times","value")))
    if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in xmat.proj"))}
  }
  if (!is.null(xmat.coltypes)) {
    check_names = all(tibble::has_name(xmat.coltypes, c("ks","value")))
    if (!all(check_names)) {stop(paste0(err.txt,"Missing columns in xmat.proj"))}
  }

  # check values
  if (!all(targets >=0)) {
    targets[targets<0] = 0
    stop(paste0(err.txt,"Negative targets!"))
  }
  if (!all(areas$value >= 0)) {stop(paste0(err.txt,"All areas must be larger or equal to zero."))}
  if (!is.null(restrictions)) {
    if (!all(restrictions$value %in% c(0,1))) {stop(paste0(err.txt,"Restrictions must be 0 or 1"))}
  }
  if (!is.null(priors)) {
    if (!all(priors$value >=0)) {stop(paste0(err.txt,"Negative priors, must be >=0"))}
  }
  # check if all targets are covered as either betas or priors
  chck.names = targets  %>% dplyr::left_join(
    betas %>% group_by(.data$lu.from,.data$lu.to) %>% dplyr::summarize(n = n(),.groups = "keep"),by = c("lu.from", "lu.to"))
  chck.names$n[is.na(chck.names$n)] = 0
  if (!is.null(priors)) {
    if (any(paste0(priors$lu.from) == paste0(priors$lu.to))) {stop(paste0(err.txt,"Priors lu.from must be unequal to lu.to."))}
    chck.names = chck.names %>%
      left_join(
        priors %>% group_by(.data$lu.from,.data$lu.to) %>% dplyr::summarize(n2 = n(),.groups = "keep"),by = c("lu.from", "lu.to"))
    chck.names$n2[is.na(chck.names$n2)] = 0
    chck.names$n = chck.names$n + chck.names$n2
  }
  if (any(chck.names$n == 0)) {stop(paste0(err.txt,"Missing betas or priors for targets!"))}
  # check if all targets are below the areas
  err.check = targets %>% group_by(.data$times) %>%
    dplyr::summarise(total = sum(.data$value))
  if (any(sum( areas$value) < err.check$total )) {stop(paste0(err.txt,"Sum of areas larger than sum of targets."))}
  # check xmat.coltypes
  if (!all(xmat.coltypes$value %in% c("static","dynamic","projected"))) {
    stop(paste0(err.txt,"All xmat.coltypes values must be either static,dynamic, or projected"))}
  # for projected columns, make sure xmat.proj is supplied
  if (any(xmat.coltypes$value == "projected")) {
    if (is.null(xmat.proj)) {stop(paste0(err.txt,"Columns are specified as projected, but xmat.proj missing."))}
    chck.xmat = expand.grid(times = unique(targets$times),
                            ks = dplyr::filter(xmat.coltypes,.data$value == "projected")$ks) %>%
      left_join(
        xmat.proj %>% group_by(.data$times,.data$ks) %>% summarize(n = n(),.groups = "keep"),by = c("times", "ks")
      )
    if (any(is.na(chck.xmat$n))) {stop(paste0(err.txt,"xmat.proj must provide values for all times and projected ks."))}
  }
  if (any(xmat.coltypes$value == "dynamic") && is.null(xmat.dyn.fun)) {
    stop(paste0(err.txt,"Dynamic columns specified but missing xmat.dyn.fun for update."))}

  # check completeness
  # betas: Check if we have all ks
  ks = unique(xmat$ks)
  if (!all(ks %in% betas$ks)) {stop(paste0(err.txt,"Missing variables in betas (reference xmat)!"))}
  # betas: Check if we have all lu.from
  lu.from = unique(targets$lu.from)
  if (!all(lu.from %in% betas$lu.from)) {stop(paste0(err.txt,"Missing lu.from in betas (reference targets)!"))}
  # areas: check all lu present
  lu.from <- unique(targets$lu.from)
  if (!all(lu.from %in% areas$lu.from)) {stop(paste0(err.txt,"Missing lu.from in areas (reference targets)!"))}
  # xmat: Check if we have all ns
  ns = unique(areas$ns)
  if (!all(ns %in% xmat$ns)) {stop(paste0(err.txt,"Missing pixels in xmat (reference areas)!"))}
  # xmat: Check if we have all combinations
  expanded = xmat %>% tidyr::expand(.data$ks,.data$ns)
  if (nrow(expanded) != nrow(xmat)) {stop(paste0(err.txt,"Missing variables for pixels."))}
  if (!is.null(priors)) {
    # priors: Check if we have all ns
    ns = unique(areas$ns)
    if (!all(ns %in% priors$ns)) {stop(paste0(err.txt,"Missing pixels in priors (reference areas)!"))}
  }
  if (!is.null(restrictions)) {
    # restrictions: Check if we have all ns
    ns = unique(areas$ns)
    if (!all(ns %in% restrictions$ns)) {stop(paste0(err.txt,"Missing pixels in restrictions (reference areas)!"))}
  }
  if (!is.null(xmat.proj)) {
    # xmat.proj: Check if we have all ns
    ns = unique(areas$ns)
    if (!all(ns %in% xmat.proj$ns)) {stop(paste0(err.txt,"Missing pixels in xmat.proj (reference areas)!"))}
    # xmat.proj: Check if we have all combinations
    expanded = xmat.proj %>% tidyr::expand(.data$times,.data$ks,.data$ns)
    if (nrow(expanded) != nrow(xmat.proj)) {stop(paste0(err.txt,"Missing variables in xmat.projfor pixels."))}
  }
}
