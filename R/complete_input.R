#' Complete input targets
#'
#' @param targets Dataframe of targets, must have columns lu.to and value
#'
#' @return Completed dataframe with lu.from, times columns and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_targets = function(targets) {
  if (!tibble::has_name(targets,"lu.from")) {
    targets = cbind(lu.from = PLCHOLD_LU,targets)
  } else {
    if (any(targets$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  if (!tibble::has_name(targets,"times")) {
    targets = cbind(times = PLCHOLD_T,targets)
  } else {
    if (any(targets$times == PLCHOLD_T)) stop(paste0("The times ",PLCHOLD_T," is reserved, use another label."))
  }
  # Add all combinations
  targets = targets %>%
    dplyr::right_join(targets %>%
                        tidyr::expand(.data$lu.from,.data$lu.to,.data$times)) %>%
    tidyr::replace_na(list(value = 0))
  return(targets)
}

#' Complete input areas
#'
#' @param areas Dataframe of areas, must have columns ns and value
#'
#' @return Completed dataframe with lu.from and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_areas = function(areas) {
  if (!tibble::has_name(areas,"lu.from")) {
    areas = cbind(lu.from = PLCHOLD_LU,areas)
  } else {
    if (any(areas$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  # Add all combinations
  areas = areas %>%
    dplyr::right_join(areas %>%
                        tidyr::expand(.data$lu.from,.data$ns)) %>%
    tidyr::replace_na(list(value = 0))
  return(areas)
}

#' Complete input xmat
#'
#' @param xmat Dataframe of xmat, must have columns ns, ks and value
#'
#' @return Completed dataframe
#'
#' Internal function. Placeholder in case of needed additional completions.
#' @keywords internal
complete_xmat = function(xmat) {
  return(xmat)
}

#' Complete input betas
#'
#' @param betas Dataframe of betas, must have columns ks, lu.from and lu.to and value
#'
#' @return Completed dataframe, with added lu.from
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_betas = function(betas) {
  if (!tibble::has_name(betas,"lu.from")) {
    betas = cbind(lu.from = PLCHOLD_LU,betas)
  } else {
    if (any(betas$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  return(betas)
}

#' Complete input priors
#'
#' @param priors Dataframe of priors, must have columns ns, lu.to and value
#'
#' @return Completed dataframe with lu.from and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_priors = function(priors) {
  if (!tibble::has_name(priors,"lu.from")) {
    priors = cbind(lu.from = PLCHOLD_LU,priors)
  }  else {
    if (any(priors$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  # Add all combinations
  priors = priors %>%
    dplyr::right_join(priors %>%
                        tidyr::expand(.data$lu.from,.data$ns)%>%
                        tidyr::expand(.data$lu.to,.data$ns)) %>%
    tidyr::replace_na(list(value = 0))
  return(priors)
}

#' Complete input restrictions
#'
#' @param restrictions Dataframe of restrictions, must have columns ns, lu.to and value
#'
#' @return Completed dataframe with lu.from and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_restrictions = function(restrictions) {
  if (!tibble::has_name(restrictions,"lu.from")) {
    restrictions = cbind(lu.from = PLCHOLD_LU,restrictions)
  } else {
    if (any(restrictions$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  # Add all combinations
  restrictions = restrictions %>%
    dplyr::right_join(restrictions %>%
                        tidyr::expand(.data$lu.from,.data$ns)%>%
                        tidyr::expand(.data$lu.to,.data$ns)) %>%
    tidyr::replace_na(list(value = 0))
  return(restrictions)
}

#' Complete input xmat.proj
#'
#' @param xmat.proj Dataframe of xmat.proj, must have columns ns, ks and value
#'
#' @return Completed dataframe with times and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_xmat.proj = function(xmat.proj) {
  if (!tibble::has_name(xmat.proj,"times")) {
    xmat.proj = cbind(times = PLCHOLD_T,xmat.proj)
  } else {
    if (any(xmat.proj$times == PLCHOLD_T)) stop(paste0("The times ",PLCHOLD_T," is reserved, use another label."))
  }
  return(xmat.proj)
}
