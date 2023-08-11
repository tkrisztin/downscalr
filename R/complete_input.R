#' Complete input targets
#'
#' @param targets Dataframe of targets, must have columns lu.to and value
#'
#' @return Completed dataframe with lu.from, times columns and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_targets = function(targets) {
  lu.from = lu.to = ns = NULL
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
                        tidyr::expand(nesting(lu.from,lu.to),.data$times),
                      by = c("times", "lu.from", "lu.to")) %>%
    filter(as.character(.data$lu.from) != as.character(.data$lu.to)) %>%
    tidyr::replace_na(list(value = 0))
  targets = dplyr::arrange(targets,.data$lu.from,.data$lu.to,.data$times)
  return(targets)
}

#' Complete input targets for population downscaling
#'
#' @param targets Dataframe of targets, must have columns pop.type and value
#'
#' @return Completed dataframe with pop.type, times columns and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_targets_pop = function(targets) {
  pop.type = ns = times = NULL
  if (!tibble::has_name(targets,"times")) {
    targets = cbind(times = PLCHOLD_T,targets)
  } else {
    if (any(targets$times == PLCHOLD_T)) stop(paste0("The times ",PLCHOLD_T," is reserved, use another label."))
  }
  # Add all combinations
  targets = targets %>%
    dplyr::right_join(targets %>%
                        tidyr::expand(pop.type,times),
                      by = c("times", "pop.type"))  %>%
    tidyr::replace_na(list(value = 0))
  targets = dplyr::arrange(targets,pop.type,times)
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
                        tidyr::expand(.data$lu.from,.data$ns),by = c("ns", "lu.from")) %>%
    tidyr::replace_na(list(value = 0))
  areas = dplyr::arrange(areas,.data$lu.from,.data$ns)
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
  xmat = dplyr::arrange(xmat,.data$ks,.data$ns) %>%
    dplyr::right_join(xmat %>%
                        tidyr::expand(.data$ns,.data$ks),
                      by = c("ns", "ks")) %>%
    tidyr::replace_na(list(value = 0)) %>%
    dplyr::mutate(ns = as.character(.data$ns), ks = as.character(.data$ks))

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
  betas = dplyr::arrange(betas,.data$lu.from,.data$lu.to,.data$ks)
  return(betas)
}

#' Complete input betas for population downscaling
#'
#' @param betas Dataframe of betas, must have columns ks, lu.from and lu.to and value
#'
#' @return Completed dataframe, with added lu.from
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_betas_pop = function(betas) {
  if (!tibble::has_name(betas,"pop.type")) {
    betas = cbind(pop.type = PLCHOLD_POPT,betas)
  } else {
    if (any(betas$pop.type == PLCHOLD_POPT)) stop(paste0("The pop.type ",PLCHOLD_POPT," is reserved, use another name."))
  }
  betas = dplyr::arrange(betas,.data$pop.type,.data$ks)
  return(betas)
}

#' Complete input priors
#'
#' @param priors Dataframe of priors, must have columns ns, lu.to and value
#' @param xmat Dataframe of xmat, must have columns ns, ks and value
#' @param xmat Dataframe of targets, must have columns lu.from and times
#'
#' @return Completed dataframe with times, lu.from, lu.to, weight and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_priors = function(priors,xmat,targets) {
  # lu.from & lu.to defined to fool the package checker with dplyr namebindings
  #   (.data$ does not work in nested function)
  lu.from = lu.to = ns =  NULL

  if (!tibble::has_name(priors,"lu.from")) {
    priors = cbind(lu.from = PLCHOLD_LU,priors)
  }  else {
    if (any(priors$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  if (tibble::has_name(priors,"weight")) {
    if (any(priors$weight < 0) || any(priors$weight >1 )) {
      stop("All prior weights must be between 0 and 1.")
    }
  } else {
    priors = cbind(priors,weight = 1)
  }
  if (!tibble::has_name(priors,"times")) {
    priors = priors %>%
      right_join(priors %>%
                   expand(times = unique(targets$times),lu.from,lu.to,ns),
                 by= c("ns","lu.from","lu.to"))
  }
  #Add all combinations
  priors = priors %>%  dplyr::right_join(
    priors %>% right_join(dplyr::select(xmat,ns) %>% distinct(),by= c("ns")) %>%
    tidyr::expand(.data$ns,nesting(lu.from,lu.to))  %>% filter(!is.na(lu.from) & !is.na(lu.to)),
                            by= c("ns", "lu.from", "lu.to")) %>%
    filter(lu.from != lu.to) %>%
    tidyr::replace_na(list(value = 0, weight = 0))
  priors = dplyr::arrange(priors,.data$lu.from,.data$lu.to,.data$ns)
  return(priors)
}

#' Complete input restrictions
#'
#' @param restrictions Dataframe of restrictions, must have columns ns, lu.to and value
#' @param xmat Dataframe of xmat, must have columns ns, ks and value
#'
#' @return Completed dataframe with lu.from and all combinations
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_restrictions = function(restrictions,xmat) {
  if (!tibble::has_name(restrictions,"lu.from")) {
    restrictions = cbind(lu.from = PLCHOLD_LU,restrictions)
  } else {
    if (any(restrictions$lu.from == PLCHOLD_LU)) stop(paste0("The lu.from ",PLCHOLD_LU," is reserved, use another name."))
  }
  # Add all combinations
  # lu.from & lu.to defined to fool the package checker with dplyr namebindings
  #   (.data$ does not work in nested function)
  lu.from = lu.to = ns = NULL
  restrictions = restrictions %>%  dplyr::right_join(
    restrictions %>% right_join(select(xmat,ns) %>% distinct(),by= c("ns")) %>%
      tidyr::expand(.data$ns,nesting(lu.from,lu.to))  %>% filter(!is.na(lu.from) & !is.na(lu.to))
    ,by= c("ns", "lu.from", "lu.to")) %>%
    tidyr::replace_na(list(value = 0))
  restrictions = dplyr::arrange(restrictions,.data$lu.from,.data$lu.to,.data$ns)
  return(restrictions)
}

#' Complete input xmat.coltypes
#'
#' @param xmat.coltypes Dataframe of xmat.coltypes, must have columns ks and value;
#' if it is NULL, this will be created as all static
#' deprecated: can also be a vector
#'
#' @return Completed dataframe with ks and values
#'
#' Internal function. Adds missing columns and completes potential sparse dataframes.
#' @keywords internal
complete_xmat.coltypes = function(xmat.coltypes,xmat) {
  if (is.null(xmat.coltypes)) {
    xmat.coltypes = data.frame(ks = unique(xmat$ks),
                               value = "static")
  }
  # Handle legacy xmat.coltypes with warning
  if (is.vector(xmat.coltypes)) {
    warning("Depreciated: xmat.coltypes should be a dataframe with columns ks and value, not a vector.\n
            Attempting to convert to appropriate coltype.\n
            This will be removed in future versions.")
    if (is.null(names(xmat.coltypes))) stop("No names in xmat.coltypes")
    ks = unique(xmat$ks)
    xmat.coltypes2 = data.frame(ks = ks,value = xmat.coltypes)
    xmat.coltypes2$value[xmat.coltypes2$ks %in% names(xmat.coltypes)] = xmat.coltypes
    xmat.coltypes = xmat.coltypes2
  }
  return(xmat.coltypes)
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
  xmat.proj = dplyr::arrange(xmat.proj,.data$times,.data$ks,.data$ns)
  return(xmat.proj)
}

#' Check for targets to area mismatch over timesteps
#'
#' @param targets Dataframe of targets, must have columns lu.to and value
#' @param areas Dataframe of areas, must have columns ns and value#' @
#'
#' Internal function. Checks if targets are fullfillable over time before downscaling.
#' @keywords internal
#'
target_area_check <- function(targets, areas){

  temp_curr.lu_levels = temp_lu.from.targets = temp_lu.to.targets = lu.from = value = times = value_net = lu.to = value.lu.to = value.lu.from = . = NULL

  timesteps <- base::sort(base::unique(targets$times))

  jjj <- timesteps[1]
  for(jjj in timesteps){

    if(jjj==timesteps[1]) {
      temp_curr.lu_levels <- areas %>%
        dplyr::group_by(lu.from) %>%
        dplyr::summarize(value=sum(value))
    }

    temp_lu.from.targets <- targets %>%
      dplyr::filter(times == c(jjj)) %>%
      dplyr::group_by(lu.from) %>%
      dplyr::summarize(value.lu.from=sum(value))

    temp_curr.lu_levels <- temp_curr.lu_levels %>%
      dplyr::left_join(temp_lu.from.targets, by=c("lu.from")) %>%
      base::replace(is.na(.),0) %>%
      dplyr::mutate(value=value-value.lu.from) %>%
      dplyr::select(!value.lu.from)

    temp_lu.to.targets <- targets %>%
      dplyr::filter(times == c(jjj), value!=0) %>%
      dplyr::left_join(temp_curr.lu_levels %>%
                         dplyr::rename(value_net="value"), by="lu.from") %>%
      dplyr::group_by(lu.from) %>%
      dplyr::mutate(value=ifelse(value_net<0, value+(value_net/n()), value)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(lu.to) %>%
      dplyr::summarize(value.lu.to=sum(value))

    if(any(temp_curr.lu_levels$value<0)) {
      base::cat(paste0("Total area to target mismatch in timestep ", jjj," in the following class(es):"),
              knitr::kable(temp_curr.lu_levels %>%
                             dplyr::filter(value<0) %>%
                             dplyr::rename(class="lu.from")), sep="\n")
    }

    temp_curr.lu_levels <- temp_curr.lu_levels %>%
      dplyr::mutate(value=ifelse(value<0,0,value)) %>%
      dplyr::left_join(temp_lu.to.targets, by=c("lu.from"="lu.to")) %>%
      base::replace(is.na(.),0) %>%
      dplyr::mutate(value=value+value.lu.to) %>%
      dplyr::select(!value.lu.to)

  }

}


