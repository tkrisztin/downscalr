#' Areas summing up function over land-use changes
#'
#' @param res Result from downscale
#' @param curr.areas Dataframe of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix dataframe
#' @param xmat.proj Projected x matrix
#'
#' @return curr.areas A
#' @export areas.sum_to
areas.sum_to = function(res, curr.areas, priors, xmat, xmat.proj) {
  curr.areas2 = res$out.res %>%
    group_by(.data$ns,.data$lu.to) %>%
    summarize(value = sum(.data$value),.groups = "keep") %>%
    rename("lu.from" = "lu.to")
  # correct for small numerical mistakes
  if (min(curr.areas2$value) < 0 && min(curr.areas2$value) > -10^-10) {
    curr.areas2$value[curr.areas2$value < 0] = 0
  }
  return(curr.areas2)
}

#' Areas identity function over land-uses
#'
#' @param res Result from downscale
#' @param curr.areas Dataframe of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix dataframe
#' @param xmat.proj Projected x matrix
#'
#' @return curr.areas A
#' @export areas.identity
areas.identity = function(res, curr.areas, priors, xmat, xmat.proj) {
  curr.areas = res$out.res
  # correct for small numerical mistakes
  if (min(curr.areas$value) < 0 && min(curr.areas$value) > -10^-10) {
    curr.areas$value[curr.areas$value < 0] = 0
  }
  return(curr.areas)
}
