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
#'
#' @examples
#' ## An example
areas.sum_to = function(res, curr.areas, priors, xmat, xmat.proj) {
  curr.areas = res$out.res %>%
    group_by(.data$ns,.data$lu.to) %>%
    summarize(value = sum(.data$value),.groups = "keep") %>%
    rename("lu.from" = "lu.to")
  # correct for small numerical mistakes
  if (min(curr.areas$value) < 0 && min(curr.areas$value) > -10^-10) {
    curr.areas$value[curr.areas$value < 0] = 0
  }
  return(curr.areas)
}
