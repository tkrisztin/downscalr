#' xmat summing up function over land-use changes
#'
#' @param res Result from downscale
#' @param curr.areas Dataframe of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix dataframe
#' @param xmat.proj Projected x matrix
#'
#' @return curr.xmat Dynamically updated xmat columns
#' @export xmat.sum_to
xmat.sum_to = function(res, curr.areas, priors, xmat, xmat.proj) {
  curr.xmat = res$out.res %>%
    group_by(.data$ns,.data$lu.to) %>%
    summarize(value = sum(.data$value),.groups = "keep") %>%
    rename("ks" = "lu.to")
  # correct for small numerical mistakes
  if (min(curr.xmat$value) < 0 && min(curr.xmat$value) > -10^-10) {
    curr.xmat$value[curr.xmat$value < 0] = 0
  }
  return(curr.xmat)
}


#' xmat identity function
#'
#' This function returns the current projections from res
#'
#' @param res Result from downscale
#' @param curr.areas Dataframe of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix dataframe
#' @param xmat.proj Projected x matrix
#'
#' @return curr.xmat Dynamically updated xmat columns
#' @export xmat.identity
xmat.identity = function(res, curr.areas, priors, xmat, xmat.proj) {
  curr.xmat = res$out.res
  # correct for small numerical mistakes
  if (min(curr.xmat$value) < 0 && min(curr.xmat$value) > -10^-10) {
    curr.xmat$value[curr.xmat$value < 0] = 0
  }
  return(curr.xmat)
}

