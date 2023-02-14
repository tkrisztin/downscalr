


#' Convert land-use transition targets to gross land-use change transition targets
#'
#' @param targets.from A dataframe with columns lu and value (has to be numeric and >=0)
#' @param targets.to A dataframe with columns lu and value (has to be numeric and >=0)
#' @param restrictions A dataframe with columns lu.from (string), lu.to (string). Each row
#'  of the dataframe specifies which transitions will be set to zero in the returned
#'  land-use change transitions.
#' @param keep_areas Either "from", "to" or "both" (default). If \code{sum(targets.from$value)} is not
#' equal to \code{sum(targets.to$value)} this parameter specifies to which value the output should sum to.
#' If "both" is given and the sums differ, an artificial land use class \code{NODATA} is created.
#'
#' @details
#' When given two land-use targets, this function creates potential gross land-use
#' changes between them using Fienberg rebalancing. The transitions are forced towards
#' keeping as much targets in the same class as possible.
#'
#' The total set of land use classes is the combination of unique \code{lu} values from
#' \code{targets.from} and \code{targets.to}.  If \code{keep_areas = "both"} the
#' \code{NODATA} land use class is also created. Land use transitions
#' are allocated to the \code{NODATA} class if the sum of \code{targets.from} is not
#' equal to \code{targets.to}.
#'
#' @return A data.frame with columns lu.from, lu.to and value.
#' @export
#'
#' @examples
#'  targets.from = data.frame(lu = c("crop","grass"),value = c(10,5))
#'  targets.to = data.frame(lu = c("crop","grass","forest"),value = c(3,5,7))
#'  res = LU_to_LUC(targets.from = targets.from,targets.to = targets.to)
LU_to_LUC = function(targets.from, targets.to, restrictions = NULL,
                     keep_areas = "both") {
  #Error check of targets.from and targets.to
  # columns exists
  # column lu has strings only
  # column value is all numeric and >=0
  lu_classes = unique(c(targets.from$lu,targets.to$lu))

  if (keep_areas == "from") {
    targets.to$value = targets.to$value / sum(targets.to$value) * sum(targets.from$value)
  } else if (keep_areas == "to") {
    targets.from$value = targets.from$value / sum(targets.from$value) * sum(targets.to$value)
  } else if (keep_areas == "both") {
    if (sum(targets.from$value) != sum(targets.to$value)) {
      lu_classes = c(lu_classes,"NODATA")
    }
  } else {stop("keep_areas has to be ('from','to','both').")  }

  n = length(lu_classes)
  matA = diag(n) + 10^-8; matA = matA / rowSums(matA)
  targets.to = data.frame(lu = lu_classes) %>% left_join(targets.to,by = join_by(lu)) %>%
    replace_na(list(value = 0))
  targets.from = data.frame(lu = lu_classes) %>% left_join(targets.from,by = join_by(lu)) %>%
    replace_na(list(value = 0))

  res = fienberg(start_mat=matA,target_from=targets.from$value,target_to=targets.to$value)

  output = data.frame(lu.from = lu_classes,res$start_mat) %>%
    pivot_longer(cols = -1,names_to = "lu.to")

  return(output)
}
