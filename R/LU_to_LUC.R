


#' Convert land-use transition targets to gross land-use change transition targets
#'
#' @param targets.from A dataframe with columns lu and value (has to be numeric and >=0)
#' @param targets.to A dataframe with columns lu and value (has to be numeric and >=0)
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
LU_to_LUC = function(targets.from, targets.to, keep_areas = "both") {
  min_cutoff = 10^-8

  #Error check of targets.from and targets.to
  check_names = all(tibble::has_name(targets.from, c("lu","value")))
  if (!all(check_names)) {stop("Missing columns in targets.from")}
  check_names = all(tibble::has_name(targets.to, c("lu","value")))
  if (!all(check_names)) {stop("Missing columns in targets.to")}
  if (!all(targets.from$value >= 0)) {stop("All targets.from values must be larger or equal to zero.")}
  if (!all(targets.to$value >= 0)) {stop("All targets.to values must be larger or equal to zero.")}
  lu_classes = unique(
    c(as.character(targets.from$lu),
      as.character(targets.to$lu)))

  if (keep_areas == "from") {
    targets.to$value = targets.to$value / sum(targets.to$value) * sum(targets.from$value)
  } else if (keep_areas == "to") {
    targets.from$value = targets.from$value / sum(targets.from$value) * sum(targets.to$value)
  } else if (keep_areas == "both") {
    if (sum(targets.from$value) != sum(targets.to$value)) {
      lu_classes = c(lu_classes,"NODATA")
      if (sum(targets.from$value) > sum(targets.to$value)) {
        targets.to = bind_rows(targets.to,
                               data.frame(lu = "NODATA",
                                          value = sum(targets.from$value) - sum(targets.to$value)))
      } else {
        targets.from = bind_rows(targets.from,
                                 data.frame(lu = "NODATA",
                                            value = sum(targets.to$value) - sum(targets.from$value)))
      }
    }
  } else {stop("keep_areas has to be ('from','to','both').")  }

  n = length(lu_classes)
  lu = value = NULL # for code check and dplyr
  targets.to = data.frame(lu = lu_classes) %>% left_join(targets.to,by = c(("lu"))) %>%
    replace_na(list(value = 0))
  targets.from =
    data.frame(lu = lu_classes) %>% left_join(targets.from,by = c("lu")) %>%
    replace_na(list(value = 0))

  matA = diag(c(targets.from$value)) + min_cutoff
  matA = diag(n) + min_cutoff; matA = matA / rowSums(matA)
  colnames(matA) = rownames(matA) = lu_classes
  if (any(targets.from$value == 0)) {
    matA = matA[-which(targets.from$value==0),]
    targets.from = targets.from[-which(targets.from$value==0),]
  }
  if (any(targets.to$value == 0)) {
    matA = matA[,-which(targets.to$value==0)]
    targets.to = targets.to[-which(targets.to$value==0),]
  }

  res = fienberg(start_mat=matA,target_from=targets.from$value,target_to=targets.to$value)

  output = data.frame(lu.from = targets.from$lu,
                      targets.from$value * res$start_mat) %>%
    pivot_longer(cols = -1,names_to = "lu.to") %>%
    filter(value > min_cutoff)

  return(output)
}
