








#' Fienberg rebalancing for two vectors with a given starting matrix
#'
#' @param start_mat An \eqn{n} by \eqn{n} matrix with non-negative elements
#' and row sums equal to 1.
#' @param target_from An \eqn{n} by 1 vector with values >=0.
#' @param target_to An \eqn{n} by 1 vector with values >=0. The sum of
#' \code{target_from} must be equal to the sum of \code{target_to}.
#'
#' Implements a specific version of the iterative proportional fitting algorithm.
#' Given an \eqn{n} by \eqn{n} matrix \code{start_mat}, the algorithm rebalances it
#' so that the \code{target_from * start_mat = target_to}. This is done by an iterative
#' fitting algorithm.
#'
#' @return The rebalanced \eqn{n} by \eqn{n} matrix \code{start_mat}.
#' @export
#'
#' @examples
#' set.seed(123)
#' n = 10
#' vec0 = runif(n) * 10
#' vec1 = runif(n); vec1 = vec1 / sum(vec1) * sum(vec0)
#' matA = diag(n) + 10^-8; matA = matA / rowSums(matA)
#'
#' res = fienberg(start_mat=matA,target_from=vec0,target_to=vec1)
fienberg = function(start_mat,target_from,target_to) {
  S = start_mat
  n = nrow(start_mat)
  k = ncol(start_mat)
  curr.optimVal = Inf
  Tol.optim = 10^-8
  curr.gain = Inf
  Tol.gain = 10^-10
  while (curr.gain > Tol.gain && curr.optimVal > Tol.optim) {
    S = start_mat * matrix(target_to / colSums(start_mat * target_from),n,k,byrow = T)
    start_mat = S * 1/rowSums(S)
    ttemp =sqrt(mean((target_to - colSums(start_mat * target_from))^2))
    curr.gain = curr.optimVal - ttemp
    curr.optimVal = ttemp
  }
  return(list(start_mat = start_mat,optimVal = curr.optimVal,gain = curr.gain))
}
