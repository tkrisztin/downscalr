#' Multinomial logit mean calculation with restrictions and cutoff
#'
#' @param x Bias correction parameters
#' @param mu Means of each class
#' @param areas Vector of areas
#' @param restrictions Optional restrictions. Binary, if 1, the log-odds are set to zero
#' @param cutoff Optional, defaults to zero. If set, all log-odds below cutoff are set to zero
#'
#' @return Returns the vector of log-odds times the areas
#'
#' @keywords internal
mu.mnl = function(x,mu,areas,restrictions = NULL,cutoff = 0) {
  mu = mu * matrix(x,nrow(mu),ncol(mu),byrow = TRUE)
  if (!is.null(restrictions)) {
    mu[restrictions == 1] = 0
  }
  mu[mu<cutoff] = 0
  mu = mu / rowSums(cbind(1,mu)) * areas
  return(mu)
}

#' Squared differences optimization function with multinomial logit
#'
#' @param x Bias correction parameters
#' @param mu Means of each class
#' @param areas Vector of areas
#' @param targets Targets to match
#' @param restrictions Optional restrictions. Binary, if 1, the log-odds are set to zero
#' @param cutoff Optional, defaults to zero. If set, all log-odds below cutoff are set to zero
#'
#' @return Squared differences of optimized values and targets
#'
#' @keywords internal
sqr_diff.mnl = function(x,mu,areas,targets,restrictions = NULL,cutoff = 0) {
  x.mu = x[1:length(targets)]
  mu = mu.mnl(x.mu,mu,areas,restrictions,cutoff)
  return( sum((colSums(mu) - targets)^2) )
}

#' Provides analytical gradients for the squared differences optimization function with multinomial logit
#'
#' @inheritParams sqr_diff.mnl
#'
#' @return Vector of length \eqn{J} containing the gradients for each value of \code{x}.
#'
#' Does not work for \code{cutoff} \eqn{> 0} or with restrictions. Generates a warning if these are supplied.
#'
#' @keywords internal
grad_sqr_diff.mnl = function(x,mu,areas,targets,restrictions = NULL,cutoff = 0) {
  if (!is.null(restrictions) || cutoff > 0) {
    stop("Gradient optimization not yet implemented for restrictions or cutoff > 0. \n
            Consider switching to a gradient free optimiser (see the nloptr documentation for details).")
  }
  n = nrow(mu); J = length(x)
  mu = cbind(1,mu * matrix(exp(x) ,n,J,byrow = TRUE))
  theta = mu / c(rowSums(mu)); theta = theta[,-1,drop=FALSE]
  diff_targets = targets - colSums(areas * theta)

  grad_x = rep(0,J)
  for (jj in 1:J) {
    own_d = sum( areas * theta[,jj] * (1 - theta[,jj])) * diff_targets[jj]
    other_d = 0
    for (jj2 in c(1:J)[-jj]) {
      other_d = other_d + sum(areas * theta[,jj] * theta[,jj2]) * diff_targets[jj2]
    }
    grad_x[jj] = -2*own_d + 2*other_d
  }
  return(grad_x)
}
