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
  mu = mu / rowSums(cbind(1,mu)) * areas
  if (!is.null(restrictions)) {
    mu[restrictions == 1] = 0
  }
  mu[mu<cutoff] = 0
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
sqr_diff.mnl = function(x,mu,areas,targets,restrictions,cutoff = 0) {
  x.mu = x[1:length(targets)]
  mu = mu.mnl(x.mu,mu,areas,restrictions,cutoff)
  return( sum((colSums(mu) - targets)^2) )
}
