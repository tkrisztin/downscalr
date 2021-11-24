#' Bayesian logit model with Pólya Gamma prior with MCMC
#'
#' @param X An n by k matrix of explanatory variables
#' @param Y An n by p matrix of dependent variables
#' @param baseline Baseline class for estimation. Parameters will be set to zero. Defaults to the p-th column.
#' @param ntot Total number of MCMC draws
#' @param nburn Burn-in draws for MCMC
#' @param A0 Prior for slope coefficients
#'
#'@details Performs MCMC estimation of a multinomial logit model as in Polson et al. (2012)
#'
#' @return A list containing
#' * \code{postb} A k x p x (ntot - nburn dimensions) array containing posterior draws of the slope coefficients.
#' * \code{marginal_fx} A k x p x (ntot - nburn dimensions) array containing posterior draws of marginal effects.
#' * \code{X, Y, baseline} The matrices of explanatory and dependent variables, as defined above and the baseline class.
#'
#' @examples
#'n = 100
#'p = 3
#'k = 2
#'X = cbind(1,matrix(rnorm(n*(k-1),0,2),n,k-1))
#'BETA <- matrix(sample(c(-3:3),k*p,replace=TRUE)  , k, p)
#'BETA[,p] = 0
#'Y = exp(X %*% BETA) / rowSums(exp(X %*% BETA))
#'res1 = mnlogit(X,Y)
#'print(BETA)
#'print(apply(res1$postb,c(1,2),mean))
mnlogit <- function(X, Y, baseline = ncol(Y),
                    ntot = 1000, nburn = 500, A0 = 10^4) {
  ### let us assign prior values
  # beta mean and variance are proper, but with high variance (= low prior information)
  n <- nrow(X)
  k <- ncol(X)
  p <- ncol(Y)
  pp <- (1:p)[-baseline]
  nn <- matrix(1, nrow(X), ncol(Y))

  beta_prior_mean <- matrix(0, k, p)
  beta_prior_var <- diag(k) * A0

  ### set-up the gibbs sampler
  # total number of draws
  niter <- ntot
  # retain only S2 and discard S1, with S = S1 + S2
  ndiscard <- nburn
  nretain <- niter - ndiscard
  # save the posterior draws here
  postb <- array(0, c(k, p, nretain))
  posty <- array(0, c(n, p, nretain))
  postom <- array(0, c(n, p, nretain))

  # starting values (won't matter after sufficient draws)
  # curr.beta = mvrnorm(1,beta_prior_mean,beta_prior_var)
  # curr.beta = solve(crossprod(X)) %*% crossprod(X,Y)
  curr.beta <- matrix(0, ncol = p, nrow = k)
  curr.beta[, p] <- 0
  curr.xb <- X %*% curr.beta
  curr.om <- matrix(0, n, p)

  # pre-calculate some terms for faster draws
  beta_prior_var_inv <- solve(beta_prior_var)
  kappa <- Y - nn / 2

  ### Gibbs sampling
  pb <- utils::txtProgressBar(min = 0, max = niter, style = 3)
  for (iter in 1:niter) {
    for (j in pp) {
      A <- rowSums(exp(X %*% curr.beta[, -j]))
      c.j <- log(A)
      eta.j <- X %*% curr.beta[, j] - c.j

      # sample omega
      curr.om[, j] <- BayesLogit::rpg(n, nn[, j], eta.j)

      # draw beta
      V <- solve(beta_prior_var_inv + t(X) %*% (X * curr.om[, j]))
      b <- V %*% (beta_prior_var_inv %*% beta_prior_mean[, j] + t(X) %*% (kappa[, j] + c.j * curr.om[, j]))
      curr.beta[, j] <- MASS::mvrnorm(1, b, V)
    }

    # we are past the burn-in, save the draws
    if (iter > ndiscard) {
      s <- iter - ndiscard
      postb[, , s] <- curr.beta
      curr.xb <- X %*% curr.beta
      posty[, , s] <- exp(curr.xb) / rowSums(exp(curr.xb))
      postom[, , s] <- curr.om
    }
    utils::setTxtProgressBar(pb,iter)
  }
  close(pb)
  ### calculate posterior mean of beta and sigma
  beta_mean_hat <- apply(postb, c(1, 2), mean)
  y_mean_hat <- apply(posty, c(1, 2), mean)

  marginal_fx <- array(0, c(k, p, nretain))
  dimnames(marginal_fx)[[1]] <- colnames(X)
  dimnames(marginal_fx)[[2]] <- colnames(Y)

  meanXs <- apply(X, c(2), mean)

  jjj <- 1
  for (jjj in 1:nretain) {
    ppp <- 1

    MU <- X %*% postb[, , jjj]
    pr <- exp(MU) / rowSums(exp(MU))

    for (ppp in 1:p) {
      bbb <- matrix(1, n, k) %*% diag(postb[, ppp, jjj])
      pr_bbb <- bbb
      for (kk in 1:k) {
        pr_bbb[, kk] <- rowSums(pr %*% diag(postb[kk, , jjj]))
      }

      partial1 <- pr[, ppp] * (bbb - pr_bbb)


      marginal_fx[, ppp, jjj] <- apply(partial1, c(2), mean)
    }
  }
  results <- list(
    postb = postb,
    marginal_fx = marginal_fx,
    X = X, Y = Y, baseline = baseline
  )
  return(results)
}
