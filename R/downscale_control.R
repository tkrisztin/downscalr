#' Set options for downscaling
#'
#' This function can be called to specify additional options to the downscalR solver.
#' The user can change the solver function, specify solver behavior, stopping conditions, and
#' change how exogeneous and econometric prior specifications are mixed.
#'
#' @param solve_fun The name of the downscaling function to use. Has to be a
#' valid R function. Defaults to \code{"solve_biascorr"}.
#' @param algorithm Solver algorithm (see the \code{\link[nloptr]{nloptr}} package for documentation and more detail).
#' @param xtol_rel Relative tolerance of solver.
#' @param xtol_abs Absolute solver tolerance.
#' @param maxeval Maximum evaluation of solver.
#' @param MAX_EXP Numerical cutoff for MNL function.
#' @param cutoff Optional cutoff to avoid MNL values close to zero.
#' @param err.txt Error text for caller identification (used for debugging)
#' @param max_diff If difference to targets is larger, redo the estimation (helps to avoid convergence errors)
#' @param redo Maximum number of repeats
#' @param prior_weights Scalar between 0 and 1; if both betas and priors are specified this value gives the weight of the priors.
#'
#' @return List with default options for bias correction solver
#'
#' @details
#' A user specified solve_fun function has to be able to process the following input parameters:
#' \code{targets}, \code{areas}, \code{xmat}, \code{betas}, \code{priors}, \code{restrictions}, and
#' \code{options} (the output of this function). The output should be a list containing two parameters:
#' a dataframe \code{out.res} containing the downscaled output and a list \code{out.solver} detailing
#' the solver results.
#'
#' @export downscale_control
#'
#' @examples
#' opts1 = downscale_control()
downscale_control = function(solve_fun = "solve_biascorr",algorithm = "NLOPT_LN_SBPLX",
                                  xtol_rel = 1.0e-20,xtol_abs = 1.0e-20,maxeval = 1600,
                                  MAX_EXP = log(.Machine$double.xmax),cutoff = 0,
                                  redo = 2,max_diff = 10^-8,err.txt = "",prior_weights = 0) {
  if (!solve_fun %in% c("solve_biascorr")) {stop("solve_fun not correctly specified.")}
  if (prior_weights<0 || prior_weights >1) {stop("prior_weights must be >=0 and <=1.")}
  return(list(solve_fun = solve_fun,algorithm = algorithm,xtol_rel = xtol_rel,xtol_abs = xtol_abs,
              maxeval = maxeval,MAX_EXP = MAX_EXP,cutoff = cutoff,redo = redo,
              max_diff = max_diff,err.txt = err.txt, prior_weights = prior_weights
  ))
}
